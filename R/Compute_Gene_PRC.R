#' Compute PR-AUC and Optionally Plot Precision-Recall Curve for a Single Gene
#'
#' @description
#' Evaluates the discriminatory power of a single gene in separating a
#' user-defined positive cell group from the rest, using the Area Under the
#' Precision-Recall Curve (PR-AUC), reported here as Average Precision (AP).
#' Unlike ROC-AUC, which is insensitive to class balance, PR-AUC summarises the
#' trade-off between precision and recall and is therefore better suited to
#' highly imbalanced settings—such as identifying a rare cell type among
#' thousands of cells—where the positive class is the minority.
#' Several scoring strategies and optional cell subsetting are provided to
#' handle dropout noise and sparse expression typical of single-cell data.
#' If \code{plot = TRUE} (the default), a publication-ready precision-recall
#' curve is automatically displayed via \pkg{ggplot2}.
#'
#' @param seurat_obj A Seurat object containing single-cell expression data.
#' @param gene A single character string specifying the gene to evaluate. Must
#'   be present in the rownames of the chosen assay.
#' @param group_col Character string giving the name of a column in
#'   \code{seurat_obj@meta.data} that defines cell groups (e.g.,
#'   \code{"seurat_clusters"}).
#' @param group_label The value in \code{group_col} that defines the positive
#'   class (e.g., \code{1} for cluster 1). All other cells are treated as
#'   negatives.
#' @param assay Character string specifying which assay to use. Default is
#'   \code{"RNA"}.
#' @param method Scoring method: \code{"raw"} (raw expression, optionally
#'   truncated by \code{min_expression}) or \code{"rank"} (expression ranks,
#'   which is robust to dropout and does not require \code{min_expression}).
#' @param min_expression Numeric threshold for expression truncation (only
#'   used when \code{method = "raw"}). Values below this threshold are set to
#'   zero. Default is \code{NULL} (no truncation).
#' @param keep_expression_above Optional numeric threshold. If set, only cells
#'   with expression greater than this value are retained for PR-AUC
#'   computation. \strong{Note:} this changes the population over which PR-AUC
#'   is evaluated, shifting the interpretation from "overall discrimination
#'   among all cells" to "discrimination among cells that express the gene
#'   above this level". Because the no-skill baseline is the prevalence of the
#'   positive class, subsetting also changes the baseline. Use with caution and
#'   always compare with the result without subsetting. Default is
#'   \code{NULL} (all cells retained).
#' @param plot Logical indicating whether to create and automatically display
#'   a ggplot2 precision-recall curve. Default is \code{TRUE}.
#' @param plot_title Character string used as the title of the PR plot.
#'   Default is \code{"Gene Precision-Recall Curve | SlimR"}.
#' @param line_color Colour of the PR curve. Default is \code{"darkmagenta"}.
#' @param line_size Numeric value for the thickness of the PR curve line.
#'   Default is \code{1}.
#'
#' @return Invisibly, a list with the following elements:
#' \describe{
#'   \item{PR_AUC}{Numeric value (between 0 and 1) of the area under the
#'     precision-recall curve, computed as Average Precision.}
#'   \item{Baseline}{Numeric value of the no-skill baseline, equal to the
#'     prevalence of the positive class. A PR-AUC above this value indicates
#'     better-than-random discrimination; a value close to the baseline
#'     indicates no discriminatory power.}
#'   \item{pr_data}{A data.frame with columns \code{recall} and
#'     \code{precision} that can be used for custom plotting.}
#'   \item{predictions}{Numeric vector of the (possibly truncated) expression
#'     values used as prediction scores.}
#'   \item{labels}{Logical vector indicating whether each cell belongs to the
#'     positive class (\code{TRUE}) or not (\code{FALSE}).}
#'   \item{pr_plot}{If \code{plot = TRUE}, a \pkg{ggplot2} object displaying
#'     the precision-recall curve; otherwise \code{NULL}.}
#' }
#'
#' @details
#' The difference between ROC-AUC and PR-AUC is important when interpreting
#' the results:
#' \itemize{
#'   \item \strong{ROC-AUC} plots the True Positive Rate against the False
#'     Positive Rate. It is insensitive to class imbalance and has a fixed
#'     no-skill baseline of 0.5.
#'   \item \strong{PR-AUC} plots Precision (Positive Predictive Value) against
#'     Recall (Sensitivity). It focuses on the performance on the positive
#'     class and has a no-skill baseline equal to the prevalence of the
#'     positive class. For rare cell populations this baseline is close to
#'     zero, so PR-AUC provides a more informative and more stringent measure
#'     of discriminatory power than ROC-AUC.
#' }
#'
#' PR-AUC is computed here as Average Precision (AP), i.e. the mean of the
#' precision values obtained at each positive cell, weighted by the
#' corresponding increase in recall:
#' \deqn{AP = \sum_{k \in \text{positives}} P(k) \, \Delta R(k)}
#' where the sum runs over all cells in the positive class. This definition is
#' equivalent to \code{sklearn.metrics.average_precision_score} and avoids the
#' systematic underestimation that a naive trapezoidal integration would
#' introduce on the non-monotonic PR curve. The no-skill \code{Baseline} is
#' returned alongside \code{PR_AUC} because the two must be interpreted
#' together: a PR-AUC of, say, 0.55 is informative when the baseline is 0.02
#' but nearly uninformative when the baseline is 0.50.
#'
#' @export
#'
#' @family Section_5_Other_Functions_Provided
#'
#' @importFrom Seurat DefaultAssay FetchData
#' @importFrom stats sd
#' @importFrom ggplot2 ggplot aes geom_line geom_hline labs theme_minimal
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous expansion
#' @importFrom ggplot2 element_line element_text element_rect theme
#' @importFrom scales number_format
#'
#' @examples
#' \dontrun{
#' # Default: raw expression, plot displayed automatically
#' res <- Compute_Gene_PRC(sce, "DCLK1", "seurat_clusters", 12)
#'
#' # Rank-based (robust to dropout)
#' res_rank <- Compute_Gene_PRC(sce, "DCLK1", "seurat_clusters", 12,
#'                              method = "rank")
#'
#' # Keep only cells expressing above 0.5 (exploratory)
#' res_sub <- Compute_Gene_PRC(sce, "DCLK1", "seurat_clusters", 12,
#'                             keep_expression_above = 0.5)
#' }
Compute_Gene_PRC <- function(
    seurat_obj,
    gene,
    group_col,
    group_label,
    assay = "RNA",
    method = c("raw", "rank"),
    min_expression = NULL,
    keep_expression_above = NULL,
    plot = TRUE,
    plot_title = "Gene Precision-Recall Curve | SlimR",
    line_color = "darkmagenta",
    line_size = 1
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("'seurat_obj' must be a Seurat object.")
  }
  if (!is.character(gene) || length(gene) != 1L) {
    stop("'gene' must be a single character string.")
  }
  if (!group_col %in% colnames(seurat_obj@meta.data)) {
    stop("'group_col' must be a column in seurat_obj@meta.data.")
  }
  if (!gene %in% rownames(seurat_obj[[assay]])) {
    stop("Gene '", gene, "' not found in the '", assay, "' assay.")
  }
  method <- match.arg(method)
  if (!is.null(min_expression) &&
      (!is.numeric(min_expression) || min_expression < 0)) {
    stop("'min_expression' must be NULL or a non-negative numeric value.")
  }
  if (!is.null(keep_expression_above) &&
      (!is.numeric(keep_expression_above) || keep_expression_above < 0)) {
    stop("'keep_expression_above' must be NULL or a non-negative numeric value.")
  }
  if (!is.logical(plot) || length(plot) != 1L) {
    stop("'plot' must be a single logical value.")
  }

  Seurat::DefaultAssay(seurat_obj) <- assay
  expr_data <- Seurat::FetchData(seurat_obj, vars = gene)
  expr_vec <- expr_data[[gene]]

  meta_col <- seurat_obj@meta.data[[group_col]]
  labels <- meta_col == group_label

  if (length(unique(labels)) < 2L) {
    stop("After binarization by 'group_label', only one class is present. ",
         "PR-AUC cannot be computed.")
  }

  if (!is.null(keep_expression_above)) {
    keep_idx <- expr_vec > keep_expression_above
    if (sum(keep_idx) == 0) {
      stop("No cells with expression above 'keep_expression_above' = ",
           keep_expression_above, ". Cannot compute PR-AUC.")
    }
    expr_vec <- expr_vec[keep_idx]
    labels <- labels[keep_idx]
    if (length(unique(labels)) < 2L) {
      stop("After subsetting by 'keep_expression_above', only one class remains. ",
           "Cannot compute PR-AUC.")
    }
  }

  if (method == "raw") {
    if (!is.null(min_expression)) {
      expr_vec[expr_vec < min_expression] <- 0
    }
    scores <- expr_vec
    method_label <- if (!is.null(min_expression)) {
      paste0("raw (min_expression = ", min_expression, ")")
    } else {
      "raw"
    }
  } else {
    scores <- rank(expr_vec, ties.method = "average") / length(expr_vec)
    method_label <- "rank"
  }

  if (!is.null(keep_expression_above)) {
    method_label <- paste0(method_label,
                           " | subset: expr > ", keep_expression_above,
                           " (n = ", length(scores), ")")
  }

  if (stats::sd(scores, na.rm = TRUE) < .Machine$double.eps) {
    warning("Scores for gene '", gene, "' have near-zero variance. ",
            "PR-AUC may be unreliable.")
  }

  pr_res <- .fastPR(scores, labels)
  pr_df <- .compute_pr_data(scores, labels)

  out <- list(
    PR_AUC      = pr_res$pr_auc,
    Baseline    = pr_res$baseline,
    pr_data     = pr_df,
    predictions = scores,
    labels      = labels,
    pr_plot     = NULL
  )

  if (plot) {
    subtitle_text <- paste0("PR-AUC = ", round(pr_res$pr_auc, 4),
                            "  |  baseline = ", round(pr_res$baseline, 4),
                            "  |  method: ", method_label)
    out$pr_plot <- ggplot2::ggplot(
      pr_df,
      ggplot2::aes(x = .data$recall, y = .data$precision)
    ) +
      ggplot2::geom_hline(yintercept = pr_res$baseline,
                          linetype = "dashed", color = "grey40") +
      ggplot2::geom_line(color = line_color, linewidth = line_size) +
      ggplot2::scale_x_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        labels = scales::number_format(accuracy = 0.1),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        labels = scales::number_format(accuracy = 0.1),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::labs(
        title    = plot_title,
        subtitle = subtitle_text,
        x        = "Recall (Sensitivity)",
        y        = "Precision (Positive Predictive Value)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(size = 14, face = "bold",
                                                 hjust = 0.5),
        plot.subtitle    = ggplot2::element_text(size = 11, hjust = 0.5),
        axis.title       = ggplot2::element_text(size = 11),
        axis.text        = ggplot2::element_text(size = 9),
        panel.grid.major = ggplot2::element_line(color = "grey90"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border     = ggplot2::element_rect(color = "black", fill = NA,
                                                 linewidth = 0.5)
      )

    print(out$pr_plot)
  }

  invisible(out)
}

# Internal: Average Precision (PR-AUC) and no-skill baseline.
# AP = mean of precision values at each positive cell, weighted by the
# corresponding increase in recall. Because recall increases by exactly
# 1 / n_pos at each positive, this simplifies to sum(precision[pos]) / n_pos.
.fastPR <- function(predictions, labels) {
  ord <- order(predictions, decreasing = TRUE)
  labels <- labels[ord]

  n_pos <- sum(labels)
  n_neg <- sum(!labels)

  if (n_pos == 0L) {
    stop("No positive samples available for PR-AUC computation.")
  }

  tp <- cumsum(labels)
  fp <- cumsum(!labels)

  precision <- tp / (tp + fp)

  pr_auc   <- sum(precision[labels]) / n_pos
  baseline <- n_pos / (n_pos + n_neg)

  list(pr_auc = pr_auc, baseline = baseline)
}

# Internal: PR curve data for plotting. The (recall = 0, precision = 1)
# anchor is prepended so that the curve starts at the top-left corner, as
# is conventional for precision-recall plots.
.compute_pr_data <- function(predictions, labels) {
  ord <- order(predictions, decreasing = TRUE)
  labels <- labels[ord]

  n_pos <- sum(labels)

  tp <- cumsum(labels)
  fp <- cumsum(!labels)

  precision <- tp / (tp + fp)
  recall    <- tp / n_pos

  precision <- c(1, precision)
  recall    <- c(0, recall)

  data.frame(recall = recall, precision = precision,
             stringsAsFactors = FALSE)
}