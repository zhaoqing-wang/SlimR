#' Compute AUC and Optionally Plot ROC Curve for a Single Gene
#'
#' @description
#' Evaluates the discriminatory power of a single gene in separating a
#' user-defined positive cell group from the rest, using the Area Under the
#' Receiver Operating Characteristic curve (AUC). Several scoring strategies
#' and optional cell subsetting are provided to handle dropout noise and
#' sparse expression typical of single-cell data.
#' If \code{plot = TRUE} (the default), a publication‑ready ROC curve is
#' automatically displayed via \pkg{ggplot2}.
#'
#' @param seurat_obj A Seurat object containing single-cell expression data.
#' @param gene A single character string specifying the gene to evaluate. Must
#'   be present in the rownames of the chosen assay.
#' @param group_col Character string giving the name of a column in
#'   `seurat_obj@meta.data` that defines cell groups (e.g., `"seurat_clusters"`).
#' @param group_label The value in `group_col` that defines the positive class
#'   (e.g., `1` for cluster 1). All other cells are treated as negatives.
#' @param assay Character string specifying which assay to use. Default is
#'   `"RNA"`.
#' @param method Scoring method: \code{"raw"} (raw expression, optionally
#'   truncated by \code{min_expression}) or \code{"rank"} (expression ranks,
#'   which is robust to dropout and does not require \code{min_expression}).
#' @param min_expression Numeric threshold for expression truncation (only
#'   used when \code{method = "raw"}). Values below this threshold are set to
#'   zero. Default is \code{NULL} (no truncation).
#' @param keep_expression_above Optional numeric threshold. If set, only cells
#'   with expression greater than this value are retained for AUC computation.
#'   **Note:** this changes the population over which AUC is evaluated, shifting
#'   the interpretation from "overall discrimination among all cells" to
#'   "discrimination among cells that express the gene above this level". Use
#'   with caution and always compare with the result without subsetting.
#'   Default is \code{NULL} (all cells retained).
#' @param plot Logical indicating whether to create and automatically display
#'   a ggplot2 ROC curve. Default is \code{TRUE}.
#' @param plot_title Character string used as the title of the ROC plot.
#'   Default is \code{"ROC Curve | SlimR"}.
#' @param line_color Colour of the ROC curve. Default is \code{"navy"}.
#' @param line_size Numeric value for the thickness of the ROC curve line.
#'   Default is \code{1}.
#'
#' @return Invisibly, a list with the following elements:
#' \describe{
#'   \item{AUC}{Numeric value (between 0 and 1) of the area under the ROC curve.}
#'   \item{roc_data}{A data.frame with columns \code{fpr} (False Positive Rate)
#'     and \code{tpr} (True Positive Rate) that can be used for custom plotting.}
#'   \item{predictions}{Numeric vector of the (possibly truncated) expression
#'     values used as prediction scores.}
#'   \item{labels}{Logical vector indicating whether each cell belongs to the
#'     positive class (\code{TRUE}) or not (\code{FALSE}).}
#'   \item{roc_plot}{If \code{plot = TRUE}, a **ggplot** object displaying the
#'     ROC curve; otherwise \code{NULL}.}
#' }
#'
#' @export
#'
#' @family Section_5_Other_Functions_Provided
#'
#' @importFrom Seurat DefaultAssay FetchData
#' @importFrom stats sd
#' @importFrom ggplot2 ggplot aes geom_line geom_abline labs theme_minimal
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous expansion
#' @importFrom ggplot2 element_line element_text element_rect theme
#' @importFrom scales number_format
#' @importFrom utils head tail
#'
#' @examples
#' \dontrun{
#' # Default: raw expression, plot displayed automatically
#' res <- Compute_Gene_AUC_ROC(sce, "BMX", "seurat_clusters", 12)
#'
#' # Rank-based (robust to dropout)
#' res_rank <- Compute_Gene_AUC_ROC(sce, "BMX", "seurat_clusters", 12,
#'                                  method = "rank")
#'
#' # Keep only cells expressing above 0.5 (exploratory)
#' res_sub <- Compute_Gene_AUC_ROC(sce, "BMX", "seurat_clusters", 12,
#'                                 keep_expression_above = 0.5)
#' }
Compute_Gene_AUC_ROC <- function(
    seurat_obj,
    gene,
    group_col,
    group_label,
    assay = "RNA",
    method = c("raw", "rank"),
    min_expression = NULL,
    keep_expression_above = NULL,
    plot = TRUE,
    plot_title = "ROC Curve | SlimR",
    line_color = "navy",
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
  if (!is.null(min_expression) && (!is.numeric(min_expression) || min_expression < 0)) {
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
         "AUC cannot be computed.")
  }

  if (!is.null(keep_expression_above)) {
    keep_idx <- expr_vec > keep_expression_above
    if (sum(keep_idx) == 0) {
      stop("No cells with expression above 'keep_expression_above' = ",
           keep_expression_above, ". Cannot compute AUC.")
    }
    expr_vec <- expr_vec[keep_idx]
    labels <- labels[keep_idx]
    if (length(unique(labels)) < 2L) {
      stop("After subsetting by 'keep_expression_above', only one class remains. ",
           "Cannot compute AUC.")
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
    warning("Scores for gene '", gene, "' have near-zero variance. AUC may be unreliable.")
  }

  auc_val <- .fastAUC(scores, labels)
  roc_df <- .compute_roc_data(scores, labels)

  out <- list(
    AUC         = auc_val,
    roc_data    = roc_df,
    predictions = scores,
    labels      = labels,
    roc_plot    = NULL
  )

  if (plot) {
    subtitle_text <- paste0("AUC = ", round(auc_val, 4),
                            "  |  method: ", method_label)
    out$roc_plot <- ggplot2::ggplot(
      roc_df,
      ggplot2::aes(x = .data$fpr, y = .data$tpr)
    ) +
      ggplot2::geom_line(color = line_color, linewidth = line_size) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", color = "grey40") +
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
        x        = "False Positive Rate (1 - Specificity)",
        y        = "True Positive Rate (Sensitivity)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle    = ggplot2::element_text(size = 11, hjust = 0.5),
        axis.title       = ggplot2::element_text(size = 11),
        axis.text        = ggplot2::element_text(size = 9),
        panel.grid.major = ggplot2::element_line(color = "grey90"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border     = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
    
    print(out$roc_plot)
  }

  invisible(out)
}

.fastAUC <- function(predictions, labels) {
  ord <- order(predictions, decreasing = TRUE)
  labels <- labels[ord]
  predictions <- predictions[ord]

  tpr <- cumsum(labels) / sum(labels)
  fpr <- cumsum(!labels) / sum(!labels)

  tpr <- c(0, tpr, 1)
  fpr <- c(0, fpr, 1)

  auc <- sum(diff(fpr) * (utils::head(tpr, -1) + utils::tail(tpr, -1)) / 2)
  return(auc)
}

.compute_roc_data <- function(predictions, labels) {
  ord <- order(predictions, decreasing = TRUE)
  labels <- labels[ord]

  n_pos <- sum(labels)
  n_neg <- sum(!labels)

  tpr <- cumsum(labels) / n_pos
  fpr <- cumsum(!labels) / n_neg

  tpr <- c(0, tpr, 1)
  fpr <- c(0, fpr, 1)

  data.frame(fpr = fpr, tpr = tpr, stringsAsFactors = FALSE)
}