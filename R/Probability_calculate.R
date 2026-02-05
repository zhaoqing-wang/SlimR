#' Calculate gene set expression and infer probabilities with control datasets (Use in package)
#'
#' @param object Enter a Seurat object.
#' @param features Enter one or a set of markers.
#' @param assay Enter the assay used by the Seurat object, such as "RNA". Default
#'     parameters use "assay = NULL".
#' @param cluster_col Enter the meta.data column in the Seurat object to be
#'     annotated, such as "seurat_cluster". Default parameters use "cluster_col = NULL".
#' @param min_expression The min_expression parameter defines a threshold value to
#'     determine whether a cell's expression of a feature is considered "expressed"
#'     or not. It is used to filter out low-expression cells that may contribute
#'     noise to the analysis. Default parameters use "min_expression = 0.1".
#' @param specificity_weight The specificity_weight parameter controls how much the
#'     expression variability (standard deviation) of a feature within a cluster
#'     contributes to its "specificity score." It amplifies or suppresses the impact
#'     of variability in the final score calculation.Default parameters use
#'     "specificity_weight = 3".
#'
#' @returns Average expression of genes in the input "Seurat" object given
#'     "cluster_col" and given "features".
#'
#' @family Section_1_Functions_Use_in_Package
#'
#' @importFrom Seurat DefaultAssay DefaultAssay<- CellsByIdentities FetchData
#' @importFrom stats sd weighted.mean
#'
calculate_probability <- function(
    object,
    features,
    assay = NULL,
    cluster_col = NULL,
    min_expression = 0.1,
    specificity_weight = 3
) {
  assay <- if (is.null(assay)) DefaultAssay(object) else assay
  DefaultAssay(object) <- assay

  cells <- unlist(CellsByIdentities(object = object))
  data.features <- FetchData(object = object, vars = features, cells = cells)

  if (!is.null(cluster_col)) {
    data.features$id <- object@meta.data[cells, cluster_col, drop = TRUE]
  } else {
    data.features$id <- Idents(object = object)[cells]
  }

  id.levels <- levels(factor(data.features$id))
  data.features$id <- factor(data.features$id, levels = id.levels)

  cluster_stats <- lapply(unique(data.features$id), function(ident) {
    cluster_cells <- data.features$id == ident
    data.use <- data.features[cluster_cells, features, drop = FALSE]

    avg.exp <- apply(data.use, 2, function(x) mean(expm1(x)))
    sd.exp <- apply(data.use, 2, sd)
    frac.expressing <- apply(data.use, 2, function(x) mean(x > min_expression))

    specificity_score <- avg.exp * frac.expressing *
      (1 + specificity_weight * (sd.exp / (mean(sd.exp) + 1e-6)))

    return(list(
      avg = avg.exp,
      frac = frac.expressing,
      specificity = specificity_score
    ))
  })

  expr_matrix <- do.call(rbind, lapply(cluster_stats, function(x) x$avg))
  rownames(expr_matrix) <- unique(data.features$id)

  frac_matrix <- do.call(rbind, lapply(cluster_stats, function(x) x$frac))
  rownames(frac_matrix) <- unique(data.features$id)

  specificity_matrix <- do.call(rbind, lapply(cluster_stats, function(x) x$specificity))
  rownames(specificity_matrix) <- unique(data.features$id)

  normalize_column <- function(x) {
    # Handle NA values
    x[is.na(x)] <- 0
    # Check if all values are the same
    range_diff <- diff(range(x, na.rm = TRUE))
    if (is.na(range_diff) || range_diff == 0) return(rep(0, length(x)))
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  normalized_matrix <- apply(specificity_matrix, 2, normalize_column)

  gene_weights <- apply(specificity_matrix, 2, function(x) {
    m <- mean(x, na.rm = TRUE)
    if (is.na(m) || m == 0) 0 else sd(x, na.rm = TRUE) / m
  })

  cluster_scores <- apply(normalized_matrix, 1, function(x) {
    weighted.mean(x, w = gene_weights, na.rm = TRUE)
  })

  names(cluster_scores) <- rownames(specificity_matrix)

  return(list(cluster_expr = expr_matrix,
              cluster_frac = frac_matrix,
              cluster_scores = cluster_scores))
}
