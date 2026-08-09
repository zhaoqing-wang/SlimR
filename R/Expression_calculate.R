#' Counts average expression of gene set (Use in package)
#'
#' @param object Enter a Seurat object.
#' @param features Enter one or a set of markers.
#' @param assay Enter the assay used by the Seurat object, such as "RNA". Default
#'     parameters use "assay = NULL".
#' @param cluster_col Enter the meta.data column in the Seurat object to be
#'     annotated, such as "seurat_cluster". Default parameters use "cluster_col = NULL".
#' @param colour_low Color for lowest expression level. (default = "white")
#' @param colour_high Color for highest expression level. (default = "black")
#'
#' @returns Average expression genes and relatied informations in the input "Seurat" object
#'     given "cluster_col" and given "features".
#'
#' @family Section_1_Functions_Use_in_Package
#'
#' @importFrom Seurat DefaultAssay DefaultAssay<- CellsByIdentities FetchData
#' @importFrom dplyr group_by summarise left_join
#'
calculate_expression <- function(
    object,
    features,
    assay = NULL,
    cluster_col = NULL,
    colour_low = "white",
    colour_high = "navy") {

  if (!is.null(cluster_col) && !(cluster_col %in% colnames(object@meta.data))) {
    stop("cluster_col not found in meta.data")
  }

  if (is.null(assay)) {
  assay <- DefaultAssay(object)
  } else if (!(assay %in% names(object@assays))) {
    stop(paste0("Assay not found in Seurat object: ", assay))
  }
  DefaultAssay(object) <- assay

  cells <- unlist(CellsByIdentities(object = object, cells = colnames(object[[assay]])))

  data.features <- FetchData(object = object, vars = features, cells = cells)

  if (!is.null(cluster_col)) {
    data.features$id <- object@meta.data[cells, cluster_col, drop = TRUE]
  } else {
    data.features$id <- Idents(object = object)[cells]
  }

  id.levels <- levels(factor(data.features$id))
  data.features$id <- factor(data.features$id, levels = id.levels)

  cluster_expr_list <- lapply(unique(data.features$id), function(ident) {
    data.use <- data.features[data.features$id == ident, features, drop = FALSE]
    avg.exp <- apply(data.use, 2, function(x) mean(expm1(x)))
    return(avg.exp)
  })

  expr_matrix <- do.call(rbind, cluster_expr_list)
  rownames(expr_matrix) <- unique(data.features$id)

  expr_df <- as.data.frame(expr_matrix) %>%
    tibble::rownames_to_column("cluster") %>%
    tidyr::pivot_longer(
      cols = -cluster,
      names_to = "gene",
      values_to = "mean_expression"
    )

  expr_df$cluster <- factor(expr_df$cluster, levels = id.levels)

  cluster_avg_exp <- expr_df %>%
    group_by(cluster) %>%
    summarise(Avg_exp = mean(mean_expression), .groups = "drop")

  expr_df <- left_join(expr_df, cluster_avg_exp, by = "cluster")

  return(expr_df)
}
