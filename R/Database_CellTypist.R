#' List of cell type markers from the CellTypist organ atlas
#'
#' A dataset containing marker genes for human cell types across 12 organs,
#' derived from the CellTypist online resource.
#'
#' @format A named list with 12 elements, each corresponding to an organ:
#'   Blood, Bone_marrow, Heart, Hippocampus, Intestine, Kidney, Liver,
#'   Lung, Lymph_node, Pancreas, Skeletal_muscle, Spleen. Each organ element
#'   is a list of data frames (one per cell type), containing marker genes
#'   and associated statistics. In total, 399 cell types are represented.
#'
#' @details This list comprises markers for 399 cell types across 12 human organs,
#'   obtained from the CellTypist organ atlas (\url{https://www.celltypist.org/organs}).
#'   Marker genes were generated using a standard Scanpy workflow: log1p‑normalised
#'   data, `sc.tl.rank_genes_groups` (Wilcoxon test), adjusted p‑value < 0.01 and
#'   log2 fold‑change > 0, ranked by log fold‑change and cut to the top 100 genes
#'   per cluster. The resulting Excel file was imported with `Read_excel_markers`.
#'   Please cite:
#'   Xu et al. (2023) \doi{10.1016/j.cell.2023.11.026} and
#'   Domínguez Conde et al. (2022) \doi{10.1126/science.abl5197}.
#'
#' @family Section_0_Database
#'
#' @source \doi{10.1016/j.cell.2023.11.026}
#' @source \doi{10.1126/science.abl5197}
"CellTypist"