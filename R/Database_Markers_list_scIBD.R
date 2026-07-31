#' List of cell type markers from scIBD
#'
#' A dataset containing marker genes for human intestinal cell types
#' derived from scIBD.
#'
#' @format A named list of 101 data frames, each representing one cell type
#'   and containing marker genes along with associated statistics.
#'
#' @details This list comprises markers for 101 human intestine cell types,
#'   obtained from the scIBD study (Nie et al. 2023).
#'   Marker genes were generated using a standard Scanpy workflow:
#'   log1p‑normalised data, `sc.tl.rank_genes_groups` (Wilcoxon test),
#'   adjusted p‑value < 0.01 and log2 fold‑change > 0, ranked by log
#'   fold‑change and cut to the top 100 genes per cluster (equivalent to
#'   `sort_by = "logFC"` and `gene_filter = 100`). The resulting Excel file
#'   was imported with `Read_excel_markers`.
#'
#' @family Section_0_Database
#'
#' @source \doi{10.1038/s43588-023-00464-9}
"Markers_list_scIBD"