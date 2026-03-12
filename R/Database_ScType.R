#' ScType raw dataset
#'
#' The original ScType marker database before processing.
#'
#' @format A tibble with 5 columns:
#' \describe{
#'   \item{tissueType}{Tissue type}
#'   \item{cellName}{Full cell type name}
#'   \item{geneSymbolmore1}{Comma-separated positive marker genes}
#'   \item{geneSymbolmore2}{Comma-separated negative marker genes (not used in processing)}
#'   \item{shortName}{Abbreviated cell type name}
#' }
#'
#' @family Section_0_Database
#'
#' @source \url{https://github.com/IanevskiAleksandr/sc-type}
"ScType_raw"

#' ScType dataset
#'
#' A processed long-format dataset containing marker genes for different cell
#' types from the ScType database. Each row represents one marker gene for a
#' given tissue type and cell type.
#'
#' @format A tibble with 3 columns:
#' \describe{
#'   \item{tissue_type}{Tissue type (e.g., "Immune system", "Brain", "Liver")}
#'   \item{cell_name}{Cell type name, formatted as "cellName(shortName)" when a
#'       short name is available, or "cellName" otherwise}
#'   \item{marker}{Gene symbol of the marker}
#' }
#'
#' @details This dataset is used to filter and create a standardized marker list.
#'     The dataset can be filtered based on tissue type and cell name to generate
#'     a list of marker genes for specific cell types using
#'     \code{\link{Markers_filter_ScType}}.
#'
#' @family Section_0_Database
#'
#' @source \url{https://github.com/IanevskiAleksandr/sc-type}
"ScType"

#' ScType metadata table
#'
#' A list of frequency tables summarizing the ScType database, useful for
#' exploring available tissue types and cell types before filtering.
#'
#' @format A list with 2 elements:
#' \describe{
#'   \item{tissue_type}{Frequency table of tissue types}
#'   \item{cell_name}{Frequency table of cell type names}
#' }
#'
#' @family Section_0_Database
#'
#' @source \url{https://github.com/IanevskiAleksandr/sc-type}
"ScType_table"
