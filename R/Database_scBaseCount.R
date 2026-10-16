#' scBaseCount dataset
#'
#' A processed cross-species marker gene resource derived from scBaseCount.
#' The dataset contains four species-specific gene set lists: human, mouse,
#' Drosophila, and zebrafish. Each species-specific element is a named list of
#' marker gene sets, where names correspond to cell types and values are
#' character vectors of marker genes.
#'
#' @format A named list with 4 elements:
#' \describe{
#'   \item{human_gene_sets}{A named list of human marker gene sets. Names are
#'       cell types; each element is a character vector of marker genes.}
#'   \item{mouse_gene_sets}{A named list of mouse marker gene sets. Names are
#'       cell types; each element is a character vector of marker genes.}
#'   \item{fly_gene_sets}{A named list of Drosophila marker gene sets. Names are
#'       cell types; each element is a character vector of marker genes.}
#'   \item{zebrafish_gene_sets}{A named list of zebrafish marker gene sets.
#'       Names are cell types; each element is a character vector of marker
#'       genes.}
#' }
#'
#' @details The scBaseCount resource was built from more than 500 million cells
#'   across four species. The species-specific gene set lists can be extracted
#'   and used directly as standardized marker lists, for example:
#'   \code{Markers_list_scBaseCount_human <- scBaseCount$human_gene_sets},
#'   \code{Markers_list_scBaseCount_mouse <- scBaseCount$mouse_gene_sets},
#'   \code{Markers_list_scBaseCount_fly <- scBaseCount$fly_gene_sets}, and
#'   \code{Markers_list_scBaseCount_zebrafish <- scBaseCount$zebrafish_gene_sets}.
#'   These lists are compatible with SlimR annotation functions.
#'
#'   When using this database, please cite Youngblut et al. (2026)
#'   \doi{10.1016/j.cell.2026.08.025}.
#'
#' @family Section_0_Database
#'
#' @source \doi{10.1016/j.cell.2026.08.025}
#'
#' @keywords datasets
"scBaseCount"