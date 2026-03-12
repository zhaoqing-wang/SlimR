#' Create Marker_list from the ScType database
#'
#' @param df Standardized ScType database. It is read as \code{data(ScType)}
#'     in the SlimR library.
#' @param tissue_type Tissue type information in the ScType database. The input
#'     can be retrieved by \code{"ScType_table"}. For more information, please
#'     refer to \url{https://github.com/IanevskiAleksandr/sc-type}.
#' @param cell_name Cell type name information in the ScType database. The input
#'     can be retrieved by \code{"ScType_table"}. For more information, please
#'     refer to \url{https://github.com/IanevskiAleksandr/sc-type}.
#'
#' @returns The standardized "Marker_list" in the SlimR package
#' @export
#' @family Section_2_Standardized_Markers_List
#'
#' @examples
#' ScType <- SlimR::ScType
#' Markers_list_ScType <- Markers_filter_ScType(
#'     ScType,
#'     tissue_type = "Immune system",
#'     cell_name = NULL
#'     )
#'
Markers_filter_ScType <- function(df, tissue_type = NULL, cell_name = NULL) {

  required_columns <- c("tissue_type", "cell_name", "marker")

  if (!all(required_columns %in% colnames(df))) {
    stop("Data frame missing necessary columns! Please make sure to include the following: ",
         paste(required_columns, collapse = ", "))
  }

  filtered_df <- df

  if (!is.null(tissue_type)) {
    filtered_df <- filtered_df[filtered_df$tissue_type %in% tissue_type, ]
  }

  if (!is.null(cell_name)) {
    filtered_df <- filtered_df[filtered_df$cell_name %in% cell_name, ]
  }

  if (nrow(filtered_df) == 0) {
    stop("Filter result is empty, check criteria")
  }

  result_df <- filtered_df[, c("cell_name", "marker")]
  result_df$marker <- toupper(result_df$marker)
  result_df <- unique(result_df)

  split_result <- split(result_df, result_df$cell_name)
  clean_result <- lapply(split_result, function(sub_df) sub_df[, -1, drop = FALSE])

  return(clean_result)
}
