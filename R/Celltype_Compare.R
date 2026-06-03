#' Compare cell type labels across two single-cell datasets after aligning cell barcodes
#'
#' This function automatically aligns cell barcodes between two Seurat objects using a
#' variety of normalization transformations, then cross-tabulates a cell type label column
#' (from the first object) against a grouping column (from the second object). It returns
#' count tables, proportion tables, a dominant mapping, and a heatmap.
#'
#' @param sce_label A Seurat object containing the cell type label column.
#' @param sce A Seurat object containing the grouping column.
#' @param label_col Character. Name of the metadata column in `sce_label` that stores
#'   cell type labels (e.g., "Sub_cell_type").
#' @param group_col Character. Name of the metadata column in `sce` that stores
#'   grouping information (e.g., "SCT_snn_res.0.3").
#' @param barcode_col Optional character. Name of a metadata column in both objects
#'   that contains the cell barcode identifiers. If `NULL`, the function uses
#'   `colnames(sce_label)` and `colnames(sce)`.
#' @param color_low Character. Color for low proportion values in the heatmap.
#'   Default: "grey70".
#' @param color_high Character. Color for high proportion values in the heatmap.
#'   Default: "navy".
#' @param show_plot Logical. If `TRUE` (default), the heatmap is automatically
#'   displayed when the function is called. Set to `FALSE` to suppress
#'   automatic plotting (e.g., in non‑interactive environments).
#'
#' @return A list with five components:
#'   \item{count_table}{A data frame (wide format) with rows = unique
#'     `label_col` values and columns = unique `group_col` values; cell values are
#'     raw counts of shared cells.}
#'   \item{prop_table}{Same shape as `count_table`; each cell shows the proportion
#'     of cells within a `group_col` column (column-wise sum = 1).}
#'   \item{main_to_sub}{A data frame mapping each `group_col` value to the most
#'     frequent `label_col` value among shared cells.}
#'   \item{plot}{A ggplot2 heatmap object visualizing the proportion table.}
#'   \item{match_info}{A tibble with columns `label_transform`, `sce_transform`,
#'     `shared_n` – the transformations used to align barcodes and the number of
#'     shared cells after alignment.}
#'
#' @details
#' **Cell barcode alignment**:
#' The function automatically tries a set of normalization functions on the cell
#' identifiers (either from `barcode_col` or from column names) to maximise the
#' number of shared barcodes between the two objects. Transformations include:
#' `identity`, `drop_numeric_suffix` (removes e.g., "-1-2"), `drop_suffix` (removes
#' "-1"), and several prefix removals. The transformation pair yielding the highest
#' number of shared identifiers is selected.
#'
#' **Proportion calculation**:
#' Proportions are computed **within each `group_col` level** (column-wise),
#' i.e. for each group, the sum of proportions across all cell types equals 1.
#'
#' **Plot**:
#' The heatmap uses `ggplot2::geom_tile()` with a fixed coordinate ratio and a
#' colour gradient from `color_low` to `color_high`.
#'
#' @examples
#' \dontrun{
#' # Basic usage with two Seurat objects and default barcode alignment
#' result <- Celltype_Compare(
#'   sce_label = seurat_obj1,
#'   sce = seurat_obj2,
#'   label_col = "cell_type",
#'   group_col = "cluster"
#' )
#'
#' # Access the proportion table
#' head(result$prop_table)
#'
#' # View the dominant mapping
#' print(result$main_to_sub)
#'
#' # Display the heatmap
#' print(result$plot)
#'
#' # Use a custom barcode column
#' result2 <- Celltype_Compare(
#'   sce_label = seurat_obj1,
#'   sce = seurat_obj2,
#'   label_col = "cell_type",
#'   group_col = "cluster",
#'   barcode_col = "cell_barcode"
#' )
#' }
#'
#' @export
#' @family Section_5_Other_Functions_Provided
#' 
Celltype_Compare <- function(
    sce_label,
    sce,
    label_col = NULL,
    group_col = NULL,
    barcode_col = NULL,
    color_low = "grey70",
    color_high = "navy",
    show_plot = TRUE
) {
  # Align by shared cell barcodes, with automatic normalization
  if (is.null(barcode_col)) {
    label_cells <- colnames(sce_label)
    sce_cells <- colnames(sce)
  } else {
    if (!barcode_col %in% colnames(sce_label@meta.data)) {
      stop(paste0("Column not found in sce_label meta.data: ", barcode_col))
    }
    if (!barcode_col %in% colnames(sce@meta.data)) {
      stop(paste0("Column not found in sce meta.data: ", barcode_col))
    }
    label_cells <- sce_label@meta.data[[barcode_col]]
    sce_cells <- sce@meta.data[[barcode_col]]
  }

  normalize_candidates <- list(
    identity = function(x) x,
    drop_numeric_suffix = function(x) sub("-\\d+(?:-\\d+)?$", "", x),
    drop_suffix = function(x) sub("-\\d+$", "", x),
    drop_prefix_underscore = function(x) sub("^.*_", "", x),
    drop_prefix_dash = function(x) sub("^.*-", "", x),
    drop_prefix_hash = function(x) sub("^.*#", "", x),
    drop_prefix_colon = function(x) sub("^.*:", "", x),
    drop_prefix_pipe = function(x) sub("^.*\\|", "", x)
  )

  choose_best_match <- function(label_cells, sce_cells, transforms) {
    results <- list()
    for (label_name in names(transforms)) {
      label_norm <- transforms[[label_name]](label_cells)
      label_map <- dplyr::tibble(cell = label_cells, norm = label_norm) |>
        dplyr::distinct(norm, .keep_all = TRUE)
      for (sce_name in names(transforms)) {
        sce_norm <- transforms[[sce_name]](sce_cells)
        sce_map <- dplyr::tibble(cell = sce_cells, norm = sce_norm) |>
          dplyr::distinct(norm, .keep_all = TRUE)
        shared_norm <- intersect(label_map$norm, sce_map$norm)
        results[[length(results) + 1]] <- dplyr::tibble(
          label_transform = label_name,
          sce_transform = sce_name,
          shared_n = length(shared_norm),
          label_map = list(label_map),
          sce_map = list(sce_map),
          shared_norm = list(shared_norm)
        )
      }
    }

    dplyr::bind_rows(results) |>
      dplyr::arrange(dplyr::desc(.data$shared_n), .data$label_transform, .data$sce_transform) |>
      dplyr::slice_head(n = 1)
  }

  match_info <- choose_best_match(label_cells, sce_cells, normalize_candidates)

  if (match_info$shared_n == 0) {
    stop("No shared cell barcodes found between sce_label and sce after normalization.")
  }

  label_map <- match_info$label_map[[1]]
  sce_map <- match_info$sce_map[[1]]
  shared_norm <- match_info$shared_norm[[1]]

  shared_cells_label <- label_map |>
    dplyr::filter(norm %in% shared_norm)
  shared_cells_sce <- sce_map |>
    dplyr::filter(norm %in% shared_norm)

  shared_cells <- dplyr::tibble(
    norm = shared_norm
  ) |>
    dplyr::left_join(shared_cells_label, by = "norm") |>
    dplyr::left_join(shared_cells_sce, by = "norm", suffix = c("_label", "_sce"))

  label_meta <- sce_label@meta.data[shared_cells$cell_label, , drop = FALSE]
  sce_meta <- sce@meta.data[shared_cells$cell_sce, , drop = FALSE]

  if (!label_col %in% colnames(label_meta)) {
    stop(paste0("Column not found in sce_label meta.data: ", label_col))
  }
  if (!group_col %in% colnames(sce_meta)) {
    stop(paste0("Column not found in sce meta.data: ", group_col))
  }

  df <- dplyr::tibble(
    cell = shared_cells$cell_label,
    sub_cell_type = label_meta[[label_col]],
    main_group = sce_meta[[group_col]]
  ) |>
    dplyr::filter(!is.na(sub_cell_type), !is.na(main_group))

  # Comparison table (counts)
  comparison_table <- df |>
    dplyr::count(sub_cell_type, main_group, name = "n") |>
    dplyr::arrange(dplyr::desc(.data$n))

  # Proportion table (within each main_group column)
  proportion_table <- comparison_table |>
    dplyr::group_by(.data$main_group) |>
    dplyr::mutate(prop = .data$n / sum(.data$n)) |>
    dplyr::ungroup()

  # Wide count and proportion tables (rows=sub_cell_type, cols=main_group)
  count_table <- comparison_table |>
    tidyr::pivot_wider(
    names_from = main_group,
    values_from = n,
    values_fill = 0
  ) |>
    dplyr::arrange(.data$sub_cell_type)

  # Column-wise proportion (robust to zero columns)
  prop_table <- count_table
  col_sums <- colSums(prop_table[, -1, drop = FALSE])
  # Avoid division by zero: if a column sum is zero, keep zeros
  prop_table[, -1] <- sweep(prop_table[, -1, drop = FALSE], 2, col_sums, FUN = "/")
  prop_table[, -1][is.na(prop_table[, -1])] <- 0

  # Main group to dominant sub_cell_type mapping
  main_to_sub <- comparison_table |>
    dplyr::group_by(.data$main_group) |>
    dplyr::slice_max(order_by = .data$n, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$main_group)

  # Plot: heatmap of proportions
  comparison_plot <- ggplot2::ggplot(
    proportion_table,
    ggplot2::aes(x = main_group, y = sub_cell_type, fill = prop)
  ) +
    ggplot2::geom_tile(color = "grey70", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(
      low = color_low,
      high = color_high
    ) +
    ggplot2::labs(
      x = group_col,
      y = label_col,
      fill = "Proportion",
      title = "Celltype mapping heatmap | SlimR"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.4)
    )

  if (show_plot) {
    print(comparison_plot)
  }

  list(
    count_table = count_table,
    prop_table = prop_table,
    main_to_sub = main_to_sub,
    plot = comparison_plot,
    match_info = match_info |>
      dplyr::select(.data$label_transform, .data$sce_transform, .data$shared_n)
  )
}