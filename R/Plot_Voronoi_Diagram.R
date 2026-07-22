#' Generate a Weighted Voronoi Treemap for Hierarchical Cell Types
#'
#' @description
#' Creates a weighted Voronoi treemap displaying the hierarchical composition
#' of single-cell data. The area of each polygon is proportional to the
#' number of cells in the corresponding sub-type, and the polygons are
#' grouped by the main cell type. Each sub-type is filled with a distinct
#' colour from the **ArchR**-inspired *stallion* palette (or a user-supplied
#' custom palette).
#'
#' **Note on colour handling:**  
#' The function uses `WeightedTreemaps::voronoiTreemap()` to compute polygon
#' coordinates, but **replaces the original drawing routine** with a custom
#' `ggplot2`-based renderer. This is necessary because the upstream
#' `add_color()` function treats discrete palettes as continuous gradients
#' and does not respect the order of polygon cells, making it impossible
#' to guarantee exact colour correspondence with functions like
#' `Plot_Hierarchy_Proportion` or `Seurat::DimPlot`.  
#' The custom drawing code extracts polygon vertices and uses
#' `scale_fill_manual()` with a named colour vector, ensuring that each
#' cell type receives exactly its intended colour regardless of internal
#' cell order.
#'
#' @section Colour palette:
#' The default colour palette is derived from the **ArchR** package
#' (\href{https://www.archrproject.com/}{ArchR project site},
#' \href{https://github.com/GreenleafLab/ArchR}{GitHub}).
#' If no custom palette is supplied, the function calls the internal
#' `paletteDiscrete()` function, which replicates the *stallion* palette
#' from ArchR. When the number of categories exceeds the palette size,
#' colours are interpolated. ArchR is described in:
#' Granja JM, Corces MR et al. ArchR is a scalable software package for
#' integrative single-cell chromatin accessibility analysis.
#' *Nature Genetics*, 2021. \doi{10.1038/s41588-021-00790-6}
#'
#' @param seurat_obj A Seurat object.
#' @param Main_cell_types Character string naming a column in
#'   \code{seurat_obj@meta.data} that holds the top-level (main) cell-type labels.
#' @param Cell_types Character string naming a column in
#'   \code{seurat_obj@meta.data} that holds the sub-type labels.
#' @param label_type Character string specifying the type of label to display
#'   inside each polygon. One of \code{"both"}, \code{"count"}, \code{"percentage"},
#'   or \code{"none"}. Default is \code{"both"}.
#' @param shape Character string specifying the shape of the polygons.
#'   \code{"rounded_rect"} (default) or \code{"circle"}.
#' @param seed Integer seed for reproducibility of the Voronoi layout.
#'   Default is \code{1}.
#' @param col_Cell_types Optional named or unnamed character vector of colours
#'   for the \code{Cell_types} categories. If provided, it overrides the
#'   default ArchR-derived palette.
#' @param label_size Numeric value controlling the size of the polygon labels.
#'   Default is \code{3}.
#' @param label_color Colour of the polygon labels. Default is \code{"black"}.
#' @param label_fontface Font face for polygon labels (\code{"plain"},
#'   \code{"italic"}, \code{"bold"}, etc.). Default is \code{"bold"}.
#' @param border_color Colour of the border lines separating polygons.
#'   Default is \code{"grey90"} (a very light grey).
#' @param subtype_border_lwd Line width for borders between sub-type polygons.
#'   Default is \code{0.15}.
#' @param main_border_lwd Line width for borders between main cell types.
#'   Default is \code{0.35}.
#' @param outer_border_lwd Line width for the outer frame of the entire plot.
#'   The frame is drawn as a tight convex hull around the outermost polygons,
#'   perfectly following the natural outline of the treemap. Default is \code{0.4}.
#' @param legend Logical indicating whether to display a colour legend.
#'   Default is \code{TRUE}.
#' @param legend_position Position of the legend when \code{legend = TRUE}.
#'   Default is \code{"right"}.
#'
#' @return
#' Invisibly, a list with \code{voronoi_treemap} (the raw treemap object) and
#' \code{plot} (the \pkg{ggplot2} object).
#'
#' @export
#' 
#' @family Section_5_Other_Functions_Provided
#'
#' @importFrom ggplot2 ggplot aes geom_polygon geom_text
#' @importFrom ggplot2 scale_fill_manual coord_equal theme_void
#' @importFrom ggplot2 element_text margin
#' @importFrom grDevices chull
#' @importFrom stats aggregate ave
#'
#' @examples
#' \dontrun{
#' # Basic rounded-rectangle treemap with both count and percentage labels
#' Plot_Voronoi_diagram(
#'   seurat_obj      = sce,
#'   Main_cell_types = "compartment",
#'   Cell_types      = "cell_type"
#' )
#'
#' # Circular polygons, show only percentages, custom colours
#' Plot_Voronoi_diagram(
#'   sce,
#'   "compartment",
#'   "cell_type",
#'   label_type     = "percentage",
#'   shape          = "circle",
#'   seed           = 42,
#'   col_Cell_types = c(
#'     "Enterocyte"   = "#1f78b4",
#'     "Goblet"       = "#33a02c",
#'     "Stem"         = "#e31a1c",
#'     "Paneth"       = "#ff7f00"
#'   ),
#'   legend_position = "right"
#' )
#'
#' # Only sub-type names, no legend
#' Plot_Voronoi_diagram(
#'   sce,
#'   "compartment",
#'   "cell_type",
#'   label_type = "none",
#'   legend     = FALSE
#' )
#' }
Plot_Voronoi_diagram <- function(
    seurat_obj,
    Main_cell_types,
    Cell_types,
    label_type = c("both", "count", "percentage", "none"),
    shape = c("rounded_rect", "circle"),
    seed = 1,
    col_Cell_types = NULL,
    label_size = 3.5,
    label_color = "black",
    label_fontface = "bold",
    border_color = "grey90",
    subtype_border_lwd = 0.5,
    main_border_lwd = 1.5,
    outer_border_lwd = 2.0,
    legend = TRUE,
    legend_position = "right"
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("'seurat_obj' must be a Seurat object.")
  }
  if (!requireNamespace("WeightedTreemaps", quietly = TRUE)) {
    stop(
      "Package 'WeightedTreemaps' is required but not installed.\n",
      "Install it with: devtools::install_github(\"m-jahn/WeightedTreemaps\")"
    )
  }
  stopifnot(is.character(Main_cell_types), length(Main_cell_types) == 1L)
  stopifnot(is.character(Cell_types), length(Cell_types) == 1L)
  if (!Main_cell_types %in% colnames(seurat_obj@meta.data)) {
    stop("'Main_cell_types' not found in seurat_obj@meta.data.")
  }
  if (!Cell_types %in% colnames(seurat_obj@meta.data)) {
    stop("'Cell_types' not found in seurat_obj@meta.data.")
  }

  label_type <- match.arg(label_type)
  shape <- match.arg(shape)

  meta <- seurat_obj@meta.data
  valid_char <- function(x) !is.na(x) & x != ""

  get_full_cats <- function(meta_col) {
    if (is.factor(meta_col)) {
      return(levels(meta_col))
    } else {
      vec <- as.character(meta_col)
      return(sort(unique(vec[valid_char(vec)])))
    }
  }

  full_main_cats <- get_full_cats(meta[[Main_cell_types]])
  full_cell_cats <- get_full_cats(meta[[Cell_types]])

  main_vec <- as.character(meta[[Main_cell_types]])
  cell_vec <- as.character(meta[[Cell_types]])

  keep <- valid_char(main_vec) & valid_char(cell_vec)
  main_vec <- main_vec[keep]
  cell_vec <- cell_vec[keep]

  present_main <- full_main_cats[full_main_cats %in% unique(main_vec)]
  present_cell <- full_cell_cats[full_cell_cats %in% unique(cell_vec)]

  tab <- as.data.frame(table(
    Main = factor(main_vec, levels = present_main),
    Cell = factor(cell_vec, levels = present_cell)
  ))
  tab <- tab[tab$Freq > 0, ]
  colnames(tab) <- c("Main_cell_type", "Cell_type", "Cell_number")
  tab$Percentage <- round(tab$Cell_number / sum(tab$Cell_number) * 100, 2)
  tab$Cell_type <- as.character(tab$Cell_type)
  tab$Main_cell_type <- as.character(tab$Main_cell_type)

  tab <- tab[order(factor(tab$Cell_type, levels = present_cell)), ]
  rownames(tab) <- NULL

  default_pal <- function(n, names = NULL) {
    if (requireNamespace("scales", quietly = TRUE)) {
      cols <- scales::hue_pal()(n)
    } else {
      cols <- grDevices::rainbow(n)
    }
    if (!is.null(names)) setNames(cols, names) else cols
  }

  get_palette <- function(cats, user_cols) {
    if (!is.null(user_cols)) {
      if (!is.null(names(user_cols))) {
        pal <- setNames(rep(NA_character_, length(cats)), cats)
        common <- intersect(cats, names(user_cols))
        pal[common] <- user_cols[common]
        missing <- cats[is.na(pal)]
        if (length(missing) > 0) {
          message("Some categories lack user colours, auto-generating...")
          pal[missing] <- default_pal(length(missing))
        }
        return(pal)
      } else {
        if (length(user_cols) < length(cats)) {
          user_cols <- rep(user_cols, length.out = length(cats))
        }
        return(setNames(user_cols[seq_along(cats)], cats))
      }
    }

    paletteDiscrete(cats)
  }

  pal_full <- get_palette(full_cell_cats, col_Cell_types)
  pal_cell <- pal_full[present_cell]

  tab$Label <- switch(
    label_type,
    both       = paste0(tab$Cell_type, "\n", tab$Cell_number, "\n", tab$Percentage, "%"),
    count      = paste0(tab$Cell_type, "\n", tab$Cell_number),
    percentage = paste0(tab$Cell_type, "\n", tab$Percentage, "%"),
    none       = tab$Cell_type
  )

  tm <- WeightedTreemaps::voronoiTreemap(
    data      = tab,
    levels    = c("Main_cell_type", "Label"),
    cell_size = "Cell_number",
    shape     = shape,
    seed      = seed
  )

  main_labels <- present_main
  sub_labels  <- tab$Label

  extract_poly <- function(cell, name) {
    mat <- cell$poly[[1]]
    if (is.null(mat) || nrow(mat) == 0) return(NULL)
    data.frame(
      x = mat[, 1],
      y = mat[, 2],
      name = name,
      stringsAsFactors = FALSE
    )
  }

  main_poly_list <- lapply(tm@cells, function(cell) {
    if (cell$name %in% main_labels) {
      extract_poly(cell, cell$name)
    } else NULL
  })
  main_poly_df <- do.call(rbind, main_poly_list)

  sub_poly_list <- lapply(tm@cells, function(cell) {
    if (cell$name %in% sub_labels) {
      extract_poly(cell, cell$name)
    } else NULL
  })
  sub_poly_df <- do.call(rbind, sub_poly_list)

  if (is.null(sub_poly_df) || nrow(sub_poly_df) == 0) {
    stop("No sub-type polygons could be extracted from the treemap.")
  }

  sub_poly_df$Cell_type <- sub("\n.*$", "", sub_poly_df$name)

  centres <- do.call(rbind, lapply(tm@cells, function(cell) {
    if (cell$name %in% sub_labels) {
      mat <- cell$poly[[1]]
      if (is.null(mat) || nrow(mat) == 0) return(NULL)
      data.frame(
        x = mean(mat[, 1]),
        y = mean(mat[, 2]),
        name = cell$name,
        Cell_type = sub("\n.*$", "", cell$name),
        stringsAsFactors = FALSE
      )
    } else NULL
  }))

  # Compute convex hull for outer frame
  if (!is.null(main_poly_df) && nrow(main_poly_df) > 2) {
    hull_idx <- grDevices::chull(main_poly_df$x, main_poly_df$y)
    hull_df <- main_poly_df[hull_idx, ]
  } else {
    hull_idx <- grDevices::chull(sub_poly_df$x, sub_poly_df$y)
    hull_df <- sub_poly_df[hull_idx, ]
  }

  p <- ggplot2::ggplot() +
    # Sub-type fill (no border)
    ggplot2::geom_polygon(
      data = sub_poly_df,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$name, fill = .data$Cell_type),
      colour = NA
    ) +
    # Sub-type thin border
    ggplot2::geom_polygon(
      data = sub_poly_df,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$name),
      fill = NA, colour = border_color, linewidth = subtype_border_lwd
    ) +
    # Main thick border (drawn on top)
    ggplot2::geom_polygon(
      data = main_poly_df,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$name),
      fill = NA, colour = border_color, linewidth = main_border_lwd
    ) +
    # Outer convex hull border
    ggplot2::geom_polygon(
      data = hull_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = NA, colour = border_color, linewidth = outer_border_lwd
    ) +
    ggplot2::geom_text(
      data = centres,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$name),
      size = label_size, colour = label_color, fontface = label_fontface,
      lineheight = 0.9
    ) +
    ggplot2::scale_fill_manual(
      values = pal_cell,
      breaks = present_cell,
      name = "Cell type"
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = if (legend) legend_position else "none",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )

  print(p)

  invisible(list(
    voronoi_treemap = tm,
    plot = p
  ))
}