#' Generate a Weighted Voronoi Treemap for Hierarchical Cell Types
#'
#' @description
#' Creates a weighted Voronoi treemap displaying the hierarchical composition
#' of single‑cell data. The area of each polygon is proportional to the
#' number of cells in the corresponding sub‑type, and the polygons are
#' grouped by the main cell type. Each sub‑type is filled with a distinct
#' colour matching the ArchR or custom palette.
#'
#' @param seurat_obj A Seurat object.
#' @param Main_cell_types Character string naming a column in
#'   \code{seurat_obj@meta.data} that holds the top‑level (main) cell‑type labels.
#' @param Cell_types Character string naming a column in
#'   \code{seurat_obj@meta.data} that holds the sub‑type labels.
#' @param label_type Character string specifying the type of label to display
#'   inside each polygon. One of \code{"both"}, \code{"count"}, \code{"percentage"},
#'   or \code{"none"}. Default is \code{"both"}.
#' @param shape Character string specifying the shape of the polygons.
#'   \code{"rounded_rect"} (default) or \code{"circle"}.
#' @param seed Integer seed for reproducibility of the Voronoi layout. Default is \code{1}.
#' @param col_Cell_types Optional named or unnamed character vector of colours
#'   for the \code{Cell_types} categories.
#' @param label_size Numeric value controlling the size of the polygon labels. Default is \code{3}.
#' @param label_color Colour of the polygon labels. Default is \code{"black"}.
#' @param legend Logical indicating whether to display a colour legend. Default is \code{TRUE}.
#' @param legend_position Position of the legend when \code{legend = TRUE}. Default is \code{"bottom"}.
#'
#' @return
#' Invisibly, a list with \code{voronoi_treemap} and \code{plot}.
#'
#' @export
#' 
#' @family Section_5_Other_Functions_Provided
#'
#' @examples
#' \dontrun{
#' # Basic rounded‑rectangle treemap with both count and percentage labels
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
#' # Only sub‑type names, no legend
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
    label_size = 3,
    label_color = "black",
    legend = TRUE,
    legend_position = "bottom"
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
  main_vec <- as.character(meta[[Main_cell_types]])
  cell_vec <- as.character(meta[[Cell_types]])

  valid_char <- function(x) !is.na(x) & x != ""
  keep <- valid_char(main_vec) & valid_char(cell_vec)
  main_vec <- main_vec[keep]
  cell_vec <- cell_vec[keep]

  tab <- as.data.frame(table(Main = main_vec, Cell = cell_vec))
  tab <- tab[tab$Freq > 0, ]
  colnames(tab) <- c("Main_cell_type", "Cell_type", "Cell_number")
  tab$Percentage <- round(tab$Cell_number / sum(tab$Cell_number) * 100, 2)

  default_pal <- function(n, names = NULL) {
    if (requireNamespace("scales", quietly = TRUE)) {
      cols <- scales::hue_pal()(n)
    } else {
      cols <- grDevices::rainbow(n)
    }
    if (!is.null(names)) setNames(cols, names) else cols
  }

  cell_cats <- sort(unique(tab$Cell_type))

  get_palette <- function(cats, user_cols) {
    if (!is.null(user_cols)) {
      if (!is.null(names(user_cols))) {
        pal <- setNames(rep(NA_character_, length(cats)), cats)
        common <- intersect(cats, names(user_cols))
        pal[common] <- user_cols[common]
        missing <- cats[is.na(pal)]
        if (length(missing) > 0) {
          message("Some categories lack user colours, auto-generating for: ",
                  paste(missing, collapse = ", "))
          pal[missing] <- default_pal(length(missing))
        }
        return(pal)
      } else {
        if (length(user_cols) < length(cats)) {
          warning("'user_cols' shorter than categories; colours will be recycled.")
          user_cols <- rep(user_cols, length.out = length(cats))
        }
        return(setNames(user_cols[seq_along(cats)], cats))
      }
    }

    if (requireNamespace("ArchR", quietly = TRUE)) {
      cols <- tryCatch(ArchR::paletteDiscrete(values = cats), error = function(e) NULL)
      if (!is.null(cols)) return(cols)
    }

    message("ArchR is not available; falling back to default ggplot2 palette.")
    default_pal(length(cats), cats)
  }

  pal_cell <- get_palette(cell_cats, col_Cell_types)

  tab$Label <- switch(
    label_type,
    both       = paste0(tab$Cell_type, "\n", tab$Cell_number, "\n", tab$Percentage, "%"),
    count      = paste0(tab$Cell_type, "\n", tab$Cell_number),
    percentage = paste0(tab$Cell_type, "\n", tab$Percentage, "%"),
    none       = tab$Cell_type
  )

  pal_label <- pal_cell[tab$Cell_type]
  names(pal_label) <- tab$Label

  tm <- WeightedTreemaps::voronoiTreemap(
    data      = tab,
    levels    = c("Main_cell_type", "Label"),
    cell_size = "Cell_number",
    shape     = shape,
    seed      = seed
  )

  p <- WeightedTreemaps::drawTreemap(
    treemap         = tm,
    color_palette   = pal_label,
    color_level     = 2,
    label_level     = 2,
    label_size      = label_size,
    label_color     = label_color,
    legend          = legend,
    legend_position = legend_position
  )

  print(p)

  invisible(list(
    voronoi_treemap = tm,
    plot            = p
  ))
}
