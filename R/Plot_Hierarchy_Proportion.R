#' Plot a Hierarchical Cell Type Tree with Optional Group Proportion Heatmap
#'
#' @description
#' Creates a publication‑quality composite figure with a hierarchical tree
#' panel and an optional proportion heatmap.  Leaf labels are connected to the
#' terminal nodes by short vertical sticks and placed at a fixed distance
#' below the tree, guaranteeing a clean separation and perfect top alignment
#' of all cell‑type names.  Short vertical ticks on top of the heatmap mark
#' each column, aligning with the labels above.  The tree height is independent
#' of label length; the heatmap sits tightly beneath the labels.
#'
#' Colours can be supplied for each annotation level.  When colours are
#' missing, the function tries to obtain a perceptually uniform palette via
#' \code{\link[ArchR:paletteDiscrete]{ArchR::paletteDiscrete()}}; if
#' \pkg{ArchR} is unavailable, it falls back to the default \pkg{ggplot2}
#' hue palette and prints an informative message.
#'
#' @param seurat_obj A Seurat object.
#' @param Main_cell_types Character string naming a column in
#'   \code{seurat_obj@meta.data} that holds the top‑level cell‑type labels.
#' @param col_Main_cell_types Optional named or unnamed character vector of
#'   colours for the \code{Main_cell_types} categories.  If \code{NULL},
#'   colours are generated automatically.
#' @param Cell_types Character string naming a column with second‑level
#'   (sub‑type) labels.  Use \code{NULL} if no second level exists.
#' @param col_Cell_types Optional colours for the \code{Cell_types} categories.
#' @param Sub_cell_types Character string naming a column with third‑level
#'   (sub‑sub‑type) labels.  Use \code{NULL} if absent.  \code{Cell_types}
#'   must be supplied when \code{Sub_cell_types} is used.
#' @param col_Sub_cell_types Optional colours for the \code{Sub_cell_types}
#'   categories.
#' @param proportion Logical scalar.  If \code{TRUE} (the default), a
#'   proportion heatmap is drawn below the tree; if \code{FALSE}, only the
#'   tree with labels is returned.
#' @param Groups Character string naming a \code{meta.data} column that
#'   defines sample groups for the proportion heatmap.  Required when
#'   \code{proportion = TRUE}.
#' @param show_labels Logical scalar. If \code{TRUE} (the default), text
#'   labels are placed next to the non‑leaf Main and Cell level nodes.
#' @param low_col Character string for the lowest proportion colour in the
#'   heatmap. Default is \code{"white"}.
#' @param high_col Character string for the highest proportion colour in the
#'   heatmap. Default is \code{"red"}.
#'
#' @return
#' Invisibly, a list with the following components:
#' \describe{
#'   \item{tree_plot}{A \pkg{ggplot2} object with the hierarchical tree,
#'     leaf sticks, and labels.}
#'   \item{prop_plot}{If \code{proportion = TRUE}, a \pkg{ggplot2} heatmap with
#'     top ticks; otherwise \code{NULL}.}
#'   \item{combined_plot}{If \code{proportion = TRUE} and \pkg{patchwork} is
#'     installed, a combined plot with title \emph{Hierarchical Proportion Plot | SlimR};
#'     otherwise \code{NULL}.}
#' }
#' The function automatically prints the \code{combined_plot} (or the
#' \code{tree_plot} when no proportion is requested) to the active graphics
#' device.
#'
#' @export
#'
#' @family Section_5_Other_Functions_Provided
#'
#' @importFrom stats ave aggregate
#' @importFrom grDevices rainbow
#' @importFrom ggplot2 ggplot aes geom_segment geom_point geom_text geom_rect
#' @importFrom ggplot2 scale_size_area scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 expansion theme_void theme element_text element_blank margin
#' @importFrom ggplot2 coord_cartesian guide_legend geom_tile scale_fill_gradient
#' @importFrom ggplot2 element_rect
#' @importFrom scales hue_pal
#'
#' @examples
#' \dontrun{
#' # Full three-level hierarchy with proportion heatmap
#' res <- Plot_Hierarchy_Proportion(
#'   seurat_obj          = sce,
#'   Main_cell_types     = "Main_cell_type",
#'   Cell_types          = "Cell_type",
#'   Sub_cell_types      = "Sub_cell_type",
#'   proportion          = TRUE,
#'   Groups              = "orig.ident"
#' )
#'
#' # Only the tree, user-specified colours
#' res <- Plot_Hierarchy_Proportion(
#'   sce,
#'   "Main_cell_type",
#'   col_Main_cell_types = c(Immune = "red", Stromal = "blue", Epithelial = "green"),
#'   "Cell_type",
#'   proportion = FALSE
#' )
#' }
Plot_Hierarchy_Proportion <- function(
    seurat_obj,
    Main_cell_types,
    col_Main_cell_types = NULL,
    Cell_types = NULL,
    col_Cell_types = NULL,
    Sub_cell_types = NULL,
    col_Sub_cell_types = NULL,
    proportion = TRUE,
    Groups = NULL,
    show_labels = TRUE,
    low_col = "white",
    high_col = "red"
) {
  # --- argument checks -------------------------------------------------------
  if (!inherits(seurat_obj, "Seurat")) {
    stop("'seurat_obj' must be a Seurat object.")
  }
  stopifnot(is.character(Main_cell_types), length(Main_cell_types) == 1L)
  if (!Main_cell_types %in% colnames(seurat_obj@meta.data)) {
    stop("'Main_cell_types' not found in seurat_obj@meta.data.")
  }
  stopifnot(is.logical(proportion), length(proportion) == 1L)
  stopifnot(is.logical(show_labels), length(show_labels) == 1L)

  has_Cell <- !is.null(Cell_types)
  if (has_Cell) {
    stopifnot(is.character(Cell_types), length(Cell_types) == 1L)
    if (!Cell_types %in% colnames(seurat_obj@meta.data)) {
      stop("'Cell_types' not found in seurat_obj@meta.data.")
    }
  }
  has_Sub <- !is.null(Sub_cell_types)
  if (has_Sub) {
    stopifnot(is.character(Sub_cell_types), length(Sub_cell_types) == 1L)
    if (!has_Cell) {
      stop("'Sub_cell_types' requires 'Cell_types' to be provided as well.")
    }
    if (!Sub_cell_types %in% colnames(seurat_obj@meta.data)) {
      stop("'Sub_cell_types' not found in seurat_obj@meta.data.")
    }
  }

  if (proportion) {
    if (is.null(Groups)) stop("'Groups' is required when proportion = TRUE.")
    if (!Groups %in% colnames(seurat_obj@meta.data)) {
      stop("'Groups' not found in seurat_obj@meta.data.")
    }
  }

  # --- extract and clean meta data -------------------------------------------
  meta <- seurat_obj@meta.data
  main_vec <- as.character(meta[[Main_cell_types]])
  if (has_Cell) {
    cell_vec <- as.character(meta[[Cell_types]])
    cell_vec[is.na(cell_vec)] <- ""
  }
  if (has_Sub) {
    sub_vec <- as.character(meta[[Sub_cell_types]])
    sub_vec[is.na(sub_vec)] <- ""
  }

  valid_char <- function(x) !is.na(x) & x != ""

  # --- determine terminal (leaf) labels and leaf ID for each cell ------------
  leaf_id <- rep(NA_character_, nrow(meta))
  if (has_Sub) {
    idx_sub <- valid_char(sub_vec)
    leaf_id[idx_sub] <- paste0("Sub:", sub_vec[idx_sub])
    idx_cell <- valid_char(cell_vec) & is.na(leaf_id)
    leaf_id[idx_cell] <- paste0("Cell:", cell_vec[idx_cell])
    idx_main <- is.na(leaf_id)
    leaf_id[idx_main] <- paste0("Main:", main_vec[idx_main])
  } else if (has_Cell) {
    idx_cell <- valid_char(cell_vec)
    leaf_id[idx_cell] <- paste0("Cell:", cell_vec[idx_cell])
    idx_main <- is.na(leaf_id)
    leaf_id[idx_main] <- paste0("Main:", main_vec[idx_main])
  } else {
    leaf_id[] <- paste0("Main:", main_vec)
  }

  # --- build node table ------------------------------------------------------
  main_cats <- unique(main_vec[valid_char(main_vec)])
  if (has_Cell) cell_cats <- unique(cell_vec[valid_char(cell_vec)]) else cell_cats <- character(0)
  if (has_Sub)  sub_cats  <- unique(sub_vec[valid_char(sub_vec)])  else sub_cats  <- character(0)

  if (anyDuplicated(main_cats)) stop("Main cell type labels must be unique.")
  if (has_Cell && anyDuplicated(cell_cats)) stop("Cell type labels must be unique within this level.")
  if (has_Sub  && anyDuplicated(sub_cats))  stop("Sub cell type labels must be unique within this level.")

  edges <- list()
  if (has_Cell) {
    sel <- valid_char(cell_vec)
    edges[["Main_Cell"]] <- unique(data.frame(
      from = main_vec[sel],
      to   = cell_vec[sel],
      stringsAsFactors = FALSE
    ))
  }
  if (has_Sub) {
    sel <- valid_char(sub_vec)
    edges[["Cell_Sub"]] <- unique(data.frame(
      from = cell_vec[sel],
      to   = sub_vec[sel],
      stringsAsFactors = FALSE
    ))
  }

  node_list <- list()
  for (cat in main_cats) {
    node_id <- paste0("Main:", cat)
    children <- if (has_Cell) {
      child_labs <- edges[["Main_Cell"]]$to[edges[["Main_Cell"]]$from == cat]
      paste0("Cell:", child_labs)
    } else character(0)
    node_list[[node_id]] <- list(
      id       = node_id,
      label    = cat,
      layer    = "Main",
      is_leaf  = length(children) == 0,
      children = children,
      n        = sum(main_vec == cat, na.rm = TRUE)
    )
  }
  if (has_Cell) {
    for (cat in cell_cats) {
      node_id <- paste0("Cell:", cat)
      children <- if (has_Sub) {
        child_labs <- edges[["Cell_Sub"]]$to[edges[["Cell_Sub"]]$from == cat]
        paste0("Sub:", child_labs)
      } else character(0)
      node_list[[node_id]] <- list(
        id       = node_id,
        label    = cat,
        layer    = "Cell",
        is_leaf  = length(children) == 0,
        children = children,
        n        = sum(cell_vec == cat, na.rm = TRUE)
      )
    }
  }
  if (has_Sub) {
    for (cat in sub_cats) {
      node_id <- paste0("Sub:", cat)
      node_list[[node_id]] <- list(
        id       = node_id,
        label    = cat,
        layer    = "Sub",
        is_leaf  = TRUE,
        children = character(0),
        n        = sum(sub_vec == cat, na.rm = TRUE)
      )
    }
  }

  # --- y coordinates (equal vertical spacing) --------------------------------
  if (has_Sub) {
    y_levels <- c(Main = 2.4, Cell = 1.2, Sub = 0)
  } else if (has_Cell) {
    y_levels <- c(Main = 1.2, Cell = 0)
  } else {
    y_levels <- c(Main = 0)
  }

  for (nm in names(node_list)) {
    nd <- node_list[[nm]]
    node_list[[nm]]$y <- y_levels[nd$layer]
  }

  # --- x coordinates (deterministic layout) ----------------------------------
  x_counter <- 1
  main_order <- sort(main_cats)
  for (mcat in main_order) {
    main_id <- paste0("Main:", mcat)
    main_nd <- node_list[[main_id]]
    if (main_nd$is_leaf) {
      node_list[[main_id]]$x <- x_counter
      x_counter <- x_counter + 1
    } else {
      child_ids <- main_nd$children
      cell_labs <- sapply(node_list[child_ids], `[[`, "label")
      child_ids <- child_ids[order(cell_labs)]
      for (cid in child_ids) {
        cell_nd <- node_list[[cid]]
        if (cell_nd$is_leaf) {
          node_list[[cid]]$x <- x_counter
          x_counter <- x_counter + 1
        } else {
          sub_ids <- cell_nd$children
          sub_labs <- sapply(node_list[sub_ids], `[[`, "label")
          sub_ids <- sub_ids[order(sub_labs)]
          for (sid in sub_ids) {
            node_list[[sid]]$x <- x_counter
            x_counter <- x_counter + 1
          }
          node_list[[cid]]$x <- mean(sapply(sub_ids, function(s) node_list[[s]]$x))
        }
      }
      node_list[[main_id]]$x <- mean(sapply(child_ids, function(c) node_list[[c]]$x))
    }
  }

  leaf_nodes <- node_list[sapply(node_list, `[[`, "is_leaf")]
  leaf_order <- order(sapply(leaf_nodes, `[[`, "x"))
  leaf_ids_ordered <- names(leaf_nodes)[leaf_order]
  leaf_labels_ordered <- sapply(leaf_nodes, `[[`, "label")[leaf_order]
  leaf_x_ordered <- sapply(leaf_nodes, `[[`, "x")[leaf_order]
  N_leaves <- length(leaf_nodes)

  # --- colours ---------------------------------------------------------------
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
      cols <- tryCatch(ArchR::paletteDiscrete(cats), error = function(e) NULL)
      if (!is.null(cols)) return(setNames(cols, cats))
    }
    message("ArchR is not available; falling back to default ggplot2 palette. ",
            "Install ArchR for a wider colour range.")
    default_pal(length(cats), cats)
  }

  default_pal <- function(n, names = NULL) {
    if (requireNamespace("scales", quietly = TRUE)) {
      cols <- scales::hue_pal()(n)
    } else {
      cols <- grDevices::rainbow(n)
    }
    if (!is.null(names)) setNames(cols, names) else cols
  }

  pal_main <- get_palette(main_cats, col_Main_cell_types)
  pal_cell <- if (has_Cell) get_palette(cell_cats, col_Cell_types) else NULL
  pal_sub  <- if (has_Sub)  get_palette(sub_cats, col_Sub_cell_types)  else NULL

  for (nm in names(node_list)) {
    nd <- node_list[[nm]]
    if (nd$layer == "Main") {
      node_list[[nm]]$fill <- pal_main[nd$label]
    } else if (nd$layer == "Cell") {
      node_list[[nm]]$fill <- pal_cell[nd$label]
    } else {
      node_list[[nm]]$fill <- pal_sub[nd$label]
    }
  }

  # --- prepare nodes & three‑segment edges -----------------------------------
  nodes_df <- do.call(rbind, lapply(node_list, function(nd) {
    data.frame(x = nd$x, y = nd$y, n = nd$n, fill = nd$fill,
               label = nd$label, layer = nd$layer, is_leaf = nd$is_leaf,
               stringsAsFactors = FALSE)
  }))

  edge_segments <- list()
  for (e_type in names(edges)) {
    edf <- edges[[e_type]]
    for (i in seq_len(nrow(edf))) {
      from_lab <- edf$from[i]
      to_lab   <- edf$to[i]
      if (e_type == "Main_Cell") {
        from_id <- paste0("Main:", from_lab)
        to_id   <- paste0("Cell:", to_lab)
      } else {
        from_id <- paste0("Cell:", from_lab)
        to_id   <- paste0("Sub:", to_lab)
      }
      from_nd <- node_list[[from_id]]
      to_nd   <- node_list[[to_id]]
      if (is.null(from_nd) || is.null(to_nd)) next
      xp <- from_nd$x; yp <- from_nd$y
      xc <- to_nd$x;   yc <- to_nd$y
      ymid <- yp - (yp - yc) / 2
      edge_segments[[length(edge_segments) + 1]] <- data.frame(
        x = xp, xend = xp, y = yp, yend = ymid,
        group = paste0("edge_v1_", e_type, "_", i)
      )
      if (abs(xp - xc) > 1e-6) {
        edge_segments[[length(edge_segments) + 1]] <- data.frame(
          x = xp, xend = xc, y = ymid, yend = ymid,
          group = paste0("edge_h_", e_type, "_", i)
        )
      }
      edge_segments[[length(edge_segments) + 1]] <- data.frame(
        x = xc, xend = xc, y = ymid, yend = yc,
        group = paste0("edge_v2_", e_type, "_", i)
      )
    }
  }
  edges_seg_df <- do.call(rbind, edge_segments)

  # --- unified x limits ------------------------------------------------------
  xlim_use <- c(0.5, N_leaves + 0.5)

  # --- 1. Tree panel with sticks and aligned labels --------------------------
  max_y <- if (has_Sub) 2.4 else if (has_Cell) 1.2 else 0

  # Fixed stick length and label offset
  stick_len <- 0.10
  label_offset <- 0.02   # tiny gap after stick before text

  # Short stick from leaf node (y=0) down
  leaf_stick_df <- data.frame(
    x    = leaf_x_ordered,
    xend = leaf_x_ordered,
    y    = 0,
    yend = -stick_len,
    stringsAsFactors = FALSE
  )
  # Labels placed just below stick end, top-aligned
  leaf_label_df <- data.frame(
    x     = leaf_x_ordered,
    y     = -(stick_len + label_offset),
    label = leaf_labels_ordered,
    stringsAsFactors = FALSE
  )

  # Dynamic lower limit: account for the longest label rotated 90°
  max_nchar <- max(nchar(leaf_labels_ordered))
  char_height <- 0.09   # data units per character at size=3 when rotated
  label_height <- max_nchar * char_height
  y_lower <- -(stick_len + label_offset + label_height + 0.05)

  p_tree <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges_seg_df,
      ggplot2::aes(x = .data$x, y = .data$y,
                   xend = .data$xend, yend = .data$yend,
                   group = .data$group),
      colour = "black", linewidth = 0.3
    ) +
    ggplot2::geom_point(
      data = nodes_df,
      ggplot2::aes(x = .data$x, y = .data$y, size = .data$n, fill = I(.data$fill)),
      shape = 21, colour = "grey20", stroke = 0.3
    ) +
    # Leaf sticks
    ggplot2::geom_segment(
      data = leaf_stick_df,
      ggplot2::aes(x = .data$x, xend = .data$xend,
                   y = .data$y, yend = .data$yend),
      colour = "black", linewidth = 0.3
    ) +
    # Leaf labels (top-aligned, rotated)
    ggplot2::geom_text(
      data = leaf_label_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      angle = 90, vjust = 0, hjust = 0.5, size = 3
    ) +
    ggplot2::scale_size_area(
      max_size = 10,
      guide = ggplot2::guide_legend(title = "Cell count")
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = xlim_use, ylim = c(y_lower, max_y + 0.4),
                             clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      plot.margin = ggplot2::margin(5, 5, 0, 5)
    )

  # Non-leaf labels (if requested)
  if (show_labels) {
    label_nodes <- nodes_df[!nodes_df$is_leaf & nodes_df$layer %in% c("Main", "Cell"), ]
    if (nrow(label_nodes) > 0) {
      p_tree <- p_tree +
        ggplot2::geom_text(
          data = label_nodes,
          ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
          nudge_y = 0.2, size = 3, fontface = "italic"
        )
    }
  }

  # --- 2. Heatmap panel (with top ticks) ------------------------------------
  p_prop <- NULL
  if (proportion) {
    grp_vec <- as.character(meta[[Groups]])
    keep <- valid_char(grp_vec)
    grp_vec <- grp_vec[keep]
    leaf_id_sub <- leaf_id[keep]
    tab <- stats::aggregate(
      list(Freq = rep(1, length(grp_vec))),
      list(Group = grp_vec,
           End   = factor(leaf_id_sub, levels = leaf_ids_ordered)),
      FUN = length
    )
    tab$Proportion <- stats::ave(tab$Freq, tab$Group, FUN = function(x) x / sum(x))
    tab$x_num <- as.numeric(tab$End)

    group_levels <- sort(unique(tab$Group))
    tab$y_num <- as.numeric(factor(tab$Group, levels = group_levels))
    n_groups <- length(group_levels)

    # Column ticks on top of heatmap
    heatmap_stick_df <- data.frame(
      x    = leaf_x_ordered,
      xend = leaf_x_ordered,
      y    = n_groups + 0.5,
      yend = n_groups + 0.85,
      stringsAsFactors = FALSE
    )
    border_df <- data.frame(
      xmin = 0.5,
      xmax = N_leaves + 0.5,
      ymin = 0.5,
      ymax = n_groups + 0.5
    )

    p_prop <- ggplot2::ggplot() +
      ggplot2::geom_tile(
        data = tab,
        ggplot2::aes(x = .data$x_num, y = .data$y_num, fill = .data$Proportion),
        width = 0.8, height = 0.8, colour = "white", linewidth = 0.3
      ) +
      ggplot2::geom_rect(
        data = border_df,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                     ymin = .data$ymin, ymax = .data$ymax),
        fill = NA, colour = "black", linewidth = 0.5
      ) +
      ggplot2::geom_segment(
        data = heatmap_stick_df,
        ggplot2::aes(x = .data$x, xend = .data$xend,
                     y = .data$y, yend = .data$yend),
        colour = "black", linewidth = 0.5
      ) +
      ggplot2::scale_fill_gradient(
        low = low_col, high = high_col,
        limits = c(0, 1), name = "Proportion"
      ) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(
        breaks = seq_len(n_groups),
        labels = group_levels,
        expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(xlim = xlim_use, ylim = c(0.5, n_groups + 0.9),
                               clip = "off") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 9, face = "bold"),
        axis.title = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(0, 5, 5, 5)
      )
  }

  # --- Combine panels --------------------------------------------------------
  combined <- NULL
  if (proportion && !is.null(p_prop)) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      tree_h <- if (has_Sub) 2.5 else if (has_Cell) 1.8 else 1.2
      combined <- p_tree / p_prop +
        patchwork::plot_layout(heights = c(tree_h, 1)) +
        patchwork::plot_annotation(
          title = "Hierarchical Proportion Plot | SlimR",
          theme = ggplot2::theme(plot.title = ggplot2::element_text(
            hjust = 0.5, face = "bold", size = 14))
        )
    } else {
      message("Package 'patchwork' is not installed. Returning separate plots.")
    }
  }

  # --- Output ----------------------------------------------------------------
  if (!is.null(combined)) {
    print(combined)
  } else {
    p_tree_titled <- p_tree +
      ggplot2::labs(title = "Hierarchical Proportion Plot | SlimR") +
      ggplot2::theme(plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 14))
    print(p_tree_titled)
    if (!is.null(p_prop)) {
      message("Displaying proportion heatmap separately.")
      print(p_prop)
    }
  }

  invisible(list(
    tree_plot      = p_tree,
    prop_plot      = p_prop,
    combined_plot  = combined
  ))
}