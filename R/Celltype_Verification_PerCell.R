#' Verify per-cell annotations with marker expression dotplot
#'
#' @description This function verifies per-cell SlimR annotations by generating a dotplot
#'     showing marker gene expression across predicted cell types.
#'
#' @param seurat_obj A Seurat object with per-cell annotations.
#' @param SlimR_percell_result A list from Celltype_Calculate_PerCell() containing
#'     Expression_list with marker genes per cell type.
#' @param assay Assay to use. Default: "RNA".
#' @param gene_number Number of top genes to show per cell type. Default: 5.
#' @param annotation_col Column in meta.data with cell type annotations. 
#'     Default: "Cell_type_PerCell_SlimR".
#' @param colour_low Color for lowest expression. Default: "white".
#' @param colour_high Color for highest expression. Default: "navy".
#' @param min_cells Minimum number of cells required for a cell type to be included
#'     in the plot. Default: 10.
#'
#' @return A ggplot object showing marker gene expression dotplot.
#'
#' @export
#' @family Section_3_Automated_Annotation
#'
#' @importFrom Seurat DotPlot Idents FetchData FindMarkers DefaultAssay
#' @importFrom dplyr top_n pull
#' @importFrom ggplot2 theme_bw element_blank element_text labs scale_color_gradientn ggtitle
#'
#' @examples
#' \dontrun{
#' # After running Celltype_Calculate_PerCell and Celltype_Annotation_PerCell
#' dotplot <- Celltype_Verification_PerCell(
#'     seurat_obj = sce,
#'     SlimR_percell_result = result,
#'     gene_number = 5,
#'     annotation_col = "Cell_type_PerCell_SlimR"
#' )
#' print(dotplot)
#' }
#'
Celltype_Verification_PerCell <- function(
    seurat_obj,
    SlimR_percell_result,
    assay = "RNA",
    gene_number = 5,
    colour_low = "white",
    colour_high = "navy",
    annotation_col = "Cell_type_PerCell_SlimR",
    min_cells = 10
) {
  
  # ============================================================================
  # Validate inputs
  # ============================================================================
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  
  if (!is.list(SlimR_percell_result)) {
    stop("SlimR_percell_result must be a list")
  }
  
  if (!"Expression_list" %in% names(SlimR_percell_result)) {
    stop("Expression_list not found in SlimR_percell_result")
  }
  
  if (!is.numeric(gene_number) || gene_number < 1) {
    stop("gene_number must be a positive integer")
  }
  
  if (!(annotation_col %in% colnames(seurat_obj@meta.data))) {
    stop(paste0(annotation_col, " not found in seurat_obj meta.data. ",
                "Run Celltype_Annotation_PerCell() first."))
  }
  
  Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[annotation_col]]
  assay <- if (is.null(assay)) Seurat::DefaultAssay(seurat_obj) else assay
  Seurat::DefaultAssay(seurat_obj) <- assay
  
  # ============================================================================
  # Get cell types and filter by minimum cell count
  # ============================================================================
  cell_type_counts <- table(seurat_obj@meta.data[[annotation_col]])
  valid_cell_types <- names(cell_type_counts[cell_type_counts >= min_cells])
  valid_cell_types <- valid_cell_types[valid_cell_types != "Unassigned"]
  
  if (length(valid_cell_types) == 0) {
    stop("No valid cell types with at least ", min_cells, " cells found")
  }
  
  message(sprintf("SlimR PerCell Verification: Verifying %d cell types with >= %d cells",
                  length(valid_cell_types), min_cells))
  
  # ============================================================================
  # Select top markers for each cell type
  # ============================================================================
  feature_list <- list()
  
  for (i in seq_along(valid_cell_types)) {
    cell_type <- valid_cell_types[i]
    message(sprintf("[%d/%d] Selecting markers for: %s", i, length(valid_cell_types), cell_type))
    
    # First try to use markers from SlimR_percell_result$Expression_list
    if (cell_type %in% names(SlimR_percell_result$Expression_list)) {
      
      expr_mat <- SlimR_percell_result$Expression_list[[cell_type]]
      prop_mat <- SlimR_percell_result$Proportion_list[[cell_type]]
      
      if (!is.null(expr_mat) && !is.null(prop_mat)) {
        markers <- colnames(expr_mat)
        
        if (length(markers) > 0) {
          # Recalculate expression stats for ranking
          ct_cells <- seurat_obj@meta.data[[annotation_col]] == cell_type
          other_cells <- seurat_obj@meta.data[[annotation_col]] != cell_type & 
                        seurat_obj@meta.data[[annotation_col]] != "Unassigned"
          
          if (sum(ct_cells) > 0 && sum(other_cells) > 0) {
            ct_expr <- Seurat::FetchData(seurat_obj, vars = markers, cells = colnames(seurat_obj)[ct_cells])
            other_expr <- Seurat::FetchData(seurat_obj, vars = markers, cells = colnames(seurat_obj)[other_cells])
            
            # Calculate mean expression
            ct_mean <- apply(ct_expr, 2, function(x) mean(expm1(x)))
            other_mean <- apply(other_expr, 2, function(x) mean(expm1(x)))
            
            # Calculate proportion expressing
            ct_prop <- apply(ct_expr, 2, function(x) mean(x > 0.1))
            
            # Calculate log2FC
            other_mean[other_mean == 0] <- 1e-5
            log2fc <- log2(ct_mean / other_mean)
            
            # Score: log2FC * proportion
            scores <- log2fc * ct_prop
            scores[is.nan(scores) | is.infinite(scores)] <- 0
            
            # Select top genes and filter out NAs
            top_genes <- names(sort(scores, decreasing = TRUE))[1:min(gene_number, length(scores))]
            top_genes <- top_genes[!is.na(top_genes)]
            feature_list[[cell_type]] <- top_genes
            
            next
          }
        }
      }
    }
    
    # Fallback: use FindMarkers
    message(sprintf("  Using FindMarkers for %s", cell_type))
    
    tryCatch({
      markers <- Seurat::FindMarkers(
        seurat_obj, 
        ident.1 = cell_type, 
        only.pos = TRUE,
        min.pct = 0.1,
        logfc.threshold = 0.25
      )
      
      if (nrow(markers) > 0) {
        markers$gene <- rownames(markers)
        markers$score <- markers$avg_log2FC * markers$pct.1
        
        top_markers <- markers %>%
          dplyr::top_n(gene_number, score) %>%
          dplyr::pull(gene)
        
        feature_list[[cell_type]] <- top_markers
      }
    }, error = function(e) {
      warning(sprintf("Could not find markers for %s: %s", cell_type, e$message))
    })
  }
  
  # ============================================================================
  # Generate dotplot
  # ============================================================================
  all_features <- unique(unlist(feature_list))
  # Remove any NA values that might have slipped through
  all_features <- all_features[!is.na(all_features)]
  
  if (length(all_features) == 0) {
    stop("No valid features found for verification")
  }
  
  message(sprintf("\nGenerating dotplot with %d features across %d cell types",
                  length(all_features), length(valid_cell_types)))
  
  # Filter Seurat object to valid cell types only
  cells_to_plot <- seurat_obj@meta.data[[annotation_col]] %in% valid_cell_types
  seurat_subset <- seurat_obj[, cells_to_plot]
  Seurat::Idents(seurat_subset) <- seurat_subset@meta.data[[annotation_col]]
  
  # Suppress Seurat warnings about scaling (expected with minimal cell types)
  dotplot <- suppressWarnings(
    Seurat::DotPlot(
      object = seurat_subset,
      features = all_features,
      assay = assay,
      group.by = annotation_col
    )
  ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = "Per-Cell Annotation Verification | SlimR"
    ) +
    ggplot2::scale_color_gradientn(
      colours = c(colour_low, colour_high),
      values = seq(0, 1, length.out = 2)
    )
  
  message("SlimR PerCell Verification: Dotplot generated successfully")
  
  return(dotplot)
}
