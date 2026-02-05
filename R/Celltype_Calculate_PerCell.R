#' Per-cell annotation using marker expression and optional UMAP spatial smoothing
#'
#' Unlike cluster-based annotation, this function assigns cell type labels to each
#' individual cell based on marker gene expression profiles. Optionally uses UMAP
#' coordinates to smooth predictions via k-nearest neighbor voting.
#'
#' @param seurat_obj Seurat object with normalized expression data.
#' @param gene_list A standardized marker list (same format as Celltype_Calculate).
#' @param species "Human" or "Mouse" for gene name formatting.
#' @param assay Assay to use (default: "RNA").
#' @param method Scoring method: "AUCell" (rank-based), "mean" (average expression),
#'   or "weighted" (expression * detection weighted). Default: "weighted".
#' @param min_expression Minimum expression threshold for detection. Default: 0.1.
#' @param use_umap_smoothing Logical. If TRUE, apply k-NN smoothing using UMAP
#'   coordinates to improve annotation consistency. Default: FALSE.
#' @param umap_reduction Name of UMAP reduction in Seurat object. Default: "umap".
#' @param k_neighbors Number of neighbors for UMAP smoothing. Default: 15.
#' @param smoothing_weight Weight for neighbor votes vs cell's own score (0-1).
#'   Higher values give more weight to neighbors. Default: 0.3.
#' @param min_score Minimum score threshold to assign a cell type. Cells below
#'   this threshold are labeled "Unassigned". Default: "auto" which adaptively
#'   sets the threshold based on number of cell types (1.5 / n_celltypes).
#'   Set to a numeric value (e.g., 0.1) to use a fixed threshold.
#' @param min_confidence Minimum confidence threshold. Cells with confidence below
#'   this value are labeled "Unassigned". Confidence is calculated as the ratio
#'   of max score to second-highest score. Default: 1.2 (max must be 20% higher
#'   than second). Set to 1.0 to disable confidence filtering.
#' @param return_scores If TRUE, return full score matrix. Default: FALSE.
#' @param ncores Number of cores for parallel processing. Default: 1.
#' @param chunk_size Number of cells to process per chunk (memory optimization).
#'   Default: 5000.
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @returns A list containing:
#' \itemize{
#'   \item Cell_annotations: Data frame with Cell_barcode, Predicted_cell_type, Max_score, Confidence
#'   \item Cell_confidence: Numeric vector of confidence scores per cell
#'   \item Summary: Summary table of cell type counts and percentages
#'   \item Expression_list: List of mean expression matrices per cell type (for verification)
#'   \item Proportion_list: List of detection proportion matrices per cell type
#'   \item Prediction_results: Summary data frame with per-cell-type statistics
#'   \item Probability_matrix: Full cell × cell_type probability matrix (normalized)
#'   \item Raw_score_matrix: Full cell × cell_type raw score matrix (before normalization)
#'   \item Parameters: List of parameters used including adaptive thresholds
#'   \item Cell_scores: (if return_scores=TRUE) Same as Probability_matrix
#' }
#'
#' @details
#' ## Scoring Methods
#'
#' **"weighted" (recommended)**: Combines normalized expression with detection rate.
#' For each cell and cell type: score = mean(expr_i * weight_i) where weight_i

#' is derived from the marker's specificity across the dataset.
#'
#' **"mean"**: Simple average of normalized marker expression. Fast but less
#' discriminative for overlapping marker sets.
#'
#' **"AUCell"**: Rank-based scoring similar to AUCell package. For each cell,
#' genes are ranked by expression, and the score is the proportion of marker
#' genes in the top X% of expressed genes. Robust to technical variation.
#'
#' ## UMAP Smoothing
#'
#' When `use_umap_smoothing = TRUE`, the function:
#' 1. Computes initial per-cell scores
#' 2. Finds k nearest neighbors in UMAP space for each cell
#' 3. Smooths scores by weighted averaging with neighbors
#' 4. Re-assigns cell types based on smoothed scores
#'
#' This helps reduce noise and improve consistency of annotations within
#' spatially coherent regions.
#'
#' @export
#' @family Section_3_Automated_Annotation
#'
#' @importFrom Seurat DefaultAssay FetchData Embeddings
#' @importFrom stats dist sd
#'
#' @examples
#' \dontrun{
#' # Basic per-cell annotation
#' result <- Celltype_Calculate_PerCell(
#'     seurat_obj = sce,
#'     gene_list = Markers_list,
#'     species = "Human",
#'     method = "weighted"
#' )
#'
#' # Add annotations to Seurat object
#' sce$Cell_type_PerCell <- result$Cell_annotations$Predicted_cell_type
#'
#' # With UMAP smoothing for more consistent annotations
#' result_smooth <- Celltype_Calculate_PerCell(
#'     seurat_obj = sce,
#'     gene_list = Markers_list,
#'     species = "Human",
#'     use_umap_smoothing = TRUE,
#'     k_neighbors = 20,
#'     smoothing_weight = 0.3
#' )
#' }
#'
Celltype_Calculate_PerCell <- function(
    seurat_obj,
    gene_list,
    species,
    assay = "RNA",
    method = c("weighted", "mean", "AUCell"),
    min_expression = 0.1,
    use_umap_smoothing = FALSE,
    umap_reduction = "umap",
    k_neighbors = 15,
    smoothing_weight = 0.3,
    min_score = "auto",
    min_confidence = 1.2,
    return_scores = FALSE,
    ncores = 1,
    chunk_size = 5000,
    verbose = TRUE
) {
  
  # ============================================================================
  # Input validation
  # ============================================================================
  method <- match.arg(method)
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object!")
  }
  if (!is.list(gene_list)) {
    stop("gene_list must be a list of data.frames!")
  }
  if (!species %in% c("Human", "Mouse")) {
    stop("species must be 'Human' or 'Mouse'")
  }
  if (use_umap_smoothing && !umap_reduction %in% names(seurat_obj@reductions)) {
    stop(paste0("UMAP reduction '", umap_reduction, "' not found. ",
                "Run Seurat::RunUMAP() first or set use_umap_smoothing = FALSE."))
  }
  if (smoothing_weight < 0 || smoothing_weight > 1) {
    stop("smoothing_weight must be between 0 and 1")
  }
  if (!is.character(min_score) && !is.numeric(min_score)) {
    stop("min_score must be 'auto' or a numeric value")
  }
  if (is.numeric(min_confidence) && min_confidence < 1) {
    stop("min_confidence must be >= 1.0 (ratio of max to second-max score)")
  }
  
  assay <- if (is.null(assay)) Seurat::DefaultAssay(seurat_obj) else assay
  Seurat::DefaultAssay(seurat_obj) <- assay
  
  # ============================================================================
  # Prepare marker gene sets
  # ============================================================================
  if (verbose) message("SlimR PerCell: Preparing marker gene sets...")
  
  cell_types <- names(gene_list)
  marker_sets <- list()
  all_markers <- character(0)
  
  for (cell_type in cell_types) {
    current_df <- gene_list[[cell_type]]
    if (ncol(current_df) < 1) next
    
    genes <- current_df[[1]]
    genes_processed <- if (species == "Human") {
      toupper(genes)
    } else {
      paste0(toupper(substr(genes, 1, 1)), tolower(substr(genes, 2, nchar(genes))))
    }
    
    # Filter to genes present in the assay
    valid_genes <- genes_processed[genes_processed %in% rownames(seurat_obj[[assay]])]
    valid_genes <- unique(valid_genes)
    
    if (length(valid_genes) > 0) {
      marker_sets[[cell_type]] <- valid_genes
      all_markers <- union(all_markers, valid_genes)
    }
  }
  
  if (length(marker_sets) == 0) {
    stop("No valid marker genes found in the Seurat object!")
  }
  
  n_celltypes <- length(marker_sets)
  n_cells <- ncol(seurat_obj)
  
  if (verbose) {
    message(sprintf("SlimR PerCell: %d cell types with valid markers, %d total cells",
                    n_celltypes, n_cells))
  }
  
  # ============================================================================
  # Extract expression data (chunked for memory efficiency)
  # ============================================================================
  if (verbose) message("SlimR PerCell: Extracting expression data...")
  
  # Get expression matrix for all marker genes
  expr_data <- Seurat::FetchData(seurat_obj, vars = all_markers, cells = colnames(seurat_obj))
  expr_matrix <- as.matrix(expr_data)
  
  # ============================================================================
  # Compute per-cell scores for each cell type
  # ============================================================================
  if (verbose) message(sprintf("SlimR PerCell: Computing scores using '%s' method...", method))
  
  score_matrix <- matrix(0, nrow = n_cells, ncol = n_celltypes)
  rownames(score_matrix) <- colnames(seurat_obj)
  colnames(score_matrix) <- names(marker_sets)
  
  if (method == "weighted") {
    # Weighted scoring: combines expression level with marker specificity
    score_matrix <- .compute_weighted_scores(
      expr_matrix = expr_matrix,
      marker_sets = marker_sets,
      min_expression = min_expression
    )
    
  } else if (method == "mean") {
    # Simple mean expression
    for (ct in names(marker_sets)) {
      markers <- marker_sets[[ct]]
      if (length(markers) == 1) {
        score_matrix[, ct] <- expr_matrix[, markers]
      } else {
        score_matrix[, ct] <- rowMeans(expr_matrix[, markers, drop = FALSE], na.rm = TRUE)
      }
    }
    
  } else if (method == "AUCell") {
    # Rank-based scoring (AUCell-like)
    score_matrix <- .compute_aucell_scores(
      expr_matrix = expr_matrix,
      marker_sets = marker_sets,
      top_percent = 0.05  # Top 5% of genes
    )
  }
  
  # Store raw scores before normalization
  raw_score_matrix <- score_matrix
  
  # Normalize scores per cell (so they sum to 1, like probabilities)
  row_sums <- rowSums(score_matrix)
  row_sums[row_sums == 0] <- 1  # Avoid division by zero
  score_matrix_norm <- score_matrix / row_sums
  
  # ============================================================================
  # Optional: UMAP-based spatial smoothing
  # ============================================================================
  if (use_umap_smoothing) {
    if (verbose) message(sprintf("SlimR PerCell: Applying UMAP smoothing (k=%d)...", k_neighbors))
    
    score_matrix_norm <- .apply_umap_smoothing(
      seurat_obj = seurat_obj,
      score_matrix = score_matrix_norm,
      umap_reduction = umap_reduction,
      k_neighbors = k_neighbors,
      smoothing_weight = smoothing_weight,
      chunk_size = chunk_size,
      verbose = verbose
    )
  }
  
  # ============================================================================
  # Calculate adaptive thresholds
  # ============================================================================
  
  # Adaptive min_score: based on number of cell types
  # With normalized scores summing to 1, uniform distribution gives 1/n_celltypes
  # We use 1.5x the uniform score as default threshold
  uniform_score <- 1 / n_celltypes
  if (is.character(min_score) && min_score == "auto") {
    effective_min_score <- 1.5 * uniform_score
    if (verbose) {
      message(sprintf("SlimR PerCell: Using adaptive min_score = %.4f (1.5 / %d cell types)",
                      effective_min_score, n_celltypes))
    }
  } else {
    effective_min_score <- as.numeric(min_score)
  }
  
  # ============================================================================
  # Assign cell types
  # ============================================================================
  if (verbose) message("SlimR PerCell: Assigning cell types...")
  
  # Compute max scores and indices, handling rows with all NA values
  max_scores <- apply(score_matrix_norm, 1, function(x) {
    # Check if all values are NA
    if (all(is.na(x))) return(0)
    m <- max(x, na.rm = TRUE)
    # If max returns -Inf (shouldn't happen after NA check), convert to 0
    if (is.infinite(m)) return(0)
    return(m)
  })
  
  max_indices <- apply(score_matrix_norm, 1, function(x) {
    idx <- which.max(x)
    # If which.max returns integer(0) (all NA), return NA_integer_
    if (length(idx) == 0) return(NA_integer_)
    return(idx)
  })
  predicted_types <- colnames(score_matrix_norm)[max_indices]
  
  # Mark cells with NA scores as Unassigned
  predicted_types[is.na(predicted_types)] <- "Unassigned"
  
  # Calculate confidence as RATIO of max to second-max (more robust than difference)
  # A ratio > 1.2 means the top choice is at least 20% better than second
  confidence <- apply(score_matrix_norm, 1, function(x) {
    sorted_x <- sort(x, decreasing = TRUE, na.last = TRUE)
    # Handle edge cases
    if (length(sorted_x) == 0 || is.na(sorted_x[1]) || sorted_x[1] == 0) return(0)
    if (length(sorted_x) == 1) return(Inf)
    if (sorted_x[2] == 0) return(Inf)
    # Confidence as ratio of top two scores
    sorted_x[1] / sorted_x[2]
  })
  
  # Also compute difference-based confidence for backward compatibility
  confidence_diff <- apply(score_matrix_norm, 1, function(x) {
    sorted_x <- sort(x, decreasing = TRUE, na.last = TRUE)
    if (length(sorted_x) < 2 || is.na(sorted_x[1])) return(0)
    sorted_x[1] - sorted_x[2]
  })
  
  # Apply thresholds
  # 1. Minimum score threshold
  predicted_types[max_scores < effective_min_score] <- "Unassigned"
  
  # 2. Minimum confidence threshold (ratio-based)
  if (!is.null(min_confidence) && min_confidence > 1) {
    low_confidence <- confidence < min_confidence & predicted_types != "Unassigned"
    if (sum(low_confidence) > 0 && verbose) {
      message(sprintf("  %d cells marked Unassigned due to low confidence (ratio < %.2f)",
                      sum(low_confidence), min_confidence))
    }
    predicted_types[low_confidence] <- "Unassigned"
  }
  
  # ============================================================================
  # Prepare output (compatible with Celltype_Annotation and Celltype_Verification)
  # ============================================================================
  cell_annotations <- data.frame(
    Cell_barcode = colnames(seurat_obj),
    Predicted_cell_type = predicted_types,
    Max_score = max_scores,
    Confidence = confidence,
    Confidence_diff = confidence_diff,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  # Summary statistics
  # Handle edge case where predicted_types might be empty (0 cells)
  tbl <- table(predicted_types)
  if (length(tbl) == 0) {
    summary_table <- data.frame(
      Cell_type = character(0),
      Count = integer(0),
      Percentage = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    summary_table <- as.data.frame.table(tbl, stringsAsFactors = FALSE)
    colnames(summary_table) <- c("Cell_type", "Count")
    summary_table$Percentage <- round(100 * summary_table$Count / n_cells, 2)
    summary_table <- summary_table[order(-summary_table$Count), ]
  }
  
  if (verbose) {
    message("\nSlimR PerCell: Annotation Summary:")
    message(sprintf("  Total cells: %d", n_cells))
    message(sprintf("  Assigned: %d (%.1f%%)", 
                    sum(predicted_types != "Unassigned"),
                    100 * sum(predicted_types != "Unassigned") / n_cells))
    message(sprintf("  Unassigned: %d (%.1f%%)",
                    sum(predicted_types == "Unassigned"),
                    100 * sum(predicted_types == "Unassigned") / n_cells))
  }
  
  # ============================================================================
  # Create Expression_list for compatibility with Celltype_Verification
  # ============================================================================
  expression_list <- list()
  proportion_list <- list()
  
  for (ct in names(marker_sets)) {
    markers <- marker_sets[[ct]]
    ct_cells <- which(predicted_types == ct)
    
    if (length(ct_cells) == 0) next
    
    # Calculate mean expression and proportion for this cell type's cells
    ct_expr_matrix <- expr_matrix[ct_cells, markers, drop = FALSE]
    
    # Mean expression (expm1 to reverse log1p transformation)
    mean_expr <- apply(ct_expr_matrix, 2, function(x) mean(expm1(x)))
    
    # Proportion expressing
    prop_expr <- apply(ct_expr_matrix, 2, function(x) mean(x > min_expression))
    
    # Store as single-row data frames (cell type name as row)
    expression_list[[ct]] <- matrix(mean_expr, nrow = 1, dimnames = list(ct, markers))
    proportion_list[[ct]] <- matrix(prop_expr, nrow = 1, dimnames = list(ct, markers))
  }
  
  # ============================================================================
  # Create pseudo Prediction_results for compatibility
  # Note: This is per-cell, so we create a summary table by cell type
  # ============================================================================
  prediction_results <- data.frame(
    Cell_type = names(marker_sets),
    Cell_count = sapply(names(marker_sets), function(ct) sum(predicted_types == ct, na.rm = TRUE)),
    Mean_score = sapply(names(marker_sets), function(ct) {
      ct_cells <- predicted_types == ct
      if (sum(ct_cells, na.rm = TRUE) == 0) return(0)
      mean(max_scores[ct_cells], na.rm = TRUE)
    }),
    Mean_confidence = sapply(names(marker_sets), function(ct) {
      ct_cells <- predicted_types == ct
      if (sum(ct_cells, na.rm = TRUE) == 0) return(0)
      mean(confidence[ct_cells], na.rm = TRUE)
    }),
    stringsAsFactors = FALSE
  )
  
  # Store parameters used for reproducibility
  params_used <- list(
    method = method,
    min_expression = min_expression,
    min_score_input = min_score,
    min_score_effective = effective_min_score,
    min_confidence = min_confidence,
    n_celltypes = n_celltypes,
    uniform_score = uniform_score,
    use_umap_smoothing = use_umap_smoothing,
    k_neighbors = k_neighbors,
    smoothing_weight = smoothing_weight
  )
  
  # Build return list
  result <- list(
    Cell_annotations = cell_annotations,
    Cell_confidence = confidence,
    Summary = summary_table,
    Expression_list = expression_list,
    Proportion_list = proportion_list,
    Prediction_results = prediction_results,
    Probability_matrix = score_matrix_norm,
    Raw_score_matrix = raw_score_matrix,
    Parameters = params_used
  )
  
  if (return_scores) {
    result$Cell_scores <- score_matrix_norm
  }
  
  return(result)
}


# ==============================================================================
# Internal helper functions
# ==============================================================================

#' Compute weighted scores for per-cell annotation
#' 
#' This function uses an improved weighting scheme that considers:
#' 1. Expression level (log-normalized)
#' 2. Detection rate (binary: above min_expression threshold)
#' 3. Marker specificity (how unique is this marker to this cell type)
#' 4. Expression variability (CV-based: more variable genes are more discriminative)
#' 
#' @keywords internal
.compute_weighted_scores <- function(expr_matrix, marker_sets, min_expression) {
  n_cells <- nrow(expr_matrix)
  n_celltypes <- length(marker_sets)
  
  score_matrix <- matrix(0, nrow = n_cells, ncol = n_celltypes)
  colnames(score_matrix) <- names(marker_sets)
  
  # Pre-compute global statistics for weighting
  all_markers <- unique(unlist(marker_sets))
  
  # ============================================================================
  # Compute marker specificity weights
  # A marker is more specific if it appears in fewer cell type marker sets
  # ============================================================================
  marker_counts <- table(unlist(marker_sets))
  specificity_weights <- 1 / as.numeric(marker_counts[all_markers])
  specificity_weights[is.na(specificity_weights)] <- 1
  names(specificity_weights) <- all_markers
  
  # ============================================================================
  # Detection rate per gene (used for IDF-like weighting)
  # Genes expressed in fewer cells are more discriminative
  # ============================================================================
  detection_rates <- colMeans(expr_matrix[, all_markers, drop = FALSE] > min_expression)
  # IDF-like weight: log(1 / detection_rate), capped to avoid Inf
  idf_weights <- log1p(1 / pmax(detection_rates, 0.01))
  idf_weights <- idf_weights / max(idf_weights, na.rm = TRUE)  # Normalize to [0, 1]
  names(idf_weights) <- all_markers
  
  # ============================================================================
  # Expression variability (CV) - genes with higher CV are more discriminative
  # ============================================================================
  gene_means <- colMeans(expr_matrix[, all_markers, drop = FALSE])
  gene_sds <- apply(expr_matrix[, all_markers, drop = FALSE], 2, sd)
  gene_cv <- ifelse(gene_means > 0, gene_sds / gene_means, 0)
  
  # Normalize CV to [0.5, 1.5] range for weighting
  cv_range <- max(gene_cv, na.rm = TRUE) - min(gene_cv, na.rm = TRUE)
  if (is.na(cv_range) || cv_range == 0) {
    cv_weights <- rep(1, length(gene_cv))
  } else {
    cv_weights <- 0.5 + (gene_cv - min(gene_cv, na.rm = TRUE)) / cv_range
  }
  names(cv_weights) <- all_markers
  
  # ============================================================================
  # Combine weights: specificity * IDF * CV
  # ============================================================================
  combined_weights <- specificity_weights * (0.5 + 0.5 * idf_weights) * cv_weights
  
  for (i in seq_along(marker_sets)) {
    ct <- names(marker_sets)[i]
    markers <- marker_sets[[ct]]
    
    if (length(markers) == 0) next
    
    # Get expression for this cell type's markers
    ct_expr <- expr_matrix[, markers, drop = FALSE]
    
    # Get combined weights for these markers
    ct_weights <- combined_weights[markers]
    
    # Detection mask: only count genes above threshold
    ct_detected <- ct_expr > min_expression
    
    # Combined score: expression * detection * combined_weight
    # Using expression directly (not squared) to avoid over-penalizing low expressers
    weighted_expr <- sweep(ct_expr * ct_detected, 2, ct_weights, "*")
    
    if (length(markers) == 1) {
      score_matrix[, ct] <- weighted_expr
    } else {
      # Use weighted mean (weighted by marker weights)
      score_matrix[, ct] <- rowSums(weighted_expr) / sum(ct_weights)
    }
  }
  
  return(score_matrix)
}


#' Compute AUCell-like rank-based scores
#' 
#' Uses a ranking approach similar to AUCell: for each cell, genes are ranked by 
#' expression, and the score is based on where marker genes fall in that ranking.
#' This method is robust to batch effects and technical variation.
#' 
#' Key improvement: Uses recovery curve area under curve (AUC) calculation
#' rather than simple proportion, giving partial credit to markers ranked
#' just outside the top threshold.
#' 
#' @keywords internal
.compute_aucell_scores <- function(expr_matrix, marker_sets, top_percent = 0.05) {
  n_cells <- nrow(expr_matrix)
  n_genes <- ncol(expr_matrix)
  n_celltypes <- length(marker_sets)
  
  # Number of genes in top X% - adaptive based on marker set sizes
  max_markers <- max(sapply(marker_sets, length))
  # Ensure n_top is at least 2x the largest marker set
  n_top <- max(max_markers * 2, round(n_genes * top_percent))
  n_top <- min(n_top, n_genes)  # Cap at total genes
  
  score_matrix <- matrix(0, nrow = n_cells, ncol = n_celltypes)
  colnames(score_matrix) <- names(marker_sets)
  
  # For each cell, rank genes by expression
  # Using a chunked approach for memory efficiency
  chunk_size <- 1000
  n_chunks <- ceiling(n_cells / chunk_size)
  
  for (chunk in seq_len(n_chunks)) {
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, n_cells)
    
    chunk_expr <- expr_matrix[start_idx:end_idx, , drop = FALSE]
    
    # Rank genes per cell (higher expression = lower rank number = better)
    chunk_ranks <- t(apply(chunk_expr, 1, function(x) rank(-x, ties.method = "average")))
    
    for (ct in names(marker_sets)) {
      markers <- marker_sets[[ct]]
      marker_cols <- which(colnames(expr_matrix) %in% markers)
      
      if (length(marker_cols) == 0) next
      
      n_markers <- length(marker_cols)
      marker_ranks <- chunk_ranks[, marker_cols, drop = FALSE]
      
      # Improved AUC calculation:
      # 1. Count markers in top n_top (binary contribution)
      # 2. Add weighted contribution based on rank (markers ranked higher get more credit)
      
      if (is.matrix(marker_ranks)) {
        # Binary: proportion of markers in top n_top
        in_top <- rowMeans(marker_ranks <= n_top)
        
        # Weighted: average of (1 - rank/n_genes) for markers in top n_top
        # This gives partial credit based on how highly ranked the marker is
        rank_scores <- 1 - marker_ranks / n_genes
        rank_scores[marker_ranks > n_top] <- 0  # Zero out markers not in top
        weighted_score <- rowMeans(rank_scores)
        
        # Combine: 70% binary, 30% rank-weighted
        score_matrix[start_idx:end_idx, ct] <- 0.7 * in_top + 0.3 * weighted_score
      } else {
        in_top <- as.numeric(marker_ranks <= n_top)
        rank_score <- ifelse(marker_ranks <= n_top, 1 - marker_ranks / n_genes, 0)
        score_matrix[start_idx:end_idx, ct] <- 0.7 * in_top + 0.3 * rank_score
      }
    }
  }
  
  return(score_matrix)
}


#' Apply UMAP-based spatial smoothing to scores
#' @keywords internal
.apply_umap_smoothing <- function(seurat_obj, score_matrix, umap_reduction, 
                                   k_neighbors, smoothing_weight, chunk_size, verbose) {
  
  # Get UMAP coordinates
  umap_coords <- Seurat::Embeddings(seurat_obj, reduction = umap_reduction)
  n_cells <- nrow(umap_coords)
  
  if (nrow(umap_coords) != nrow(score_matrix)) {
    stop("UMAP coordinates and score matrix have different numbers of cells!")
  }
  
  # Try to use RANN for fast approximate nearest neighbors if available
  use_rann <- requireNamespace("RANN", quietly = TRUE)
  
  if (use_rann) {
    if (verbose) message("  Using RANN for fast k-NN computation...")
    
    # RANN provides O(n log n) k-NN search
    nn_result <- RANN::nn2(umap_coords, umap_coords, k = k_neighbors + 1)
    
    # nn_result$nn.idx: matrix of neighbor indices (includes self as first column)
    # nn_result$nn.dists: matrix of distances
    
    # Remove self (first column)
    neighbor_idx <- nn_result$nn.idx[, -1, drop = FALSE]
    neighbor_dists <- nn_result$nn.dists[, -1, drop = FALSE]
    
    # Compute smoothed scores using vectorized operations
    smoothed_scores <- matrix(0, nrow = n_cells, ncol = ncol(score_matrix))
    colnames(smoothed_scores) <- colnames(score_matrix)
    
    # Distance weights (inverse distance)
    dist_weights <- 1 / (neighbor_dists + 1e-6)
    dist_weights <- dist_weights / rowSums(dist_weights)
    
    # For each cell type, compute weighted neighbor average
    for (ct_idx in seq_len(ncol(score_matrix))) {
      ct_scores <- score_matrix[, ct_idx]
      
      # Get neighbor scores for this cell type
      neighbor_ct_scores <- matrix(ct_scores[neighbor_idx], nrow = n_cells, ncol = k_neighbors)
      
      # Weighted average
      neighbor_avg <- rowSums(neighbor_ct_scores * dist_weights)
      
      # Blend with original
      smoothed_scores[, ct_idx] <- (1 - smoothing_weight) * ct_scores + 
        smoothing_weight * neighbor_avg
    }
    
  } else {
    # Fallback: chunked distance computation (slower but works without RANN)
    if (verbose) {
      message("  RANN not available, using chunked computation (install RANN for faster processing)...")
    }
    
    smoothed_scores <- score_matrix
    n_chunks <- ceiling(n_cells / chunk_size)
    
    for (chunk in seq_len(n_chunks)) {
      if (verbose && n_chunks > 1) {
        message(sprintf("  Processing chunk %d/%d...", chunk, n_chunks))
      }
      
      start_idx <- (chunk - 1) * chunk_size + 1
      end_idx <- min(chunk * chunk_size, n_cells)
      chunk_indices <- start_idx:end_idx
      chunk_size_actual <- length(chunk_indices)
      
      # Compute distance matrix for this chunk
      chunk_umap <- umap_coords[chunk_indices, , drop = FALSE]
      
      # Vectorized distance computation for chunk
      # d[i,j] = sqrt(sum((chunk_umap[i,] - umap_coords[j,])^2))
      dist_matrix <- as.matrix(dist(rbind(chunk_umap, umap_coords)))[1:chunk_size_actual, (chunk_size_actual + 1):(chunk_size_actual + n_cells)]
      
      for (i in seq_along(chunk_indices)) {
        cell_idx <- chunk_indices[i]
        dists <- dist_matrix[i, ]
        
        # Exclude self (distance = 0)
        dists[cell_idx] <- Inf
        
        # Find k nearest neighbors
        neighbors <- order(dists)[1:k_neighbors]
        
        # Get neighbor scores
        neighbor_scores <- score_matrix[neighbors, , drop = FALSE]
        
        # Distance-weighted average
        neighbor_dists <- dists[neighbors]
        dist_weights <- 1 / (neighbor_dists + 1e-6)
        dist_weights <- dist_weights / sum(dist_weights)
        
        neighbor_avg <- colSums(sweep(neighbor_scores, 1, dist_weights, "*"))
        
        # Blend
        smoothed_scores[cell_idx, ] <- (1 - smoothing_weight) * score_matrix[cell_idx, ] +
          smoothing_weight * neighbor_avg
      }
    }
  }
  
  # Re-normalize after smoothing
  row_sums <- rowSums(smoothed_scores)
  row_sums[row_sums == 0] <- 1
  smoothed_scores <- smoothed_scores / row_sums
  
  return(smoothed_scores)
}
