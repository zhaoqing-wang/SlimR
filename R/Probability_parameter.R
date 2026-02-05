#' Adaptive Parameter Tuning for Single-Cell Data Annotation in SlimR
#'
#' This function automatically determines optimal min_expression, specificity_weight,
#' and threshold parameters for single-cell data analysis based on dataset characteristics
#' using adaptive algorithms derived from empirical analysis of single-cell datasets.
#'
#' @param seurat_obj A Seurat object containing single-cell data
#' @param features Character vector of feature names (genes) to analyze. If NULL,
#'        will use highly variable features from the Seurat object.
#' @param assay Name of assay to use (default: default assay)
#' @param cluster_col Column name in metadata containing cluster information
#' @param n_celltypes Expected number of cell types in marker database (default: 50).
#'        Used for threshold recommendation calculation.
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item min_expression: Recommended expression threshold
#'   \item specificity_weight: Recommended specificity weight
#'   \item threshold: Recommended probability threshold for candidate selection
#'   \item dataset_features: Extracted dataset characteristics
#'   \item parameter_rationale: Explanation of parameter choices
#' }
#'
#' @export
#' @family Section_3_Automated_Annotation
#' 
#' @importFrom stats dist median sd aggregate quantile pnorm
#'
#' @examples
#' \dontrun{
#' SlimR_params <- Parameter_Calculate(
#'   seurat_obj = sce,
#'   features = c("CD3E", "CD4", "CD8A"),
#'   assay = "RNA",
#'   cluster_col = "seurat_clusters",
#'   n_celltypes = 98,
#'   verbose = TRUE
#'   )
#' }
#'
Parameter_Calculate <- function(
    seurat_obj,
    features = NULL,
    assay = NULL,
    cluster_col = NULL,
    n_celltypes = 50,
    verbose = TRUE
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input object must be a Seurat object")
  }
  
  assay <- if (is.null(assay)) Seurat::DefaultAssay(seurat_obj) else assay
  Seurat::DefaultAssay(seurat_obj) <- assay
  
  if (is.null(features) || length(features) == 0) {
    if (verbose) message("SlimR parameter calculate: No features provided, using variable features.")
    features <- Seurat::VariableFeatures(seurat_obj)
    if (length(features) == 0) {
      features <- head(rownames(seurat_obj[[assay]]), 2000)
    }
    features <- head(features, 500)
  }
  
  valid_features <- features[features %in% rownames(seurat_obj[[assay]])]
  if (length(valid_features) < 3) {
    warning("Fewer than 3 valid features found, results may be unreliable")
  }
  
  if (verbose) message("SlimR parameter calculate: Extracting dataset features from ", length(valid_features), " genes.")
  
  dataset_features <- extract_dataset_features(seurat_obj, valid_features, assay, cluster_col)
  
  if (verbose) message("SlimR parameter calculate: Computing adaptive parameters.")
  
  optimal_params <- compute_adaptive_parameters(dataset_features, n_celltypes)
  
  if (verbose) {
    message("SlimR parameter calculate: Parameter recommendation: ")
    message("  min_expression: ", optimal_params$min_expression)
    message("  specificity_weight: ", round(optimal_params$specificity_weight, 2))
    message("  threshold: ", optimal_params$threshold)
    message("  Rationale: ", optimal_params$rationale)
  }
  
  result <- list(
    min_expression = optimal_params$min_expression,
    specificity_weight = optimal_params$specificity_weight,
    threshold = optimal_params$threshold,
    dataset_features = dataset_features,
    parameter_rationale = optimal_params$rationale
  )
  
  return(result)
}

#' Extract Dataset Characteristics for Adaptive Parameter Calculation (Use in package)
#'
#' Computes various statistical features from single-cell data that are used
#' as input for the parameter prediction model.
#'
#' @param seurat_obj Seurat object
#' @param features Features to analyze
#' @param assay Assay name
#' @param cluster_col Cluster column name
#'
#' @return List of dataset characteristics including expression statistics,
#'         variability measures, and cluster properties
#'
#' @family Section_1_Functions_Use_in_Package
#' 
#' @importFrom stats dist median sd aggregate quantile var
#' 
extract_dataset_features <- function(seurat_obj, features, assay = NULL, cluster_col = NULL) {
  assay <- if (is.null(assay)) Seurat::DefaultAssay(seurat_obj) else assay
  Seurat::DefaultAssay(seurat_obj) <- assay
  
  cells <- unlist(Seurat::CellsByIdentities(object = seurat_obj))
  data.features <- Seurat::FetchData(object = seurat_obj, vars = features, cells = cells)
  
  # Assign cluster identities
  if (!is.null(cluster_col)) {
    data.features$id <- seurat_obj@meta.data[cells, cluster_col, drop = TRUE]
  } else {
    data.features$id <- Seurat::Idents(object = seurat_obj)[cells]
  }
  
  features <- setdiff(features, "id")
  expr_matrix <- as.matrix(data.features[, features])
  
  # Compute gene expression statistics
  expression_stats <- sapply(features, function(gene) {
    expr <- data.features[[gene]]
    c(
      mean_expr = mean(expr),
      sd_expr = stats::sd(expr),
      zero_frac = mean(expr == 0),
      median_expr = stats::median(expr),
      cv_expr = stats::sd(expr) / (mean(expr) + 1e-6)
    )
  })
  
  # Compute expression quantiles for non-zero values
  nonzero_expr <- expr_matrix[expr_matrix > 0]
  expr_quantiles <- stats::quantile(nonzero_expr, 
                                     probs = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.75, 0.9), 
                                     na.rm = TRUE)
  
  # Compute cluster-level statistics for threshold estimation
  cluster_ids <- unique(data.features$id)
  n_clusters <- length(cluster_ids)
  
  # Calculate within-cluster and between-cluster variance
  cluster_means <- stats::aggregate(expr_matrix, by = list(cluster = data.features$id), mean)
  cluster_matrix <- as.matrix(cluster_means[, -1])
  
  # Between-cluster variance (signal)
  global_mean <- colMeans(expr_matrix)
  between_var <- mean(apply(cluster_matrix, 2, function(x) stats::var(x)))
  

  # Within-cluster variance (noise) 
  within_var <- mean(sapply(cluster_ids, function(cid) {
    cluster_expr <- expr_matrix[data.features$id == cid, , drop = FALSE]
    mean(apply(cluster_expr, 2, stats::var))
  }))
  
  # Signal-to-noise ratio
  snr <- if (!is.na(within_var) && within_var > 0) between_var / within_var else between_var
  
  # Calculate expression dynamic range
  dynamic_range <- log10(max(nonzero_expr) / (min(nonzero_expr) + 1e-10) + 1)
  
  # Gene detection rate per cluster
  detection_rates <- sapply(cluster_ids, function(cid) {
    cluster_expr <- expr_matrix[data.features$id == cid, , drop = FALSE]
    mean(colMeans(cluster_expr > 0))
  })
  
  dataset_features <- list(
    n_genes = length(features),
    n_cells = nrow(data.features),
    n_clusters = n_clusters,
    
    # Expression statistics
    global_mean_expression = mean(expr_matrix),
    global_zero_fraction = mean(expr_matrix == 0),
    expression_sparsity = mean(apply(expr_matrix, 2, function(x) mean(x == 0))),
    
    # Gene variability
    mean_gene_cv = mean(expression_stats["cv_expr", ], na.rm = TRUE),
    sd_gene_cv = stats::sd(expression_stats["cv_expr", ], na.rm = TRUE),
    median_gene_cv = stats::median(expression_stats["cv_expr", ], na.rm = TRUE),
    
    # Cluster metrics
    cluster_variability = calculate_cluster_variability(data.features, features),
    between_cluster_var = between_var,
    within_cluster_var = within_var,
    signal_to_noise = snr,
    
    # Distribution characteristics
    expression_skewness = calculate_expression_skewness(expr_matrix),
    dynamic_range = dynamic_range,
    batch_effect_score = estimate_batch_effect(seurat_obj, assay),
    
    # Detection rates
    mean_detection_rate = mean(detection_rates),
    sd_detection_rate = stats::sd(detection_rates),
    
    # Detailed quantiles
    expression_quantiles = expr_quantiles,
    
    # Per-gene statistics summary
    gene_mean_distribution = stats::quantile(expression_stats["mean_expr", ], 
                                              probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE),
    gene_zero_distribution = stats::quantile(expression_stats["zero_frac", ], 
                                              probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
  )
  
  return(dataset_features)
}

#' Calculate Cluster Variability (Use in package)
#'
#' Measures the degree of separation between different cell clusters
#' based on expression patterns.
#'
#' @param data.features Data frame containing expression data and cluster labels
#' @param features Feature names to include in analysis
#'
#' @return Numeric value representing cluster separation strength
#'
#' @family Section_1_Functions_Use_in_Package
#' 
#' @importFrom stats dist aggregate
#' 
calculate_cluster_variability <- function(data.features, features) {
  cluster_means <- stats::aggregate(data.features[, features], 
                            by = list(cluster = data.features$id), 
                            mean)
  
  cluster_matrix <- as.matrix(cluster_means[, -1])
  
  if (nrow(cluster_matrix) > 1) {
    # Compute mean distance between cluster centroids
    dist_matrix <- stats::dist(cluster_matrix)
    variability <- mean(as.matrix(dist_matrix))
  } else {
    variability <- 0  # Only one cluster
  }
  
  return(variability)
}

#' Calculate Expression Distribution Skewness (Use in package)
#'
#' Computes the average skewness of gene expression distributions
#' across all features.
#'
#' @param expression_matrix Matrix of expression values
#'
#' @return Mean absolute skewness across all genes
#'
#' @family Section_1_Functions_Use_in_Package
#' 
calculate_expression_skewness <- function(expression_matrix) {
  skew_vals <- apply(expression_matrix, 2, function(x) {
    sd_val <- stats::sd(x, na.rm = TRUE)
    if (is.na(sd_val) || sd_val == 0) return(0)
    m <- mean(x, na.rm = TRUE)
    mean((x - m)^3, na.rm = TRUE) / (sd_val^3)  # Fisher-Pearson coefficient of skewness
  })
  return(mean(abs(skew_vals), na.rm = TRUE))
}

#' Estimate Batch Effect Strength (Use in package)
#'
#' Roughly estimates the potential impact of batch effects
#' using available metadata.
#'
#' @param seurat_obj Seurat object
#' @param assay Assay name
#'
#' @return Batch effect score (0 indicates no detectable batch effect)
#'
#' @family Section_1_Functions_Use_in_Package
#' 
estimate_batch_effect <- function(seurat_obj, assay) {
  if ("batch" %in% colnames(seurat_obj@meta.data)) {
    batch_groups <- unique(seurat_obj@meta.data$batch)
    if (length(batch_groups) > 1) {
      return(length(batch_groups) * 0.1)
    }
  }
  return(0)
}

#' Compute Adaptive Parameters Based on Dataset Features (Use in package)
#'
#' Calculates optimal min_expression, specificity_weight, and threshold parameters
#' using continuous adaptive algorithms based on dataset characteristics.
#'
#' @param dataset_features List of dataset characteristics from extract_dataset_features()
#' @param n_celltypes Expected number of cell types in marker database
#'
#' @return List containing min_expression, specificity_weight, threshold, and rationale
#'
#' @family Section_1_Functions_Use_in_Package
#'
#' @importFrom stats quantile pnorm
#'
compute_adaptive_parameters <- function(dataset_features, n_celltypes = 50) {
  
  # Extract features
  sparsity <- dataset_features$global_zero_fraction
  mean_expr <- dataset_features$global_mean_expression
  cv <- dataset_features$mean_gene_cv
  median_cv <- dataset_features$median_gene_cv
  cluster_var <- dataset_features$cluster_variability
  n_clusters <- dataset_features$n_clusters
  skewness <- dataset_features$expression_skewness
  snr <- dataset_features$signal_to_noise
  dynamic_range <- dataset_features$dynamic_range
  detection_rate <- dataset_features$mean_detection_rate
  expr_quantiles <- dataset_features$expression_quantiles
  gene_mean_dist <- dataset_features$gene_mean_distribution
  
  rationale_parts <- character(0)
  
  # ============================================================================

  # PART 1: min_expression - Continuous Adaptive Calculation

  # ============================================================================
  
  # Base value from expression quantiles (use 10th percentile of non-zero expression)
  if (!is.null(expr_quantiles) && length(expr_quantiles) >= 3) {
    # Use interpolation between 5th and 15th percentile based on sparsity
    q05 <- expr_quantiles["5%"]
    q10 <- expr_quantiles["10%"]
    q15 <- expr_quantiles["15%"]
    q20 <- expr_quantiles["20%"]
    
    # Higher sparsity -> use lower quantile; lower sparsity -> use higher quantile
    sparsity_weight <- (sparsity - 0.5) / 0.5  # Normalize to [-1, 1] range
    sparsity_weight <- max(-1, min(1, sparsity_weight), na.rm = TRUE)
    
    if (is.na(sparsity_weight) || sparsity_weight > 0) {
      # High sparsity: interpolate between q05 and q10
      base_min_expr <- q05 + (q10 - q05) * (1 - ifelse(is.na(sparsity_weight), 0, sparsity_weight))
    } else {
      # Low sparsity: interpolate between q10 and q20
      base_min_expr <- q10 + (q20 - q10) * (-sparsity_weight)
    }
  } else {
    # Fallback: use mean expression scaled by sparsity
    base_min_expr <- mean_expr * (0.3 - 0.2 * sparsity)
  }
  
  # Ensure positive base value
  base_min_expr <- max(0.01, base_min_expr)
  
  # Adjustment 1: Skewness correction (continuous)
  # High skewness indicates long-tail distribution, need lower threshold
  skewness_factor <- 1 - 0.05 * log1p(skewness)  # Logarithmic dampening
  skewness_factor <- max(0.7, min(1.2, skewness_factor))
  
  # Adjustment 2: Dynamic range correction
  # Large dynamic range suggests need for relative rather than absolute threshold
  if (!is.null(dynamic_range) && !is.na(dynamic_range) && dynamic_range > 3) {
    range_factor <- 1 - 0.08 * (dynamic_range - 3)
    range_factor <- max(0.6, range_factor)
  } else {
    range_factor <- 1
  }
  
  # Adjustment 3: Detection rate consideration
  # Low detection rate across clusters needs lower threshold
  if (!is.null(detection_rate)) {
    detection_factor <- 0.7 + 0.6 * detection_rate  # Range: [0.7, 1.3]
  } else {
    detection_factor <- 1
  }
  
  # Adjustment 4: Gene-level mean distribution
  # If most genes have low mean expression, lower the threshold
  if (!is.null(gene_mean_dist)) {
    median_gene_mean <- gene_mean_dist["50%"]
    gene_factor <- sqrt(median_gene_mean / (mean_expr + 1e-6))
    gene_factor <- max(0.5, min(1.5, gene_factor))
  } else {
    gene_factor <- 1
  }
  
  # Combine all factors
  min_expression <- base_min_expr * skewness_factor * range_factor * detection_factor * gene_factor
  
  # Final bounds with finer granularity
  min_expression <- max(0.01, min(1.0, min_expression))
  
  # Round to 3 significant figures for cleaner output
  min_expression <- signif(min_expression, 3)
  
  rationale_parts <- c(rationale_parts, 
                       sprintf("min_expr base=%.3f (Q%.0f%%)", base_min_expr, 
                               ifelse(sparsity > 0.5, 5 + 5*(1-sparsity_weight), 10 - 10*sparsity_weight)))
  
  # ============================================================================
  # PART 2: specificity_weight - Continuous Adaptive Calculation
  # ============================================================================
  
  # Base weight from signal-to-noise ratio
  if (!is.null(snr) && snr > 0) {
    # Higher SNR -> lower weight needed (clusters already well-separated)
    # Use inverse relationship with logarithmic scaling
    base_weight <- 6 - 2 * log1p(snr)
    base_weight <- max(1, min(6, base_weight))
  } else {
    # Fallback to cluster variability
    base_weight <- 6 - log1p(cluster_var)
    base_weight <- max(1, min(6, base_weight))
  }
  
  # Adjustment 1: CV-based fine-tuning (continuous)
  # Higher CV means genes are more variable -> can rely more on expression patterns
  cv_factor <- 1.3 - 0.15 * cv
  cv_factor <- max(0.7, min(1.3, cv_factor))
  
  # Adjustment 2: Number of clusters scaling
  # More clusters -> slightly higher weight to improve discrimination
  cluster_factor <- 1 + 0.01 * (n_clusters - 10)
  cluster_factor <- max(0.85, min(1.2, cluster_factor))
  
  # Adjustment 3: Sparsity interaction
  # Very high sparsity datasets need different weighting strategy
  if (sparsity > 0.9) {
    sparsity_adj <- 1 + 0.5 * (sparsity - 0.9) / 0.1
  } else if (sparsity < 0.5) {
    sparsity_adj <- 1 - 0.2 * (0.5 - sparsity) / 0.5
  } else {
    sparsity_adj <- 1
  }
  
  # Combine factors
  specificity_weight <- base_weight * cv_factor * cluster_factor * sparsity_adj
  
  # Final bounds
  specificity_weight <- max(0.5, min(8, specificity_weight))
  
  # Round to 2 decimal places
  specificity_weight <- round(specificity_weight, 2)
  
  rationale_parts <- c(rationale_parts,
                       sprintf("weight base=%.1f (SNR=%.2f)", base_weight, 
                               ifelse(!is.null(snr), snr, cluster_var)))
  
  # ============================================================================
  # PART 3: threshold - Adaptive Calculation for Candidate Selection
  # ============================================================================
  
  # Goal: Select threshold that balances:

  # - Computational efficiency (fewer candidates = faster)
  # - Annotation accuracy (more candidates = better coverage)
  # Range: [0.55, 0.85], where:
  #   - Lower threshold (0.55-0.65): More candidates, slower but thorough
  #   - Medium threshold (0.65-0.75): Balanced
  #   - Higher threshold (0.75-0.85): Fewer candidates, faster but may miss
  
  # Factor 1: Cluster quality (SNR-based)
  # Poor separation -> need more candidates to ensure coverage
  if (!is.null(snr) && snr > 0) {
    # SNR typically ranges 0.1-5; transform to contribution
    snr_contrib <- 0.1 * tanh(snr)  # Saturates around ±0.1
  } else {
    snr_contrib <- -0.05  # Default: slightly lower threshold
  }
  
  # Factor 2: Number of cell types in database
  # More cell types -> more potential confusion -> lower threshold to verify
  celltype_contrib <- -0.1 * log10(n_celltypes / 50)  # 50 as baseline
  celltype_contrib <- max(-0.1, min(0.1, celltype_contrib))
  
  # Factor 3: Number of clusters
  # Many clusters -> heterogeneous sample -> may need lower threshold
  cluster_contrib <- -0.05 * log(n_clusters / 20 + 1)
  cluster_contrib <- max(-0.08, min(0.05, cluster_contrib))
  
  # Factor 4: Expression characteristics
  # High sparsity makes probability estimates less reliable -> lower threshold
  sparsity_contrib <- -0.1 * (sparsity - 0.7)  # 0.7 as typical value
  sparsity_contrib <- max(-0.08, min(0.08, sparsity_contrib))
  
  # Factor 5: Gene variability
  # High CV -> markers more distinctive -> can use higher threshold
  cv_contrib <- 0.03 * (cv - 1.5)
  cv_contrib <- max(-0.05, min(0.08, cv_contrib))
  
  # Combine: start from 0.70 as neutral baseline
  threshold <- 0.70 + snr_contrib + celltype_contrib + cluster_contrib + 
               sparsity_contrib + cv_contrib
  
  # Bound to [0.55, 0.85]
  threshold <- max(0.55, min(0.85, threshold))
  
  # Round to 2 decimal places
  threshold <- round(threshold, 2)
  
  # Estimate expected candidates for reporting
  # Assuming roughly uniform probability distribution, top (1-threshold) fraction are candidates
  est_candidates <- round(n_celltypes * (1 - threshold))
  est_candidates <- max(1, est_candidates)
  
  rationale_parts <- c(rationale_parts,
                       sprintf("threshold=%.2f (~%d candidates/cluster)", threshold, est_candidates))
  
  # ============================================================================
  # Compile rationale
  # ============================================================================
  
  rationale <- paste(rationale_parts, collapse = "; ")
  
  return(list(
    min_expression = min_expression,
    specificity_weight = specificity_weight,
    threshold = threshold,
    rationale = rationale
  ))
}
