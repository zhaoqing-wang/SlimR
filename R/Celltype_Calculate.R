#' Uses "marker_list" to calculate probability, prediction results, AUC and generate heatmap for cell annotation
#'
#' @param seurat_obj Enter the Seurat object with annotation columns such as
#'     "seurat_cluster" in meta.data to be annotated.
#' @param gene_list A list of cells and corresponding gene controls, the name of
#'     the list is cell type, and the first column of the list corresponds to markers.
#'     Lists can be generated using functions such as "Markers_filter_Cellmarker2 ()",
#'     "Markers_filter_PanglaoDB ()", "read_excel_markers ()", "read_seurat_markers ()", etc.
#' @param species This parameter selects the species "Human" or "Mouse" for standard
#'     gene format correction of markers entered by "Marker_list".
#' @param cluster_col Enter annotation columns such as "seurat_cluster" in meta.data
#'     of the Seurat object to be annotated. Default parameters use "cluster_col =
#'     'seurat_clusters'".
#' @param assay Enter the assay used by the Seurat object, such as "RNA". Default
#'     parameters use "assay = 'RNA'".
#' @param min_expression The min_expression parameter defines a threshold value to
#'     determine whether a cell's expression of a feature is considered "expressed"
#'     or not. It is used to filter out low-expression cells that may contribute
#'     noise to the analysis. Default parameters use "min_expression = 0.1".
#' @param specificity_weight The specificity_weight parameter controls how much the
#'     expression variability (standard deviation) of a feature within a cluster
#'     contributes to its "specificity score." It amplifies or suppresses the impact
#'     of variability in the final score calculation.Default parameters use
#'     "specificity_weight = 3".
#' @param threshold This parameter refers to the normalized similarity between the
#'     "alternative cell type" and the "predicted cell type" in the returned results.
#'     (the default parameter is 0.6)
#' @param compute_AUC Logical indicating whether to calculate AUC values for predicted
#'     cell types. AUC measures how well the marker genes distinguish the cluster from
#'     others. When TRUE, adds an AUC column to the prediction results. (default: TRUE)
#' @param plot_AUC The logic indicates whether to draw an AUC curve for the predicted cell
#'     type. When TRUE, add an AUC_plot to result. (default: TRUE)
#' @param AUC_correction Logical value controlling AUC-based correction. (default = FALSE)
#'     When set to TRUE:
#'     1.Computes AUC values for candidate cell types. (probability > threshold)
#'     2.Selects the cell type with the highest AUC as the final predicted type.
#'     3.Records the selected type's AUC value in the "AUC" column.
#' @param colour_low Color for lowest probability level in Heatmap visualization of
#'     probability matrix. (default = "navy")
#' @param colour_high Color for highest probability level Heatmap visualization of
#'     probability matrix. (default = "firebrick3")
#'
#' @returns A list containing:
#' \itemize{
#'   \item Expression_list: List of expression matrices for each cell type
#'   \item Proportion_list: List of proportion of expression for each cell type
#'   \item Expression_scores_matrix: Matrix of expression scores
#'   \item Probability_matrix: Matrix of normalized probabilities
#'   \item Prediction_results: Data frame with cluster annotations including:
#'     \itemize{
#'       \item cluster_col: Cluster identifier
#'       \item Predicted_cell_type: Primary predicted cell type
#'       \item AUC: Area Under the Curve value (when compute_AUC = TRUE)
#'       \item Alternative_cell_types: Semi-colon separated alternative cell types
#'     }
#'   \item Heatmap_plot: Heatmap visualization of probability matrix
#'   \item AUC_plot: AUC visualization of Predicted cell type
#'   \item AUC_list: The resulting list of AUC values calculated for genes in alternative cell types above the approximate threshold
#' }
#'
#' @export
#' @family Section_3_Automated_Annotation
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom utils tail
#' @importFrom stats runif
#' @importFrom ggplot2 ggplot aes geom_line geom_abline scale_color_manual
#' @importFrom ggplot2 theme_minimal labs theme element_text element_blank
#' @importFrom ggplot2 guide_legend guides scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 element_line expand_limits
#' @importFrom ggplot2 .data
#'
#' @examples
#' \dontrun{
#' SlimR_anno_result <- Celltype_Calculate(seurat_obj = sce,
#'     gene_list = Markers_list,
#'     species = "Human",
#'     cluster_col = "seurat_clusters",
#'     assay = "RNA",
#'     min_expression = 0.1,
#'     specificity_weight = 3,
#'     threshold = 0.6,
#'     compute_AUC = TRUE,
#'     plot_AUC = TRUE,
#'     AUC_correction = FALSE,
#'     colour_low = "navy",
#'     colour_high = "firebrick3"
#'     )
#'     }
#'
Celltype_Calculate <- function(
    seurat_obj,
    gene_list,
    species,
    cluster_col = "seurat_clusters",
    assay = "RNA",
    min_expression = 0.1,
    specificity_weight = 3,
    threshold = 0.6,
    compute_AUC = TRUE,
    plot_AUC = TRUE,
    AUC_correction = FALSE,
    colour_low = "navy",
    colour_high = "firebrick3"
) {
  required_packages <- c("ggplot2", "patchwork", "dplyr", "scales", "tidyr", "gridExtra", "gtable", "grid", "pheatmap")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Please install the required package: %s", pkg))
    }
    library(pkg, character.only = TRUE)
  }

  if (plot_AUC) compute_AUC <- TRUE
  if (AUC_correction) compute_AUC <- TRUE

  if (!inherits(seurat_obj, "Seurat")) stop("Input object must be a Seurat object!")
  if (!is.list(gene_list)) stop("Gene list must be a list of data.frames!")
  if (species != "Human" && species != "Mouse") stop("species must be 'Human' or 'Mouse'")

  assay <- if (is.null(assay)) DefaultAssay(seurat_obj) else assay

  cluster_scores_list <- list()
  cluster_mean_list <- list()
  cluster_frac_list <- list()
  valid_genes_list <- list()

  cell_types <- names(gene_list)
  total <- length(cell_types)
  cycles <- 0

  message(paste0("SlimR calculate: The input 'Markers_list' has ",total," cell types to be calculated."))

  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    message(paste0("\n","[", i, "/", total, "] Processing cell type: ", cell_type))

    current_df <- gene_list[[cell_type]]

    if (ncol(current_df) < 1) {
      warning(paste("Skipping", cell_type, ": Requires at least a gene column"))
      next
    }

    genes <- current_df[[1]]
    genes_processed <- if (species == "Human") {
      toupper(genes)
    } else {
      paste0(toupper(substr(genes, 1, 1)), tolower(substr(genes, 2, nchar(genes))))
    }

    valid_idx <- genes_processed %in% rownames(seurat_obj[[assay]])
    if (sum(valid_idx) == 0) {
      warning(paste("No valid genes for", cell_type))
      next
    }

    valid_data <- data.frame(
      original = genes[valid_idx],
      processed = genes_processed[valid_idx],
      stringsAsFactors = FALSE
    ) %>% distinct(processed, .keep_all = TRUE)

    gene_order_processed <- valid_data$processed
    gene_order_original <- valid_data$original

    prob_expression <- calculate_probability(object = seurat_obj,
                                             cluster_col = cluster_col,
                                             assay = assay,
                                             features = gene_order_processed,
                                             min_expression = min_expression,
                                             specificity_weight = specificity_weight)
    cluster_scores_list[[cell_type]] <- prob_expression$cluster_scores
    cluster_mean_list[[cell_type]] <- prob_expression$cluster_expr
    cluster_frac_list[[cell_type]] <- prob_expression$cluster_frac

    valid_genes_list[[cell_type]] <- unique(colnames(prob_expression$cluster_expr))

    message(paste0("[", i, "/", total, "] ",cell_type," characteristic genes expression calculated."))
    cycles <- cycles + 1
  }
  message(paste0("\n","SlimR calculate: Out of the ",total," cell types in 'Markers_list', ",cycles," cell types have been calculated. You can see the reason for not calculating cell types by 'warnings()'."))

  expr_list <- cluster_mean_list
  frac_list <- cluster_frac_list
  scores_matrix <- do.call(rbind, cluster_scores_list)

  normalize_row <- function(x) {
  x[is.na(x)] <- 0
  
  if (diff(range(x)) == 0) {
    return(rep(0, length(x)))
  }
  
  (x - min(x)) / (max(x) - min(x))
  }

  normalize_matrix <- apply(scores_matrix, 2, normalize_row)
  result_matrix <- t(normalize_matrix)

  # Check if all values are the same (would cause pheatmap to fail)
  if (length(unique(as.vector(result_matrix))) == 1) {
    # Add small random noise to prevent identical values
    result_matrix <- result_matrix + matrix(runif(length(result_matrix), -1e-10, 1e-10), 
                                             nrow = nrow(result_matrix))
  }

  p <- pheatmap::pheatmap(result_matrix,
                          main = "Cell annotation heatmap | SlimR",
                          color = colorRampPalette(c(colour_low, "white", colour_high))(50),
                          fontsize = 12,
                          cluster_rows = T,
                          cluster_cols = T,
                          legend_breaks = c(0,1),
                          legend_labels = c("Low probability","High probability"))

  generate_prediction_table <- function(df, threshold = threshold) {
    clusters <- rownames(df)
    predicted_cell_types <- vector("character", length = length(clusters))
    alternative_cell_types <- vector("character", length = length(clusters))
    candidate_types_list <- list()

    for (i in seq_along(clusters)) {
      cluster <- clusters[i]
      row_values <- as.numeric(unlist(df[i, ]))
      cell_types <- names(df[i, ])
      max_index <- which.max(row_values)
      predicted <- cell_types[max_index]
      candidate_types <- cell_types[row_values > threshold]
      candidate_types_list[[cluster]] <- candidate_types

      alt <- candidate_types[candidate_types != predicted]
      alt_str <- if (length(alt) > 0) paste(alt, collapse = "; ") else NA_character_
      predicted_cell_types[i] <- predicted
      alternative_cell_types[i] <- alt_str
    }
    result_df <- data.frame(
      cluster_col = clusters,
      Predicted_cell_type = predicted_cell_types,
      Alternative_cell_types = alternative_cell_types,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    attr(result_df, "candidate_types") <- candidate_types_list
    return(result_df)
  }

  scores_matrix <- as.data.frame(t(scores_matrix))
  probability_matrix <- as.data.frame(result_matrix)
  prediction_results <- generate_prediction_table(probability_matrix, threshold = threshold)
  candidate_types_list <- attr(prediction_results, "candidate_types")

  fastAUC <- function(predictions, labels) {
    ord <- order(predictions, decreasing = TRUE)
    labels <- labels[ord]
    predictions <- predictions[ord]

    tpr <- cumsum(labels) / sum(labels)
    fpr <- cumsum(!labels) / sum(!labels)

    tpr <- c(0, tpr, 1)
    fpr <- c(0, fpr, 1)

    auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
    return(auc)
  }

  compute_max_gene_auc <- function(expr_data, labels) {
    if (ncol(expr_data) == 0) {
      return(list(max_auc = NA, individual_aucs = numeric(0), best_gene_expr = NULL, 
                  gene_names = character(0), sorted_indices = integer(0)))
    }
    
    n_genes <- ncol(expr_data)
    individual_aucs <- numeric(n_genes)
    gene_names <- colnames(expr_data)
    
    for (i in seq_len(n_genes)) {
      gene_expr <- expr_data[, i]
      
      gene_sd <- sd(gene_expr, na.rm = TRUE)
      if (is.na(gene_sd) || gene_sd < 1e-6) {
        individual_aucs[i] <- NA
        next
      }
      
      if (length(unique(labels)) >= 2) {
        individual_aucs[i] <- fastAUC(gene_expr, labels)
      } else {
        individual_aucs[i] <- NA
      }
    }
    
    valid_idx <- !is.na(individual_aucs)
    
    if (sum(valid_idx) == 0) {
      return(list(max_auc = NA, individual_aucs = individual_aucs, best_gene_expr = NULL,
                  gene_names = gene_names, sorted_indices = integer(0)))
    }
    
    sorted_indices <- order(individual_aucs, decreasing = TRUE, na.last = TRUE)
    max_auc <- individual_aucs[sorted_indices[1]]
    best_gene_expr <- expr_data[, sorted_indices[1]]
    
    return(list(
      max_auc = max_auc,
      individual_aucs = individual_aucs,
      best_gene_expr = best_gene_expr,
      gene_names = gene_names,
      sorted_indices = sorted_indices
    ))
  }

  compute_roc_data <- function(predictions, labels) {
    ord <- order(predictions, decreasing = TRUE)
    labels <- labels[ord]
    predictions <- predictions[ord]

    n_pos <- sum(labels)
    n_neg <- sum(!labels)
    tpr <- cumsum(labels) / n_pos
    fpr <- cumsum(!labels) / n_neg

    tpr <- c(0, tpr, 1)
    fpr <- c(0, fpr, 1)

    data.frame(fpr = fpr, tpr = tpr)
  }

  if (compute_AUC) {

    Seurat::DefaultAssay(seurat_obj) <- assay
    
    auc_list_storage <- list()

    if (AUC_correction) {
      message(paste0("\n","SlimR AUC correction: Performing AUC correction for all candidate cell types (threshold > ",threshold,")."))
      new_predicted <- character(nrow(prediction_results))
      new_aucs <- numeric(nrow(prediction_results))
      new_alt_list <- character(nrow(prediction_results))

      for (i in seq_len(nrow(prediction_results))) {
        cluster_id <- prediction_results$cluster_col[i]
        candidate_types <- candidate_types_list[[cluster_id]]

        if (length(candidate_types) == 0) {
          new_predicted[i] <- NA
          new_aucs[i] <- NA
          new_alt_list[i] <- NA
          next
        }

        auc_results_list <- list()
        
        for (j in seq_along(candidate_types)) {
          cell_type <- candidate_types[j]
          features <- valid_genes_list[[cell_type]]

          all_cells <- colnames(seurat_obj)
          expr_data <- FetchData(seurat_obj, vars = features, cells = all_cells)
          
          labels <- seurat_obj@meta.data[all_cells, cluster_col] == cluster_id
          if (length(unique(labels)) < 2) {
            warning(paste("Skipping AUC for cluster", cluster_id, "and cell type", cell_type, ": Only one class present"))
            next
          } else {
            auc_results_list[[cell_type]] <- compute_max_gene_auc(expr_data, labels)
            
            if (!cell_type %in% names(auc_list_storage)) {
              auc_list_storage[[cell_type]] <- matrix(NA, 
                                                       nrow = length(unique(seurat_obj@meta.data[[cluster_col]])),
                                                       ncol = length(features))
              rownames(auc_list_storage[[cell_type]]) <- unique(seurat_obj@meta.data[[cluster_col]])
              colnames(auc_list_storage[[cell_type]]) <- features
            }
            
            auc_list_storage[[cell_type]][cluster_id, ] <- auc_results_list[[cell_type]]$individual_aucs
          }
        }

        if (length(auc_results_list) == 0) {
          new_predicted[i] <- candidate_types[1]
          new_aucs[i] <- NA
          new_alt_list[i] <- NA
          next
        }
        
        ranked_types <- names(auc_results_list)
        type_scores <- matrix(NA, nrow = length(ranked_types), 
                             ncol = max(sapply(auc_results_list, function(x) length(x$individual_aucs))))
        rownames(type_scores) <- ranked_types
        
        for (ct in ranked_types) {
          sorted_aucs <- auc_results_list[[ct]]$individual_aucs[auc_results_list[[ct]]$sorted_indices]
          type_scores[ct, 1:length(sorted_aucs)] <- sorted_aucs
        }
        
        final_order <- ranked_types
        for (col_idx in 1:ncol(type_scores)) {
          col_vals <- type_scores[final_order, col_idx]
          if (all(is.na(col_vals))) break
          final_order <- final_order[order(col_vals, decreasing = TRUE, na.last = TRUE)]
          if (length(unique(col_vals[!is.na(col_vals)])) == 1) next
        }
        
        best_type <- final_order[1]
        best_auc <- auc_results_list[[best_type]]$max_auc
        
        new_predicted[i] <- best_type
        new_aucs[i] <- best_auc
        
        alt_types <- final_order[-1]
        alt_strs <- character(0)
        for (k in seq_along(alt_types)) {
          alt_auc <- auc_results_list[[alt_types[k]]]$max_auc
          alt_strs[k] <- paste0(alt_types[k], " (", round(alt_auc, digits = 7), ")")
        }
        new_alt_list[i] <- paste(alt_strs, collapse = " ; ")
      }

      prediction_results$Predicted_cell_type <- new_predicted
      prediction_results$AUC <- new_aucs
      prediction_results$Alternative_cell_types <- new_alt_list

      message(paste0("SlimR AUC correction: The predicted cell types were corrected by AUC values."))

    } else {
      message("\n","SlimR AUC compute: Calculating AUC values for predicted cell type.")
      auc_values <- numeric(nrow(prediction_results))

      for (i in seq_len(nrow(prediction_results))) {
        cluster_id <- prediction_results$cluster_col[i]
        cell_type <- prediction_results$Predicted_cell_type[i]

        if (!cell_type %in% names(valid_genes_list)) {
          warning(paste("Skipping AUC for cluster", cluster_id, ": No valid genes for", cell_type))
          auc_values[i] <- NA
          next
        }
        features <- valid_genes_list[[cell_type]]

        all_cells <- colnames(seurat_obj)
        expr_data <- FetchData(seurat_obj, vars = features, cells = all_cells)
        
        labels <- seurat_obj@meta.data[all_cells, cluster_col] == cluster_id
        if (length(unique(labels)) < 2) {
          warning(paste("Skipping AUC for cluster", cluster_id, ": Only one class present"))
          auc_values[i] <- NA
        } else {
          auc_result <- compute_max_gene_auc(expr_data, labels)
          auc_values[i] <- auc_result$max_auc
          
          if (!cell_type %in% names(auc_list_storage)) {
            auc_list_storage[[cell_type]] <- matrix(NA, 
                                                     nrow = length(unique(seurat_obj@meta.data[[cluster_col]])),
                                                     ncol = length(features))
            rownames(auc_list_storage[[cell_type]]) <- unique(seurat_obj@meta.data[[cluster_col]])
            colnames(auc_list_storage[[cell_type]]) <- features
          }
          
          auc_list_storage[[cell_type]][cluster_id, ] <- auc_result$individual_aucs
        }
      }
      prediction_results$AUC <- auc_values

      message("SlimR AUC compute: AUC values for predicting cell types were calculated.")
    }

  } else {
    prediction_results$AUC <- NA
  }

  prediction_results <- prediction_results[, c("cluster_col", "Predicted_cell_type", "AUC", "Alternative_cell_types")]

  heatmap_plot <- p

  auc_plot <- NULL
  if (plot_AUC && compute_AUC) {

    Seurat::DefaultAssay(seurat_obj) <- assay

    message(paste0("\n","SlimR AUC plot: Generating combined AUC plot for predicted cell types."))

    roc_data_list <- list()

    for (i in seq_len(nrow(prediction_results))) {
      cluster_id <- prediction_results$cluster_col[i]
      cell_type <- prediction_results$Predicted_cell_type[i]
      auc_value <- prediction_results$AUC[i]

      if (is.na(cell_type)) next
      if (!cell_type %in% names(valid_genes_list)) next

      features <- valid_genes_list[[cell_type]]
      all_cells <- colnames(seurat_obj)

      expr_data <- FetchData(seurat_obj, vars = features, cells = all_cells)
      
      labels <- seurat_obj@meta.data[all_cells, cluster_col] == cluster_id
      
      if (length(unique(labels)) < 2) {
        warning(paste("Skipping AUC plot for cluster", cluster_id, ": Insufficient classes"))
        next
      }
      
      auc_result <- compute_max_gene_auc(expr_data, labels)
      
      if (is.null(auc_result$best_gene_expr) || is.na(auc_result$max_auc)) {
        warning(paste("Skipping AUC plot for cluster", cluster_id, ": No valid genes"))
        next
      }
      
      cell_scores <- auc_result$best_gene_expr

      roc_data <- compute_roc_data(cell_scores, labels)

      curve_label <- sprintf("%s - %s - %.3f", cluster_id, cell_type, auc_value)

      roc_data_list[[curve_label]] <- roc_data
    }

    if (length(roc_data_list) > 0) {
      combined_df <- do.call(rbind, lapply(names(roc_data_list), function(label) {
        data.frame(
          FPR = roc_data_list[[label]]$fpr,
          TPR = roc_data_list[[label]]$tpr,
          label = label
        )
      }))

      color_count <- length(unique(combined_df$label))
      colors <- scales::hue_pal()(color_count)

      auc_plot <- ggplot(combined_df, aes(.data$FPR, y = .data$TPR, color = .data$label)) +
        geom_line(size = 1) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
        labs(
          title = "ROC Curves for Predicted Cell Types | SlimR",
          x = "False Positive Rate (1 - Specificity)",
          y = "True Positive Rate (Sensitivity)",
          color = "Cluster - Cell Type - AUC"
        ) +
        scale_color_manual(values = colors) +
        scale_x_continuous(
          limits = c(0, 1),
          breaks = seq(0, 1, by = 0.2),
          labels = scales::number_format(accuracy = 0.1)
        ) +
        scale_y_continuous(
          limits = c(0, 1),
          breaks = seq(0, 1, by = 0.2),
          labels = scales::number_format(accuracy = 0.1)
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.ticks = element_line(color = "black"),
          axis.line = element_line(color = "black"),
          legend.position = "right",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        guides(color = guide_legend(override.aes = list(size = 2))) +
        expand_limits(x = 0, y = 0)
    } else {
      warning("No valid AUC data available for plotting.")
    }
    message(paste0("SlimR AUC plot: AUC graphs for the predicted cell types have been generated."))
  }

  return_list <- list(
    Expression_list = expr_list,
    Proportion_list = frac_list,
    Expression_scores_matrix = scores_matrix,
    Probability_matrix = probability_matrix,
    Prediction_results = prediction_results,
    Heatmap_plot = heatmap_plot
  )

  if (!is.null(auc_plot)) {
    return_list$AUC_plot <- auc_plot
  }
  
  if (compute_AUC) {
    if (exists("auc_list_storage")) {
      return_list$AUC_list <- auc_list_storage
    } else {
      return_list$AUC_list <- list()
    }
  }

  return(return_list)
}
