# Helper functions and test data creation for SlimR tests

#' Create a small test Seurat object
#' @keywords internal
create_test_seurat <- function(n_cells = 100, n_genes = 50, n_clusters = 3) {
  # Suppress Seurat messages
  suppressWarnings(suppressMessages({
    # Create expression matrix
    set.seed(42)
    expr_mat <- matrix(
      rnorm(n_cells * n_genes, mean = 1, sd = 0.5),
      nrow = n_genes,
      ncol = n_cells
    )
    
    # Make some genes more expressed in specific clusters
    cluster_genes <- split(1:n_genes, rep(1:n_clusters, length.out = n_genes))
    clusters <- rep(1:n_clusters, length.out = n_cells)
    
    for (i in 1:n_clusters) {
      cluster_cells <- which(clusters == i)
      cluster_specific_genes <- cluster_genes[[i]]
      expr_mat[cluster_specific_genes, cluster_cells] <- 
        expr_mat[cluster_specific_genes, cluster_cells] + 2
    }
    
    rownames(expr_mat) <- paste0("GENE", 1:n_genes)
    colnames(expr_mat) <- paste0("CELL", 1:n_cells)
    
    # Create Seurat object
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat package required for tests")
    }
    
    sce <- Seurat::CreateSeuratObject(counts = expr_mat, project = "test")
    sce <- Seurat::NormalizeData(sce, verbose = FALSE)
    
    # Find variable features first (required for PCA)
    sce <- Seurat::FindVariableFeatures(sce, selection.method = "vst", 
                                        nfeatures = min(2000, n_genes),
                                        verbose = FALSE)
    
    # Scale only variable features to avoid issues
    sce <- Seurat::ScaleData(sce, features = Seurat::VariableFeatures(sce), verbose = FALSE)
    
    # Add cluster information
    sce$seurat_clusters <- factor(clusters)
    Seurat::Idents(sce) <- sce$seurat_clusters
    
    # Add PCA and UMAP (for per-cell tests)
    var_features <- Seurat::VariableFeatures(sce)
    if (length(var_features) > 0) {
      npcs_use <- min(10, n_genes - 1, length(var_features))
      sce <- Seurat::RunPCA(sce, features = var_features, verbose = FALSE, npcs = npcs_use)
      sce <- Seurat::RunUMAP(sce, dims = 1:npcs_use, verbose = FALSE)
    }
    
    return(sce)
  }))
}

#' Create a test marker list
#' @keywords internal
create_test_markers <- function(n_genes = 50, n_types = 3) {
  genes_per_type <- n_genes %/% n_types
  marker_list <- list()
  
  for (i in 1:n_types) {
    start_gene <- (i - 1) * genes_per_type + 1
    end_gene <- min(i * genes_per_type, n_genes)
    
    marker_list[[paste0("Type_", LETTERS[i])]] <- data.frame(
      marker = paste0("GENE", start_gene:end_gene),
      stringsAsFactors = FALSE
    )
  }
  
  return(marker_list)
}

#' Skip tests if Seurat is not available
#' @keywords internal
skip_if_no_seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    testthat::skip("Seurat not available")
  }
}

#' Skip tests if running in limited environment
#' @keywords internal
skip_if_limited <- function() {
  # Skip if CRAN checks or limited memory
  if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
    testthat::skip("Skipping on CRAN")
  }
}
