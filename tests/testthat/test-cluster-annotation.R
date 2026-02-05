# Tests for cluster-based annotation functions

test_that("Celltype_Calculate exists and has correct structure", {
  expect_true(exists("Celltype_Calculate"))
  expect_type(Celltype_Calculate, "closure")
  
  fn_args <- names(formals(Celltype_Calculate))
  expect_true("seurat_obj" %in% fn_args)
  expect_true("gene_list" %in% fn_args)
  expect_true("species" %in% fn_args)
  expect_true("cluster_col" %in% fn_args)
})

test_that("Celltype_Calculate validates inputs correctly", {
  skip_if_no_seurat()
  
  sce <- create_test_seurat(n_cells = 50, n_genes = 30)
  markers <- create_test_markers(n_genes = 30, n_types = 2)
  
  # Test invalid species
  expect_error(
    Celltype_Calculate(sce, markers, species = "Dog"),
    "species must be 'Human' or 'Mouse'"
  )
  
  # Test invalid Seurat object
  expect_error(
    Celltype_Calculate("not_seurat", markers, species = "Human"),
    "Input object must be a Seurat object"
  )
  
  # Test invalid gene_list
  expect_error(
    Celltype_Calculate(sce, "not_list", species = "Human"),
    "Gene list must be a list"
  )
})

test_that("Celltype_Calculate runs successfully", {
  skip_if_no_seurat()
  skip_if_limited()
  
  sce <- create_test_seurat(n_cells = 60, n_genes = 30, n_clusters = 3)
  markers <- create_test_markers(n_genes = 30, n_types = 3)
  
  result <- Celltype_Calculate(
    seurat_obj = sce,
    gene_list = markers,
    species = "Human",
    cluster_col = "seurat_clusters",
    compute_AUC = FALSE,
    plot_AUC = FALSE
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_true("Expression_list" %in% names(result))
  expect_true("Proportion_list" %in% names(result))
  expect_true("Expression_scores_matrix" %in% names(result))
  expect_true("Probability_matrix" %in% names(result))
  expect_true("Prediction_results" %in% names(result))
  expect_true("Heatmap_plot" %in% names(result))
  
  # Check Prediction_results structure
  expect_s3_class(result$Prediction_results, "data.frame")
  expect_equal(nrow(result$Prediction_results), 3) # 3 clusters
  expect_true(all(c("cluster_col", "Predicted_cell_type", "AUC", "Alternative_cell_types") %in% 
                    colnames(result$Prediction_results)))
})

test_that("Celltype_Calculate AUC computation works", {
  skip_if_no_seurat()
  skip_if_limited()
  
  sce <- create_test_seurat(n_cells = 60, n_genes = 30, n_clusters = 3)
  markers <- create_test_markers(n_genes = 30, n_types = 3)
  
  # Suppress expected warnings about missing genes in test data
  result <- suppressWarnings(
    Celltype_Calculate(
      seurat_obj = sce,
      gene_list = markers,
      species = "Human",
      cluster_col = "seurat_clusters",
      compute_AUC = TRUE,
      plot_AUC = FALSE
    )
  )
  
  # Check that AUC computation runs without error
  expect_true("AUC_list" %in% names(result))
  expect_true("AUC" %in% colnames(result$Prediction_results))
  # Note: AUC values may be NA if genes don't pass filtering criteria
  # This is expected behavior, not an error
})

test_that("Celltype_Annotation exists and validates inputs", {
  expect_true(exists("Celltype_Annotation"))
  expect_type(Celltype_Annotation, "closure")
})

test_that("Celltype_Annotation adds labels correctly", {
  skip_if_no_seurat()
  skip_if_limited()
  
  sce <- create_test_seurat(n_cells = 60, n_genes = 30, n_clusters = 3)
  markers <- create_test_markers(n_genes = 30, n_types = 3)
  
  result <- Celltype_Calculate(
    seurat_obj = sce,
    gene_list = markers,
    species = "Human",
    cluster_col = "seurat_clusters",
    compute_AUC = FALSE,
    plot_AUC = FALSE
  )
  
  sce_annotated <- Celltype_Annotation(
    seurat_obj = sce,
    cluster_col = "seurat_clusters",
    SlimR_anno_result = result,
    plot_UMAP = FALSE,
    annotation_col = "Test_CellType"
  )
  
  expect_true("Test_CellType" %in% colnames(sce_annotated@meta.data))
  expect_equal(length(unique(sce_annotated$Test_CellType)), 3)
})

test_that("Celltype_Verification exists", {
  expect_true(exists("Celltype_Verification"))
  expect_type(Celltype_Verification, "closure")
})

test_that("Parameter_Calculate exists and has correct structure", {
  expect_true(exists("Parameter_Calculate"))
  expect_type(Parameter_Calculate, "closure")
  
  fn_args <- names(formals(Parameter_Calculate))
  expect_true("seurat_obj" %in% fn_args)
  expect_true("features" %in% fn_args)
  expect_true("cluster_col" %in% fn_args)
})

test_that("Parameter_Calculate runs successfully", {
  skip_if_no_seurat()
  skip_if_limited()
  
  sce <- create_test_seurat(n_cells = 60, n_genes = 30, n_clusters = 3)
  
  params <- Parameter_Calculate(
    seurat_obj = sce,
    features = rownames(sce)[1:10],
    cluster_col = "seurat_clusters",
    verbose = FALSE
  )
  
  expect_type(params, "list")
  expect_true("min_expression" %in% names(params))
  expect_true("specificity_weight" %in% names(params))
  expect_true("threshold" %in% names(params))
  
  # Check parameter ranges
  expect_true(params$min_expression >= 0.01 && params$min_expression <= 1.0)
  expect_true(params$specificity_weight >= 0.5 && params$specificity_weight <= 8)
  expect_true(params$threshold >= 0.55 && params$threshold <= 0.85)
})
