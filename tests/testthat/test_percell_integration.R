# Integration tests for per-cell annotation workflow
# Tests the complete workflow from calculation to verification

test_that("Complete per-cell workflow runs successfully", {
  skip_if_no_seurat()
  skip_if_limited()
  
  # Step 1: Create test data
  sce <- create_test_seurat(n_cells = 60, n_genes = 30, n_clusters = 3)
  markers <- create_test_markers(n_genes = 30, n_types = 3)
  
  # Step 2: Calculate per-cell annotations
  result <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = markers,
    species = "Human",
    method = "weighted",
    min_score = 0.1,  # Use lower threshold for test data
    min_confidence = 1.0,  # Disable confidence filtering
    verbose = FALSE
  )
  
  expect_type(result, "list")
  expect_true("Cell_annotations" %in% names(result))
  
  # Step 3: Annotate Seurat object
  sce <- Celltype_Annotation_PerCell(
    seurat_obj = sce,
    SlimR_percell_result = result,
    plot_UMAP = FALSE,
    annotation_col = "Cell_type_Test"
  )
  
  expect_true("Cell_type_Test" %in% colnames(sce@meta.data))
  
  # Step 4: Verify annotations (only if there are enough assigned cells)
  n_assigned <- sum(sce@meta.data$Cell_type_Test != "Unassigned")
  
  if (n_assigned >= 5) {
    dotplot <- Celltype_Verification_PerCell(
      seurat_obj = sce,
      SlimR_percell_result = result,
      annotation_col = "Cell_type_Test",
      min_cells = 5
    )
    expect_s3_class(dotplot, "ggplot")
  } else {
    # Verification should gracefully fail with informative error
    expect_error(
      Celltype_Verification_PerCell(
        seurat_obj = sce,
        SlimR_percell_result = result,
        annotation_col = "Cell_type_Test",
        min_cells = 5
      ),
      "No valid cell types"
    )
  }
})

test_that("Per-cell and cluster-based workflows produce compatible outputs", {
  skip_if_no_seurat()
  skip_if_limited()
  
  sce <- create_test_seurat(n_cells = 60, n_genes = 30, n_clusters = 3)
  markers <- create_test_markers(n_genes = 30, n_types = 3)
  
  # Cluster-based
  cluster_result <- Celltype_Calculate(
    seurat_obj = sce,
    gene_list = markers,
    species = "Human",
    cluster_col = "seurat_clusters",
    compute_AUC = FALSE,
    plot_AUC = FALSE
  )
  
  # Per-cell
  percell_result <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = markers,
    species = "Human",
    verbose = FALSE
  )
  
  # Both should have Expression_list
  expect_true("Expression_list" %in% names(cluster_result))
  expect_true("Expression_list" %in% names(percell_result))
  
  # Both should have compatible structures
  expect_type(cluster_result$Expression_list, "list")
  expect_type(percell_result$Expression_list, "list")
})

test_that("Per-cell functions are properly exported", {
  # Check functions are in NAMESPACE
  expect_true(exists("Celltype_Calculate_PerCell"))
  expect_true(exists("Celltype_Annotation_PerCell"))
  expect_true(exists("Celltype_Verification_PerCell"))
  
  # Check they are functions
  expect_type(Celltype_Calculate_PerCell, "closure")
  expect_type(Celltype_Annotation_PerCell, "closure")
  expect_type(Celltype_Verification_PerCell, "closure")
})
