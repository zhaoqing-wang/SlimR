# Tests for per-cell annotation functions

test_that("Celltype_Calculate_PerCell exists and has correct structure", {
  expect_true(exists("Celltype_Calculate_PerCell"))
  expect_type(Celltype_Calculate_PerCell, "closure")
  
  # Check function has expected parameters
  fn_args <- names(formals(Celltype_Calculate_PerCell))
  expect_true("seurat_obj" %in% fn_args)
  expect_true("gene_list" %in% fn_args)
  expect_true("species" %in% fn_args)
  expect_true("method" %in% fn_args)
  expect_true("use_umap_smoothing" %in% fn_args)
})

test_that("Celltype_Calculate_PerCell validates inputs correctly", {
  skip_if_no_seurat()
  
  sce <- create_test_seurat(n_cells = 50, n_genes = 30)
  markers <- create_test_markers(n_genes = 30, n_types = 2)
  
  # Test invalid species
  expect_error(
    Celltype_Calculate_PerCell(sce, markers, species = "Dog"),
    "species must be 'Human' or 'Mouse'"
  )
  
  # Test invalid input type
  expect_error(
    Celltype_Calculate_PerCell("not_seurat", markers, species = "Human"),
    "Input must be a Seurat object"
  )
  
  # Test invalid gene_list
  expect_error(
    Celltype_Calculate_PerCell(sce, "not_list", species = "Human"),
    "gene_list must be a list"
  )
  
  # Test invalid smoothing_weight
  expect_error(
    Celltype_Calculate_PerCell(sce, markers, species = "Human", smoothing_weight = 1.5),
    "smoothing_weight must be between 0 and 1"
  )
})

test_that("Celltype_Calculate_PerCell runs with weighted method", {
  skip_if_no_seurat()
  skip_if_limited()
  
  sce <- create_test_seurat(n_cells = 50, n_genes = 30)
  markers <- create_test_markers(n_genes = 30, n_types = 2)
  
  result <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = markers,
    species = "Human",
    method = "weighted",
    verbose = FALSE
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_true("Cell_annotations" %in% names(result))
  expect_true("Cell_confidence" %in% names(result))
  expect_true("Summary" %in% names(result))
  expect_true("Expression_list" %in% names(result))
  expect_true("Proportion_list" %in% names(result))
  
  # Check Cell_annotations structure
  expect_s3_class(result$Cell_annotations, "data.frame")
  expect_equal(nrow(result$Cell_annotations), 50)
  expect_true(all(c("Cell_barcode", "Predicted_cell_type", "Max_score", "Confidence") %in% 
                    colnames(result$Cell_annotations)))
  
  # Check all cells are annotated
  expect_equal(length(result$Cell_confidence), 50)
  
  # Check predictions are valid
  expect_true(all(result$Cell_annotations$Predicted_cell_type %in% 
                    c(names(markers), "Unassigned")))
})

test_that("Celltype_Calculate_PerCell works with mean method", {
  skip_if_no_seurat()
  skip_if_limited()
  
  sce <- create_test_seurat(n_cells = 50, n_genes = 30)
  markers <- create_test_markers(n_genes = 30, n_types = 2)
  
  result <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = markers,
    species = "Human",
    method = "mean",
    verbose = FALSE
  )
  
  expect_type(result, "list")
  expect_true("Cell_annotations" %in% names(result))
  expect_equal(nrow(result$Cell_annotations), 50)
})

test_that("Celltype_Calculate_PerCell works with AUCell method", {
  skip_if_no_seurat()
  skip_if_limited()
  
  sce <- create_test_seurat(n_cells = 50, n_genes = 30)
  markers <- create_test_markers(n_genes = 30, n_types = 2)
  
  result <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = markers,
    species = "Human",
    method = "AUCell",
    verbose = FALSE
  )
  
  expect_type(result, "list")
  expect_true("Cell_annotations" %in% names(result))
  expect_equal(nrow(result$Cell_annotations), 50)
})

test_that("Celltype_Calculate_PerCell UMAP smoothing works without RANN", {
  skip_if_no_seurat()
  skip_if_limited()
  
  # Unload RANN if loaded to test fallback
  if ("RANN" %in% loadedNamespaces()) {
    try(unloadNamespace("RANN"), silent = TRUE)
  }
  
  sce <- create_test_seurat(n_cells = 50, n_genes = 30)
  markers <- create_test_markers(n_genes = 30, n_types = 2)
  
  result <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = markers,
    species = "Human",
    use_umap_smoothing = TRUE,
    k_neighbors = 10,
    smoothing_weight = 0.3,
    chunk_size = 25,
    verbose = FALSE
  )
  
  expect_type(result, "list")
  expect_equal(nrow(result$Cell_annotations), 50)
})

test_that("Celltype_Calculate_PerCell UMAP smoothing requires UMAP", {
  skip_if_no_seurat()
  
  # Create Seurat without UMAP
  sce <- create_test_seurat(n_cells = 30, n_genes = 20)
  sce@reductions$umap <- NULL
  
  markers <- create_test_markers(n_genes = 20, n_types = 2)
  
  expect_error(
    Celltype_Calculate_PerCell(
      seurat_obj = sce,
      gene_list = markers,
      species = "Human",
      use_umap_smoothing = TRUE
    ),
    "UMAP reduction.*not found"
  )
})

test_that("Celltype_Calculate_PerCell handles Mouse species correctly", {
  skip_if_no_seurat()
  skip_if_limited()
  
  sce <- create_test_seurat(n_cells = 50, n_genes = 30)
  # Change gene names to mouse format
  rownames(sce) <- paste0("Gene", 1:30)
  
  # Create new markers with correct number of genes
  markers <- list(
    "Type_A" = data.frame(marker = paste0("Gene", 1:10), stringsAsFactors = FALSE),
    "Type_B" = data.frame(marker = paste0("Gene", 11:20), stringsAsFactors = FALSE)
  )
  
  result <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = markers,
    species = "Mouse",
    verbose = FALSE
  )
  
  expect_type(result, "list")
  expect_equal(nrow(result$Cell_annotations), 50)
})

test_that("Celltype_Calculate_PerCell return_scores parameter works", {
  skip_if_no_seurat()
  skip_if_limited()
  
  sce <- create_test_seurat(n_cells = 50, n_genes = 30)
  markers <- create_test_markers(n_genes = 30, n_types = 2)
  
  result_with_scores <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = markers,
    species = "Human",
    return_scores = TRUE,
    verbose = FALSE
  )
  
  result_without_scores <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = markers,
    species = "Human",
    return_scores = FALSE,
    verbose = FALSE
  )
  
  expect_true("Cell_scores" %in% names(result_with_scores))
  expect_false("Cell_scores" %in% names(result_without_scores))
  
  # Check score matrix dimensions
  expect_equal(nrow(result_with_scores$Cell_scores), 50)
  expect_equal(ncol(result_with_scores$Cell_scores), 2) # 2 cell types
})

test_that("Celltype_Calculate_PerCell handles empty marker sets gracefully", {
  skip_if_no_seurat()
  
  sce <- create_test_seurat(n_cells = 30, n_genes = 20)
  
  # Create markers with non-existent genes
  markers <- list(
    "Type_A" = data.frame(marker = c("NONEXISTENT1", "NONEXISTENT2")),
    "Type_B" = data.frame(marker = c("NONEXISTENT3", "NONEXISTENT4"))
  )
  
  expect_error(
    Celltype_Calculate_PerCell(
      seurat_obj = sce,
      gene_list = markers,
      species = "Human",
      verbose = FALSE
    ),
    "No valid marker genes found"
  )
})

test_that("Celltype_Annotation_PerCell exists and validates inputs", {
  expect_true(exists("Celltype_Annotation_PerCell"))
  expect_type(Celltype_Annotation_PerCell, "closure")
})

test_that("Celltype_Annotation_PerCell adds metadata correctly", {
  skip_if_no_seurat()
  skip_if_limited()
  
  sce <- create_test_seurat(n_cells = 50, n_genes = 30)
  markers <- create_test_markers(n_genes = 30, n_types = 2)
  
  result <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = markers,
    species = "Human",
    verbose = FALSE
  )
  
  sce_annotated <- Celltype_Annotation_PerCell(
    seurat_obj = sce,
    SlimR_percell_result = result,
    plot_UMAP = FALSE,
    annotation_col = "Test_Annotation"
  )
  
  # Check metadata columns added
  expect_true("Test_Annotation" %in% colnames(sce_annotated@meta.data))
  expect_true("Test_Annotation_score" %in% colnames(sce_annotated@meta.data))
  expect_true("Test_Annotation_confidence" %in% colnames(sce_annotated@meta.data))
  
  # Check values match
  expect_equal(
    sce_annotated$Test_Annotation,
    result$Cell_annotations$Predicted_cell_type
  )
})

test_that("Celltype_Verification_PerCell exists and validates inputs", {
  expect_true(exists("Celltype_Verification_PerCell"))
  expect_type(Celltype_Verification_PerCell, "closure")
})

test_that("Celltype_Verification_PerCell creates dotplot", {
  skip_if_no_seurat()
  skip_if_limited()
  
  sce <- create_test_seurat(n_cells = 50, n_genes = 30)
  markers <- create_test_markers(n_genes = 30, n_types = 2)
  
  result <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = markers,
    species = "Human",
    verbose = FALSE
  )
  
  sce <- Celltype_Annotation_PerCell(
    seurat_obj = sce,
    SlimR_percell_result = result,
    plot_UMAP = FALSE,
    annotation_col = "Test_Annotation"
  )
  
  dotplot <- Celltype_Verification_PerCell(
    seurat_obj = sce,
    SlimR_percell_result = result,
    annotation_col = "Test_Annotation",
    min_cells = 5
  )
  
  expect_s3_class(dotplot, "ggplot")
})
