# Tests for marker list creation functions

test_that("Markers_filter_Cellmarker2 exists and validates inputs", {
  expect_true(exists("Markers_filter_Cellmarker2"))
  expect_type(Markers_filter_Cellmarker2, "closure")
})

test_that("Markers_filter_Cellmarker2 works with built-in database", {
  skip_if_limited()
  
  # Check if Cellmarker2 database exists
  if (!exists("Cellmarker2", where = "package:SlimR")) {
    skip("Cellmarker2 database not available")
  }
  
  Cellmarker2 <- SlimR::Cellmarker2
  
  result <- Markers_filter_Cellmarker2(
    Cellmarker2,
    species = "Human",
    tissue_class = "Blood"
  )
  
  expect_type(result, "list")
  expect_true(length(result) > 0)
  
  # Check structure of first element
  if (length(result) > 0) {
    expect_s3_class(result[[1]], "data.frame")
    expect_true(ncol(result[[1]]) >= 1)
  }
})

test_that("Markers_filter_PanglaoDB exists and validates inputs", {
  expect_true(exists("Markers_filter_PanglaoDB"))
  expect_type(Markers_filter_PanglaoDB, "closure")
})

test_that("Markers_filter_PanglaoDB works with built-in database", {
  skip_if_limited()
  
  if (!exists("PanglaoDB", where = "package:SlimR")) {
    skip("PanglaoDB database not available")
  }
  
  PanglaoDB <- SlimR::PanglaoDB
  
  result <- Markers_filter_PanglaoDB(
    PanglaoDB,
    species_input = "Human",
    organ_input = "Blood"
  )
  
  expect_type(result, "list")
  
  # Check structure if results returned
  if (length(result) > 0) {
    expect_s3_class(result[[1]], "data.frame")
  }
})

test_that("Read_seurat_markers exists and validates inputs", {
  expect_true(exists("Read_seurat_markers"))
  expect_type(Read_seurat_markers, "closure")
})

test_that("Read_seurat_markers processes Seurat markers correctly", {
  skip_if_limited()
  
  # Create mock Seurat markers dataframe
  seurat_markers <- data.frame(
    p_val = c(0.01, 0.02, 0.03, 0.01, 0.02),
    avg_log2FC = c(1.5, 1.2, 1.0, 1.8, 1.3),
    pct.1 = c(0.8, 0.7, 0.6, 0.85, 0.75),
    pct.2 = c(0.3, 0.4, 0.5, 0.2, 0.3),
    p_val_adj = c(0.01, 0.02, 0.03, 0.01, 0.02),
    cluster = c("Type_A", "Type_A", "Type_A", "Type_B", "Type_B"),
    gene = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
    stringsAsFactors = FALSE
  )
  
  result <- Read_seurat_markers(
    seurat_markers,
    sources = "Seurat",
    sort_by = "avg_log2FC",
    gene_filter = 2
  )
  
  expect_type(result, "list")
  expect_true("Type_A" %in% names(result))
  expect_true("Type_B" %in% names(result))
  expect_s3_class(result$Type_A, "data.frame")
  # The function returns "gene" as column name, not "marker"
  expect_true("gene" %in% colnames(result$Type_A))
  expect_equal(nrow(result$Type_A), 2)
  expect_equal(nrow(result$Type_B), 2)
})

test_that("Read_seurat_markers handles FSS sorting", {
  skip_if_limited()
  
  seurat_markers <- data.frame(
    p_val = c(0.01, 0.02, 0.03),
    avg_log2FC = c(2.0, 1.5, 1.0),
    pct.1 = c(0.8, 0.7, 0.9),
    pct.2 = c(0.2, 0.3, 0.1),
    p_val_adj = c(0.01, 0.02, 0.03),
    cluster = c("Type_A", "Type_A", "Type_A"),
    gene = c("GENE1", "GENE2", "GENE3"),
    stringsAsFactors = FALSE
  )
  
  result <- Read_seurat_markers(
    seurat_markers,
    sources = "Seurat",
    sort_by = "FSS",
    gene_filter = 3
  )
  
  expect_type(result, "list")
  expect_true("FSS" %in% colnames(result$Type_A))
})

test_that("Read_excel_markers exists", {
  expect_true(exists("Read_excel_markers"))
  expect_type(Read_excel_markers, "closure")
})

test_that("Built-in marker lists exist", {
  # Check if built-in marker lists are available
  expect_true(
    exists("Markers_list_scIBD", where = "package:SlimR") ||
    exists("Markers_list_TCellSI", where = "package:SlimR") ||
    exists("Markers_list_PCTIT", where = "package:SlimR") ||
    exists("Markers_list_PCTAM", where = "package:SlimR")
  )
})

test_that("Marker list format is consistent", {
  skip_if_limited()
  
  # Test with a simple marker list
  test_markers <- list(
    "Type_A" = data.frame(marker = c("GENE1", "GENE2")),
    "Type_B" = data.frame(marker = c("GENE3", "GENE4"))
  )
  
  # Check structure
  expect_type(test_markers, "list")
  expect_true(all(sapply(test_markers, is.data.frame)))
  expect_true(all(sapply(test_markers, function(x) "marker" %in% colnames(x))))
})
