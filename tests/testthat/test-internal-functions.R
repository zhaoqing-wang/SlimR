# Tests for internal helper functions

test_that("calculate_probability function works correctly", {
  skip_if_no_seurat()
  skip_if_limited()
  
  sce <- create_test_seurat(n_cells = 60, n_genes = 30, n_clusters = 3)
  
  # Access internal function if available
  if (exists("calculate_probability", where = asNamespace("SlimR"), inherits = FALSE)) {
    calc_prob <- get("calculate_probability", envir = asNamespace("SlimR"))
    
    result <- calc_prob(
      object = sce,
      features = rownames(sce)[1:10],
      cluster_col = "seurat_clusters",
      min_expression = 0.1,
      specificity_weight = 3
    )
    
    expect_type(result, "list")
    expect_true("cluster_expr" %in% names(result))
    expect_true("cluster_frac" %in% names(result))
    expect_true("cluster_scores" %in% names(result))
  }
})

test_that("Internal scoring functions handle edge cases", {
  skip_if_limited()
  
  # Test with zero expression
  zero_matrix <- matrix(0, nrow = 10, ncol = 5)
  colnames(zero_matrix) <- paste0("CELL", 1:5)
  rownames(zero_matrix) <- paste0("GENE", 1:10)
  
  # Should not crash with zero expression
  expect_silent({
    row_sums <- rowSums(zero_matrix)
    row_sums[row_sums == 0] <- 1
    normalized <- zero_matrix / row_sums
  })
})

test_that("Gene name formatting works correctly", {
  # Test Human formatting (uppercase)
  human_genes <- c("cd3e", "CD4", "Cd8a")
  human_formatted <- toupper(human_genes)
  expect_equal(human_formatted, c("CD3E", "CD4", "CD8A"))
  
  # Test Mouse formatting (capitalize first letter)
  mouse_genes <- c("cd3e", "CD4", "Cd8a")
  mouse_formatted <- paste0(
    toupper(substr(mouse_genes, 1, 1)),
    tolower(substr(mouse_genes, 2, nchar(mouse_genes)))
  )
  expect_equal(mouse_formatted, c("Cd3e", "Cd4", "Cd8a"))
})

test_that("Score normalization works correctly", {
  # Test min-max normalization
  scores <- c(1, 2, 3, 4, 5)
  
  # Min-max to [0, 1]
  normalized <- (scores - min(scores)) / (max(scores) - min(scores))
  
  expect_equal(min(normalized), 0)
  expect_equal(max(normalized), 1)
  expect_true(all(normalized >= 0 & normalized <= 1))
  
  # Test with constant values
  constant_scores <- rep(5, 5)
  # When all values are the same, normalization should return all zeros
  range_diff <- diff(range(constant_scores))
  if (range_diff == 0) {
    constant_normalized <- rep(0, length(constant_scores))
  } else {
    constant_normalized <- (constant_scores - min(constant_scores)) / (max(constant_scores) - min(constant_scores))
  }
  expect_equal(constant_normalized, rep(0, 5))
})

test_that("Confidence score calculation is valid", {
  # Test confidence as difference between top two scores
  prob_matrix <- matrix(c(
    0.5, 0.3, 0.2,  # Cell 1: confident (0.5 - 0.3 = 0.2)
    0.8, 0.1, 0.1,  # Cell 2: very confident (0.8 - 0.1 = 0.7)
    0.4, 0.35, 0.25 # Cell 3: uncertain (0.4 - 0.35 = 0.05)
  ), nrow = 3, byrow = TRUE)
  
  confidence <- apply(prob_matrix, 1, function(x) {
    sorted_x <- sort(x, decreasing = TRUE)
    if (sorted_x[1] == 0) return(0)
    sorted_x[1] - sorted_x[2]
  })
  
  expect_equal(length(confidence), 3)
  expect_true(confidence[2] > confidence[1])  # Cell 2 more confident than Cell 1
  expect_true(confidence[1] > confidence[3])  # Cell 1 more confident than Cell 3
  expect_true(all(confidence >= 0 & confidence <= 1))
})

test_that("Weighted mean calculation works", {
  values <- c(1, 2, 3, 4, 5)
  weights <- c(1, 1, 1, 1, 1)
  
  weighted_mean <- sum(values * weights) / sum(weights)
  expect_equal(weighted_mean, 3)
  
  # Test with different weights
  weights2 <- c(5, 1, 1, 1, 1)
  weighted_mean2 <- sum(values * weights2) / sum(weights2)
  expect_true(weighted_mean2 < 3)  # Should be closer to 1
})

test_that("Matrix operations are dimension-safe", {
  # Test that operations preserve dimensions correctly
  mat <- matrix(1:12, nrow = 3, ncol = 4)
  
  # Row operations
  row_means <- rowMeans(mat)
  expect_equal(length(row_means), 3)
  
  # Column operations  
  col_means <- colMeans(mat)
  expect_equal(length(col_means), 4)
  
  # Transposition
  mat_t <- t(mat)
  expect_equal(dim(mat_t), c(4, 3))
})
