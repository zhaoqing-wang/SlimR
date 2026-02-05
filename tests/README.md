# SlimR Tests

This directory contains the test suite for the SlimR package using the `testthat` framework.

## Test Structure

```
tests/
├── testthat.R                           # Main test runner
├── testthat/
│   ├── helper-test_data.R              # Helper functions and test data creation
│   ├── test-percell-annotation.R       # Tests for per-cell annotation functions
│   ├── test-cluster-annotation.R       # Tests for cluster-based annotation functions
│   ├── test-marker-functions.R         # Tests for marker list functions
│   ├── test_percell_integration.R      # Integration tests for complete workflows
│   └── test-internal-functions.R       # Tests for internal helper functions
```

## Running Tests

### Run all tests
```r
devtools::test()
```

### Run specific test file
```r
devtools::test(filter = "percell-annotation")
devtools::test(filter = "cluster-annotation")
devtools::test(filter = "marker-functions")
```

### Run tests with coverage
```r
covr::package_coverage()
```

### Run during package check
```r
devtools::check()
```

## Test Categories

### 1. Per-Cell Annotation Tests (`test-percell-annotation.R`)

Tests for the new per-cell annotation system:
- `Celltype_Calculate_PerCell()` function validation
- Three scoring methods (weighted, mean, AUCell)
- UMAP spatial smoothing (with and without RANN)
- Input validation and error handling
- `Celltype_Annotation_PerCell()` metadata addition
- `Celltype_Verification_PerCell()` dotplot generation

**Key Test Cases:**
- Basic per-cell calculation with weighted method
- Method comparison (weighted, mean, AUCell)
- UMAP smoothing with and without RANN package
- Species-specific gene formatting (Human/Mouse)
- Empty marker set handling
- Return scores parameter

### 2. Cluster-Based Annotation Tests (`test-cluster-annotation.R`)

Tests for traditional cluster-based annotation:
- `Celltype_Calculate()` function validation
- AUC calculation and correction
- `Celltype_Annotation()` label assignment
- `Celltype_Verification()` dotplot generation
- `Parameter_Calculate()` adaptive ML parameter tuning

**Key Test Cases:**
- Basic cluster annotation workflow
- AUC computation and plotting
- Parameter validation
- Adaptive parameter calculation
- Cluster-to-cell mapping

### 3. Marker Functions Tests (`test-marker-functions.R`)

Tests for marker list creation and management:
- `Markers_filter_Cellmarker2()` database filtering
- `Markers_filter_PanglaoDB()` database filtering
- `Read_seurat_markers()` Seurat marker processing
- `Read_excel_markers()` Excel file reading
- Built-in marker list availability

**Key Test Cases:**
- Database filtering by species/tissue
- Seurat marker dataframe processing
- FSS (Feature Significance Score) calculation
- Marker list format validation

### 4. Integration Tests (`test_percell_integration.R`)

Complete workflow tests:
- End-to-end per-cell annotation workflow
- Compatibility between cluster-based and per-cell methods
- Function export verification

**Key Test Cases:**
- Calculate → Annotate → Verify workflow
- Per-cell vs cluster-based output compatibility
- NAMESPACE exports

### 5. Internal Functions Tests (`test-internal-functions.R`)

Tests for internal helper functions:
- `calculate_probability()` score calculation
- Gene name formatting (Human/Mouse)
- Score normalization (min-max)
- Confidence score calculation
- Matrix operations

**Key Test Cases:**
- Zero expression handling
- Gene name case conversion
- Min-max normalization edge cases
- Weighted mean calculations
- Dimension preservation

## Test Helpers

### Helper Functions (`helper-test_data.R`)

**`create_test_seurat(n_cells, n_genes, n_clusters)`**
- Creates a small test Seurat object
- Includes normalized data, PCA, and UMAP
- Cluster-specific gene expression patterns

**`create_test_markers(n_genes, n_types)`**
- Creates a standardized marker list
- Compatible with test Seurat objects

**`skip_if_no_seurat()`**
- Skips test if Seurat is not available

**`skip_if_limited()`**
- Skips test in limited environments (CRAN checks)

## Test Data

Test objects are created dynamically using helper functions rather than storing static data. This ensures:
- Tests run with current package versions
- Reduced package size
- Reproducible test data (seed = 42)

## Expected Test Behavior

### All Tests Pass
All tests should pass without errors or warnings.

### Expected Skips
- Tests requiring Seurat may be skipped if not installed
- Large dataset tests may be skipped on CRAN
- UMAP smoothing tests without RANN will use fallback

### Test Coverage

**Target Coverage:** >80% of code

**Priority Areas:**
1. Main user-facing functions (100% coverage goal)
2. Per-cell annotation workflow (>90%)
3. Cluster-based annotation workflow (>90%)
4. Input validation (>95%)
5. Error handling (>85%)

## Adding New Tests

When adding new functionality:

1. **Create test file** named `test-<feature>.R`
2. **Use descriptive test names** with `test_that("description", { ... })`
3. **Include edge cases** and error conditions
4. **Use helper functions** for test data creation
5. **Add skip conditions** for optional dependencies
6. **Document test purpose** in comments

### Test Template

```r
test_that("Function does what it should", {
  skip_if_no_seurat()
  skip_if_limited()
  
  # Setup
  test_data <- create_test_seurat()
  
  # Execute
  result <- some_function(test_data)
  
  # Assert
  expect_type(result, "list")
  expect_true("expected_element" %in% names(result))
  expect_equal(nrow(result$data), expected_value)
})
```

## Continuous Integration

Tests are run automatically on:
- Local development: `devtools::test()`
- Package check: `R CMD check`
- GitHub Actions (if configured)
- Before CRAN submission

## Debugging Failed Tests

### Run single test interactively
```r
testthat::test_file("tests/testthat/test-percell-annotation.R")
```

### Run with browser
```r
# Add browser() to test code
test_that("test name", {
  data <- create_test_seurat()
  browser()  # Pause here
  result <- function_under_test(data)
  expect_true(...)
})
```

### Check test data
```r
# Manually create test data to inspect
source("tests/testthat/helper-test_data.R")
sce <- create_test_seurat()
str(sce)
```

## Dependencies for Testing

**Required:**
- testthat (>= 3.0.0)

**Suggested:**
- Seurat (for annotation tests)
- RANN (for UMAP smoothing tests)
- covr (for coverage reports)

## Notes

- Tests use `set.seed(42)` for reproducibility
- Seurat output is suppressed with `verbose = FALSE`
- UMAP plots are not generated (`plot_UMAP = FALSE`) in tests
- Tests are designed to run quickly (<5 minutes total)

## Contact

For test-related issues:
- File an issue: https://github.com/zhaoqing-wang/SlimR/issues
- Email: zhaoqingwang@mail.sdu.edu.cn
