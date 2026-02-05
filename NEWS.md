# SlimR 1.1.1 (2026-02-05)

## New Features

-   **Per-Cell Annotation System**: Introduced cell-level annotation as a complement to cluster-based annotation, enabling finer-grained cell type identification
    -   Three scoring methods available: weighted (default), mean, and AUCell
    -   Optional UMAP spatial smoothing for noise reduction and improved annotation consistency
    -   Optimized for large datasets with vectorized operations and memory-efficient processing

## Enhancements

-   **Improved Parameter Calculation**: Enhanced adaptive machine learning algorithms for automatic parameter optimization
-   **Better Error Handling**: Added comprehensive NA value checks and validation in `compute_adaptive_parameters()` to prevent runtime errors
-   **Performance Optimizations**: Implemented RANN package support for 10-100× faster k-NN computation in spatial smoothing

## Bug Fixes

-   Fixed NA handling issues in parameter calculation functions that could cause conditional statement errors
-   Resolved metadata naming inconsistencies in per-cell annotation workflows
-   Corrected AUC computation edge cases with filtered marker genes

## Testing

-   Added comprehensive test suite with 147 tests covering core functionality
-   Implemented robust test cases for NA handling, parameter validation, and annotation workflows
-   All tests pass with zero failures and zero warnings

## Documentation

-   Reorganized README with clearer structure distinguishing cluster-based and per-cell annotation approaches
-   Added detailed workflow guides and method comparison tables
-   Enhanced function documentation with usage examples and parameter descriptions

## Dependencies

-   Added RANN to Suggests for optional fast k-NN computation (recommended for large datasets)

# SlimR 1.1.0 (2026-01-20)

## Improvements

-   Optimized AUC calculation in `Celltype_Calculate()` to use individual gene AUCs for more robust predictions
-   Enhanced adaptive machine learning in `Parameter_Calculate()` for better model generalization
-   Added threshold parameter prediction to `Parameter_Calculate()` function
-   Updated documentation and README structure

# SlimR 1.0.9 (2025-12-18)

## New Databases

-   Integrated PCTIT database: Pan-cancer T cell markers from Zheng et al. (2021) <doi:10.1126/science.abe6474>
-   Integrated PCTAM database: Pan-cancer macrophage markers from Ma et al. (2022) <doi:10.1016/j.it.2022.04.008>

## Enhancements

-   Added `has_colnames` parameter to `Read_excel_markers()` for Excel files without headers
-   Updated documentation

# SlimR 1.0.8 (2025-10-08)

## New Features

-   Machine learning-based parameter recognition using Random Forest, Gradient Boosting, SVM, and Ensemble Learning

## Improvements

-   Optimized data filtering for `Markers_list_scIBD` using `sort_by = "logFC"` and `gene_filter = 20`
-   Enhanced FSS calculation in `Read_seurat_markers()` for presto sources
-   Improved console output in `Celltype_Verification()`

# SlimR 1.0.7 (2025-08-19)

## New Features

-   Added `Celltype_Verification()` for validation dotplot generation
-   Custom color parameters (`colour_low`, `colour_high`) for all plotting functions

## Improvements

-   Enhanced `Read_seurat_markers()` with presto compatibility and FSS calculation
-   Standardized function naming (renamed multiple functions for consistency)
-   Improved console message system
-   CRAN compliance updates

## Bug Fixes

-   Resolved various user-reported issues

# SlimR 1.0.6 (2025-08-06)

## New Features

-   Integrated scIBD human intestine reference database
-   AUC calculation and visualization in `Celltype_Calculate()`
-   AUC-based prediction correction

## Improvements

-   Streamlined output formatting
-   CRAN compliance updates

## Bug Fixes

-   Fixed critical bugs in prediction pipeline

# SlimR 1.0.5 (2025-08-05)

## New Features

-   Added TCellSI T-cell reference database
-   Introduced `Celltype_Calculate()` for automated scoring
-   Added `Celltype_Annotation()` for end-to-end annotation

## Improvements

-   Enhanced message output system
-   CRAN compliance updates

## Bug Fixes

-   Resolved multiple code errors

# SlimR 1.0.4 (2025-07-30)

## Improvements

-   Optimized `Celltype_annotation_Heatmap()` performance
-   Enhanced probability calculation in `calculate_probability()`
-   Changed license from GPL-3 to MIT
-   CRAN compliance updates

# SlimR 1.0.3 (2025-07-28)

## Improvements

-   Updated `Celltype_annotation_Heatmap()` to use `calculate_probability()`
-   CRAN compliance updates

# SlimR 1.0.1 (2025-07-19)

## Changes

-   Renamed `Celltype_annotation_Bar()` to `Celltype_annotation_Box()` with improved visualization

# SlimR 1.0.0 (2025-07-07)

-   Initial CRAN release with core annotation framework and visualization functions
