## Version 1.1.3 (2026-03-16)

*   **Bug Fixes**
    *   **AUCell scoring now uses full-transcriptome ranking.** Previously, genes were ranked only among marker genes rather than the entire transcriptome, producing biologically incorrect percentile thresholds. The corrected implementation ranks all genes in the assay for each cell (chunked for memory efficiency), consistent with the original AUCell methodology.
    *   Removed improper `library()` calls from `Celltype_Calculate()` to comply with R packaging best practices. All dependencies are now accessed exclusively through NAMESPACE imports, preventing namespace pollution.
    *   Fixed ggplot2 deprecation warning: replaced `geom_line(size = ...)` with `geom_line(linewidth = ...)` in ROC curve plotting.

## Version 1.1.2 (2026-03-12)

*   **New Features**
    *   Integrated the ScType marker database for tissue-aware cell type annotation.
    *   Added `Markers_filter_ScType()` for filtering ScType markers by tissue type and cell name.
    *   Added `ScType`, `ScType_raw`, and `ScType_table` data objects.

## Version 1.1.1 (2026-02-05)

This release introduces a major new feature for fine-grained cell type identification alongside significant improvements to accuracy, performance, and usability.

*   **New Feature: Per-Cell Annotation System**
    *   Added a new workflow for annotating individual cells, complementing the existing cluster-based approach.
    *   Offers three scoring methods: `weighted` (default), `mean`, and `AUCell`.
    *   Includes optional UMAP-based spatial smoothing to improve annotation consistency.

*   **Major Improvements to Per-Cell Annotation**
    *   Introduced an adaptive `min_score = "auto"` threshold that scales with the number of cell types, preventing excessive "Unassigned" labels in large marker sets (e.g., 30+ types).
    *   Added a new `min_confidence` parameter for ratio-based filtering, providing more robust cell type discrimination than simple score differences.
    *   Enhanced the `weighted` scoring method with new marker specificity and IDF-like weights.
    *   Improved the `AUCell` method with adaptive `n_top` calculation and a combined binary/rank-weighted scoring strategy.
    *   Annotation results now include a `Raw_score_matrix` and a `Parameters` list for better reproducibility and debugging.

*   **Performance & Stability**
    *   Optimized for large datasets with vectorized operations and memory-efficient processing.
    *   Integrated the `RANN` package for optional 10-100x faster k-NN computations in spatial smoothing.
    *   Added comprehensive validation and error handling in core functions like `compute_adaptive_parameters()`.

*   **Testing & Documentation**
    *   Added a comprehensive test suite with 147 tests covering core functionality, NA handling, and workflows.
    *   Reorganized documentation and README with clearer distinctions between cluster-based and per-cell annotation.

## Version 1.1.0 (2026-01-20)

*   **Improvements**
    *   Optimized AUC calculation in `Celltype_Calculate()` to use individual gene AUCs for more robust predictions.
    *   Enhanced the adaptive machine learning algorithm in `Parameter_Calculate()` for better model generalization.
    *   Extended `Parameter_Calculate()` to include threshold parameter prediction.
    *   Updated general documentation and README structure.

## Version 1.0.9 (2025-12-18)

*   **New Features**
    *   Integrated two new pan-cancer immune cell reference databases:
        *   **PCTIT:** Pan-cancer T cell markers.
        *   **PCTAM:** Pan-cancer macrophage markers.

*   **Improvements**
    *   Added a `has_colnames` parameter to `Read_excel_markers()` to support reading Excel files without column headers.
    *   Updated documentation.

## Version 1.0.8 (2025-10-08)

*   **New Features**
    *   Implemented machine learning-based parameter recognition using Random Forest, Gradient Boosting, SVM, and an Ensemble learner.

*   **Improvements**
    *   Optimized data filtering for the `Markers_list_scIBD` database.
    *   Enhanced FSS (Fraction of Samples Significant) calculation in `Read_seurat_markers()` for outputs from the `presto` package.
    *   Improved console output formatting in `Celltype_Verification()`.

## Version 1.0.7 (2025-08-19)

*   **New Features**
    *   Added the `Celltype_Verification()` function for generating validation dot plots.
    *   Introduced custom color parameters (`colour_low`, `colour_high`) for all plotting functions.

*   **Improvements**
    *   Enhanced `Read_seurat_markers()` with better compatibility for `FindMarkers` results from the `presto` package.
    *   Standardized function names across the package for consistency.
    *   Improved the internal messaging system for clearer user feedback.
    *   General updates for CRAN compliance.

*   **Bug Fixes**
    *   Resolved various user-reported issues.

## Version 1.0.6 (2025-08-06)

*   **New Features**
    *   Integrated the scIBD human intestine cell reference database.
    *   Added AUC (Area Under the Curve) calculation and visualization to the `Celltype_Calculate()` function.
    *   Implemented AUC-based prediction correction.

*   **Improvements**
    *   Streamlined the formatting of function outputs.
    *   General updates for CRAN compliance.

*   **Bug Fixes**
    *   Fixed critical bugs in the core prediction pipeline.

## Version 1.0.5 (2025-08-05)

*   **New Features**
    *   Added the TCellSI T-cell reference database.
    *   Introduced the `Celltype_Calculate()` function for automated cluster scoring.
    *   Introduced the `Celltype_Annotation()` function for an end-to-end annotation workflow.

*   **Improvements**
    *   Enhanced the console message output system.
    *   General updates for CRAN compliance.

*   **Bug Fixes**
    *   Resolved multiple code errors.

## Version 1.0.4 (2025-07-30)

*   **Improvements**
    *   Optimized the performance of the `Celltype_annotation_Heatmap()` function.
    *   Enhanced the probability calculation in the helper function `calculate_probability()`.
    *   Changed the package license from GPL-3 to MIT.
    *   General updates for CRAN compliance.

## Version 1.0.3 (2025-07-28)

*   **Improvements**
    *   Updated `Celltype_annotation_Heatmap()` to use the new `calculate_probability()` function.
    *   General updates for CRAN compliance.

## Version 1.0.1 (2025-07-19)

*   **Changes**
    *   Renamed `Celltype_annotation_Bar()` to `Celltype_annotation_Box()` and improved its visualization output.

## Version 1.0.0 (2025-07-07)

*   Initial release on CRAN.
*   Provides the core framework for cluster-based cell type annotation and visualization.