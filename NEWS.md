# SlimR 1.1.1 (2026-02-05)

## Major New Features

-   **Added Per-Cell Annotation System**: Introduced three new functions for individual cell-level annotation as an alternative to cluster-based annotation:
    -   `Celltype_Calculate_PerCell()`: Main per-cell annotation engine with three scoring methods ("weighted", "mean", "AUCell") and optional UMAP spatial smoothing via k-NN.
    -   `Celltype_Annotation_PerCell()`: Integrates per-cell predictions into Seurat objects with confidence scores.
    -   `Celltype_Verification_PerCell()`: Generates validation dotplots for per-cell annotations.
-   Per-cell annotation provides finer-grained resolution than cluster-based annotation and is particularly useful for heterogeneous clusters, rare cell types, and continuous cell states.
-   UMAP spatial smoothing leverages spatial context to reduce noise and improve annotation consistency. Uses RANN package for O(n log n) k-NN computation when available.
-   Performance optimizations include vectorized operations, chunked processing for memory efficiency, and configurable parameters for large datasets.

## Improvements

-   Restructured README.md with separate sections for Cluster-Based Annotation (3.2) and Per-Cell Annotation (3.3) within the Automated Annotation Workflow.
-   Added comprehensive documentation for per-cell annotation workflow including usage examples, parameter descriptions, and method comparisons.
-   Enhanced Overview section in README.md to clearly describe both annotation approaches.

## Dependencies

-   Added RANN to Suggests for optional fast k-NN computation in UMAP spatial smoothing (10-100× speedup for large datasets).

## Documentation

-   Modified and optimized the README and NEWS files with new per-cell annotation sections.
-   Added detailed workflow documentation in `inst/doc/` including:
    -   `PerCell_Annotation_Guide.md`: Comprehensive user guide
    -   `PerCell_Summary.md`: Quick reference
    -   `Workflow_Diagram.md`: Visual comparison of annotation methods
    -   `PERCELL_ADDITIONS.md`: Technical implementation details

# SlimR 1.1.0 (2026-01-20)
-   This version optimizes the AUC calculation in `Celltype_Calculate ()` to be more robust by test all gene AUCs instead of average gene expression AUCs in the previous version.
-   The machine learning algorithm in the `Parameter_Calculate()` function is optimized to be adaptive machine learning to improve the generalization ability of the model.
-   Add prediction to the `Parameter_Calculate()` function for the threshold parameter used by the `Celltype_Calculate ()` function.
-   Modify and optimize the README and NEWS files.

# SlimR 1.0.9 (2025-12-18)
-   This version incorporates the T Cell Markers database PCTIT, with the data sourced from the article "Pan-cancer single cell landscape of tumor-infiltrating T cells". The reference literature is: L. Zheng et al. (2021) <doi:10.1126/science.abe6474>.
-   This version incorporates the Macrophage Markers database PCTAM, with the data sourced from the article "Macrophage diversity in cancer revisited in the era of single-cell omics". The reference literature is: Ruo-Yu Ma et al. (2022) <doi:10.1016/j.it.2022.04.008>.
-   Added a 'has_colnames' parameter to 'Read_excel_markers()' function to support reading Excel files without column headers by automatically naming the first column as "Markers".
-   Modify and optimize the README and NEWS files.

# SlimR 1.0.8 (2025-10-08)

-   This version adds the function of machine learning (e.g., 'Random Forest', 'Gradient Boosting', 'Support Vector Machine', 'Ensemble Learning') for cell types probability calculation parameter recognition.
-   Optimize the data filter mode of "Markers_list_scIBD" in the package, and filter through `sort_by = "logFC"` and `gene_filter = 20` parameter.
-   Adjust the calculation process of the 'FSS' value in the `read_seurat_markers()` function when 'resources' is set to 'presto'.
-   Optimize the prompt output during the execution of the `Celltype_Verification()` function.
-   Modify and optimize the README and NEWS files.

# SlimR 1.0.7 (2025-08-19)

-   Added new function `Celltype_Verification()` for predicted cell types validation and generated the validation dotplot.
-   Optimize the function 'Read_seurat_markers()', which is compatible with the 'presto::wilcoxauc()' source tag, and the 'Feature Significance Score' (FSS, which is the product of `log2FC` and `Expression ratio`) can be calculated and sorted accordingly.
-   Add custom color parameters `colour_low` and `colour_high` to all plotting output functions.
-   Renamed `Celltype_annotation_Dotplot()` to `Celltype_Annotation_Features()`, `Celltype_annotation_Box()` to `Celltype_Annotation_Combined()`, `read_seurat_markers()` to `Read_seurat_markers()`, `read_excel_markers()` to `Read_excel_markers()` for unified function naming structure.
-   Enhanced README with detailed process descriptions.
-   Optimized message output system for cleaner console feedback.
-   Resolved various code bugs reported by users.
-   Modified codebase to meet CRAN standards and policies.

# SlimR 1.0.6 (2025-08-06)

-   Integrated "scIBD" human intestine reference database.
-   Added AUC calculation and visualization to `Celltype_Calculate()`.
-   Implemented AUC-based prediction correction in cell typing.
-   Streamlined code output formatting.
-   Fixed critical bugs in the prediction pipeline.
-   Modified code to meet CRAN submission requirements.

# SlimR 1.0.5 (2025-08-05)

-   Added "TCellSI" T-cell reference database.
-   Introduced `Celltype_Calculate()` for automated scoring.
-   Added `Celltype_Annotation()` for end-to-end cell typing.
-   Improved message output system.
-   Resolved multiple code errors.
-   Modified code to meet CRAN standards and policies.

# SlimR 1.0.4 (2025-07-30)

-   Optimized `Celltype_annotation_Heatmap()` performance.
-   Enhanced probability calculation in `calculate_probability()`.
-   Modified code to meet CRAN submission requirements.
-   Change the License type from "GPL-3" to "MIT".

# SlimR 1.0.3 (2025-07-28)

-   Replaced `calculate_mean_expression()` with `calculate_probability()` in `Celltype_annotation_Heatmap()`.
-   Modified code to meet CRAN standards and policies.

# SlimR 1.0.1 (2025-07-19)

-   Changed `Celltype_annotation_Bar()` to `Celltype_annotation_Box()` with improved visualization capabilities.

# SlimR 1.0.0 (2025-07-07)

-   Initial release of the SlimR package with the core cell type annotation framework and basic visualization functions.