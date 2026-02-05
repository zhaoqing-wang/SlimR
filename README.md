# SlimR: Adaptive Machine Learning-Powered, Context-Matching Tool for Single-Cell and Spatial Transcriptomics Annotation

[![CRAN Package Version](https://img.shields.io/cran/v/SlimR?label=CRAN)](https://cran.r-project.org/package=SlimR) [![CRAN License](https://img.shields.io/cran/l/SlimR?label=License&color=green)](https://cran.r-project.org/package=SlimR) [![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/SlimR)](https://cran.r-project.org/package=SlimR) [![GitHub Package Version](https://img.shields.io/github/r-package/v/zhaoqing-wang/SlimR?label=GitHub&color=green)](https://github.com/zhaoqing-wang/SlimR/releases) [![GitHub Maintainer](https://img.shields.io/badge/Maintainer-Zhaoqing_Wang-blue)](https://github.com/zhaoqing-wang)

## Overview

<img src="docs/Sticker.png" alt="Sticker" width="233.28" height="270" align="right"/>

SlimR is an R package for annotating single-cell and spatial transcriptomics datasets. It creates a unified marker list (`Markers_list`) from multiple sources: built-in curated databases (`Cellmarker2`, `PanglaoDB`, `scIBD`, `TCellSI`, `PCTIT`, `PCTAM`), Seurat objects with cell labels, or user-provided Excel tables.

SlimR first uses adaptive machine learning for parameter optimization, and then offers two automated annotation approaches: "cluster-based" and "per-cell". Cluster-based annotation assigns one label per cluster, expression-based probability calculation, and AUC validation. Per-cell annotation assigns labels to individual cells using three scoring methods, with optional UMAP spatial smoothing, making it ideal for heterogeneous clusters and rare cell types. The package also supports semi-automated workflows with heatmaps, feature plots, and combined visualizations for manual annotation.

## Table of Contents

1.  [Preparation](#1-preparation)
    -   [1.1 Installation](#11-installation)
    -   [1.2 Loading SlimR](#12-loading-slimr)
    -   [1.3 Prepare Seurat Object](#13-prepare-seurat-object)
    -   [1.4 Dependencies (Alternative)](#14-dependencies-alternative)
2.  [Standardized Markers_list Input](#2-standardized-markers_list-input)
    -   [2.1 From Cellmarker2 Database](#21-from-cellmarker2-database)
    -   [2.2 From PanglaoDB Database](#22-from-panglaodb-database)
    -   [2.3 From Seurat Objects](#23-from-seurat-objects)
    -   [2.4 From Excel Tables](#24-from-excel-tables)
    -   [2.5 Example: From Article scIBD](#25-example-from-article-scibd)
    -   [2.6 Example: From Tool TCellSI](#26-example-from-tool-tcellsi)
    -   [2.7 Example: From Atlas of Pan Cancer T Cells](#27-example-from-atlas-of-pan-cancer-t-cells)
    -   [2.8 Example: From Review of Pan Cancer Macrophages](#28-example-from-review-of-pan-cancer-macrophages)
3.  [Automated Annotation Workflow](#3-automated-annotation-workflow)
    -   [3.1 Calculate Parameter](#31-calculate-parameter)
    -   [3.2 Cluster-Based Annotation](#32-cluster-based-annotation)
        -   [3.2.1 Calculate Cell Types](#321-calculate-cell-types)
        -   [3.2.2 Annotate Cell Types](#322-annotate-cell-types)
        -   [3.2.3 Verify Cell Types](#323-verify-cell-types)
    -   [3.3 Per-Cell Annotation](#33-per-cell-annotation)
        -   [3.3.1 Calculate Per-Cell Types](#331-calculate-per-cell-types)
        -   [3.3.2 Annotate Per-Cell Types](#332-annotate-per-cell-types)
        -   [3.3.3 Verify Per-Cell Types](#333-verify-per-cell-types)
4.  [Semi-Automated Annotation Workflow](#4-semi-automated-annotation-workflow)
    -   [4.1 Annotation Heat Map](#41-annotation-heat-map)
    -   [4.2 Annotation Feature Plots](#42-annotation-feature-plots)
    -   [4.3 Annotation Combined Plots](#43-annotation-combined-plots)
5.  [Other Functions Provided by SlimR](#5-other-functions-provided-by-slimr)
6.  [Conclusion](#6-conclusion)

------------------------------------------------------------------------

## 1. Preparation

### 1.1 Installation

**Option One: CRAN** [![CRAN Version](https://img.shields.io/cran/v/SlimR?label=CRAN)](https://cran.r-project.org/package=SlimR)

Install the stable version from CRAN (recommended when CRAN and GitHub versions match):

``` r
install.packages("SlimR")
```

*Note: If you encounter version mismatches, try adjusting the CRAN mirror to `Global (CDN)` or use `BiocManager::install("SlimR")`.*

**Option Two: GitHub** [![GitHub R package version](https://img.shields.io/github/r-package/v/zhaoqing-wang/SlimR?label=GitHub&color=green)](https://github.com/zhaoqing-wang/SlimR/releases)

Install the development version from GitHub (recommended when GitHub version is newer):

``` r
devtools::install_github("zhaoqing-wang/SlimR")
```

*Note: If the function doesn't work, please run `install.packages('devtools')` first.*

### 1.2 Loading SlimR

Load the package in your R environment:

``` r
library(SlimR)
```

### 1.3 Prepare Seurat Object

For Seurat objects with multiple layers, join layers before annotation:

``` r
# Example: Join layers in the RNA assay
sce@assays$RNA <- SeuratObject::JoinLayers(sce@assays$RNA)
```

**Important: Ensure your Seurat object has completed standard preprocessing (normalization, scaling, clustering) and batch effect correction for accurate annotation.**

*Tip: Use the `clustree` package to determine optimal clustering resolution.*

### 1.4 Dependencies (Optional)

SlimR requires R (≥ 3.5) and the following packages: `cowplot`, `dplyr`, `ggplot2`, `patchwork`, `pheatmap`, `readxl`, `scales`, `Seurat`, `tidyr`, `tools`. If installation fails, manually install missing dependencies:

``` r
# Install dependencies if needed:
install.packages(c("cowplot", "dplyr", "ggplot2", "patchwork", 
                   "pheatmap", "readxl", "scales", "Seurat", 
                   "tidyr", "tools"))
```

**Optional dependency for Per-Cell Annotation:**

For faster UMAP spatial smoothing in per-cell annotation (10-100× speedup), install the RANN package:

``` r
# Optional: Install RANN for fast k-NN computation
install.packages("RANN")
```

*Note: RANN is optional. Per-cell annotation works without it but uses a slower fallback method for UMAP smoothing.*

## 2. Standardized Markers_list Input

SlimR uses a standardized list format where:
- List names = cell types (required)
- First column = marker genes (required)  
- Additional columns = metrics (optional)

### 2.1 From Cellmarker2 Database

Cellmarker2 is a comprehensive database of cell types and markers across multiple species and tissues.

**Reference:** *Hu et al. (2023) <doi:10.1093/nar/gkac947>*

#### 2.1.1 Load Database:

``` r
Cellmarker2 <- SlimR::Cellmarker2
```

#### 2.1.2 Optional Metadata Exploration:

``` r
Cellmarker2_table <- SlimR::Cellmarker2_table
View(Cellmarker2_table)
```

#### 2.1.3 Generate `Markers_list`:

``` r
Markers_list_Cellmarker2 <- Markers_filter_Cellmarker2(
  Cellmarker2,
  species = "Human",
  tissue_class = "Intestine",
  tissue_type = NULL,
  cancer_type = NULL,
  cell_type = NULL
)
```

**Important: Specify at least `species` and `tissue_class` for accurate annotations.**

*The resulting `Markers_list` can be used in automated (Section 3) and semi-automated (Section 4) workflows. [Jump to Section 3](#3-automated-annotation-workflow)*

### 2.2 From PanglaoDB Database

PanglaoDB is a database of cell types and markers across multiple species and organs.

**Reference:** *Franzén et al. (2019) <doi:10.1093/database/baz046>*

#### 2.2.1 Load Database:

``` r
PanglaoDB <- SlimR::PanglaoDB
```

#### 2.2.2 Optional Metadata Exploration:

``` r
PanglaoDB_table <- SlimR::PanglaoDB_table
View(PanglaoDB_table)
```

#### 2.2.3 Generate `Markers_list`:

``` r
Markers_list_panglaoDB <- Markers_filter_PanglaoDB(
  PanglaoDB,
  species_input = 'Human',
  organ_input = 'GI tract'
)
```

**Important: Select the `species_input` and `organ_input` parameters to ensure the accuracy of the annotation.**

*Link: Output `Markers_list` usable in sections 3.1, 4.1, 4.2, 4.3, and 5.2. [Click to section 3 automated annotation workflow.](#3-automated-annotation-workflow)*

### 2.3 From Seurat Objects

#### 2.3.1 Identify Markers and Generate `Markers_list`

Generate `Markers_list` from Seurat differential expression results:

``` r
seurat_markers <- Seurat::FindAllMarkers(
    object = sce,
    group.by = "Cell_type",
    only.pos = TRUE)

Markers_list_Seurat <- Read_seurat_markers(seurat_markers,
    sources = "Seurat",
    sort_by = "FSS",
    gene_filter = 20
    )
```

*Tip: Use `sort_by = "FSS"` for Feature Significance Score (log2FC × Expression ratio) ranking, or `sort_by = "avg_log2FC"` for fold-change ranking.*

#### 2.3.2 Use `presto` for Speed (Optional)

For large datasets, `presto::wilcoxauc()` provides ~10× faster computation with slight accuracy trade-offs:

``` r
seurat_markers <- dplyr::filter(
    presto::wilcoxauc(
      X = sce,
      group_by = "Cell_type",
      seurat_assay = "RNA"
      ),
    padj < 0.05, logFC > 0.5
    )

Markers_list_Seurat <- Read_seurat_markers(seurat_markers,
    sources = "presto",
    sort_by = "FSS",
    gene_filter = 20
    )
```

**Important: Install presto with `devtools::install_github('immunogenomics/presto')` if needed.**

*Tip: Use `sort_by = "logFC"` or `sort_by = "FSS"` for marker ranking.*

*The resulting `Markers_list` can be used in automated and semi-automated workflows. [Jump to Section 3](#3-automated-annotation-workflow)*

### 2.4 From Excel Tables

**Format Requirements**:

-   Each sheet name = cell type (essential)

-   First row = column headers (essential)

-   First column = markers (essential)

-   Subsequent columns = metrics (can be omitted)

``` r
Markers_list_Excel <- Read_excel_markers("D:/Laboratory/Marker_load.xlsx")
```

**Important: If your Excel file lacks column headers, set `has_colnames = FALSE` in `Read_excel_markers()`.**

*The resulting `Markers_list` can be used in automated and semi-automated workflows. [Jump to Section 3](#3-automated-annotation-workflow)*

### 2.5 Example: From Article scIBD

scIBD: The Human Intestinal Cell Database (Inflammatory Bowel Disease).

Reference: *Nie et al. (2023) <doi:10.1038/s43588-023-00464-9>*.

``` r
Markers_list_scIBD <- SlimR::Markers_list_scIBD
```

**Important: This is for human intestinal annotation only. The input Seurat object was ensured to be of a human intestinal type to ensure the accuracy of the labeling.**

*Note: The `Markers_list_scIBD` was generated using section 2.3.2 and the parameters `sort_by = "logFC"` and `gene_filter = 20` were set.*

*Link: Output `Markers_list` usable in sections 3.1, 4.1, 4.2, 4.3, and 5.3. [Click to section 3 automated annotation workflow.](#3-automated-annotation-workflow)*

### 2.6 Example: From Tool TCellSI

TCellSI: A database of T cell markers of different subtypes.

Reference: *Yang et al. (2024) <doi:10.1002/imt2.231>*.

``` r
Markers_list_TCellSI <- SlimR::Markers_list_TCellSI
```

**Important: This is only used for annotation of T cell subsets. It was ensured that the input Seurat subjects were T cell subsets to ensure the accuracy of labeling.**

*Note: The `Markers_list_TCellSI` was generated using section 2.4.*

*Link: Output `Markers_list` usable in sections 3.1, 4.1, 4.2, 4.3, and 5.4. [Click to section 3 automated annotation workflow.](#3-automated-annotation-workflow)*

### 2.7 Example: From Atlas of Pan Cancer T Cells

PCTIT: List of T cell subtype markers in the article "Pan-cancer single cell landscape of tumor-infiltrating T cells".

Reference: *L. Zheng et al. (2021) <doi:10.1126/science.abe6474>*.

``` r
Markers_list_PCTIT <- SlimR::Markers_list_PCTIT
```

**Important: This is only used for annotation of T cell subsets. It was ensured that the input Seurat subjects were T cell subsets to ensure the accuracy of labeling.**

*Note: The `Markers_list_PCTIT` was generated using section 2.4.*

*Link: Output `Markers_list` usable in sections 3.1, 4.1, 4.2, 4.3, and 5.4. [Click to section 3 automated annotation workflow.](#3-automated-annotation-workflow)*

### 2.8 Example: From Review of Pan Cancer Macrophages

PCTAM: List of Macrophage subtype markers in the article "Macrophage diversity in cancer revisited in the era of single-cell omics".

Reference: *Ruo-Yu Ma et al. (2022) <doi:10.1016/j.it.2022.04.008>*.

``` r
Markers_list_PCTAM <- SlimR::Markers_list_PCTAM
```

**Important: This is only used for annotation of Macrophage subsets. It was ensured that the input Seurat subjects were Macrophage subsets to ensure the accuracy of labeling.**

*Note: The `Markers_list_PCTAM` was generated using section 2.4.*

*Link: Output `Markers_list` usable in sections 3.1, 4.1, 4.2, 4.3, and 5.4. [Click to section 3 automated annotation workflow.](#3-automated-annotation-workflow)*

## 3. Automated Annotation Workflow

SlimR provides two automated annotation approaches: **Cluster-Based Annotation** (Section 3.2) and **Per-Cell Annotation** (Section 3.3). Both workflows share the same parameter calculation step (Section 3.1) and use the same standardized `Markers_list` format.

**Comparison of Annotation Approaches:**

| Feature | Cluster-Based | Per-Cell |
|---------|---------------|----------|
| **Annotation Unit** | Cluster (all cells in cluster get same label) | Individual cell |
| **Speed** | Fast (~10-30s for 50k cells) | Slower (~2-3min for 50k cells) |
| **Memory** | Low (~800MB for 50k cells) | Higher (~2-2.5GB for 50k cells) |
| **Resolution** | Coarse (cluster-level) | Fine (cell-level) |
| **Best For** | Homogeneous, well-separated clusters | Mixed clusters, rare cell types, continuous states |
| **Confidence Scores** | Cluster-level | Cell-level |
| **Spatial Context** | Not used | Optional (UMAP smoothing) |

**Recommendation:** Start with cluster-based annotation for initial exploration. Use per-cell annotation when clusters contain mixed populations or when finer resolution is needed.

### 3.1 Calculate Parameters

SlimR uses adaptive machine learning to automatically determine optimal `min_expression`, `specificity_weight`, and `threshold` parameters for cell type probability calculation.

``` r
# Basic usage uses default genes
SlimR_params <- Parameter_Calculate(
  seurat_obj = sce,
  features = c("CD3E", "CD4", "CD8A"),
  assay = "RNA",
  cluster_col = "seurat_clusters",
  verbose = TRUE
  )
 
 # Use with custom method: use the genes corresponding to a specific cell type in 'Markers_list' as input
 SlimR_params <- Parameter_Calculate(
  seurat_obj = sce,
  features = unique(Markers_list_Cellmarker2$`B cell`$marker),
  assay = "RNA",
  cluster_col = "seurat_clusters",
  verbose = TRUE
  )
```

**Important: This step is optional. Skip to Section 3.2 to use default parameters.**

### 3.2 Cluster-Based Annotation

Cluster-based annotation assigns a single cell type label to all cells within each cluster. This approach is computationally efficient and works well when clusters are homogeneous.

#### 3.2.1 Calculate Cell Types (Core)

Calculate cell type probabilities, generate predictions with optional AUC validation, and create heatmaps and ROC curves:

``` r
SlimR_anno_result <- Celltype_Calculate(seurat_obj = sce,
    gene_list = Markers_list,
    species = "Human",
    cluster_col = "seurat_clusters",
    assay = "RNA",
    min_expression = 0.1,
    specificity_weight = 3,
    threshold = 0.6,
    compute_AUC = TRUE,
    plot_AUC = TRUE,
    AUC_correction = FALSE,
    colour_low = "navy",
    colour_high = "firebrick3"
    )
```

*Tip: If you ran `Parameter_Calculate()` in Section 3.1, use the calculated parameters:*
```r
min_expression = SlimR_params$min_expression,
specificity_weight = SlimR_params$specificity_weight,
threshold = SlimR_params$threshold
```

**Important: Use the same `cluster_col` value in `Celltype_Calculate()` and `Celltype_Annotation()` to avoid mismatches.**

*Note: `AUC_correction = TRUE` increases runtime by ~40% but improves prediction accuracy. Lower `threshold` values check more alternative cell types, increasing computation time.*

**Error Handling:** If you see "duplicate 'row.names' are not allowed", fix with:

``` r
rownames(sce) <- base::make.unique(rownames(sce))
```

**View Heatmap (Optional)**

Check annotation probabilities for clusters and cell types:

``` r
print(SlimR_anno_result$Heatmap_plot)
```

*Tip: If the heatmap doesn't display, load pheatmap: `library(pheatmap)`*

**View Predictions (Optional)**

View predicted cell type results:

``` r
View(SlimR_anno_result$Prediction_results)
```

**View ROC Curves (Optional)**

View ROC curves and AUC values for predictions:

``` r
print(SlimR_anno_result$AUC_plot)
```

**Important: Requires `plot_AUC = TRUE` in `Celltype_Calculate()`.**

*Tip: If plots don't display, load ggplot2: `library(ggplot2)`*

**Correct Predictions (Optional)**

After reviewing predictions and AUC values, manually correct cell types:

Example 1:

``` r
# For example, cluster '15' in 'cluster_col' corresponds to cell type 'Intestinal stem cell'.
SlimR_anno_result$Prediction_results$Predicted_cell_type[
  SlimR_anno_result$Prediction_results$cluster_col == 15
] <- "Intestinal stem cell"
```

Example 2:

``` r
# For example, a predicted cell type with an AUC of 0.5 or less should be labeled 'Unknown'.
SlimR_anno_result$Prediction_results$Predicted_cell_type[
  SlimR_anno_result$Prediction_results$AUC <= 0.5
] <- "Unknown"
```

View the updated predictions:

``` r
View(SlimR_anno_result$Prediction_results)
```

**Important: When correcting, preferably use cell types from `Alternative_cell_types` column.**

#### 3.2.2 Annotate Cell Types

Transfer predicted cell types from SlimR results to Seurat object metadata:

``` r
sce <- Celltype_Annotation(seurat_obj = sce,
    cluster_col = "seurat_clusters",
    SlimR_anno_result = SlimR_anno_result,
    plot_UMAP = TRUE,
    annotation_col = "Cell_type_SlimR"
    )
```

**Important: Use matching `cluster_col` values in `Celltype_Calculate()` and `Celltype_Annotation()`, and matching `annotation_col` values in `Celltype_Annotation()` and `Celltype_Verification()`.**

#### 3.2.3 Verify Cell Types

Generate validation dotplot using Feature Significance Score (FSS = log2FC × Expression ratio) for marker ranking:

``` r
Celltype_Verification(seurat_obj = sce,
    SlimR_anno_result = SlimR_anno_result,
    gene_number = 5,
    assay = "RNA",
    colour_low = "white",
    colour_high = "navy",
    annotation_col = "Cell_type_SlimR"
    )
```

**Important: Use the same `annotation_col` value in both `Celltype_Annotation()` and `Celltype_Verification()`.**

*Note: Markers from `Expression_list` are used for cell types in `Prediction_results`; other cell types use markers from `FindMarkers()`.*

### 3.3 Per-Cell Annotation

Per-cell annotation assigns cell type labels to individual cells based on marker gene expression profiles, providing finer-grained resolution than cluster-based annotation. This approach is particularly useful when clusters contain heterogeneous populations or when cell states exist on a continuum.

**When to use Per-Cell Annotation:**

-   Clusters contain mixed cell types or transitional states
-   Need fine-grained resolution for rare cell types
-   Cell states are continuous (e.g., differentiation gradients)
-   Want to leverage spatial context via UMAP smoothing

**When to use Cluster-Based Annotation:**

-   Clusters are well-separated and homogeneous
-   Computational efficiency is critical (cluster-based is faster)
-   Dataset is very large (\>200k cells)
-   Want stable, discrete categories

#### 3.3.1 Calculate Per-Cell Types (Core)

Uses `markers_list` to calculate per-cell scores and assign cell type labels to individual cells. Three scoring methods are available: `"weighted"` (default, recommended), `"mean"` (fast baseline), and `"AUCell"` (rank-based, robust to batch effects).

``` r
SlimR_percell_result <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = Markers_list,
    species = "Human",
    assay = "RNA",
    method = "weighted",
    min_expression = 0.1,
    use_umap_smoothing = FALSE,
    umap_reduction = "umap",
    k_neighbors = 15,
    smoothing_weight = 0.3,
    min_score = 0.1,
    return_scores = FALSE,
    verbose = TRUE
    )
```

You can use the `min_expression = SlimR_params$min_expression` parameter in the function `Celltype_Calculate_PerCell()` if you have run the `Parameter_Calculate ()` function in section 3.1 above.

**Important: Per-cell annotation requires normalized data. Make sure your Seurat object has been processed with `NormalizeData()`.**

*Note: The three scoring methods differ in their approach:*

-   **"weighted"**: Combines expression level, detection rate, and marker specificity (CV-based weighting). Best for general use.
-   **"mean"**: Simple average of normalized marker expression. Fastest, good for initial exploration.
-   **"AUCell"**: Rank-based scoring using proportion of markers in top 5% expressed genes. Robust to batch effects and technical variation.

**UMAP Spatial Smoothing (Optional)**

Enable UMAP-based spatial smoothing to reduce noise and improve annotation consistency by incorporating information from spatially neighboring cells:

``` r
SlimR_percell_result <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = Markers_list,
    species = "Human",
    method = "weighted",
    use_umap_smoothing = TRUE,
    k_neighbors = 20,
    smoothing_weight = 0.3
    )
```

**Important: UMAP smoothing requires a UMAP reduction in the Seurat object. Run `RunUMAP()` first if not already computed. For faster k-NN computation, install the RANN package: `install.packages("RANN")`.**

*Note: The `k_neighbors` parameter controls how many neighboring cells to consider (recommended: 15-30). The `smoothing_weight` parameter controls the blend between a cell's own score and its neighbors' average (0-1, where 0.3 means 30% weight to neighbors). Higher values produce smoother annotations but may blur boundaries.*

**View Per-Cell Annotation Summary (Optional)**

Cell type annotation summary can be viewed with the following code:

``` r
View(SlimR_percell_result$Summary)
```

**View Per-Cell Annotations (Optional)**

Individual cell annotations with confidence scores can be viewed with the following code:

``` r
View(SlimR_percell_result$Cell_annotations)
```

#### 3.3.2 Annotate Per-Cell Types

Assigns SlimR per-cell predicted cell types information from `SlimR_percell_result$Cell_annotations$Predicted_cell_type` directly to individual cells in the Seurat object, and stores the results into `seurat_obj@meta.data$annotation_col`.

``` r
sce <- Celltype_Annotation_PerCell(
    seurat_obj = sce,
    SlimR_percell_result = SlimR_percell_result,
    plot_UMAP = TRUE,
    annotation_col = "Cell_type_PerCell_SlimR",
    plot_confidence = TRUE
    )
```

**Important: The parameter `annotation_col` in the functions `Celltype_Annotation_PerCell()` and `Celltype_Verification_PerCell()` must be strictly the same to avoid false matches.**

*Note: This function also adds `annotation_col_score` (max score per cell) and `annotation_col_confidence` (confidence score per cell) to the Seurat object's meta.data for quality control purposes.*

#### 3.3.3 Verify Per-Cell Types

Use the cell type identity information in `seurat_obj@meta.data$annotation_col` and use the 'Feature Significance Score' (FSS, product value of `log2FC` and `Expression ratio`) as the ranking basis to generate validation dotplot.

``` r
Celltype_Verification_PerCell(
    seurat_obj = sce,
    SlimR_percell_result = SlimR_percell_result,
    gene_number = 5,
    assay = "RNA",
    colour_low = "white",
    colour_high = "navy",
    annotation_col = "Cell_type_PerCell_SlimR",
    min_cells = 10
    )
```

**Important: The parameter `annotation_col` in the function `Celltype_Annotation_PerCell()` and the function `Celltype_Verification_PerCell()` must be strictly the same to avoid false matches.**

*Note: Cell types with fewer than `min_cells` cells (default: 10) are excluded from the verification plot. Cell types in `SlimR_percell_result$Expression_list` are verified using the markers information from that list; cell types not in the list are validated using markers from the function `FindMarkers()`.*

## 4. Semi-Automated Annotation Workflow

### 4.1 Annotation Heat Map

Generate a heat map to estimate the likelihood that various cell clusters exhibited similarity to control cell types:

``` r
Celltype_Annotation_Heatmap(
  seurat_obj = sce,
  gene_list = Markers_list,
  species = "Human",
  cluster_col = "seurat_cluster",
  min_expression = 0.1,
  specificity_weight = 3,
  colour_low = "navy",
  colour_high = "firebrick3"
)
```

*Note: Now this function has been incorporated into `Celltype_Calculate()`, and it is recommended to use `Celltype_Calculate()` instead.*

### 4.2 Annotation Feature Plots

Generates per-cell-type expression dot plot with metric heat map (when the metric information exists):

``` r
Celltype_Annotation_Features(
  seurat_obj = sce,
  cluster_col = "seurat_clusters",
  gene_list = Markers_list,
  gene_list_type = "Cellmarker2",
  species = "Human",
  save_path = "./SlimR/Celltype_Annotation_Features/",
  colour_low = "white",
  colour_high = "navy",
  colour_low_mertic = "white",
  colour_high_mertic = "navy"
  )
```

Each resulting combined image consists of a dot plot above and a heat map below (if metric information is present). The dot plot illustrates the relationship between the expression level and expression ratio of the cell type and its corresponding markers. Below it, a metric heat map is displayed for the corresponding markers (if metric information is available).

### 4.3 Annotation Combined Plots

Generates per-cell-type expression combined plots:

``` r
Celltype_Annotation_Combined(
  seurat_obj = sce,
  gene_list = Markers_list, 
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_Annotation_Combined/",
  colour_low = "white",
  colour_high = "navy"
)
```

Each generated combined plot shows the box plot of the expression levels of the corresponding markers for that cell type, with the colors corresponding to the average expression levels of the markers.

## 5. Other Functions Provided by SlimR

Functions in sections 5.1, 5.2, 5.3, and 5.4 have been incorporated into `Celltype_Annotation_Features()`, and it is recommended to use `Celltype_Annotation_Features()` and set corresponding parameters (for example, `gene_list_type = "Cellmarker2"`) instead. For more information, please refer to section 4.2.

#### 5.1 Annotation Feature Plots with Cellmarker2 Database

``` r
Celltype_annotation_Cellmarker2(
  seurat_obj = sce,
  gene_list = Markers_list_Cellmarker2,
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_Cellmarkers2/",
  colour_low = "white",
  colour_high = "navy",
  colour_low_mertic = "white",
  colour_high_mertic = "navy"
)
```

*Note: To call this function, set the parameter `gene_list_type = "Cellmarker2"` in the function `Celltype_Annotation_Features()`.*

#### 5.2 Annotation Feature Plots with PanglaoDB Database

``` r
Celltype_annotation_PanglaoDB(
  seurat_obj = sce,
  gene_list = Markers_list_panglaoDB,
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_PanglaoDB/",
  colour_low = "white",
  colour_high = "navy",
  colour_low_mertic = "white",
  colour_high_mertic = "navy"
)
```

*Note: To call this function, set the parameter `gene_list_type = "PanglaoDB"` in the function `Celltype_Annotation_Features()`.*

#### 5.3 Annotation Feature Plots with Seurat-Based Markers List

``` r
Celltype_annotation_Seurat(
  seurat_obj = sce,
  gene_list = Markers_list_Seurat,
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_Seurat/",
  colour_low = "white",
  colour_high = "navy",
  colour_low_mertic = "white",
  colour_high_mertic = "navy"
)
```

*Note: To call this function, set the parameter `gene_list_type = "Seurat"` in the function `Celltype_Annotation_Features()`.*

#### 5.4 Annotation Feature Plots with Excel-Based Markers List

``` r
Celltype_annotation_Excel(
  seurat_obj = sce,
  gene_list = Markers_list_Excel,
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_Excel/",
  colour_low = "white",
  colour_high = "navy",
  colour_low_mertic = "white",
  colour_high_mertic = "navy"
)
```

*Note: To call this function, set the parameter `gene_list_type = "Excel"` in the function `Celltype_Annotation_Features`. This function also works with `Markers_list` that contains either no metric information or metric information generated in other ways.*

## 6. Conclusion

Thank you for using SlimR. For questions, issues, or suggestions, please submit them in the issue section or discussion section on GitHub (suggested) or send an email (alternative):

zhaoqingwang\@mail.sdu.edu.cn

**Zhaoqing Wang**