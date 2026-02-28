# SlimR: Adaptive Machine Learning-Powered, Context-Matching Tool for Single-Cell and Spatial Transcriptomics Annotation

[![CRAN Package Version](https://img.shields.io/cran/v/SlimR?label=CRAN)](https://cran.r-project.org/package=SlimR) [![CRAN License](https://img.shields.io/cran/l/SlimR?label=License&color=green)](https://cran.r-project.org/package=SlimR) [![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/SlimR)](https://cran.r-project.org/package=SlimR) [![GitHub Package Version](https://img.shields.io/github/r-package/v/zhaoqing-wang/SlimR?label=GitHub&color=green)](https://github.com/zhaoqing-wang/SlimR/releases) [![GitHub Maintainer](https://img.shields.io/badge/Maintainer-Zhaoqing_Wang-blue)](https://github.com/zhaoqing-wang)

## Overview

<img src="docs/Sticker.png" alt="Sticker" width="200"  align="right"/>

SlimR is an R package for cell-type annotation in single-cell and spatial transcriptomics. Existing marker-based annotation methods typically rely on manually tuned thresholds and operate at a single analytical granularity, limiting their adaptability across diverse datasets. SlimR addresses these challenges through three methodological contributions: **(1)** a **context-matching** framework that standardizes heterogeneous marker sources via multi-level biological filtering; **(2)** a **dataset-adaptive parameterization** strategy that infers optimal annotation hyperparameters from intrinsic data characteristics, eliminating manual calibration; and **(3)** a **dual-granularity scoring** architecture that provides both cluster-level probabilistic assignment and per-cell resolution with manifold-aware spatial smoothing for continuous cell states. A unified Feature Significance Score ensures biologically interpretable marker ranking throughout the workflow.

## Table of Contents

1. [Preparation](#1-preparation)
    - [1.1 Installation](#11-installation)
    - [1.2 Prepare Seurat Object](#12-prepare-seurat-object)
2. [Standardized Markers_list Input](#2-standardized-markers_list-input)
    - [2.1 From Cellmarker2 Database](#21-from-cellmarker2-database)
    - [2.2 From PanglaoDB Database](#22-from-panglaodb-database)
    - [2.3 From Seurat Objects](#23-from-seurat-objects)
    - [2.4 From Excel Tables](#24-from-excel-tables)
    - [2.5 Built-in Markers Lists](#25-built-in-markers-lists)
3. [Automated Annotation Workflow](#3-automated-annotation-workflow)
    - [3.1 Calculate Parameter](#31-calculate-parameter)
    - [3.2 Cluster-Based Annotation](#32-cluster-based-annotation)
    - [3.3 Per-Cell Annotation](#33-per-cell-annotation)
4. [Semi-Automated Annotation Workflow](#4-semi-automated-annotation-workflow)
5. [Citation](#5-citation)
6. [License](#6-license)
7. [Contact](#7-contact)

---

## 1. Preparation

### 1.1 Installation

**Option One: CRAN** [![CRAN Version](https://img.shields.io/cran/v/SlimR?label=CRAN)](https://cran.r-project.org/package=SlimR)

``` r
install.packages("SlimR")
```

**Option Two: GitHub** [![GitHub R package version](https://img.shields.io/github/r-package/v/zhaoqing-wang/SlimR?label=GitHub&color=green)](https://github.com/zhaoqing-wang/SlimR/releases)

``` r
devtools::install_github("zhaoqing-wang/SlimR")
```

<details>
<summary><b>Dependencies & optional packages</b></summary>

**Required:** R (≥ 3.5), cowplot, dplyr, ggplot2, patchwork, pheatmap, readxl, scales, Seurat, tidyr, tools

``` r
install.packages(c("cowplot", "dplyr", "ggplot2", "patchwork", 
                   "pheatmap", "readxl", "scales", "Seurat", 
                   "tidyr", "tools"))
```

**Optional:** RANN (10–100× faster UMAP spatial smoothing in per-cell annotation)

``` r
install.packages("RANN")
```

</details>

### 1.2 Prepare Seurat Object

``` r
library(SlimR)

# For Seurat objects with multiple layers, join layers first
sce@assays$RNA <- SeuratObject::JoinLayers(sce@assays$RNA)
```

**Important: Ensure your Seurat object has completed standard preprocessing (normalization, scaling, clustering) and batch effect correction.**

---

## 2. Standardized Markers_list Input

SlimR uses a standardized list format: list names = cell types, first column = marker genes, additional columns = metrics (optional).

### 2.1 From Cellmarker2 Database

**Reference:** *Hu et al. (2023) [doi:10.1093/nar/gkac947](https://doi.org/10.1093/nar/gkac947)*

``` r
Cellmarker2 <- SlimR::Cellmarker2

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

<details>
<summary><b>Optional: Explore database metadata</b></summary>

``` r
Cellmarker2_table <- SlimR::Cellmarker2_table
View(Cellmarker2_table)
```

</details>

### 2.2 From PanglaoDB Database

**Reference:** *Franzén et al. (2019) [doi:10.1093/database/baz046](https://doi.org/10.1093/database/baz046)*

``` r
PanglaoDB <- SlimR::PanglaoDB

Markers_list_panglaoDB <- Markers_filter_PanglaoDB(
  PanglaoDB,
  species_input = 'Human',
  organ_input = 'GI tract'
)
```

<details>
<summary><b>Optional: Explore database metadata</b></summary>

``` r
PanglaoDB_table <- SlimR::PanglaoDB_table
View(PanglaoDB_table)
```

</details>

### 2.3 From Seurat Objects

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

*Tip: `sort_by = "FSS"` ranks by Feature Significance Score (log2FC × Expression ratio). Use `sort_by = "avg_log2FC"` for fold-change ranking.*

<details>
<summary><b>Use presto for ~10× faster marker detection</b></summary>

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

*Install presto: `devtools::install_github('immunogenomics/presto')`*

</details>

### 2.4 From Excel Tables

**Format:** Each sheet name = cell type, first row = headers, first column = markers, subsequent columns = metrics (optional).

``` r
Markers_list_Excel <- Read_excel_markers("D:/Laboratory/Marker_load.xlsx")
```

*If your Excel file lacks column headers, set `has_colnames = FALSE`.*

### 2.5 Built-in Markers Lists

SlimR includes curated marker lists for specific annotation tasks:

| List | Scope | Reference |
|------|-------|-----------|
| `Markers_list_scIBD` | Human intestinal cells (IBD) | *Nie et al. (2023) [doi:10.1038/s43588-023-00464-9](https://doi.org/10.1038/s43588-023-00464-9)* |
| `Markers_list_TCellSI` | T cell subtypes | *Yang et al. (2024) [doi:10.1002/imt2.231](https://doi.org/10.1002/imt2.231)* |
| `Markers_list_PCTIT` | Pan-cancer T cell subtypes | *L. Zheng et al. (2021) [doi:10.1126/science.abe6474](https://doi.org/10.1126/science.abe6474)* |
| `Markers_list_PCTAM` | Pan-cancer macrophage subtypes | *Ruo-Yu Ma et al. (2022) [doi:10.1016/j.it.2022.04.008](https://doi.org/10.1016/j.it.2022.04.008)* |

``` r
# Example: Load built-in markers
Markers_list <- SlimR::Markers_list_scIBD
```

**Important: Ensure your input Seurat object matches the tissue/cell type scope of the selected marker list.**

---

## 3. Automated Annotation Workflow

SlimR provides two automated approaches: **Cluster-Based** (one label per cluster, fast) and **Per-Cell** (individual cell labels, finer resolution). Both share the same parameter calculation step and `Markers_list` format.

| Feature | Cluster-Based | Per-Cell |
|---------|---------------|----------|
| **Unit** | Cluster | Individual cell |
| **Speed** | ~10–30s (50k cells) | ~2–3min (50k cells) |
| **Resolution** | Coarse | Fine |
| **Best For** | Homogeneous clusters | Mixed clusters, rare cell types |
| **Spatial Context** | Not used | Optional (UMAP smoothing) |

### 3.1 Calculate Parameter

SlimR uses adaptive machine learning to determine optimal `min_expression`, `specificity_weight`, and `threshold` parameters. **This step is optional — skip to Section 3.2 to use defaults.**

``` r
SlimR_params <- Parameter_Calculate(
  seurat_obj = sce,
  features = c("CD3E", "CD4", "CD8A"),
  assay = "RNA",
  cluster_col = "seurat_clusters",
  verbose = TRUE
  )
```

<details>
<summary><b>Custom method: use markers from a specific cell type</b></summary>

``` r
SlimR_params <- Parameter_Calculate(
  seurat_obj = sce,
  features = unique(Markers_list_Cellmarker2$`B cell`$marker),
  assay = "RNA",
  cluster_col = "seurat_clusters",
  verbose = TRUE
  )
```

</details>

### 3.2 Cluster-Based Annotation

Three steps: **Calculate → Annotate → Verify**.

**Step 1: Calculate Cell Types**

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

*If you ran `Parameter_Calculate()`, use: `min_expression = SlimR_params$min_expression`, `specificity_weight = SlimR_params$specificity_weight`, `threshold = SlimR_params$threshold`.*

<details>
<summary><b>View results & correct predictions</b></summary>

``` r
# View heatmap, predictions, and ROC curves
print(SlimR_anno_result$Heatmap_plot)
View(SlimR_anno_result$Prediction_results)
print(SlimR_anno_result$AUC_plot)   # Requires plot_AUC = TRUE

# Manually correct predictions
SlimR_anno_result$Prediction_results$Predicted_cell_type[
  SlimR_anno_result$Prediction_results$cluster_col == 15
] <- "Intestinal stem cell"

# Label low-confidence predictions as Unknown
SlimR_anno_result$Prediction_results$Predicted_cell_type[
  SlimR_anno_result$Prediction_results$AUC <= 0.5
] <- "Unknown"
```

*When correcting, preferably use cell types from the `Alternative_cell_types` column.*

</details>

**Step 2: Annotate Cell Types**

``` r
sce <- Celltype_Annotation(seurat_obj = sce,
    cluster_col = "seurat_clusters",
    SlimR_anno_result = SlimR_anno_result,
    plot_UMAP = TRUE,
    annotation_col = "Cell_type_SlimR"
    )
```

**Step 3: Verify Cell Types**

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

**Important: Use matching `cluster_col` and `annotation_col` values across all three functions.**

### 3.3 Per-Cell Annotation

Three steps: **Calculate → Annotate → Verify**. Ideal for heterogeneous clusters, rare cell types, and continuous differentiation states.

**Step 1: Calculate Per-Cell Types**

``` r
SlimR_percell_result <- Celltype_Calculate_PerCell(
    seurat_obj = sce,
    gene_list = Markers_list,
    species = "Human",
    assay = "RNA",
    method = "weighted",
    min_expression = 0.1,
    use_umap_smoothing = FALSE,
    min_score = "auto",
    min_confidence = 1.2,
    verbose = TRUE
    )
```

*Three scoring methods: `"weighted"` (default, recommended), `"mean"` (fast baseline), `"AUCell"` (rank-based, robust to batch effects).*

<details>
<summary><b>UMAP spatial smoothing & parameter tuning</b></summary>

``` r
# Enable UMAP smoothing for noise reduction
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

*Install RANN for 10–100× faster k-NN: `install.packages("RANN")`*

| Scenario | `min_score` | `min_confidence` |
|----------|-------------|------------------|
| Few cell types (<15) | `"auto"` | 1.2 (default) |
| Many cell types (>30) | `"auto"` | 1.1–1.15 |
| Strict annotation | `"auto"` | 1.3–1.5 |
| Liberal annotation | `"auto"` | 1.0 (disable) |

</details>

**Step 2: Annotate Per-Cell Types**

``` r
sce <- Celltype_Annotation_PerCell(
    seurat_obj = sce,
    SlimR_percell_result = SlimR_percell_result,
    plot_UMAP = TRUE,
    annotation_col = "Cell_type_PerCell_SlimR",
    plot_confidence = TRUE
    )
```

**Step 3: Verify Per-Cell Types**

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

**Important: Use matching `annotation_col` values in `Celltype_Annotation_PerCell()` and `Celltype_Verification_PerCell()`.**

---

## 4. Semi-Automated Annotation Workflow

For expert-guided manual annotation using visualizations:

<details>
<summary><b>4.1 Annotation Heat Map</b></summary>

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

*Note: This function is now incorporated into `Celltype_Calculate()`. Use `Celltype_Calculate()` instead for automated workflows.*

</details>

<details>
<summary><b>4.2 Annotation Feature Plots</b></summary>

Generates per-cell-type expression dot plot with metric heat map:

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

*Set `gene_list_type` to `"Cellmarker2"`, `"PanglaoDB"`, `"Seurat"`, or `"Excel"` to match your marker source.*

</details>

<details>
<summary><b>4.3 Annotation Combined Plots</b></summary>

Generates per-cell-type box plots of marker expression levels:

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

</details>

---

## 5. Citation

```
Wang Z (2026). SlimR: Adaptive Machine Learning-Powered, Context-Matching Tool for Single-Cell and Spatial Transcriptomics Annotation.
https://github.com/zhaoqing-wang/SlimR
```

## 6. License

[MIT](LICENSE)

## 7. Contact

**Author:** Zhaoqing Wang ([ORCID](https://orcid.org/0000-0001-8348-7245)) | **Email:** <zhaoqingwang@mail.sdu.edu.cn> | **Issues:** [SlimR Issues](https://github.com/zhaoqing-wang/SlimR/issues)
