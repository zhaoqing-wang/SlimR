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
    - [2.3 From ScType Database](#23-from-sctype-database)
    - [2.4 From Seurat Objects](#24-from-seurat-objects)
    - [2.5 From Excel Tables](#25-from-excel-tables)
    - [2.6 Built-in Markers Lists](#26-built-in-markers-lists)
3. [Automated Annotation Workflow](#3-automated-annotation-workflow)
    - [3.1 Calculate Parameter](#31-calculate-parameter)
    - [3.2 Cluster-Based Annotation](#32-cluster-based-annotation)
    - [3.3 Per-Cell Annotation](#33-per-cell-annotation)
4. [Semi-Automated Annotation Workflow](#4-semi-automated-annotation-workflow)
5. [Other Functions Provided](#5-other-functions-provided)
    - [5.1 Cell type mapping](#51-cell-type-mapping)
    - [5.2 Single-Gene AUC and ROC Analysis](#52-single-gene-auc-and-roc-analysis)
    - [5.3 Hierarchical Proportion Plot](#53-hierarchical-proportion-plot)
6. [Citation](#6-citation)
7. [License](#7-license)
8. [Contact](#8-contact)

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

### 2.3 From ScType Database

**Reference:** *Ianevski et al. (2022) [doi:10.1038/s41467-022-28803-w](https://doi.org/10.1038/s41467-022-28803-w)*

``` r
ScType <- SlimR::ScType

Markers_list_ScType <- Markers_filter_ScType(
  ScType,
  tissue_type = "Intestine",
  cell_name = NULL
)
```

**Important: Specify `tissue_type` for accurate annotations.**

<details>
<summary><b>Optional: Explore database metadata</b></summary>

``` r
ScType_table <- SlimR::ScType_table
View(ScType_table)
```

</details>

### 2.4 From Seurat Objects

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

### 2.5 From Excel Tables

**Format:** Each sheet name = cell type, first row = headers, first column = markers, subsequent columns = metrics (optional).

``` r
Markers_list_Excel <- Read_excel_markers("D:/Laboratory/Marker_load.xlsx")
```

*If your Excel file lacks column headers, set `has_colnames = FALSE`.*

### 2.6 Built-in Markers Lists

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
    AUC_correction = TRUE,
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

## 5. Other Functions Provided

### 5.1 Cell type mapping

Cross‑tabulate cell type labels from one Seurat object with a grouping column from another Seurat object. The function automatically aligns cell barcodes using multiple normalization strategies and returns count tables, column‑wise proportion tables, a dominant mapping, and a heatmap.

```r
result <- Celltype_Compare(
  sce_label = seurat_obj1,
  sce = seurat_obj2,
  label_col = "cell_type",
  group_col = "cluster"
)

# Access results
head(result$prop_table)   # column-wise proportions
print(result$plot)         # heatmap of proportions
result$main_to_sub         # dominant cell type per group
```

### 5.2 Single-Gene AUC and ROC Analysis

Quickly assess the discriminative power of a single gene for a user‑defined cell group. The function returns the AUC, ROC data for custom plotting, and an optional ggplot2 curve.

```r
result <- Compute_Gene_AUC_ROC(
  seurat_obj  = sce,
  gene        = "CD3D",
  group_col   = "Cell Types",
  group_label = "T cells",
  assay       = "RNA",
  method      = "rank",
  plot        = TRUE,
  line_color  = "navy",
  line_size   = 1
)

# Access results
result$AUC              # numeric AUC value
head(result$roc_data)   # data.frame with fpr and tpr
result$roc_plot         # ggplot object (when plot = TRUE)
```

<details>
<summary><b>Detailed parameter guide</b></summary>

- `method`: `"raw"` (raw expression, optionally truncated by `min_expression`) or `"rank"` (dropout‑robust rank‑based scores).  
- `min_expression`: when `method = "raw"`, values below this are set to zero.  
- `keep_expression_above`: optional threshold – keep only cells with expression above it. **Warning:** this shifts the AUC interpretation to “discrimination among expressing cells” and should be compared with the default all‑cell result.  
- `plot`, `plot_title`, `line_color`, `line_size`: control the ROC plot appearance.

</details>

### 5.3 Hierarchical Proportion Plot

Create a publication‑ready composite figure that visualises the hierarchical classification of single‑cell data from broad cell types down to fine sub‑types.  
The **upper panel** draws a layered tree diagram (bubble size ∝ cell count, parent‑child links shown as three‑segment step lines). The **lower panel** (optional) displays per‑group cell‑type proportions as a heatmap perfectly aligned with the terminal leaves.

```r
# Full three-level hierarchy with proportion heatmap
res <- Plot_Hierarchy_Proportion(
  seurat_obj        = sce,
  Main_cell_types   = "Main_type",
  Cell_types        = "Cell_type",
  Sub_cell_types    = "Sub_type",
  proportion        = TRUE,
  Groups            = "orig.ident"
)

# Access individual plot components
res$tree_plot        # ggplot object – tree including labels & short sticks
res$prop_plot        # ggplot object – proportion heatmap
res$combined_plot    # combined plot (requires patchwork)

# Custom colours, only the tree
res_tree <- Plot_Hierarchy_Proportion(
  seurat_obj          = sce,
  Main_cell_types     = "Main_type",
  col_Main_cell_types = c(Immune = "red", Stromal = "blue", Epithelial = "green"),
  Cell_types          = "Cell_type",
  proportion          = FALSE
)
```

<details>
<summary><b>Detailed parameter guide</b></summary>

- **Hierarchy levels**  
  `Main_cell_types`, `Cell_types`, `Sub_cell_types` are character strings naming columns in `seurat_obj@meta.data`.  
  Use `NULL` to omit a level. If `Sub_cell_types` is given, `Cell_types` must also be provided.  
  Category names (e.g. “T cell”) must be **unique within each level** (they can repeat across levels).

- **Partial sub‑clustering**  
  It is common that only a subset of cells receives a finer annotation (e.g., only T cells are split into subtypes). The function automatically handles this: a cell without a valid sub‑label becomes a leaf at the deepest level where it has a label. The proportion heatmap is then built from the **union of all terminal leaf labels** – so no population is lost.

- **Label placement & adaptive height**  
  Leaf labels are drawn directly below the terminal nodes inside the tree panel, rotated 90°, with short black sticks connecting nodes to labels. The tree panel’s lower limit automatically expands to accommodate the longest cell‑type name – no label is ever clipped, and the heatmap sits immediately beneath the labels.

- **Colour control**  
  `col_Main_cell_types`, `col_Cell_types`, `col_Sub_cell_types` accept named or unnamed colour vectors. When missing, the function tries to obtain a palette via `ArchR::paletteDiscrete()`. If **ArchR** is not installed, it falls back to the standard **ggplot2** hue palette and prints an informative message.

- **Proportion heatmap**  
  `proportion = TRUE` (default) adds a lower panel showing the fraction of each terminal cell type per group (column `Groups`).  
  `Groups` is required only when `proportion = TRUE`. The heatmap uses the same leaf order as the tree, has a tight black border, and uses a white‑to‑red colour gradient (customisable via `low_col` and `high_col`). Group labels are shown in **bold** on the y‑axis.

- **Non‑leaf annotations**  
  `show_labels = TRUE` (default) places italic text next to non‑leaf Main and Cell level nodes, helping identify broad categories at a glance.

- **Output**  
  The function returns a list with `tree_plot`, `prop_plot` (NULL if `proportion = FALSE`), and `combined_plot` (NULL unless **patchwork** is installed). All are **ggplot2** objects that can be further customised. The combined plot is automatically printed to the active graphics device.

</details>

## 6. Citation

```
Wang Z (2026). SlimR: Adaptive Machine Learning-Powered, Context-Matching Tool for Single-Cell and Spatial Transcriptomics Annotation.
https://github.com/zhaoqing-wang/SlimR
```

## 7. License

[MIT](https://github.com/zhaoqing-wang/SlimR/blob/main/LICENSE.md)

## 8. Contact

**Author:** Zhaoqing Wang ([ORCID](https://orcid.org/0000-0001-8348-7245)) | **Email:** <zhaoqingwang@mail.sdu.edu.cn> | **Issues:** [SlimR Issues](https://github.com/zhaoqing-wang/SlimR/issues)
