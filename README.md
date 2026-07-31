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
    - [2.4 From CellTypist Organ Atlas](#24-from-celltypist-organ-atlas)
    - [2.5 From Seurat Objects](#25-from-seurat-objects)
    - [2.6 From Scanpy (Python) Objects](#26-from-scanpy-python-objects)
    - [2.7 From Excel Tables](#27-from-excel-tables)
    - [2.8 Built-in Markers Lists](#28-built-in-markers-lists)
3. [Automated Annotation Workflow](#3-automated-annotation-workflow)
    - [3.1 Calculate Parameter](#31-calculate-parameter)
    - [3.2 Cluster-Based Annotation](#32-cluster-based-annotation)
    - [3.3 Per-Cell Annotation](#33-per-cell-annotation)
4. [Semi-Automated Annotation Workflow](#4-semi-automated-annotation-workflow)
    - [4.1 Annotation Heat Map](#41-annotation-heat-map)
    - [4.2 Annotation Feature Plots](#42-annotation-feature-plots)
    - [4.3 Annotation Combined Plots](#43-annotation-combined-plots)
5. [Other Functions Provided](#5-other-functions-provided)
    - [5.1 Cell type mapping](#51-cell-type-mapping)
    - [5.2 Single-Gene AUC and ROC Analysis](#52-single-gene-auc-and-roc-analysis)
    - [5.3 Hierarchical Proportion Plot](#53-hierarchical-proportion-plot)
    - [5.4 Weighted Voronoi Plot](#54-weighted-voronoi-plot)
    - [5.5 Built‑in Colour Palettes](#55-builtin-colour-palettes)
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

### 2.4 From CellTypist Organ Atlas

**Reference:**  
*Xu et al. (2023)* [doi:10.1016/j.cell.2023.11.026](https://doi.org/10.1016/j.cell.2023.11.026)  
*Domínguez Conde et al. (2022)* [doi:10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197)

SlimR provides a pre‑computed marker list derived from the [CellTypist organ atlas](https://www.celltypist.org/organs).  
It covers **12 human organs** (Blood, Bone_marrow, Heart, Hippocampus, Intestine, Kidney, Liver, Lung, Lymph_node, Pancreas, Skeletal_muscle, Spleen) and **399 cell types**, with markers obtained via the Scanpy workflow (log1p‑normalised data, Wilcoxon test, adjusted p‑value < 0.01, log2 fold‑change > 0, then ranked by log fold‑change; **top 100 genes** per cell type). The data have been imported using `Read_excel_markers` and are directly usable.

```r
# Load the built-in list
CellTypist <- SlimR::CellTypist

# Access markers for one organ (e.g., Intestine)
Markers_list_CellTypist <- CellTypist$Intestine

# Each organ contains a named list of data frames (one per cell type)
names(Markers_list_CellTypist)

# The data frames are pre‑sorted by log fold‑change (descending).
# To restrict to the top 20 markers for every cell type in this organ:
Markers_list_CellTypist_top20 <- lapply(Markers_list_CellTypist, function(df) head(df, 20))
```

**Key points:**
- Use `$organ_name` to extract an organ; the organ names are exactly as shown above (case‑sensitive).
- Each cell‑type data frame is already ranked by `logfoldchanges` (descending) – simply use `head(df, n)` to obtain the top *n* markers.
- The full list can be passed directly to SlimR’s annotation functions as a standard `Markers_list` object.

### 2.5 From Seurat Objects

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

**Important: To avoid long running time, for data with more than 100,000 cells, it is recommended to use scanpy for DEGs calculation (Section 2.6).**

### 2.6 From Scanpy (Python) Objects

Differential expression results from a Scanpy AnnData object can be exported to an Excel file and then loaded directly into SlimR’s standard format using `Read_excel_markers`.

<details>
<summary><b>Process Codes</b></summary>

```python
import scanpy as sc
import pandas as pd
import numpy as np
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
import re

# Load data
adata = sc.read_h5ad("adata.h5ad")

# ------------------------------------------------------------
# Ensure expression data is log1p‑normalised.
# If adata.X contains raw counts, normalise and log1p now:
#   sc.pp.normalize_total(adata, target_sum=1e4)
#   sc.pp.log1p(adata)
#
# If raw counts are in a layer (e.g., 'counts'), move them to .X first:
#   adata.X = adata.layers['counts'].copy()
#   sc.pp.normalize_total(adata, target_sum=1e4)
#   sc.pp.log1p(adata)
#
# If .X already contains log1p data, you can skip the step above.
# ------------------------------------------------------------

# Cluster column (adjust to your metadata column name)
cluster_key = "Curated_annotation"
adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")
clusters = adata.obs[cluster_key].cat.categories

# Wilcoxon test (one‑vs‑rest)
sc.tl.rank_genes_groups(adata, groupby=cluster_key,
                        method="wilcoxon", n_jobs=-1)

# Collect filtered results per cluster
de_dict = {}
for clust in clusters:
    df = sc.get.rank_genes_groups_df(adata, group=clust)
    df = df[(df["pvals_adj"] < 0.01) & (df["logfoldchanges"] > 0)]
    df = df.sort_values("logfoldchanges", ascending=False).head(100)
    df = df.rename(columns={"names": "gene"})
    # Round numeric columns for cleaner output
    for col in df.select_dtypes(include=[np.number]).columns:
        df[col] = df[col].round(4)
    de_dict[clust] = df

# Write to Excel (one sheet per cluster)
def sanitize_sheet_name(name):
    return re.sub(r'[\[\]:*?/\\]', '_', str(name))[:31]

wb = Workbook()
wb.remove(wb.active)
for clust in clusters:
    ws = wb.create_sheet(title=sanitize_sheet_name(clust))
    for row in dataframe_to_rows(de_dict[clust], index=False, header=True):
        ws.append(row)
wb.save("DEGs.xlsx")
```

</details>

**Important:**  
- Differential expression must be computed on **log1p‑normalised** data. If your `.X` still holds raw counts, normalise (e.g., `normalize_total` + `log1p`) before calling `rank_genes_groups`.  
- Adapt `groupby` to your actual annotation column (e.g., `"Cell_type"`, `"leiden"`).  
- You can adjust the significance threshold (`pvals_adj`), fold‑change direction, and number of genes (`head(100)`) to suit your analysis.

After saving the `DEGs.xlsx` file, use the `Read_excel_markers` function from Section 2.7 to import it into R.


### 2.7 From Excel Tables

**Format:** Each sheet name = cell type, first row = headers, first column = markers, subsequent columns = metrics (optional).

``` r
Markers_list_Excel <- Read_excel_markers("D:/Laboratory/Marker_load.xlsx")
```

*If your Excel file lacks column headers, set `has_colnames = FALSE`.*

### 2.8 Built-in Markers Lists

SlimR includes curated marker lists for specific annotation tasks:

| List | Scope | Reference |
|------|-------|-----------|
| `Markers_list_scIBD` | Human intestinal cells (IBD) | *Nie et al. (2023) [doi:10.1038/s43588-023-00464-9](https://doi.org/10.1038/s43588-023-00464-9)* |
| `Markers_list_TCellSI` | T cell subtypes | *Yang et al. (2024) [doi:10.1002/imt2.231](https://doi.org/10.1002/imt2.231)* |
| `Markers_list_PCTIT` | Pan-cancer T cell subtypes | *L. Zheng et al. (2021) [doi:10.1126/science.abe6474](https://doi.org/10.1126/science.abe6474)* |
| `Markers_list_PCTAM` | Pan-cancer macrophage subtypes | *Ruo-Yu Ma et al. (2022) [doi:10.1016/j.it.2022.04.008](https://doi.org/10.1016/j.it.2022.04.008)* |

``` r
# Example: Load built-in markers
Markers_list_scIBD <- SlimR::Markers_list_scIBD

# The data frames are pre‑sorted by log fold‑change (descending).
# To restrict to the top 20 markers for every cell type in this organ:
Markers_list_scIBD_top20 <- lapply(Markers_list_scIBD, function(df) head(df, 20))
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

<details>
<summary><b>Parameter descriptions</b></summary>

- **`seurat_obj`**: Seurat object containing annotation columns (e.g., `seurat_cluster`) in `meta.data`.
- **`gene_list`**: A named list of markers, where each element is a data frame with marker genes in the first column. Can be generated by `Markers_filter_Cellmarker2()`, `Markers_filter_PanglaoDB()`, `read_excel_markers()`, or `read_seurat_markers()`.
- **`species`**: `"Human"` or `"Mouse"` – used for standardising gene symbols in the marker list.
- **`cluster_col`**: Column name in `meta.data` that defines clusters (default: `"seurat_clusters"`).
- **`assay`**: Assay to use (default: `"RNA"`).
- **`min_expression`**: Threshold for considering a gene “expressed” in a cell; low‑expression cells are filtered to reduce noise (default: `0.1`).
- **`specificity_weight`**: Controls how much expression variability (standard deviation) within a cluster contributes to the specificity score; higher values amplify variability (default: `3`).
- **`threshold`**: Normalised similarity threshold between the alternative and predicted cell types; used for filtering uncertain assignments (default: `0.6`).
- **`compute_AUC`**: If `TRUE`, calculates AUC values for each predicted cell type to measure marker discriminative power (default: `TRUE`).
- **`plot_AUC`**: If `TRUE`, generates an ROC curve plot for the predicted cell types (default: `TRUE`).
- **`AUC_correction`**: If `TRUE`, uses the highest‑AUC cell type among candidates (probability > threshold) as the final prediction, and records its AUC in the `AUC` column (default: `FALSE`).
- **`colour_low`**: Colour for the lowest probability in the heatmap (default: `"navy"`).
- **`colour_high`**: Colour for the highest probability in the heatmap (default: `"firebrick3"`).

</details>

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

*If you ran `Parameter_Calculate()`, use: `min_expression = SlimR_params$min_expression`, `specificity_weight = SlimR_params$specificity_weight`, `threshold = SlimR_params$threshold`.*

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

**Please note: When performing cell-by-cell annotation, the annotation results based on cell resolution are subject to instability.**

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

### 4.1 Annotation Heat Map

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

### 4.2 Annotation Feature Plots

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

### 4.3 Annotation Combined Plots

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
# Full three-level hierarchy with proportion heatmap (default: row‑wise proportions)
res <- Plot_Hierarchy_Proportion(
  seurat_obj        = sce,
  Main_cell_types   = "Main_type",
  Cell_types        = "Cell_type",
  Sub_cell_types    = "Sub_type",
  proportion        = TRUE,
  Groups            = "orig.ident",
  low_col           = "white",
  high_col          = "navy"
)

# When plotting sub‑types of a larger population (e.g., immune subsets)
# where total group sizes differ, use adjust_by_group = TRUE
res <- Plot_Hierarchy_Proportion(
  seurat_obj        = sce,
  Main_cell_types   = "Immune_Main_type",
  Cell_types        = "Immune_Cell_type",
  Sub_cell_types    = "Immune_Sub_type",
  proportion        = TRUE,
  Groups            = "condition",
  adjust_by_group   = TRUE
)

# Access individual plot components
res$tree_plot        # ggplot object – tree including labels & short sticks
res$prop_plot        # ggplot object – proportion heatmap
res$combined_plot    # combined plot (requires patchwork)
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
  `col_Main_cell_types`, `col_Cell_types`, `col_Sub_cell_types` accept named or unnamed colour vectors. When missing, the function generates a palette using the internal `paletteDiscrete()` function, which replicates the *stallion* palette from the **ArchR** package.  No external ArchR installation is required.

- **Proportion heatmap**  
  `proportion = TRUE` (default) adds a lower panel showing the fraction of each terminal cell type per group (column `Groups`).  
  `Groups` is required only when `proportion = TRUE`. The heatmap uses the same leaf order as the tree, has a tight black border, and uses a white‑to‑red colour gradient (customisable via `low_col` and `high_col`). Group labels are shown in **bold** on the y‑axis, and, if available, the number of cells in each group is appended (e.g. “Control (1254)”).  

  - **Adjustment for unequal group sizes (`adjust_by_group`)**  
    When analysing sub‑types derived from a larger population (e.g., immune subsets) and the total number of cells in each group differs substantially, set `adjust_by_group = TRUE`. This option multiplies the row‑wise proportion by the ratio of the group’s cell count to the mean group cell count. The resulting heatmap then reflects both within‑group composition and between‑group abundance differences. The colour scale is automatically normalised across all cells and groups. For broad cell type visualisation or when group sizes are balanced, keep the default `FALSE`.

- **Non‑leaf annotations**  
  `show_labels = TRUE` (default) places italic text next to non‑leaf Main and Cell level nodes, helping identify broad categories at a glance.

- **Output**  
  The function returns a list with `tree_plot`, `prop_plot` (NULL if `proportion = FALSE`), and `combined_plot` (NULL unless **patchwork** is installed). All are **ggplot2** objects that can be further customised. The combined plot is automatically printed to the active graphics device.

</details>

### 5.4 Weighted Voronoi Plot

Generate a weighted Voronoi treemap that visualizes the hierarchical composition of single‑cell data. Polygons are grouped by the main cell type, and the area of each sub‑type polygon is proportional to its cell count. Colours follow the same palette logic as other SlimR functions, derived from ArchR but fully built into the package.

The plot is drawn using a custom `ggplot2`‑based renderer to ensure exact colour matching with `Plot_Hierarchy_Proportion` and `DimPlot`, bypassing the limited colour handling of the upstream `WeightedTreemaps` package.

```r
# Basic treemap with rounded rectangles, displaying both count and percentage
res <- Plot_Voronoi_diagram(
  seurat_obj      = sce,
  Main_cell_types = "Main_type",
  Cell_types      = "Cell_type",
  label_type      = "both",
  shape           = "rounded_rect",
  seed            = 1
)

# Access the underlying treemap object or the final ggplot
res$voronoi_treemap   # the treemap object from WeightedTreemaps
res$plot              # the ggplot object
```

<details>
<summary><b>Detailed parameter guide</b></summary>

- **Data & hierarchy**  
  `Main_cell_types` and `Cell_types` are column names in `seurat_obj@meta.data` defining the two‑level hierarchy.  
  Only cells with valid (non‑missing, non‑empty) labels in both columns are used.  
  The voronoi diagram groups cells by main type (level 1) and further splits each main type into sub‑type polygons (level 2).

- **Polygon labels (`label_type`)**  
  Controls the text displayed inside each sub‑type polygon:
  - `"both"` (default) – shows the sub‑type name, cell count, and percentage of total cells (each on a new line).
  - `"count"` – shows sub‑type name and cell count.
  - `"percentage"` – shows sub‑type name and percentage of total cells.
  - `"none"` – shows only the sub‑type name.

- **Polygon shape (`shape`)**  
  `"rounded_rect"` (default) produces rounded rectangles; `"circle"` yields circular polygons.  
  The layout is non‑deterministic but can be made reproducible via the `seed` parameter.

- **Reproducibility (`seed`)**  
  A single integer passed to the Voronoi layout algorithm. The same seed yields the same polygon arrangement across runs.

- **Colour control (`col_Cell_types`)**  
  Accepts a named or unnamed character vector of colours for the `Cell_types` categories. If `NULL`, colours are automatically generated via the internal `paletteDiscrete()` function (replicates the ArchR *stallion* palette). No external ArchR installation is needed.

- **Label appearance**  
  `label_size` controls the text size inside polygons (default `3`).  
  `label_color` sets the text colour (default `"black"`).  
  `label_fontface` controls the font face (`"plain"`, `"italic"`, `"bold"`; default `"bold"`).

- **Borders and frames**  
  The function draws three types of borders:
  - **Sub‑type borders** (`subtype_border_lwd`, default `0.15`) – thin lines between individual sub‑type polygons.
  - **Main borders** (`main_border_lwd`, default `0.35`) – thicker lines between main cell type regions, drawn on top of sub‑type borders for clear separation.
  - **Outer frame** (`outer_border_lwd`, default `0.4`) – a convex hull tightly surrounding the entire plot, following the natural outline of the treemap.
  All borders share the same `border_color` (default `"grey90"`, a very light grey). This parameter can be customised to any valid R colour.

- **Legend**  
  `legend = TRUE` (default) shows a colour legend; `legend_position` controls its placement (default `"right"`, also accepts `"left"`, `"bottom"`, `"top"`, or `"none"`).

- **Output**  
  The function invisibly returns a list with two components:
  - `voronoi_treemap`: the raw `voronoiTreemap` object from `WeightedTreemaps`, containing polygon coordinates and metadata.
  - `plot`: the final `ggplot` object, produced by a custom drawing routine that extracts polygon vertices and applies colours via `scale_fill_manual()`.
  The plot is automatically printed to the active graphics device.

- **Dependencies**  
  Requires the `WeightedTreemaps` package, available from GitHub. If not installed, an error is thrown with installation instructions.  
  The function only uses `WeightedTreemaps::voronoiTreemap()` to compute polygon layouts; all rendering is done with `ggplot2`.

</details>

### 5.5 Built‑in Colour Palettes

Since *ArchR* is not available on CRAN, SlimR incorporates its colour palettes directly (via the internal function `paletteDiscrete()`) so that users can enjoy the same publication‑quality colours without any additional installation. The palettes, including the default stallion, are hard‑coded in the package and require no external dependencies.

You can call the palette generator directly:
```r
# Display "orig.ident" using the built-in palette
col.clr <- SlimR::paletteDiscrete(values = c(names(table(sce$orig.ident))))
DimPlot(sce,
  reduction = "umap",
  group.by = "orig.ident",
  cols = col.clr,label = TRUE) + NoAxes()

# Display "cell_type" using the built-in palette
col.clr <- SlimR::paletteDiscrete(levels(sce$cell_type))
DimPlot(sce,
  reduction = "umap",
  group.by = "cell_type",
  cols = col.clr,label = TRUE) + NoAxes()
```

The function returns a named vector of hex colours, arranged horizontally according to the input vector. When the number of categories exceeds the palette size, colours are interpolated smoothly.

All SlimR plotting functions that accept `col_...` parameters automatically use this palette when no custom colours are supplied, ensuring a consistent and publication‑ready colour scheme across different types of plots.

<details>
<summary><b>Attribution</b></summary>

The palettes are derived from ArchR, a scalable software package for integrative single‑cell chromatin accessibility analysis:

- Granja JM, Corces MR et al. (2021) **ArchR is a scalable software package for integrative single‑cell chromatin accessibility analysis**. *Nature Genetics* **53**, 403–411. doi:10.1038/s41588-021-00790-6
- Project website: [https://www.archrproject.com/](https://www.archrproject.com/)
- GitHub repository: [https://github.com/GreenleafLab/ArchR](https://github.com/GreenleafLab/ArchR)

**License**  
ArchR is distributed under the **MIT License**. SlimR respects the original license by including the palette data directly and documenting its provenance.

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
