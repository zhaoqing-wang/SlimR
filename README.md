# SlimR: Adaptive Machine Learning-Powered, Context-Matching Tool for Single-Cell and Spatial Transcriptomics Annotation

[![CRAN Package Version](https://img.shields.io/cran/v/SlimR?label=CRAN)](https://cran.r-project.org/package=SlimR) [![CRAN License](https://img.shields.io/cran/l/SlimR?label=License&color=green)](https://cran.r-project.org/package=SlimR) [![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/SlimR)](https://cran.r-project.org/package=SlimR) [![GitHub Package Version](https://img.shields.io/github/r-package/v/zhaoqing-wang/SlimR?label=GitHub&color=green)](https://github.com/zhaoqing-wang/SlimR/releases) [![GitHub Maintainer](https://img.shields.io/badge/Maintainer-Zhaoqing_Wang-blue)](https://github.com/zhaoqing-wang)

## Overview

<img src="docs/Sticker.png" alt="Sticker" width="233.28" height="270" align="right"/>

SlimR is an R package designed for annotating single-cell and spatial-transcriptomics (ST) datasets. It supports the creation of a unified marker list, `Markers_list`, using sources including: the package's built-in curated species-specific cell type and marker reference databases (e.g., 'Cellmarker2', 'PanglaoDB', 'scIBD', 'TCellSI','PCTIT','PCTAM'), Seurat objects containing cell label information, or user-provided Excel tables mapping cell types to markers.

SlimR can predict calculation parameters by adaptive machine learning algorithms by `Parameter_Calculate()`, and based on `Markers_list`, calculate gene expression of different cell types and predict annotation information and calculate corresponding AUC by `Celltype_Calculate()`, and annotate it by `Celltype_Annotation()`, then verify it by `Celltype_Verification()`. At the same time, it can calculate gene expression corresponding to the cell type to generate a reference map for manual annotation (e.g., 'Heat Map', 'Feature Plots', 'Combined Plots').

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
    -   [3.2 Calculate Cell Types](#32-calculate-cell-types)
    -   [3.3 Annotate Cell Types](#33-annotate-cell-types)
    -   [3.4 Verify Cell Types](#34-verify-cell-types)
4.  [Semi-Automated Annotation Workflow](#4-semi-automated-annotation-workflow)
    -   [4.1 Annotation Heat Map](#41-annotation-heat-map)
    -   [4.2 Annotation Feature Plots](#42-annotation-feature-plots)
    -   [4.3 Annotation Combined Plots](#43-annotation-combined-plots)
5.  [Other Functions Provided by SlimR](#5-other-functions-provided-by-slimr)
6.  [Conclusion](#6-conclusion)

------------------------------------------------------------------------

## 1. Preparation

### 1.1 Installation

Option One: [![CRAN Version](https://img.shields.io/cran/v/SlimR?label=CRAN)](https://cran.r-project.org/package=SlimR)

Install SlimR directly from CRAN using: (Stable version, recommended when the version is equivalent to the GitHub package version)

``` r
install.packages("SlimR")
```

*Note: Try adjusting the CRAN image to `Global (CDN)` or use `BiocManager::install("SlimR")` if you encounter a version mismatch during installation.*

Option Two: [![GitHub R package version](https://img.shields.io/github/r-package/v/zhaoqing-wang/SlimR?label=GitHub&color=green)](https://github.com/zhaoqing-wang/SlimR/releases)

Install SlimR directly from GitHub using: (Development version, recommended when the version is higher than the CRAN package version)

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

For Seurat objects with multiple layers in the assay, please run `SeuratObject::JoinLayers()` first.

``` r
# For example, if you want to use the 'RNA' layer in the multilayered Seurat object assay.
sce@assays$RNA <- SeuratObject::JoinLayers(sce@assays$RNA)
```

**Important: To ensure accuracy of the annotation, make sure that the entered Seurat object has run the standard process and removed batch effects.**

*Note: It is recommended to use the `clustree` package to determine the appropriate resolution for the input Seurat object.*

### 1.4 Dependencies (Alternative)

SlimR requires R (≥ 3.5) and depends on the following packages: `cowplot`, `dplyr`, `ggplot2`, `patchwork`, `pheatmap`, `readxl`, `scales`, `Seurat`, `tidyr`, `tools`. If installation fails, please install missing dependencies using:

``` r
# Install dependencies if needed:
install.packages(c("cowplot", "dplyr", "ggplot2", "patchwork", 
                   "pheatmap", "readxl", "scales", "Seurat", 
                   "tidyr", "tools"))
```

## 2. Standardized Markers_list Input

SlimR requires a standardized list format for storing marker information, metrics (can be omitted), and corresponding cell types (list names = cell types (essential), first column = markers (essential), subsequent columns = metrics (can be omitted)).

### 2.1 From Cellmarker2 Database

Cellmarkers2: A database of cell types and markers covering different species and tissue types.

Reference: *Hu et al. (2023) <doi:10.1093/nar/gkac947>*.

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

**Important: Select at least the `species` and `tissue_class` parameters to ensure the accuracy of the annotation.**

*Link: Output `Markers_list` usable in sections 3.1, 4.1, 4.2, 4.3, and 5.1. [Click to section 3 automated annotation workflow.](#3-automated-annotation-workflow)*

### 2.2 From PanglaoDB Database

PanglaoDB: Database of cell types and markers covering different species and tissue types.

Reference: *Franzén et al. (2019) <doi:10.1093/database/baz046>*.

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

#### 2.3.1 Identify Markers and Generate `Markers_list`:

The standard `Markers_list` can be generated by the built-in `read_seurat_markers()` function after obtaining Markers through the `Seurat::FindAllMarkers()` function.

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

*Note: Recommend using the parameter `sort_by = "FSS"` to use the 'Feature Significance Score' (FSS, product value of `log2FC` and `Expression ratio`) or use the parameter `sort_by = "avg_log2FC"` as the ranking basis.*

#### 2.3.2 Use `presto` to Speed Up: (Alternative)

For large data sets, the `presto::wilcoxauc()` function can be used to speed up the operation. (Alternative, \~10x faster, sacrifice partial accuracy)

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

**Important: This feature depends on the `presto` package. Please first run 'devtools::install_github('immunogenomics/presto') 'if prompted to install the package.**

*Note: Recommend using the parameter `sort_by = "logFC"` or using the parameter `sort_by = "FSS"` to use the 'Feature Significance Score' (FSS, product value of `log2FC` and `Expression ratio`) as the ranking basis.*

*Link: Output `Markers_list` usable in sections 3.1, 4.1, 4.2, 4.3, and 5.3. [Click to section 3 automated annotation workflow.](#3-automated-annotation-workflow)*

### 2.4 From Excel Tables

**Format Requirements**:

-   Each sheet name = cell type (essential)

-   First row = column headers (essential)

-   First column = markers (essential)

-   Subsequent columns = metrics (can be omitted)

``` r
Markers_list_Excel <- Read_excel_markers("D:/Laboratory/Marker_load.xlsx")
```

**Important: If the "First row" is missing in the input Excel file, please set the parameter 'has_colnames=FALSE' in the function 'Read_excel_markers()'**

*Link: Output `Markers_list` usable in sections 3.1, 4.1, 4.2, 4.3, and 5.4. [Click to section 3 automated annotation workflow.](#3-automated-annotation-workflow)*

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

### 3.1 Calculate Parameter

SlimR integrates adaptive machine learning algorithms to automatically determine the optimal `min_expression` and `specificity_weight` parameters in Section 3.2 for calculating the probability of cell types.

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

**Important: This scheme is optional and can be skipped to section 3.2 for cell type probability calculation using default parameters.**

### 3.2 Calculate Cell Types

#### 3.2.1 Calculate Cell Types (Core)

Uses `markers_list` to calculate probability, prediction results, calculate corresponding AUC (optional), and generate heat map and ROC graphs (optional) for cell annotation.

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

You can use the `min_expression = SlimR_params$min_expression` and `specificity_weight = SlimR_params$specificity_weight` and `threshold = SlimR_params$threshold` parameters in the function `Celltype_Calculate()` if you have run the `Parameter_Calculate ()` function in section 3.1 above.

**Important: The parameter `cluster_col` in the function `Celltype_Calculate()` and the function `Celltype_Annotation()` must be strictly the same to avoid false matches.**

*Note: Using the parameter `AUC_correction = TRUE` takes a little longer to compute (\~20% longer than only setting parameter `plot_AUC = TRUE`; \~40% longer than only setting parameter `compute_AUC = TRUE`), but it is recommended to correct the predicted cell type this way to obtain more accurate cell type prediction results. The lower the parameter `threshold`, the more alternative cell types AUC will check, and the longer the run time will be.*

Error handling: If you encounter the error message `Error in .rowNamesDF<-: ! duplicate 'row.names' are not allowed` when running `Celltype_Calculate()`, please run `base::make.unique()` first. (Alternative)

``` r
# If you encounter the error message `Error in .rowNamesDF<-: ! duplicate 'row.names' are not allowed` when running `Celltype_Calculate()`.
rownames(sce) <- base::make.unique(rownames(sce))
```

#### 3.2.2 Plot Heat Map (Optional)

Check the annotation probability of the cell type to be annotated in the input `cluster_col` column and cell types in `Markers_list` with the following code.

``` r
print(SlimR_anno_result$Heatmap_plot)
```

*Note: If the heat map is not generated properly, please run the function `library(pheatmap)` first.*

#### 3.2.3 View Prediction Results (Optional)

Cell type information results predicted by SlimR can be viewed with the following code.

``` r
View(SlimR_anno_result$Prediction_results)
```

#### 3.2.4 Plot ROC Curve and AUC Value (Optional)

Furthermore, the ROC curve and AUC value of the corresponding `cluster_col` and predicted cell types can be viewed by the following code.

``` r
print(SlimR_anno_result$AUC_plot)
```

**Important: This feature depends on the parameter `plot_AUC = TRUE`.**

*Note: If the heat map is not generated properly, please run the function `library(ggplot2)` first.*

#### 3.2.5 Correction for Predicted Cell Types (Alternative)

After viewing the list of predicted cell types and the corresponding AUC values, the predicted cell types can be corrected with the following code.

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

After modifying the corresponding predicted cell type, the following code is used to view the updated table of predicted cell types.

``` r
View(SlimR_anno_result$Prediction_results)
```

**Important: It is strongly recommended that if you need to correct the cell type, use cell types in `SlimR_anno_result$Prediction_results$Alternative_cell_type`.**

### 3.3 Annotate Cell Types

Assigns SlimR predicted cell types information in `SlimR_anno_result$Prediction_results$Predicted_cell_type` to the Seurat object based on cluster annotations, and stores the results into `seurat_obj@meta.data$annotation_col`.

``` r
sce <- Celltype_Annotation(seurat_obj = sce,
    cluster_col = "seurat_clusters",
    SlimR_anno_result = SlimR_anno_result,
    plot_UMAP = TRUE,
    annotation_col = "Cell_type_SlimR"
    )
```

**Important: The parameter `cluster_col` in the function `Celltype_Calculate()` and the function `Celltype_Annotation()` must be strictly the same to avoid false matches. The parameter `annotation_col` in the functions `Celltype_Annotation()` and `Celltype_Verification()` must be strictly the same to avoid false matches.**

### 3.4 Verify Cell Types

Use the cell group identity information in `seurat_obj@meta.data$annotation_col` and use the 'Feature Significance Score' (FSS, product value of `log2FC` and `Expression ratio`) as the ranking basis.

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

**Important: The parameter `annotation_col` in the function `Celltype_Annotation()` and the function `Celltype_Verification()` must be strictly the same to avoid false matches.**

*Note: Cell types located in `SlimR_anno_result$Prediction_results` were verified using the markers information from `SlimR_anno_result$Expression_list`; cell types that are not in the above list are validated using the markers information from the function `FindMarkers()`.*

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