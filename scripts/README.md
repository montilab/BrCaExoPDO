## Overview
This folder contains all the scripts used for the analysis and processing tasks in this project. The scripts are organized based on their functionality, which includes custom functions, data preprocessing, data integration, and other essential computational tasks.

## Structure

The `scripts` folder is structured as follows:

- **`utils.R`**: Scripts for required packages, color vectors, and custom functions.
- **`preprocess.R`**: Scripts for cleaning and preprocessing the raw data.
- **`annotation.R`**: Scripts for annotating cells from reference data.
- **`objectQC.R`**: Scripts for quality control of processed data.
- **`integration.R`**: Scripts for integration of data by sample.
- **`finetuning`**: Scripts for validating and refining annotations.
- **`unintegrated_analyses.R`**: Scripts for performing analyses done using unintegrated data.
- **`integrated_analyses.R`**: Scripts for performing analyses done using integrated data.
- **`figures.R`**: Scripts for generating figures seen in manuscript.

## Publicly Available Data and Resources

The scripts in this repository utilize the following publicly available data and tools:

### 1. SCSubtype
- **SCSubtype**: [Nature Genetics Article](https://www.nature.com/articles/s41588-021-00911-1) | [GEO Dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078)
  - This dataset provides single-cell RNA sequencing data used during automatic cell annotation.

### 2. scVI
- **scVI**: [Setup Instructions](https://docs.scvi-tools.org/en/stable/)
  - scVI (Single-cell Variational Inference) is a tool used for probabilistic modeling of single-cell gene expression data. Follow the link for installation and setup instructions.

### 3. scIB
- **scIB**: [Setup Documentation](https://scib-metrics.readthedocs.io/en/stable/index.html)
  - scIB (Single-cell Integration Benchmarking) is a benchmarking toolkit used to assess the performance of various single-cell integration methods. Refer to the documentation for setup and usage details.

### 4. inferCNV
- **inferCNV**: [Github Repository](https://github.com/broadinstitute/inferCNV/wiki)
  - inferCNV (infer Copy Number Variations) is a tool used to analyze copy number variations in single-cell RNA-seq data. The methods documentation provides detailed information on how to implement and interpret the results.

### 5. GSEA/msigDB
- **GSEA/msigDB**: [Website](https://www.gsea-msigdb.org/gsea/index.jsp)
  - The Gene Set Enrichment Analysis (GSEA) and Molecular Signatures Database (msigDB) are used for functional enrichment analysis. Visit the website for more information and access to the tools.


## Dependencies

To ensure smooth execution of the scripts, please install the following R packages. Some packages are available on CRAN and Bioconductor, while others need to be installed directly from GitHub.

```r
install.packages(c(
  "Seurat", "ggplot2", "ggpubr", "ggbreak", "ggfortify", "mosaic", 
  "RColorBrewer", "dplyr", "data.table", "pheatmap", "tidyverse", 
  "clustree", "harmony", "anndata", "chisq.posthoc.test", 
  "patchwork", "NMF", "ggalluvial", "reticulate", "SeuratWrappers", 
  "magrittr", "sctransform", "scCustomize", "fpc", "ineq", 
  "plyr", "cowplot", "dendsort", "lmerTest", "emmeans", 
  "rstatix", "expm"
))

BiocManager::install(c(
  "EnhancedVolcano", "fgsea", "qusage", "celldex", "SingleR", 
  "ComplexHeatmap", "infercnv", "glmGamPoi", "BiocParallel", 
  "genomicInstability"
))

# Install the remotes package if you haven't already
install.packages("remotes")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
remotes::install_github("immunogenomics/presto")
remotes::install_github("jinworks/CellChat")
remotes::install_github('cole-trapnell-lab/monocle3')
remotes::install_github("cellgeni/sceasy")

```
