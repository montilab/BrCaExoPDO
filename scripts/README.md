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

## Publicly Available Resources

The scripts in this repository utilize the following publicly available tools:

### 1. Seurat
- **Seurat**: [Paper](https://www.nature.com/articles/s41587-023-01767-y) | [Website](https://satijalab.org/seurat/)
  - Seurat is a comprehensive R package designed for the analysis, visualization, and integration of single-cell data. See their website for extensive documentation and tutorials.
 
### 2. DoubletFinder
- **DoubletFinder**: [Paper](https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30073-0) | [GitHub Repository](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
  - DoubletFinder is an algorithm designed to identify doublets in single-cell data. See the GitHub repository for installation instructions and troubleshooting.
 
### 3. SingleR
- **SingleR**: [Paper](https://www.nature.com/articles/s41590-018-0276-y) | [Github Repository](https://github.com/dviraran/SingleR?tab=readme-ov-file)
  - SingleR is a tool for the automatic annotation of single-cell data. It assigns cell labels based on reference datasets, facilitating the identification of cell types.

### 4. scVI
- **scVI**: [Setup Instructions](https://docs.scvi-tools.org/en/stable/)
  - scVI (Single-cell Variational Inference) is a tool used for probabilistic modeling of single-cell gene expression data. Follow the link for installation and setup instructions.

### 5. scIB
- **scIB**: [Setup Documentation](https://scib-metrics.readthedocs.io/en/stable/index.html)
  - scIB (Single-cell Integration Benchmarking) is a benchmarking toolkit used to assess the performance of various single-cell integration methods. Refer to the documentation for setup and usage details.

### 6. inferCNV
- **inferCNV**: [Github Repository](https://github.com/broadinstitute/inferCNV/wiki)
  - inferCNV (infer Copy Number Variations) is a tool used to analyze copy number variations in single-cell data. The methods documentation provides detailed information on how to implement and interpret the results.

### 7. GSEA/msigDB
- **GSEA/msigDB**: [Website](https://www.gsea-msigdb.org/gsea/index.jsp)
  - The Gene Set Enrichment Analysis (GSEA) and Molecular Signatures Database (msigDB) are used for functional enrichment analysis. Visit the website for more information and access to the tools.

### 8. DAVID
- **DAVID**: [Website](https://david.ncifcrf.gov/)
  - The Database for Annotation, Visualization, and Integrated Discovery (DAVID) is a bioinformatics resource used for functional annotation and enrichment analysis of gene lists. The website provides access to the tools and documentation.
 
### 9. CellChat
- **CellChat**: [Paper](https://www.nature.com/articles/s41467-021-21246-9) | [GitHub Repository](https://github.com/jinworks/CellChat)
  - CellChat is a tool used for inferring and analyzing cell-cell communication networks from single-cell data. The GitHub repository provides installation instructions and usage guidelines.

### 10. Monocle3
- **Monocle3**: [Paper](https://www.nature.com/articles/s41586-019-0969-x) | [Website](http://cole-trapnell-lab.github.io/monocle3/)
  - Monocle3 is a tool for trajectory inference and pseudotime analysis of single-cell data. Visit the website for tutorials and detailed usage information.


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
