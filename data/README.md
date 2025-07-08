## Overview

This folder contains two key datasets generated from our study: the single-cell metadata (`pdo_metadata.csv`) and organoid circularity measurements (`organoid_circularity.csv`). 

## Access and Explore the Data

- **Raw Data**: The raw single-cell RNA sequencing data associated with this study can be downloaded from the Gene Expression Omnibus (GEO) database from [GSE302054](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302054).

## Datasets

### 1. `pdo_metadata.csv`

This file contains the single-cell metadata with the following characteristics for each cell in the dataset:

- **barcode**: Unique cell identifier.
- **nCount_RNA**: Total RNA molecule counts per cell.
- **nFeature_RNA**: Number of detected genes per cell.
- **nCount_SCT**: Total RNA counts per cell after SCTransform.
- **nFeature_SCT**: Number of genes per cell after SCTransform.
- **scVI_clusters**: Cell clusters from scVI integration.
- **perc_mito**: Percentage of mitochondrial gene reads.
- **perc_ribo**: Percentage of ribosomal gene reads.
- **log10genes_per_UMI**: Log10-transformed ratio of genes per UMI, indicating library complexity.
- **S_Score**: Cell cycle S phase score.
- **G2M_Score**: Cell cycle G2/M phase score.
- **phase**: Predicted cell cycle phase (G1, S, G2M).
- **treatment**: Exosome treatment condition.
- **patientID**: Patient identifier for PDO origin.
- **singleR_labels**: Cell type labels from SingleR.
- **singleR_pruned_labels**: Refined cell type labels post-SingleR pruning.
- **malignancy_score**: Malignancy likelihood score per cell, based on correlation with the top 5% most unstable CNV profiles. See [inferCNV.R](scripts/inferCNV.R).
- **cosmic_CGC**: Module score for COSMIC Cancer Gene Census gene enrichment.
- **gobp_mamm_epi_prolif**: Module score for mammary epithelial cell proliferation (Gene Ontology).
- **coarse_anno**: Coarse-level cell type annotation.
- **broad_anno**: Broad-level cell type annotation.
- **fine_anno**: Fine-level cell type annotation.

### 2. `organoid_circularity.csv`

This file contains the morphological measurements of PDOs, including circularity, area, and perimeter, obtained through image analysis. The key points from the image analysis process are as follows:

- **Imaging System**: Reflectance confocal microscopy was performed using the Live-Duo LSM 710 system (Zeiss) at the BU Microscopy Core.
- **Imaging Conditions**: Images were acquired at 37°C in a humidified atmosphere with 5% CO2, using a Plan-Apochromat 20X/0.8 objective and a 543 nm solid-state laser.
- **Image Processing**: Images were processed using FIJI software. Brightness and contrast were adjusted, and automatic thresholding was applied.
- **Morphological Measurements**: Circularity, area, and perimeter of the PDOs were calculated via particle analysis in FIJI. Single cells were excluded by applying size filters greater than 300 microns.

## SingleR Annotation Reference

The cell type annotations in `pdo_metadata.csv`, particularly the cluster labels, were informed by previous work, including the study by Wu et al. 2021:

- **SCSubtype**: [Nature Genetics Article](https://www.nature.com/articles/s41588-021-00911-1) | [GEO Dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078)

## Usage Notes

- **Privacy Considerations**: The data provided here does not include any patient-identifiable information. Ensure that any downstream analyses also comply with relevant ethical guidelines.
- **Citation**: If you use this dataset in your research, please cite [our paper](https://www.biorxiv.org/content/10.1101/2024.09.13.612950v1.abstract) as well as the Wu et al. paper for the annotations.
