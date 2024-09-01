## Overview

This folder contains the results generated from various stages of the analysis. The results are organized into three subfolders: `finetuning`, `analyses`, and `annotations`. Each subfolder contains data frames and output files from different analytical steps.

## Experimental Conditions
- **NDexo**: PDOs treated with exosomes derived from non-diabetic plasma.
- **T2Dexo**: PDOs treated with exosomes derived type 2 diabetic plasma.
- **UT**: Untreated PDOs.

## Subfolders and Contents

### 1. `finetuning`

This subfolder contains the output of the `FindAllMarkers` function from Seurat for different coarse annotations:

- **`Epi_FindAllMarkers.csv`**: Differential expression markers for epithelial cells.
- **`Stro_FindAllMarkers.csv`**: Differential expression markers for stromal cells.
- **`Imm_FindAllMarkers.csv`**: Differential expression markers for immune cells.

### 2. `analyses`

This subfolder contains differential expression results and gene module analyses for various cell types and treatment conditions:

- **`lumtum1_markers.csv`**: Differential expression markers for Luminal Tumor 1 cells compared to all other luminal-like epithelial populations.
- **`lumtum3_markers.csv`**: Differential expression markers for Luminal Tumor 3 cells compared to all other luminal-like epithelial populations.
- **`monocle_gene_modules.csv`**: Gene modules identified during trajectory analysis using Monocle3.
- **`nc_allcell_markers.csv`**: Differential expression comparing NDexo vs UT across all cells.
- **`tc_allcell_markers.csv`**: Differential expression comparing T2Dexo vs UT across all cells.
- **`tn_allcell_markers.csv`**: Differential expression comparing T2Dexo vs NDexo across all cells.
- **`nc_epi_markers.csv`**: Differential expression comparing NDexo vs UT in epithelial cells.
- **`tc_epi_markers.csv`**: Differential expression comparing T2Dexo vs UT in epithelial cells.
- **`tn_epi_markers.csv`**: Differential expression comparing T2Dexo vs NDexo in epithelial cells.
- **`nc_stro_markers.csv`**: Differential expression comparing NDexo vs UT in stromal cells.
- **`tc_stro_markers.csv`**: Differential expression comparing T2Dexo vs UT in stromal cells.
- **`tn_stro_markers.csv`**: Differential expression comparing T2Dexo vs NDexo in stromal cells.
- **`nc_imm_markers.csv`**: Differential expression comparing NDexo vs UT in immune cells.
- **`tc_imm_markers.csv`**: Differential expression comparing T2Dexo vs UT in immune cells.
- **`tn_imm_markers.csv`**: Differential expression comparing T2Dexo vs NDexo in immune cells.

### 3. `annotations`

This subfolder contains the predicted cell type annotations for each sample, as inferred using SingleR with the reference dataset from the [Wu et al. 2021](https://www.nature.com/articles/s41588-021-00911-1) paper:

- **`predicted_celltypes_swarbrick_XXX.csv`**: Predicted cell types for sample `XXX`. 

## Usage Notes

- For reproducibility and further analysis, these files can be re-analyzed from provided [scripts](../scripts) or used as input for downstream steps.
