# Unveiling the Impact of Exosomal Signaling in Obesity-Driven Diabetes using Breast Patient-Derived Organoids

## Overview

This repository contains the data, scripts, and results associated with our study on the impact of exosomal signaling within the tumor microenvironment (TME) in obesity-driven diabetes. Our research demonstrates how exosomes derived from Type 2 diabetic patient plasma influence breast tumor aggression, emphasizing pathways linked to epithelial-to-mesenchymal transition, cancer stemness, and immune evasion.

### Key Findings:
In this study, we developed patient-derived organoids (PDOs) from breast tumor resections, successfully preserving native tumor-infiltrating lymphocytes for the first time in breast PDO research. PDOs were treated with exosomes derived from type 2 diabetic and non-diabetic patient plasma to assess their impact on tumor aggressiveness and intratumoral heterogeneity. Using single-cell RNA sequencing, we measured transcriptional changes in response to exosomal signaling, providing insights into how systemic metabolic dysregulation, particularly in diabetes, contributes to a more aggressive tumor phenotype.

## Repository Contents

### 1. [Data](data)
   - Contains single-cell metadata (`pdo_metadata.csv`) and organoid circularity measurements (`organoid_circularity.csv`).
   - Raw data can be accessed via [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302054).

### 2. [Scripts](scripts)
   - Includes all scripts used for data processing, analysis, and figure generation.
   - Custom pipelines for single-cell RNA sequencing analysis, differential expression analysis, data integration, and other key computational tasks.

### 3. [Results](results)
   - Contains processed data frames, including differential expression results, gene module analyses, and cell type annotations.
   - Organized into subfolders for ease of use.

## Abstract

Women with obesity-driven type 2 diabetes (T2D) face worse breast cancer outcomes, yet metabolic status does not fully inform current standards of care. We previously identified plasma exosomes as key drivers of tumor progression; however, their effect on immune cells within the tumor microenvironment (TME) remains unclear. Using a novel patient-derived organoid (PDO) system that preserves native tumor-infiltrating lymphocytes (TILs), we show that T2D plasma exosomes induce a 13.6-fold expansion of immunosuppressive TILs relative to nondiabetic controls. This immune dysfunction may promote micrometastatic survival and resistance to checkpoint blockade, a known issue in T2D cancer patients. Tumor-intrinsic analysis revealed a 1.5-fold increase in intratumoral heterogeneity and 2.3-fold upregulation of aggressive signaling networks. These findings reveal how T2D-associated metabolic dysregulation alters tumor–immune crosstalk through previously underappreciated exosomal signaling, impairing antitumor immunity and accelerating progression. Understanding these dynamics could inform tailored therapies for this high-risk, underserved patient population.

## How to Cite

If you use any data or scripts from this repository in your research, please cite our paper as follows:
- Ennis, C.S., Seen, M., Chen, A. et al. Plasma exosomes from individuals with type 2 diabetes drive breast cancer aggression in patient-derived organoids. Commun Biol **8**, 1276 (2025). https://doi.org/10.1038/s42003-025-08663-y
