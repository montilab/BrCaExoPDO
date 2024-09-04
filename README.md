# Unveiling the Impact of Exosomal Signaling in Obesity-Driven Diabetes using Breast Patient-Derived Organoids

## Overview

This repository contains the data, scripts, and results associated with our study on the impact of exosomal signaling within the tumor microenvironment (TME) in obesity-driven diabetes. Our research demonstrates how exosomes derived from Type 2 diabetic patient plasma influence breast tumor aggression, emphasizing pathways linked to epithelial-to-mesenchymal transition, cancer stemness, and immune evasion.

### Key Findings:
In this study, we developed patient-derived organoids (PDOs) from breast tumor resections, successfully preserving native tumor-infiltrating lymphocytes for the first time in breast PDO research. PDOs were treated with exosomes derived from Type 2 diabetic and non-diabetic patient plasma to assess their impact on tumor aggressiveness and intratumoral heterogeneity. Using single-cell RNA sequencing, we measured transcriptional changes in response to exosomal signaling, providing insights into how systemic metabolic dysregulation, particularly in diabetes, contributes to a more aggressive tumor phenotype.

## Repository Contents

### 1. [Data](data)
   - Contains single-cell metadata (`pdo_metadata.csv`) and organoid circularity measurements (`organoid_circularity.csv`).
   - Raw data can be accessed via [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXX).
   - Interactive exploration available on [CellxGene](https://cellxgene.cziscience.com/collections/xxxxxx).

### 2. [Scripts](scripts)
   - Includes all scripts used for data processing, analysis, and figure generation.
   - Custom pipelines for single-cell RNA sequencing analysis, differential expression analysis, data integration, and other key computational tasks.

### 3. [Results](results)
   - Contains processed data frames, including differential expression results, gene module analyses, and cell type annotations.
   - Organized into subfolders for ease of use.

## BioRxiv Preprint

For detailed methodology and comprehensive results, please refer to our full preprint available on [bioRxiv](https://www.biorxiv.org/)

## Abstract

Women with obesity-driven diabetes are predisposed to more aggressive breast cancers. However, patient metabolic status does not currently inform clinical management. We previously identified plasma exosomes as functionally critical actors in intercellular communication and drivers of tumor progression. Here, we generated patient-derived organoids (PDOs) from breast tumor resections to model signaling within the tumor microenvironment (TME). Novel techniques and a short (1-week) culture preserved native tumor-infiltrating lymphocytes for the first time in breast tumor PDOs. After 3-day exosome treatment, we measured the impact of exosomal signaling on PDOs via single-cell RNA sequencing. Exosomes derived from Type 2 diabetic patient plasma significantly upregulated pathways associated with epithelial-to-mesenchymal transition, invasiveness, and cancer stemness, compared to non-diabetic exosome controls. Intratumoral heterogeneity and immune evasion increased in the diabetic context, consistent with enhanced tumor aggressiveness and metastatic potential of these PDOs. Our model of systemic metabolic dysregulation and perturbed transcriptional networks enhances understanding of dynamic interactions within the TME in obesity-driven diabetes and offers new insights into novel exosomal communication.

## How to Cite

If you use any data or scripts from this repository in your research, please cite our paper as follows:
- Ennis, C.S., Seen, M., Chen, A., Kang, H., Ilinski, A., Mahdaviani, K., Ko, N., Monti, S., and Denis, G.V., 2024. Plasma exosomes from individuals with type 2 diabetes drive breast cancer aggression in patient-derived organoids. (leaving blank until submitted)
