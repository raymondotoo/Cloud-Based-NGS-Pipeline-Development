# ADNI–ADRC Integrated Proteomics Analysis
 
## Overview
 
This repository contains the full analytical workflow for integrative proteomics analysis of two independent Alzheimer's disease cohorts:
 
- Alzheimer's Disease Neuroimaging Initiative (ADNI)
- Alzheimer's Disease Research Center (ADRC)
 
The project focuses on cross-cohort harmonization, systems-level network analysis, functional enrichment, module scoring, and machine learning to identify molecular signatures associated with neurodegeneration and cognitive decline.
 
All analyses are implemented in R Markdown for reproducibility and modular execution.

The original `.Rmd` notebooks are the canonical analysis scripts in this repository. Legacy `*_fixed.Rmd` variants have been moved under `archive/fixed_variants/` for reference only. The main workflow should be read and run from the original notebooks together with `project_setup.R`.
 
---
 
## Repository Structure
 
### 00_Data_preprocessing.Rmd
 
**Purpose:** Initial data cleaning and preprocessing.
 
**Key steps:**
- Raw data import (ADNI and ADRC)
- Quality control filtering
- Missingness filtering
- Normalization (if applicable)
- Metadata harmonization
- Preparation of expression matrices
 
**Output:** Cleaned datasets prepared for harmonization.
 
---
 
### 02a_ADNI_ADRC_harmonized.Rmd
 
**Purpose:** Cross-cohort harmonization.
 
**Key steps:**
- Batch correction between ADNI and ADRC
- Alignment of shared proteins
- Covariate adjustment
- Construction of harmonized expression matrix
 
**Output:** Harmonized proteomics dataset used in downstream analyses.
 
---
 
### 02b_PCA_full.Rmd
 
**Purpose:** Exploratory dimensionality reduction.
 
**Key steps:**
- Principal Component Analysis (PCA)
- Variance explained assessment
- Cohort separation visualization
- Outlier detection
 
**Output:** PCA plots and variance summaries.
 
---
 
### 03_WGCNA_running.Rmd
 
**Purpose:** Weighted Gene Co-expression Network Analysis (WGCNA).
 
**Key steps:**
- Soft-threshold power selection
- Adjacency and TOM construction
- Module detection
- Module eigengene calculation
- Module–trait correlation analysis
 
**Output:**
- Network objects
- Module assignments
- Module eigengenes
 
---
 
### 04_Functional_Analysis.Rmd
 
**Purpose:** Functional annotation of detected modules.
 
**Key steps:**
- Gene/protein ID mapping
- Gene Ontology enrichment
- Pathway enrichment analysis
- Biological interpretation of key modules
 
---
 
### 05_Module_Scoring_Selection.Rmd
 
**Purpose:** Module scoring and feature selection.
 
**Key steps:**
- Calculation of module scores
- Hub protein identification
- Feature reduction strategies
- Selection of predictive modules
 
---
 
### 06_ADNI_ADRC_ML.Rmd
 
**Purpose:** Machine learning modeling.
 
**Key steps:**
- Train/test splitting
- Model development (e.g., regularized regression, tree-based methods)
- Cross-validation
- Feature importance analysis
- Model performance evaluation
 
**Output:** Predictive performance metrics and selected biomarker signatures.
 
---
 
### 07a_Functional_Analysis_GSEA.Rmd
 
**Purpose:** Gene Set Enrichment Analysis (GSEA).
 
**Key steps:**
- Ranked protein/gene list generation
- Pre-ranked GSEA
- Enrichment score calculation
- Pathway visualization
 
---
 
### 07b_Functional_Analysis_ORA.Rmd
 
**Purpose:** Over-Representation Analysis (ORA).
 
**Key steps:**
- Significant protein/module selection
- Hypergeometric testing
- Pathway enrichment analysis
- Comparative interpretation alongside GSEA results
 
---
 
### 08_Supplementary_analysis_v3.Rmd
 
**Purpose:** Supplementary and sensitivity analyses.
 
**Includes:**
- Stratified analyses
- Robustness checks
- Additional figures and tables
- Validation analyses
 
---
 
## Analytical Workflow

1. Data preprocessing  
2. Cross-cohort harmonization  
3. Exploratory PCA  
4. Network construction (WGCNA)  
5. Functional annotation  
6. Module scoring  
7. Machine learning modeling  
8. GSEA  
9. ORA  
10. Supplementary analyses
 
Scripts are designed to be run sequentially in numeric order.

## Pathing And Outputs

- `project_setup.R` now provides portable project-root discovery and standard output directories.
- Canonical generated outputs are organized under `results/pca`, `results/wgcna`, `results/functional_analysis`, `results/functional_analysis_gsea`, `results/functional_analysis_ora`, `results/network_annotation`, `results/ml`, and `results/supplementary`.
- The Shiny app in `app.R` is intended to mirror these same notebook object contracts rather than maintain a separate analysis path.
 
---
 
## Requirements
 
- R (≥ 4.0)
- Recommended packages:
  - WGCNA
  - tidyverse
  - limma
  - clusterProfiler
  - fgsea
  - caret and/or glmnet
  - ggplot2
  - data.table
 
---
 
## Reproducibility Notes
 
- All analyses are implemented in R Markdown notebooks.
- Intermediate objects are saved as `.RData` files where necessary to reduce recomputation time.
- Ensure working directories are properly configured before execution.
- Large intermediate objects may not be version-controlled depending on repository settings.
 
---
 
## Maintainer
 
Raymond Otoo  
Ph.D. Bioinformatics


Last updated: 2026-03-22