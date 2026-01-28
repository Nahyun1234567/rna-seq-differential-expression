# RNA-seq Differential Expression Analysis

End-to-end RNA-seq differential expression analysis pipeline for identifying
disease-associated genes using Python and statistical methods.

## Overview
This project analyzes RNA-seq gene expression count data to identify
differentially expressed genes (DEGs) between disease and control samples.
The pipeline includes data preprocessing, normalization, statistical testing,
and visualization.

## Dataset
- Source: GEO (Gene Expression Omnibus)
- Data type: RNA-seq count matrix
- Samples: Disease vs. Control

## Methods
- Log normalization of raw counts
- Differential expression analysis using t-tests
- Multiple hypothesis testing correction (Benjamini–Hochberg FDR)
- Visualization with volcano plots and heatmaps

## Repository Structure
- `data/`: Raw and processed datasets
- `notebooks/`: Step-by-step analysis notebooks
- `src/`: Reusable Python modules
- `results/`: Figures and output tables
- `report/`: Summary of findings

## Skills Demonstrated
- Bioinformatics data analysis
- Statistical testing
- Data visualization
- Reproducible research workflow
