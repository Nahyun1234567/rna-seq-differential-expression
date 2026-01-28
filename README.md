# RNA-seq Differential Expression Analysis

End-to-end RNA-seq differential expression analysis pipeline for identifying
disease-associated genes using Python and statistical methods.

## Overview
This project investigates how a longevity-associated genetic variant affects gene expression patterns in hepatocellular carcinoma (HCC) cells using RNA-seq data analysis.

Recent studies have identified a rare variant of the SIRT6 gene (N308K/A313S) in centenarians, which is associated with improved DNA repair and cancer suppression.
However, the molecular mechanisms by which this variant influences cancer cell behavior remain unclear.

In this project, we analyze RNA-seq gene expression data from liver cancer cells expressing either the wild-type or the centenarian-associated SIRT6 variant, and compare them to control cells.

## Dataset
- Source: GEO (Gene Expression Omnibus)
- Data type: RNA-seq count matrix
- Samples: Disease vs. Control
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314110

### Dataset Description

Data source: NCBI Gene Expression Omnibus (GEO), accession GSE314110
Organism: Homo sapiens
Experiment type: RNA-seq (expression profiling by high-throughput sequencing)
#### Experimental Design
The dataset includes two human hepatocellular carcinoma cell lines:
- HepG2
- Huh-7

Each cell line was subjected to three experimental conditions:
1. Control: Luciferase-only expression
2. SIRT6 WT: Expression of wild-type human SIRT6
3. SIRT6 Variant: Expression of centenarian-associated SIRT6 (N308K/A313S)
Each condition includes three biological replicates, resulting in a total of 18 RNA-seq samples.

## Methods
- Log normalization of raw counts
- Differential expression analysis using t-tests
- Multiple hypothesis testing correction (Benjamini–Hochberg FDR)
- Visualization with volcano plots and heatmaps

## Goal
Identify differentially expressed genes between control and SIRT6-treated cancer cells.
Examine how the centenarian-associated SIRT6 variant may suppress genes related to cancer aggressiveness and invasion.

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
