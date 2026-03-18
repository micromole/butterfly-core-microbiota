Butterfly Core Microbiota 16S Metabarcoding Pipeline
========================
by Arne Weinhold (LMU Munich) arne.weinhold@bio.lmu.de

## What it is
This is the full processing pipeline for microbiome analysis from 16S Illumina data. It contains all code required to reproduce the analyses, figures and statistics presented in:
> Weinhold A*, Pinos A, Correa-Carmona Y, Holzmann K L, Alonso-Alonso P, Yon F, Steffan-Dewenter I, Peters M K, Brehm G, Keller A: Host-specific and common core microbiota in adult butterflies across continents, bioRxiv https://doi.org/10.64898/2026.02.27.708436

## How to
* Download all files of this repository
* Load ```R_16S_AW_Butterfly_pipeline_v01.R``` into R
* Install libraries and adapt path of working directory ```setwd("../butterfly-core-microbiota") ```
* Run pipeline

## What it does
* Loads additional functions ```R_16S_AW_functions.R```
* Imports sample data into ps object ```data.comp```
* Creates output directory for *.pdf files

## Output files
```
00_data pre-processing and filtering
01_sample taxa composition and core analysis
02_sample_div alpha and beta diversity analysis
03_sample_taxa ASV abundance
05_ggtree core tree analysis
06_alluvial network alluvial plot
07_final figures
```

## Alternative input
* Process Illumina data according to https://github.com/chiras/metabarcoding_pipeline/
* Merge Taxonomy, Community table and Metadata into ```data.comp```
