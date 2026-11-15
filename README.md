Butterfly Core Microbiota 16S Metabarcoding Pipeline
========================
by Arne Weinhold (LMU Munich) arne.weinhold@bio.lmu.de

## What it is
This repository provides the full processing pipeline for microbiome analysis based on 16S Illumina sequencing data. It contains all data and code required to reproduce sample processing, statistical analysis and figures presented in:
> Weinhold A, Pinos A, Correa-Carmona Y, Holzmann K L, Alonso-Alonso P, Yon F, Steffan-Dewenter I, Peters M K, Brehm G, Keller A: **Host-structured adult butterfly microbiomes share a common core across continents**, accepted at Nature Communications (bioRxiv https://doi.org/10.64898/2026.02.27.708436)

## How to
* Download all files of this repository
* Load ```R_16S_AW_Butterfly_pipeline_v02_Revision.R``` into R
* Install libraries and adapt path of working directory ```setwd("../butterfly-core-microbiota") ```
* Run pipeline in R

## What it does
* Loads additional functions ```R_16S_AW_functions.R```
* Imports sample ```data``` into ps object ```data.comp```
* Creates output directory for *.pdf files
* Saves intermediate and final figure into output directory 
* Run time <15 min on a standard desktop computer

## Output files
```
00_data pre-processing and filtering
01_sample taxa composition and core analysis
02_sample_div alpha diversity analysis
03_sample_div beta diversity analysis
04_sample_taxa ASV abundance, endosymbionts
05_ggtree core tree analysis
06_alluvial network alluvial plot
07_Hym_Lep bee data import
08_final figures
```

## Supplementary Information
* Run ```Supplementary_information.tex``` in LaTeX to create the ```Supplementary_information.pdf``` file.

##  Data and code availability

The original processed data and metadata under ```data``` are licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.

When using these data, please cite the associated research article:

> Weinhold A, Pinos A, Correa-Carmona Y, Holzmann K L, Alonso-Alonso P, Yon F, Steffan-Dewenter I, Peters M K, Brehm G, Keller A: **Host-structured adult butterfly microbiomes share a common core across continents**, accepted at Nature Communications 

The R scripts and analytical pipeline are licensed under the MIT License.

