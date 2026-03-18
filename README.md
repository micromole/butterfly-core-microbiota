Butterfly Core Microbiota 16S Metabarcoding Pipeline
========================
by Arne Weinhold (LMU Munich) arne.weinhold@bio.lmu.de

## What it is
This is the full processing pipeline for microbiome analysis from 16S Illumina data. It contains all code required to reproduce the analyses, figures and statistics presented in:
> Weinhold A*, Pinos A, Correa-Carmona Y, Holzmann K L, Alonso-Alonso P, Yon F, Steffan-Dewenter I, Peters M K, Brehm G, Keller A: Host-specific and common core microbiota in adult butterflies across continents, bioRxiv https://doi.org/10.64898/2026.02.27.708436

## How to
* Download all files of this repository
* Load ```R_16S_AW_Butterfly_pipeline_01.R``` into R
* Install required libraries
* Adapt path of your specific working directory ```setwd("../your_folder") ```
* Run pipeline

## What it does
* Loads additional functions ```R_16S_AW_functions.R```
* Reimports sample data into ps object ```data.comp```
* Creates output directory for *.pdf figures


## What you get
Main steps of the analysis pipeline:
* data.comp       # Raw project sample data file
* data.bacteria   # Cyano reads / plant reads / unresolved taxa removed
* data.fixed      # Low quality samples removed (low PCR / high cyano reads)
* data.prevfilter # Prevalence cutoff / low stringent filtering / rare taxa removed <0.01%
* data.pruned     # positive control / spike-in taxa removed (mock community)
* data.decontam   # decontam applied on cleared controls
* data.high       # LT2000 Low throughput samples removed
* sample.ASV      # controls removed / only samples on ASV level
* sample.species  # controls removed / only samples on genus / final dataset for analysis
* sample.filter   # optional: filter low abundant genera to simplify phylo tree
