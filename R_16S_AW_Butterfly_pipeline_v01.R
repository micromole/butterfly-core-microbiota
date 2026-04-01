# ------------------------------------------------------------- x
# Full R pipeline for microbiota metabarcoding analysis
# based on https://github.com/chiras/metabarcoding_pipeline
# by Alexander Keller (LMU) keller@bio.lmu.de
#
# Created: Do 11. Apr 12:31:38 CEST 2024
# Project: 16S_AW_Butterfly_Peru
# Marker: 16S
# Author: Arne Weinhold (LMU) arne.weinhold@bio.lmu.de
# ------------------------------------------------------------- x

# Clear workspace
rm(list = ls())

# Loading in necessary libraries
library(shades) # color saturation, load first to avoid distance masking problem
library(ggplot2)
library(tidyverse)
library(bipartite) # to sort data frames
library(microViz) # comp barplot
library(microbiome) # aggregate_taxa plot_core function
library(ggpubr) # stat_regline_equation in fig
library(RColorBrewer) # figure color
library(vegan)
library(decontam)
library(ggtree) 
library(ape) # tree as.DNAbin()
library(Biostrings) # import sequence data for alignment
library(DECIPHER) # for alignment and tree construction
library(ggtreeExtra) # function "geom_fruit" 
library(phyloseq)
library(circlize) # circos chordplot
library(ggVennDiagram)
library(ggnewscale) # new_scale_fill()
library(ggalluvial) # Alluvial plot
library(patchwork) # combine and arrange plots
library(car) # Type II ANOVA, Levene test, VIF
library(effectsize) # eta squared

# sessionInfo() # R and package versions
# search() 

## Setting working directory 
setwd("../LRZ Sync+Share/data/16S_AW_Butterfly21_Peru_pipeline")
setwd("../16S_AW_Butterfly21_Peru_pipeline")
getwd()

# Overview about main ps objects of the processing pipeline:
# data.comp       # Selected project data
# data.bacteria   # Non-bacterial reads and unresolved taxa removed
# data.fixed      # Low quality samples removed (low PCR / high cyano reads)
# data.prevfilter # Prevalence filter to remove rare taxa (<0.01%) 
# data.pruned     # positive control / spike-in taxa removed
# data.decontam   # decontam package applied on cleared controls 
# data.high       # Low throughput samples removed LT2000
# sample.species  # Samples on genus level -> final dataset for most analysis
# sample.filter   # Optional: Filter minor genera to simplify phylo tree

# Analysis pipeline:  
# 01 Composition Core / 02 Alpha and beta diversity / 03 Taxa abundance / 05 Phylo tree / 07 final figures

###  Load custom themes and functions 
source('./R_16S_AW_functions.R')

# Create output folder for plots and data
out_dir <- "plots_peru" # adjust output folder name for project
dir.create(out_dir, showWarnings = FALSE)
data_dir <- "data"      # set data directory
dir.create(data_dir, showWarnings = FALSE) # create data_dir

sink("sessionInfo.txt")
sessionInfo()
sink()

## 00 data.comp  -------
# Build from sequencing data or re-import from csv export file (line 117)
### Loading data: Taxonomy
# data.tax <- tax_table(as.matrix(read.table("taxonomy.vsearch", header=T,row.names=1,fill=T,sep=",")))
## Community table
# data.otu <- otu_table(read.table("asv_table.merge.txt"), taxa_are_rows=T)
## Sample metadata 
# data.map <- 	sample_data(read.table("samples.csv", header=T, row.names=2,  sep=";", fill=T))
# sample_names(data.map) <- gsub("-",".",sample_names(data.map)) # optional if sample names include "-")

## check metadata vs. samples in sequencing data consistency
# sample_names(data.map)%in%sample_names(data.otu)
# sample_names(data.otu)%in%sample_names(data.map)
# sample_names(data.map )[!(sample_names(data.map ) %in% sample_names(data.otu))]
# sample_names(data.otu )[!(sample_names(data.otu ) %in% sample_names(data.map))]

## check taxa names
# taxa_names(data.tax )[!(taxa_names(data.tax ) %in% taxa_names(data.otu))]
# taxa_names(data.otu )[!(taxa_names(data.otu ) %in% taxa_names(data.tax))]

## merge the three tables to a single phylseq object
# data.ps <- merge_phyloseq(data.otu, data.tax, data.map)

## given hierarchical classification options at the end, we have to propagate the taxonomy over taxonomic levels to not throw out stuff only classified to higher tax levels
# data.ps <- propagate_incomplete_taxonomy(data.ps)

## 00 data.comp select samples for projects
# unique(sample_data(data.ps)$study) # show study projects available on the chip
# data.comp <- subset_samples(data.ps, study=="Peru" | study=="Other" )  

# Cleanup pieline
# rm(data.ps)
# rm(data.otu)

## Export ps data.comp
# data.comp.df <- as(sample_data(data.comp),"data.frame") 
# write.csv(data.comp.df, file.path(data_dir, "data.comp.metadata.csv"), row.names = TRUE)
# data.comp.otu <- as.data.frame(otu_table(data.comp))
# write.csv(data.comp.otu, file.path(data_dir, "data.comp.otu.csv"), row.names = TRUE)
# data.comp.tax <- as.data.frame(tax_table(data.comp))
# write.csv(data.comp.tax, file.path(data_dir, "data.comp.tax.csv"), row.names = TRUE)


# Re-import data.comp from csv export file
re_otu <- read.csv(file.path(data_dir, "data.comp.otu.csv"), row.names = 1)
otu_tab <- otu_table(as.matrix(re_otu), taxa_are_rows = T)
re_tax <- read.csv(file.path(data_dir, "data.comp.tax.csv"), row.names = 1)
tax_tab <- tax_table(as.matrix(re_tax))
re_meta <- read.csv(file.path(data_dir, "data.comp.metadata.csv"), row.names = 1)
sample_tab <- sample_data(re_meta)

# Reconstruct phyloseq object
data.comp <- phyloseq(otu_tab, tax_tab, sample_tab)

### General data overview
data.comp
tail(tax_table(data.comp))
rank_names(data.comp)
table(tax_table(data.comp)[, "kingdom"], exclude = NULL)
table(tax_table(data.comp)[, "phylum"], exclude = NULL)
sample_variables(data.comp)

# Cleanup pipeline
rm(re_otu)
rm(otu_tab)

## 00 data.bacteria ------------------------
# Removal of non-bacterial reads: Cyanobacteria / plant organelles / unresolved taxa

sample.comp <- subset_samples(data.comp, type=="sample")

# Histogram of sample read counts, Make a data frame for the read counts of each sample
data.comp.df <- data.frame(sample_data(data.comp))
data.comp.df$LibrarySize <- sample_sums(data.comp)
data.comp.df <- data.comp.df[order(data.comp.df$LibrarySize), ]
data.comp.df$Rank <- seq_len(nrow(data.comp.df))

first.view.histogram <- ggplot(data.comp.df, aes(x = LibrarySize)) +
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  scale_x_continuous(breaks = seq(0, 100000, 5000)) +
  labs(title = "sequencing depth", x = "sample size (reads)") +
  theme_line2() 

first.view.lib <- ggplot(data.comp.df, aes(x = Rank, y = LibrarySize, color = type)) +
  geom_point() + labs(title = "sequencing depth", x = "sample rank", y = "sample size (reads)") + theme_grid()

sampling.depth <- ggplot(data.comp.df, aes(x = chip, y = LibrarySize)) +
  geom_boxplot() +
  geom_point(aes(color = type), alpha = 0.5) +
  theme_line2() +
  labs(title = "sequencing depth", y = "sample size (reads)")

# First view barplot
# filter data.comp first as there can be trouble with minor taxa names
data.comp.p <- tax_filter(data.comp, min_prevalence = 2, min_total_abundance = 50, min_sample_abundance = 10)
first.view.comp <- data.comp.p %>%
  #ps_filter(type == "sample") %>%
  comp_barplot(tax_level = "phylum", n_taxa = 15, sample_order = "asis",
               label = NULL ) +
  facet_wrap(vars(country), scales = "free") +  
  coord_flip() + ggtitle("data.comp", "Cyanobacteria and Plantae" ) +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

## filter Plant, Chloroplast, Algae, fungi and unresolved taxa (adjust according to names in taxonomy.vsearch)
(data.bacteria = subset_taxa(data.comp, kingdom=="d:Bacteria" | kingdom== "d:Archaea" )) # 7997 
(data.bacteria = subset_taxa(data.bacteria, phylum!="p:Cyanobacteria/Chloroplast")) # 7997
(data.bacteria = subset_taxa(data.bacteria, phylum!="p:Cyanobacteria"))  # 7929
(data.bacteria = subset_taxa(data.bacteria, family!="f:Mitochondria"))  # 7929
(data.bacteria = subset_taxa(data.bacteria, order!="f:Mitochondria"))  # 7929
(data.bacteria = subset_taxa(data.bacteria, genus!="d:Bacteria_spc_spc_spc_spc")) # 7380
(data.bacteria = subset_taxa(data.bacteria, kingdom!="")) # 7380
#(data.bacteria <- subset_taxa(data.bacteria, order != "p:Proteobacteria_spc" &  order != "p:Firmicutes_spc" &  order != "p:Actinobacteria_spc"))  


data.bacteria.p <- tax_filter(data.bacteria, min_prevalence = 2, min_total_abundance = 50, min_sample_abundance = 10) # simplyfy data for figure
second.view.comp <- data.bacteria.p %>%
  #ps_filter(type == "sample") %>%
  comp_barplot(tax_level = "phylum", n_taxa = 15, sample_order = "asis",
               label = NULL ) +
  facet_wrap(vars(country), scales = "free") +  
  coord_flip() + ggtitle("data.bacteria","Cyano / plants removed") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

# Check what has been removed
cyano.removed  <- prune_taxa(setdiff(names(taxa_sums(data.comp)), names(taxa_sums(data.bacteria))), data.comp)
cyano.removed.sample <- subset_samples(cyano.removed, type=="sample")


# Calculate percentage of cyano reads removed
percentage_removed_cyano_sample <- (sum(sample_sums(cyano.removed.sample)) / 
                                      sum(sample_sums(sample.comp ))) * 100

percentage_removed_cyano_total <- (sum(sample_sums(cyano.removed)) / 
                                     sum(sample_sums(data.comp))) * 100

percentage_removed_cyano_sample # 0.374%
percentage_removed_cyano_total  # 0.388%


# Show top removed taxa with genus names e.g. to identify d:Bacteria_spc_spc_spc_spc with Blast search
cyano.taxa.abundance <- data.frame(sort(taxa_sums(cyano.removed), decreasing = TRUE)[1:50])
cyano.taxa.genus  <- data.frame(
  genus = as.character(tax_table(cyano.removed)[names(sort(taxa_sums(cyano.removed), decreasing = TRUE)[1:50]), "genus"]),
  Taxa_Sum = sort(taxa_sums(cyano.removed), decreasing = TRUE)[1:50])

# Show list of removed ASVs
cyano.taxa.genus

# Make dataframe with percent removed Cyano reads per sample
cyano.removed.frame <- data.frame(
  cyano = sample_sums(cyano.removed),
  samplesum = sample_sums(data.comp),
  host_genus = sample_data(cyano.removed)$host_genus,
  type = sample_data(cyano.removed)$type,
  PCR = sample_data(cyano.removed)$PCR,
  percent_removed = (sample_sums(cyano.removed) / sample_sums(data.comp)) * 100)
# Sort the dataframe in decreasing order
cyano.removed.frame.sorted <- cyano.removed.frame[order(cyano.removed.frame$cyano), ]

# Show list of samples and percent cyano removed 
tail(cyano.removed.frame.sorted, 100)


# Mean percentage of cyano removed by genus
mean_cyanoremoved_by_genus <- cyano.removed.frame.sorted %>%
  mutate(percent_removed = as.numeric(percent_removed)) %>%
  group_by(host_genus) %>%
  summarise(cyano_percent = round(mean(percent_removed, na.rm = TRUE), 2)) %>% 
  arrange(desc(cyano_percent))

# List percentage removed per genus
as.data.frame(mean_cyanoremoved_by_genus)

# Show removed taxa
# simplify to higher tax rank to speed up figures (and safe space)
cyano.removed.p <- tax_glom(cyano.removed,taxrank="phylum")
cyano.removed.melt <- psmelt(cyano.removed.p)

cyano.removed.lowrank <- ggplot(cyano.removed.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 60)) + ggtitle("Cyano / plants removed") 

cyano.removed.plot <- ggplot(cyano.removed.melt, aes(x=Sample, y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity") + ggtitle("Cyano / plants removed") + theme(axis.text.x = element_blank()) + facet_wrap(~host_genus, scales="free_x")

# Compare data frames before and after Cyano removal  
data.bacteria.df <- as(sample_data(data.bacteria),"data.frame")
data.bacteria.df$samplesums <- sample_sums(data.bacteria)
data.bacteria.df$sumsbefore <- sample_sums(data.comp)
data.bacteria.df$cyanoremoved <- sample_sums(cyano.removed)
data.bacteria.df$percentcyano <-  ((data.bacteria.df$cyanoremoved) / (data.bacteria.df$sumsbefore)) * 100
data.bacteria.df$Shannon <- estimate_richness(data.bacteria, measures = "Shannon")$Shannon

cyano.percentremoved <- ggplot(data.bacteria.df  , aes(x=percentcyano , y=samplesums, color=subtype, size = cyanoremoved)) +
    geom_point(alpha=0.7) + ylim(0, 50000) + xlim(0, 50) + facet_wrap(~type) +
    labs(x = "percent cyano removed", y = "sample sums", title = "Cyano / plants removed") 

cyano.percentremoved.total <- ggplot(data.bacteria.df  , aes(x=percentcyano, y=  cyanoremoved, color=subtype, size = samplesums)) + 
  geom_point(alpha=0.7) + xlim(0, 50)  + facet_wrap(~type) + 
  geom_hline(yintercept = 500, color = "red", linetype = "dashed") + geom_vline(xintercept = 5, color = "blue", linetype = "dashed") +
  labs(x = "percent cyano removed", y = "cyano reads removed", title = "Cyano / plants removed") 

cyano.percent.boxplot <- ggplot(data.bacteria.df, aes(x = type, y = percentcyano, fill = type)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  
  geom_jitter(width = 0.2, shape = 21, size = 2,  alpha = 0.7) +  
  labs(x = "", y = "percent cyano removed", title = "Cyano / plants removed") 
  
data.bacteria.shannon <- ggplot(data.bacteria.df, aes(x=Shannon, y=samplesums, size = cyanoremoved, color=type, shape=country))  + geom_point(alpha=0.7)  + ggtitle(
  "data.bacteria"  )  +  ylim(0, 75000) + xlim(0, 6)


pdf(file.path(out_dir, "00_data_comp_first_view_cyano.pdf"), width=12, height=6)
first.view.histogram
first.view.lib
sampling.depth
first.view.comp
second.view.comp
cyano.removed.lowrank
cyano.removed.plot
cyano.percentremoved
cyano.percentremoved.total
cyano.percent.boxplot
data.bacteria.shannon
dev.off()


options(max.print = 4000) 
sink(file.path(out_dir, "00_data_comp_first_view_cyano.txt"))
"data.comp"
data.comp
table(sample_data(sample.comp)$country)
rank_names(data.comp)
sample_variables(data.comp)
"data.comp"
table(tax_table(data.comp)[, "kingdom"], exclude = NULL)
"data.bacteria"
data.bacteria
"data.bacteria"
table(tax_table(data.bacteria)[, "kingdom"], exclude = NULL)
"cyano.removed"
table(tax_table(cyano.removed)[, "phylum"], exclude = NULL)
"cyano taxa removed"
cyano.taxa.genus 
"cyano removed per sample sorted"
print(cyano.removed.frame.sorted)
cat("Percent of reads removed cyano total:", round(percentage_removed_cyano_total, 2), "%\n")
cat("Percent of reads removed cyano sample:", round(percentage_removed_cyano_sample, 2), "%\n")
"mean_cyanoremoved_by_genus"
as.data.frame(mean_cyanoremoved_by_genus)
sink()



### > Single sample check / Taxon Check----------------

# Inspect sample community composition for potential outlier 
bacteria.comp.plot <- data.bacteria.p %>% 
  ps_filter(host_genus == "Aglais") %>% 
  comp_barplot(tax_level = "genus", n_taxa = 30,merge_other = F, sample_order = "bray",
               label = "SAMPLE" ) +
  facet_wrap(vars(host_genus), scales = "free") +  
  coord_flip() + ggtitle("data.bacteria") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

# Inspect unidentified ASVs in a single sample for manual BLAST and classification
single.sample.check <- prune_samples(sample_names(data.comp) == "AW.B1.G067_S67", data.comp)  
sample_sums(single.sample.check) # list total sample sum

# Show genus names and taxa sums of top 20 ASVs in a sample
data.frame(genus = as.character(tax_table(single.sample.check)[names(sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:20]), "genus"]),
  Taxa_Sum = sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:20])

# Check samples for a specific ASV
# e.g. ASV341,d:Bacteria,p:uncultured_bacteria 
tag_ASV <- c("ASV13") 
singleASV <- prune_taxa(tag_ASV, data.comp)
sort(data.frame(sum = sample_sums(singleASV)), decreasing = FALSE)
tax_top(singleASV, n = 5, rank = "genus")

### > Taxon check 
# e.g. p:Cyanobacteria, p:Abditibacteriota p:Thermotogae p:Gemmatimonadetes p:Verrucomicrobia 
# check_taxa <- subset_taxa(data.bacteria, genus == "p:Proteobacteria_spc_spc_spc")
# check_taxa <- subset_taxa(data.bacteria, genus == "o:Enterobacterales_spc_spc")
# check_taxa <- subset_taxa(data.bacteria, genus == "f:Enterobacteriaceae_spc")
check_taxa <- subset_taxa(data.bacteria, phylum == "p:Deinococcus-Thermus")

# Show samples with that taxon
sort(data.frame(sum = sample_sums(check_taxa),
  host_genus = sample_data(check_taxa)$host_genus), decreasing = FALSE)

# Show top ASVs with that taxon
sort(data.frame(sum = taxa_sums(check_taxa)), decreasing = FALSE)


## 00 data.fixed  -----------
# Remove low quality samples (troublemaker):
# PCR negative samples and samples with high cyano reads e.g. >10%
# Afterwards fix taxonomy with replace_tax_prefixes and phyloseq_validate functions

data.bacteria

# Insepct PCR quality high, med, low
PCR.depth <- data.bacteria.p %>%
  #ps_filter(type == "sample") %>%
  comp_barplot(tax_level = "genus", n_taxa = 30, sample_order = "bray",
               label = "host_genus" ) +
  facet_wrap(vars(PCR), scales = "free") +  
  coord_flip() + ggtitle("data.bacteria","PCR depth") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

# Select low amplification / negative PCR samples
data.bacteria.low = subset_samples(data.bacteria, PCR=="low")

# Show sample sum of low samples
data.frame(sort(sample_sums(data.bacteria.low), decreasing = F))

# Check taxa in low PCR samples 
check_taxa <- subset_taxa(data.bacteria.low, genus == "g:Wolbachia")
sort(data.frame(sum = sample_sums(check_taxa),
                host_genus = sample_data(check_taxa)$host_genus), decreasing = FALSE)
sort(data.frame(sum = taxa_sums(check_taxa)), decreasing = FALSE)

# Ad sample sum to low PCR samples
sample_data(data.bacteria.low)$samplesums <- sample_sums(data.bacteria.low)

# PCR low samples and sample sums
PCR.low.samplesums <- data.bacteria.low %>%
  #ps_filter(host_genus == "Aglais") %>%
  comp_barplot(tax_level = "genus", n_taxa = 30, merge_other = TRUE, sample_order = "bray",
               label = "samplesums") + # SAMPLE
  facet_wrap(vars(host_genus), scales = "free") +
  coord_flip() +
  ggtitle("data.bacteria", "PCR low") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

# Select sample names with low PCR
metadata.bacteria.low <- sample_data(data.bacteria.low)

PCR.low.names <- sample_names(data.bacteria.low)[
  metadata.bacteria.low $type == "sample" &
  metadata.bacteria.low $host_genus != "Aglais" ] # Exclude Aglais to keep it in the dataset

# Define thresholds to remove samples with high Cyano reads
max_percent <- 10 # max 10%
max_reads <- 500 # max 500 reads

# Identify samples with high Cyano reads
high_cyano_read_samples <- rownames(cyano.removed.frame)[
  cyano.removed.frame$type == "sample" &   
    (cyano.removed.frame$percent_removed > max_percent | cyano.removed.frame$cyano > max_reads) ]

# Show samples with high Cyano reads
high_cyano.reads.subset <- prune_samples(sample_names(data.comp) %in% high_cyano_read_samples, data.comp)
high_cyano.reads.subset <- tax_filter(high_cyano.reads.subset , min_prevalence = 2, min_total_abundance = 10, min_sample_abundance = 10) # simplify dataset for figure
high.cyano.samples.comp <- 
  high_cyano.reads.subset %>%
  comp_barplot(tax_level = "genus", n_taxa = 30,merge_other = F, sample_order = "bray",
               label = "SAMPLE" ) +
  facet_wrap(vars(host_genus), scales = "free") +  
  coord_flip() + ggtitle("data.bacteria", "samples with > 500 cyano reads") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))


# Combine high Cyano with low PCR samples (partly overlapping)
tag_troublemaker <- union(high_cyano_read_samples, PCR.low.names) 

# Show samples troublemaker samples
troublemaker.subset <- prune_samples(sample_names(data.bacteria) %in% tag_troublemaker, data.bacteria)
troublemaker.subset <- tax_filter(troublemaker.subset , min_prevalence = 2, min_total_abundance = 10, min_sample_abundance = 10) # simplify dataset for figure
troublemaker.comp <- 
  troublemaker.subset %>%
  comp_barplot(tax_level = "genus", n_taxa = 30,merge_other = F, sample_order = "bray",
               label = "SAMPLE" ) +
  facet_wrap(vars(host_genus), scales = "free") +  
  coord_flip() + ggtitle("troublemaker.subset") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

# Make ordination to highlight troublemaker
all.PCR.low.names <- sample_names(data.bacteria.low)
all_troublemaker <- union(high_cyano_read_samples, all.PCR.low.names) 

data.bacteria.rel <- transform_sample_counts(data.bacteria.p, function(x) x/sum(x))
bacteria.PCoA <- ordinate(data.bacteria.rel, method="PCoA",distance = "bray")

pcoa_df <- plot_ordination(data.bacteria.rel, bacteria.PCoA, justDF = TRUE)
pcoa_df$SampleID <- rownames(pcoa_df)
pcoa_df$label <- ifelse(pcoa_df$SampleID %in% all_troublemaker, pcoa_df$SampleID, NA)
pcoa_df$nolabel <- ifelse(!(pcoa_df$SampleID %in% all_troublemaker), pcoa_df$SampleID, NA)

troublemaker.ordinate.plot <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = type, shape = type)) +
  geom_point(size = 4) +
  geom_label(aes(label = label), size = 3, na.rm = TRUE) +
  theme_grid() 

# Remove troublemaker from all datasets (needed later for filterframe)
data.comp.subset <- prune_samples(!(sample_names(data.comp) %in% tag_troublemaker), data.comp)
data.bacteria.subset <- prune_samples(!(sample_names(data.bacteria) %in% tag_troublemaker), data.bacteria)
cyano.removed.subset <- prune_samples(!(sample_names(cyano.removed) %in% tag_troublemaker), cyano.removed)

# Many ASVs might contain zero counts (since hyper diverse samples are now removed)
data.frame(sort(taxa_sums(data.bacteria.subset), decreasing = F)[1:1000])
# Remove zero counts from dataset (singletons are removed later)
data.bacteria.condensed = prune_taxa(taxa_sums(data.bacteria.subset )>0, data.bacteria.subset)
data.frame(sort(taxa_sums(data.bacteria.condensed), decreasing = F)[1:1000])

### Make taxa labels nice for plots 
# removes 'd:' 'p:' 'o:' in taxa names
data.bacteria.condensed <- replace_tax_prefixes(data.bacteria.condensed)

### Check the names
tail(tax_table(data.bacteria.condensed))

# Validate data with tax_fix, since replace_tax_prefix makes shorter names which can be problematic
data.bacteria.condensed <- phyloseq_validate(data.bacteria.condensed)
# tax_fix_interactive(data.bacteria.condensed) # adjust conditions with tax interactive
data.fixed <- tax_fix(data.bacteria.condensed,
                      min_length = 4,
                      unknowns = c("NA"),
                      sep = " ", anon_unique = TRUE,)

phyloseq_validate(data.fixed,
                  remove_undetected = TRUE,
                  min_tax_length = 4,
                  verbose = TRUE )


pdf(file.path(out_dir, "00_data_delete_troublemaker.pdf"), width=12, height=6)
PCR.depth
PCR.low.samplesums 
high.cyano.samples.comp
troublemaker.comp
troublemaker.ordinate.plot
dev.off()

# Cleanup pipeline
rm(data.comp.p)
rm(data.bacteria.subset)
rm(data.bacteria.condensed)
rm(data.bacteria.low)
rm(data.bacteria.rel)
rm(data.bacteria.p)
rm(PCR.depth)
rm(high.cyano.samples.comp)
rm(troublemaker.comp)


# Pre-Check which samples have unique and low abundant ASVs, optionally remove samples from which more than xx% would be removed
data.prevcheck <- tax_filter(data.fixed , min_prevalence = 3, min_total_abundance = 50, min_sample_abundance = 5 )

prevcheck.removed <-  prune_taxa(setdiff(taxa_names(data.fixed), taxa_names(data.prevcheck)), data.fixed)

prevcheck.removed.frame <- data.frame(
  prevfilter = sample_sums(prevcheck.removed),
  samplesum = sample_sums(data.fixed),
  host_genus = sample_data(prevcheck.removed)$host_genus,
  prevfilter_percent = (sample_sums(prevcheck.removed) / sample_sums(data.fixed)) * 100 )

# Sort the dataframe in decreasing order
prevcheck.removed.frame.sorted <- prevcheck.removed.frame[order(prevcheck.removed.frame$prevfilter_percent), ]

# Pre check to identify hyperdiverse samples (from which more then xx% ASVs would be removed)
print(prevcheck.removed.frame.sorted)

## 00 data.prevfilter  ---------------------
# Prevalence filter to remove rare taxa / spurious phyla / singletons
# ASVs <0.001% / genera <0.01%

data.fixed
sample.fixed <- subset_samples(data.fixed, type=="sample")

# Low stringent filtering to remove spurious phyla e.g. Acidobacteria, Tenericutes, Verrucomicrobia, candidate division WPS−1 etc.
# Show prevalence of phyla 
prev_results <- calc_prevalence(data.fixed, rank = "phylum")
prevdf1 <- subset(prev_results$taxa_table, phylum %in% get_taxa_unique(data.fixed, "phylum"))

### Prevalence Plot inspect low prevalent phyla
prevalence <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(data.fixed),color=phylum)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + # set min prevalence
  geom_vline(xintercept = 100.0, alpha = 0.5, linetype = 2) + # set min abundance
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none") + ggtitle("data.fixed",
    "example cutoff (min total abundance = 100, min prevalence = 5%)") 

prevalence_cutoff <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(data.fixed),color=phylum)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + # set min prevalence
  geom_vline(xintercept = 100.0, alpha = 0.5, linetype = 2) + # set min abundance 
  # scale_x_log10() +  
  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none") +  coord_cartesian(
    xlim = c(0, 400),
    ylim = c(0, 0.25),
    expand = TRUE,
    default = FALSE,
    clip = "on") + ggtitle("data.fixed", 
      "example cutoff (min total abundance = 100, min prevalence = 5%)"   ) 

# Adjust conditions for min abundance and min prevalence filtering
# Calculate good cutoff e.g. less than 0.01 percent of data set (0.01/100)
data.sums <- data.frame(sum = sample_sums(data.fixed))
cutoff_min_total_abundance <- ((0.01/100 ) * sum(data.sums)) # 0.01% 
cutoff_min_total_abundance # equals about 600 reads  

cutoff_percent <- cutoff_min_total_abundance*100/(sum(data.sums))  
cutoff_percent # 0.01% cutoff


# Inspect data: Most ASVs have less than 100 reads 
data.frame(sort(taxa_sums(data.fixed), decreasing = F)[1:1000])
data.fixed.prune500		= prune_taxa(taxa_sums(data.fixed)<500, data.fixed) # speed up figure
ASV_sums_df <- data.frame(ASV_reads = taxa_sums(data.fixed.prune500))

# Highlight singletons and low abundant ASVs
ASV.histogram.prefilter <- ggplot(ASV_sums_df, aes(x = ASV_reads)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2) +
  ggtitle("Low abundant ASV read distribution prefilter") + 
  coord_cartesian(ylim = c(0, 800), xlim = c(0, 500)) +
  theme(axis.title.y = element_blank())

# Calculate min_prevalence e.g. 1.5% of samples
cutoff_min_prevalence <- (1.5*nsamples(data.fixed))/100
cutoff_min_prevalence # equals about 3 samples --> set this number as min_prevalence
percent_of_samples <- (3*100)/nsamples(data.fixed)
percent_of_samples


### Low stringent filtering ASV level
# For Peru set min_prevalence = 3, min_total_abundance = 60 0.001 percent (60 reads)
data.filter.ASV <- tax_filter(data.fixed,  min_prevalence = 3, min_total_abundance = 60, min_sample_abundance = 10 )
data.frame(sort(taxa_sums(data.filter.ASV), decreasing = F)[1:500])

ASV_percentage_removed <- 100 * (1 - sum(sample_sums(data.filter.ASV)) / sum(sample_sums(data.fixed)))
ASV_percentage_removed # 0.8% removed

# Low stringent filtering genus level
data.filter.g <- aggregate_taxa(data.filter.ASV, "genus")
data.frame(sort(taxa_sums(data.filter.g))) # remove genus below  0.01 percent (600 reads)
data.filter.genus <- tax_filter(data.filter.ASV, tax_level = "genus", min_prevalence = 3, min_total_abundance = 600, min_sample_abundance = 10 )
data.filter.g <- aggregate_taxa(data.filter.genus, "genus")
data.frame(sort(taxa_sums(data.filter.g)))

genus_percentage_removed <- 100 * (1 - sum(sample_sums(data.filter.genus)) / sum(sample_sums(data.filter.ASV)))
genus_percentage_removed  # 0.39% removed

# optional: Filter on higher taxonomic ranks
data.filter.f <- aggregate_taxa(data.filter.genus, "family")
data.frame(sort(taxa_sums(data.filter.f))) 

data.filter.o <- aggregate_taxa(data.filter.genus, "order")
data.frame(sort(taxa_sums(data.filter.o))) # optional remove unresolved order

data.filter.p <- aggregate_taxa(data.filter.genus, "phylum")
data.frame(sort(taxa_sums(data.filter.p))) 

data.prevfilter <- data.filter.genus
#sort(data.frame(sum = sample_sums(data.prevfilter)), decreasing = TRUE)

data.filter.prune500		= prune_taxa(taxa_sums(data.prevfilter)<500, data.prevfilter) # speed up figure
ASV_sums_df_postfilter <- data.frame(ASV_reads = taxa_sums(data.filter.prune500))

ASV.histogram.postfilter <- ggplot(ASV_sums_df_postfilter, aes(x = ASV_reads)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2) +
  ggtitle("Low abundant ASV read distribution postfilter") + 
  coord_cartesian(ylim = c(0, 800), xlim = c(0, 500)) +
  theme(axis.title.y = element_blank())

## Check what has been filtered
prevfilter.removed  <- prune_taxa(setdiff(names(taxa_sums(data.fixed)), names(taxa_sums(data.prevfilter))), data.fixed)
prevfilter.removed.sample <- subset_samples(prevfilter.removed , type=="sample")

# Calculate percentage of reads removed (sample level and total)
percentage_removed_sample <- (sum(sample_sums(prevfilter.removed.sample)) / 
                                sum(sample_sums(sample.fixed))) * 100

percentage_removed_total <- (sum(sample_sums(prevfilter.removed)) / 
                               sum(sample_sums(data.fixed))) * 100

percentage_removed_sample # 1.08%
percentage_removed_total #  1.19%

# Show filtered ASVs and genus names
prevfilter.removed.taxa.genus <- data.frame(
  Genus = tax_table(prevfilter.removed)[names(sort(taxa_sums(prevfilter.removed), decreasing = TRUE)[1:50]), "genus"],
  Abundance = sort(taxa_sums(prevfilter.removed), decreasing = TRUE)[1:50] )

# Show list of removed ASVs
prevfilter.removed.taxa.genus

# Make dataframe with percent removed Prevfilter per sample 
prevfilter.removed.frame <- data.frame(
  prevfilter = sample_sums(prevfilter.removed),
  samplesum  = sample_sums(data.fixed),
  host_genus = sample_data(prevfilter.removed)$host_genus,
  type       = sample_data(prevfilter.removed)$type,
  PCR        = sample_data(prevfilter.removed)$PCR,
  prevfilter_percent = (sample_sums(prevfilter.removed) / sample_sums(data.fixed)) * 100)
# Sort the dataframe in decreasing order
prevfilter.removed.frame.sorted <- prevfilter.removed.frame[order(prevfilter.removed.frame$prevfilter_percent), ]

print(prevfilter.removed.frame.sorted)


# Mean percentage of prevfilter removed reads per genus
mean_prevfilter_by_genus <- prevfilter.removed.frame.sorted %>%
  mutate(prevfilter_percent = as.numeric(prevfilter_percent)) %>%
  group_by(host_genus) %>%
  summarise(prevfilter_percent = round(mean(prevfilter_percent, na.rm = TRUE), 2)) %>% 
  arrange(desc(prevfilter_percent))

# List percentage removed per genus
as.data.frame(mean_prevfilter_by_genus)

# Check what had been removed for individual sample
single.sample.check <- prune_samples(sample_names(prevfilter.removed) == "AW.B1.G032_S33", prevfilter.removed)
sample_sums(single.sample.check) # list total sample sum
data.frame(
  genus = as.character(tax_table(single.sample.check)[names(sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:20]), "genus"]),
  Taxa_Sum = sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:20])

# Show prevfilter removed taxa 
prevfilter.removed.p <- tax_glom(prevfilter.removed,taxrank="phylum") # simplify to phylum rank to speed up figures
prevfilter.removed.p <- prune_taxa(taxa_sums(prevfilter.removed.p) > 100, prevfilter.removed.p) # remove minor stuff
prevfilter.removed.melt <- psmelt(prevfilter.removed.p)

prevfilter.lowrank <- ggplot(prevfilter.removed.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 50)) + ggtitle("prevfilter removed taxa (low stringent filtering)")

prevfilter.removed.plot <- ggplot(prevfilter.removed.melt, aes(x=Sample, y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity") + ggtitle("prevfilter removed taxa (low stringent filtering)") + theme(axis.text.x = element_blank()) + facet_wrap(~host_genus, scales="free_x")


# Compare data frames before after  
data.prevfilter.df <- as(sample_data(data.prevfilter),"data.frame")
data.prevfilter.df$sumsfilter <- sample_sums(data.prevfilter)
data.prevfilter.df$sumsfixed <- sample_sums(data.fixed)
data.prevfilter.df$sumsremoved <- sample_sums(prevfilter.removed)
data.prevfilter.df$percentremoved <- ((data.prevfilter.df$sumsremoved) / (data.prevfilter.df$sumsfixed)) * 100

prevfilter.percentremoved <- ggplot(data.prevfilter.df , aes(x=percentremoved, y=sumsfilter, color=host_family, size = sumsremoved)) +
  geom_point(alpha=0.7) + facet_wrap(~type) +
 ylim(0, 80000) + xlim(0, 50) + xlab("percent removed") + ylab("read counts") 


# Show prevalence of phyla after filtering
prev_results2 <- calc_prevalence(data.prevfilter, rank = "phylum")
prevdf2 <- subset(prev_results2$taxa_table, phylum %in% get_taxa_unique(data.prevfilter, "phylum"))

prevalence2 <- ggplot(prevdf2, aes(TotalAbundance, Prevalence / nsamples(data.prevfilter),color=phylum)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.015, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept = 60.0, alpha = 0.5, linetype = 2) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none") + ggtitle("data.prevfilter" ,
       "ASV cutoff (min total abundance = 60, min prevalence = 1.5%) ~0.001 percent "                ) 

data.fixed.prune		= prune_taxa(taxa_sums(data.fixed)>20, data.fixed)
data.fixed.phyla <- data.fixed.prune %>%
  comp_barplot(tax_level = "phylum", n_taxa = 30, merge_other = FALSE, sample_order = "asis",
               label = NULL ) +
  facet_wrap(vars(country), scales = "free") +  
  coord_flip() + ggtitle("data.fixed" ) +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

data.prevfilter.phyla <- data.prevfilter %>%
  comp_barplot(tax_level = "phylum", n_taxa = 30, merge_other = FALSE, sample_order = "asis",
               label = NULL ) +
  facet_wrap(vars(country), scales = "free") +  
  coord_flip() + ggtitle("data.prevfilter") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

# Boxplot phyla per sample type
#data.prevfilter.rel <- transform_sample_counts(data.prevfilter, function(x) 100*x/sum(x))
#data.prevfilter.rel.glom <- tax_glom(data.prevfilter.rel, taxrank = 'phylum')
#df.glom <- psmelt(data.prevfilter.rel.glom) # create data frame
#boxplot.phyla.per.sample.type <- ggplot(df.glom, aes(x=phylum, y=Abundance, fill=type)) + geom_boxplot()

# Visualize rare phyla
data.minorphyla = subset_taxa(data.fixed, phylum!="Proteobacteria" & phylum!= "Actinobacteria" & phylum!= "Bacteroidetes" & phylum!= "Firmicutes" & phylum!= "Tenericutes" )
data.minorphyla.p <- tax_glom(data.minorphyla,taxrank="phylum") # compress ps object
data.minorphyla.melt <- psmelt(data.minorphyla.p)

minorphyla <- ggplot(data.minorphyla.melt, aes(x=Sample, y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity", linewidth = 5) + ggtitle("samples with rare phyla") + theme(axis.text.x = element_blank()) + facet_wrap(~host_genus, scales="free_x")

minorphyla.lowrank <- ggplot(data.minorphyla.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 30)) + ggtitle("samples with rare phyla")

pdf(file.path(out_dir, "00_data_prevfilter.pdf"), width=12, height=6)
prevalence
prevalence_cutoff
ASV.histogram.prefilter
ASV.histogram.postfilter
prevfilter.lowrank
prevfilter.removed.plot
prevfilter.percentremoved
prevalence2 
data.fixed.phyla
data.prevfilter.phyla
minorphyla
minorphyla.lowrank
dev.off()

options(max.print=4000)
sink(file.path(out_dir, "00_data_prevfilter.txt"))
"data.fixed"
data.fixed
table(tax_table(data.fixed)[, "phylum"], exclude = NULL)
"data.prevfilter"
data.prevfilter
table(tax_table(data.prevfilter)[, "phylum"], exclude = NULL)
"prevfilter.removed"
prevfilter.removed
table(tax_table(prevfilter.removed)[, "phylum"], exclude = NULL)
"percent removed per sample [%]"
print(prevfilter.removed.frame.sorted)
cat("ASV filter percent of reads removed:", round(ASV_percentage_removed, 2), "%\n")
cat("Genus filter percent of reads removed:", round(genus_percentage_removed, 2), "%\n")
cat("Percent of reads removed total:", round(percentage_removed_total, 2), "%\n")
cat("Percent of reads removed sample:", round(percentage_removed_sample, 2), "%\n")
"mean_prevfilter_by_genus"
as.data.frame(mean_prevfilter_by_genus)
"top filtered ASVs"
prevfilter.removed.taxa.genus
sink()

# Cleanup pipeline
rm(prevfilter.removed.sample)
rm(data.fixed.prune500)
rm(data.fixed.prune)
rm(data.fixed.phyla)
rm(data.minorphyla.melt)

# 00 data.pruned1  ----------------
# Remove positive controls (mock community) or spike-in control taxa from the dataset

data.prevfilter
sample.prevfilter <- subset_samples(data.prevfilter, type=="sample")

# Identify Mock taxa spillover on ASV level
data.mocks <- subset_samples(data.prevfilter, type == "positive") # select mocks
sample_names(data.mocks)
sample_sums(data.mocks)

# Select top ASVs from data.mocks
mock.top <- names(sort(taxa_sums(data.mocks), decreasing = TRUE)[1:15])
# Show genus names of ASVs
mock.top.genus <- data.frame(
  #OTU = mock.top,
  Sum = taxa_sums(data.mocks)[mock.top],
  genus = as.character(tax_table(data.mocks)[mock.top, "genus"]))

mock.top.genus # Show taxonomy information for the top ASVs

# (optional) select specific ASVs that should be kept
keep_ASV <- c("ASV5", "ASV7", "ASV21")  # Keep those ASVs as they are not contamination 
mock.select <- setdiff(mock.top, keep_ASV)
#mock.select <- mock.top # if all can be used 

# Check what has been removed
spillover.removed1 <- prune_taxa(mock.select, data.prevfilter)
spillover.removed1.sample <- subset_samples(spillover.removed1, type=="sample")

# Should only remove mock community from positive controls
spillover.melt <- psmelt(spillover.removed1)

spillover.bar  <- ggplot(spillover.melt, aes(x=Sample, y=Abundance, fill=genus)) + geom_bar(#position="fill",
  stat = "identity", linewidth = 5) + ggtitle("pos control removal") + theme(axis.text.x = element_blank()) + facet_wrap(~subtype, scales="free_x")

spillover.lowrank <- ggplot(spillover.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=genus)) + geom_bar(#position="fill",
  stat = "identity") + ggtitle("pos control removal") + coord_flip(xlim = c(0, 30))

# Show spillover on plate frame
spillover.plate <- plate_frame(spillover.removed1)

# Remove spillover of mock taxa from dataset (top mock ASVs)
data.pruned1 <- prune_taxa(!(taxa_names(data.prevfilter) %in% mock.select), data.prevfilter)

data.prevfilter # pre-filter
data.pruned1 # mock ASVs removed


# Collect mock genera names for comp barplot
mock_genera <- spillover.removed1%>% tax_top(n = 10, rank = "genus")
# Count taxa number for comp barplot
mock.number = length(unique(spillover.melt$genus))

# Check MOCK removal
spillover.prefilter.comp <- data.prevfilter %>%
  comp_barplot(tax_level = "genus", n_taxa = mock.number, merge_other = T, sample_order = "asis",
               label = NULL, tax_order = mock_genera
  ) +  facet_wrap(vars(country), scales = "free") +  
  coord_flip() + ggtitle(
    "pos control pre-filter"  )

spillover.postfilter.comp <- data.pruned1 %>%
  comp_barplot(tax_level = "genus", n_taxa = mock.number, merge_other = T, sample_order = "asis",
               label = NULL, tax_order = mock_genera
  ) +   facet_wrap(vars(country), scales = "free") +  
  coord_flip() + ggtitle(
    "pos control post-filter"  )


pdf(file.path(out_dir, "00_data_pruned1_pos_control_removal.pdf"), width=12, height=6)
spillover.bar
spillover.lowrank
spillover.plate
spillover.prefilter.comp
spillover.postfilter.comp
dev.off()


## 00 data.pruned2  -----------------
# Positive control removal fine tuning (remove remaining ASVs from mock controls)

data.pruned1

mocks.leftover <- subset_samples(data.pruned1, type == "positive")
mocks.leftover.p = prune_taxa(taxa_sums(mocks.leftover)>3, mocks.leftover) # simplify graph
#plot_bar(mocks.leftover.p,x="genus", title="mock leftover") # select groups to be removed
list_of_mockgenera <- mocks.leftover.p %>%   tax_top(n = 30, rank = "genus")

# Taxa in ZymoBIOMICS™ Microbial Community Standard:
# Bacillus subtilis, Enterococcus faecalis, Escherichia coli, Lactobacillus fermentum -> Limosilactobacillus, Listeria monocytogenes, Pseudomonas aeruginosa, Salmonella enterica, Staphylococcus aureus
# Taxa often not fully resolved: Bacillales spc, Enterobacterales spc
# Taxa in ZymoBIOMICS™ Spike-in Control I: Imtechella halotolerans, Allobacillus halotolerans

table(tax_table(data.pruned1)[, "order"], exclude = NULL)
# Optional: Remove unresolved taxa Actinobacteria spc, Firmicutes spc, Proteobacteria spc


# Collect "list_of_mockgenera" names to filter out for fine tuning. Keep taxa from the samples e.g. Gilliamella 
mock.fine <- subset_taxa(mocks.leftover,
                           genus=="Allobacillus" | 
                           genus=="Actinobacteria spc" |
                           genus=="Bacillus" |
                           genus=="Bacillales spc"  |
                           #genus=="Enterobacterales spc"  |
                           #genus=="Enterobacteriaceae spc"  |
                           genus=="Firmicutes spc"| 
                           genus=="Imtechella" |
                           genus=="Limosilactobacillus" |
                           genus=="Listeria" | 
                           genus=="Proteobacteria spc" |
                           genus=="Pseudomonas" |
                           genus=="Pseudescherichia" |
                           genus=="Salmonella" |
                           genus=="Staphylococcus"  )

mock.fine.select <- plot_bar(mock.fine,x="genus", title="Mock genera fine tuning to be removed")

# (optional) select specific ASVs that should be kept (can be decided two steps later)
mock.top2 <- names(taxa_sums(mock.fine))
keep_ASV2 <- c("ASV96", "ASV1165" , "ASV17" , "ASV25" , "ASV52" , "ASV5378" , "ASV6" ,
               "ASV267" , "ASV163" , "ASV80", "ASV50", "ASV272" , "ASV62" , "ASV298") 
mock.taxa.fine <- setdiff(mock.top2, keep_ASV2)

## List samples with spillover
spillover.removed2 = prune_taxa(mock.taxa.fine, data.pruned1)
sort(data.frame(sum = sample_sums(spillover.removed2)), decreasing = FALSE)
sort(data.frame(sum = taxa_sums(spillover.removed2)), decreasing = FALSE)

percentage_spillover2_total <- (sum(sample_sums(spillover.removed2)) / 
                                  sum(sample_sums(data.pruned1))) * 100
# total percent removed
percentage_spillover2_total  #0.36

# Inspect specific genus
check.spillover2 <- subset_taxa(spillover.removed2, genus=="Pseudomonas" )
sort(data.frame(sum = sample_sums(check.spillover2 )), decreasing = FALSE)

# (optional) check individual samples if only mock gets removed, if not, put ASV manually in keep_ASV2
single.sample.check <- prune_samples(sample_names(spillover.removed2) == "AW.B2.P62_S215", spillover.removed2) # AW.B2.P80_S234 AW.B2.P82_S237
# single.sample.check <- subset_samples(spillover.removed2, study=="Peru") # look at entire projects  
sample_sums(single.sample.check) # list total sample sum
# Show genus names and taxa sums of top ASVs in a sample
data.frame(
  genus = as.character(tax_table(single.sample.check)[names(sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:10]), "genus"]),
  Taxa_Sum = sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:10])

spillover.removed2.g <- tax_glom(spillover.removed2,taxrank="genus") # speed up figures (and safe space) 
spillover2.melt <- psmelt(spillover.removed2.g)

spillover.bar2 <- ggplot(spillover2.melt, aes(x=Sample, y=Abundance, fill=genus)) + geom_bar(#position="fill",
  stat = "identity") + ggtitle("pos control removal finetuning") + theme(axis.text.x = element_blank()) + facet_wrap(~subtype, scales="free_x")

spillover.lowrank2 <- ggplot(spillover2.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=genus)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 30))

# Show spillover2 on plate frame
spillover.plate2 <- plate_frame(spillover.removed2)

# Remove spillover of mock taxa from dataset
data.pruned2 <- prune_taxa(!(taxa_names(data.pruned1) %in% mock.taxa.fine), data.pruned1)

# merge both spillover.removed
spillover.removed <- merge_phyloseq(spillover.removed1, spillover.removed2)
sort(data.frame(sum = sample_sums(spillover.removed)), decreasing = FALSE)
spillover.removed.sample <- subset_samples(spillover.removed, type=="sample")

# Calculate percentage of spillover removed
percentage_spillover_sample <- (sum(sample_sums(spillover.removed.sample)) / 
                                  sum(sample_sums(sample.prevfilter ))) * 100

percentage_spillover_total <- (sum(sample_sums(spillover.removed)) / 
                                 sum(sample_sums(data.prevfilter))) * 100

percentage_spillover_sample # 0.450%
percentage_spillover_total  # 2.049%

# Percentage of spillover removed per sample
spillover.removed.frame <- data.frame(
  spillover = sample_sums(spillover.removed),
  samplesum = sample_sums(data.prevfilter),
  host_genus = sample_data(spillover.removed)$host_genus,
  spillover_percent = (sample_sums(spillover.removed) / sample_sums(data.prevfilter)) * 100)

# Sort the dataframe in decreasing order
spillover.removed.frame.sorted <- spillover.removed.frame[order(spillover.removed.frame$spillover_percent), ]
print(spillover.removed.frame.sorted)

data.prevfilter # pre-filter
data.pruned1 # pos control removal
data.pruned2 # pos control removal fine tuning

# Collect mock genera2 names for comp barplot
mock_genera2 <- spillover.removed2%>% tax_top(n = 20, rank = "genus")
# Count taxa number for comp barplot
mock.number2 = length(unique(spillover2.melt$genus))

spillover.prefilter.comp2 <- data.pruned1 %>%
  comp_barplot(tax_level = "genus", n_taxa = mock.number2, merge_other = TRUE, sample_order = "asis",
               label = NULL, tax_order = mock_genera2 ) +
  facet_wrap(vars(country), scales = "free") +  
  coord_flip() + ggtitle("pos control pre-filter"  )

spillover.postfilter.comp2 <- data.pruned2 %>%
  comp_barplot(tax_level = "genus", n_taxa = 4, merge_other = TRUE, sample_order = "asis",
               label = NULL, tax_order = mock_genera2 ) +
  facet_wrap(vars(country), scales = "free") +  
  coord_flip() + ggtitle("pos control post-filter"  )

pdf(file.path(out_dir, "00_data_pruned2_pos_control_removal.pdf"), width=12, height=6)
mock.fine.select
spillover.bar2
spillover.lowrank2
spillover.plate2
spillover.prefilter.comp2
spillover.postfilter.comp2
dev.off()

# Cleanup pipeline
rm(spillover.bar2)
rm(spillover.lowrank2)
rm(spillover2.melt)


sink(file.path(out_dir, "00_data_pruned2_pos_control_removal.txt"))
"mock_genera"
mock_genera
"mock_genera2"
mock_genera2
"spillover.removed.frame.sorted"
print(spillover.removed.frame.sorted)
"data.prevfilter "
data.prevfilter 
"data.pruned2"
data.pruned2 
"percentage_spillover_sample"
percentage_spillover_sample
"percentage_spillover_total"
percentage_spillover_total
sink()


## 00 data.decontam  -------------------
# Apply decontam package on cleared controls
# More as a final control step, does not filter much  

data.pruned2
sample.pruned2 <- subset_samples(data.pruned2, type=="sample")

#### decontam.neg ----
# Identify ASVs as potential contamination using negative controls 
sample_data(data.pruned2)$is.neg <- sample_data(data.pruned2)$type == "negative"  
contam.neg.prev <- isContaminant(data.pruned2, method="prevalence", neg="is.neg", threshold=0.2) # 0.2
table(contam.neg.prev$contaminant) # adjust threshold 0.1 to 0.4
head(which(contam.neg.prev$contaminant))

data.pruned2.pa <- transform_sample_counts(data.pruned2, function(abund) 1*(abund>0))
data.pruned2.pa.neg <- prune_samples(sample_data(data.pruned2.pa)$type == "negative", data.pruned2.pa)
data.pruned2.pa.sam <- prune_samples(sample_data(data.pruned2.pa)$type == "sample", data.pruned2.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(data.pruned2.pa.sam), pa.neg=taxa_sums(data.pruned2.pa.neg),
                    contaminant=contam.neg.prev$contaminant)
decontam.neg <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# create phyloseq object with contaminant ASVs removed  
data.decontam.neg <- prune_taxa(!contam.neg.prev$contaminant, data.pruned2)



# Identify ASVs removed by decontam
decontam_asvs <- setdiff(taxa_names(data.pruned2), taxa_names(data.decontam.neg))
keep_decontam <- c( "ASV3" , "ASV193" , "ASV5871"  , "ASV225" ,  "ASV250" , "ASV1445", 
                    "ASV263" ,  "ASV266" ,  "ASV268",   "ASV274" ,  "ASV333" , 
                    "ASV334"  , "ASV3919",  "ASV457"  , "ASV62" ,   "ASV834" , "ASV49"  ) 

decontam_me_neg <- setdiff(decontam_asvs, keep_decontam )
#decontam_me_neg <- (decontam_asvs)
decontam.removed.neg <- prune_taxa(decontam_me_neg, data.pruned2) 

# Show top ASVs and genus names in decontam.removed.neg
data.frame(
  genus = as.character(tax_table(decontam.removed.neg)[names(sort(taxa_sums(decontam.removed.neg), decreasing = TRUE)[1:10]), "genus"]),
  Taxa_Sum = sort(taxa_sums(decontam.removed.neg), decreasing = TRUE)[1:10])

# check individual ASVs in dataset if contamination or not
# manually put in "keep_decontam" if no contamination
tag_ASV <- c("ASV49") 
singleASV <- prune_taxa(tag_ASV, data.pruned2)
sort(data.frame(sum = sample_sums(singleASV)), decreasing = FALSE)

data.frame(sort(sample_sums(decontam.removed.neg), decreasing = T)[1:10])
sum(taxa_sums(decontam.removed.neg))

single.sample.check <- prune_samples(sample_names(decontam.removed.neg) == "AW.B1.G032_S33", decontam.removed.neg)
data.frame(
  genus = as.character(tax_table(single.sample.check)[names(sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:10]), "genus"]),
  Taxa_Sum = sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:10])

# Show to be removed taxa
decontam.removed.neg.p <- tax_glom(decontam.removed.neg,taxrank="order") # speed up figures (and safe space)
decontam.removed.neg.melt <- psmelt(decontam.removed.neg.p)

decontam.removed.neg.lowrank <- ggplot(decontam.removed.neg.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=order)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 60)) + ggtitle("decontam neg")

decontam.removed.neg.bargraph <- ggplot(decontam.removed.neg.melt, aes(x=Sample, y=Abundance, fill=order)) + geom_bar(#position="fill",
  stat = "identity", linewidth = 5) + ggtitle("decontam neg") + theme(axis.text.x = element_blank()) + facet_wrap(~host_subfamily, scales="free_x")

# Show decontam.neg on plate frame
plate.decontam.neg  <- plate_frame(decontam.removed.neg)

#### decontam pos ----
# Identify ASVs as potential contamination using positive controls
sample_data(data.pruned2)$is.pos <- sample_data(data.pruned2)$type == "positive"  
contam.pos.prev <- isContaminant(data.pruned2, method="prevalence", neg="is.pos", threshold=0.1) #0.1 threshold=0.04
table(contam.pos.prev$contaminant)
head(which(contam.pos.prev$contaminant))

data.pruned2.pa.pos <- prune_samples(sample_data(data.pruned2.pa)$type == "positive", data.pruned2.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(data.pruned2.pa.sam), pa.neg=taxa_sums(data.pruned2.pa.pos),
                    contaminant=contam.pos.prev$contaminant)
decontam.pos <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Positive Controls)") + ylab("Prevalence (True Samples)")

# create phyloseq object with contaminant ASVs removed  
data.decontam.pos <- prune_taxa(!contam.pos.prev$contaminant, data.pruned2)

# Identify ASVs removed by decontam.pos
decontam.pos_asvs <- setdiff(taxa_names(data.pruned2), taxa_names(data.decontam.pos))
keep_decontam.pos <- c( "ASV85", "ASV334" , "ASV8304",  "ASV14859" , "ASV1473" , "ASV12739" ,
                        "ASV121", "ASV194", "ASV3076", "ASV12115","ASV8279","ASV2235" )   
#decontam_me_pos <- (decontam.pos_asvs)
decontam_me_pos <- setdiff(decontam.pos_asvs, keep_decontam.pos)
decontam.removed.pos <- prune_taxa(decontam_me_pos, data.pruned2)

# Top ASVs and genus names in decontam.removed.pos
data.frame(
  genus = as.character(tax_table(decontam.removed.pos)[names(sort(taxa_sums(decontam.removed.pos), decreasing = TRUE)[1:20]), "genus"]),
  Taxa_Sum = sort(taxa_sums(decontam.removed.pos), decreasing = TRUE)[1:20])
# manually put in "keep_decontam"

# Top decontam samples 
data.frame(sort(sample_sums(decontam.removed.pos), decreasing = T)[1:10])
sum(taxa_sums(decontam.removed.pos))

single.sample.check <- prune_samples(sample_names(data.pruned2) == "AW.B1.P003_S114", decontam.removed.pos) # AW.B1.G056_S58 AW.B1.G029_S29
data.frame(
  genus = as.character(tax_table(single.sample.check)[names(sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:10]), "genus"]),
  Taxa_Sum = sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:10])

# Show to be removed taxa
decontam.removed.pos.p <- tax_glom(decontam.removed.pos,taxrank="family") # speed up figures (and safe space)
decontam.removed.pos.melt <- psmelt(decontam.removed.pos.p)

decontam.removed.pos.lowrank <- ggplot(decontam.removed.pos.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=family)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 60)) + ggtitle("decontam pos")

decontam.removed.pos.bargraph <- ggplot(decontam.removed.pos.melt, aes(x=Sample, y=Abundance, fill=family)) + geom_bar(#position="fill",
  stat = "identity", linewidth = 5) + ggtitle("decontam pos") + theme(axis.text.x = element_blank()) + facet_wrap(~host_subfamily, scales="free_x")

# Show decontam.pos on plate frame
plate.decontam.pos <- plate_frame(decontam.removed.pos)

# Combine decontam neg and pos
decontam_me_combined <- union(decontam_me_neg, decontam_me_pos)

decontam.removed <- prune_taxa(decontam_me_combined, data.pruned2)

# Remove decontam taxa from dataset data.pruned3
data.decontam <- prune_taxa(!(taxa_names(data.pruned2) %in% decontam_me_combined), data.pruned2)

data.pruned2  # pos controls removed
data.decontam # decontam removed

# Relative amount removed by decontam
data.pruned2.rel <- transform_sample_counts(data.pruned2, function(x) x / sum(x))
decontam.removed.rel <- subset_taxa(data.pruned2.rel, taxa_names(data.pruned2.rel) %in% taxa_names(decontam.removed))
decontam.removed.rel <- tax_glom(decontam.removed.rel,taxrank="order") # compress
decontam.removed.rel.melt <- psmelt(decontam.removed.rel)

decontam.removed.rel.plot <- ggplot(decontam.removed.rel.melt, aes(x=Sample, y=Abundance, fill=order)) + geom_bar(#position="fill",
  stat = "identity", linewidth = 5) + ggtitle("percent decontam removed") + theme(axis.text.x = element_blank()) + facet_wrap(~host_subfamily, scales="free_x")

# Make dataframe with percent removed decontam per sample
decontam.removed.frame <- data.frame(
  decontam = sample_sums(decontam.removed),
  samplesum = sample_sums(data.pruned2),
  host_genus = sample_data(decontam.removed)$host_genus,
  type = sample_data(decontam.removed)$type,
  decontam_percent = (sample_sums(decontam.removed) / sample_sums(data.pruned2)) * 100 )

# Sort by decontam_percent
decontam.removed.frame.sorted <- decontam.removed.frame[
  order(decontam.removed.frame$decontam_percent), ]

print(decontam.removed.frame.sorted)

# Compare data frames before after
data.decontam.df <- as(sample_data(data.decontam),"data.frame")
data.decontam.df$sump2 <- sample_sums(data.pruned2)
data.decontam.df$sump3 <- sample_sums(data.decontam)
data.decontam.df$decontam <- sample_sums(decontam.removed)
data.decontam.df$percentdecontam <- ((data.decontam.df$decontam) / (data.decontam.df$sump2)) * 100

decontam.percentremoved <- ggplot(data.decontam.df, aes(x=percentdecontam, y=sump3, shape=type, color=host_family, size = decontam)) +
  geom_point(alpha=0.7) + ylim(0, 50000) + xlim(0, 75) + xlab("percent removed") + ylab("read counts") + facet_wrap(~type)

## Check what has been decontamed
decontam.removed.sample <- subset_samples(decontam.removed, type=="sample")

# Divide each count by the total count and multiply by 100 to get the percentage
percentage_removed_decontam_sample <- (sum(sample_sums(decontam.removed.sample)) / 
                                         sum(sample_sums(sample.pruned2))) * 100
percentage_removed_decontam_total <- (sum(sample_sums(decontam.removed)) / 
                                     sum(sample_sums(data.pruned2))) * 100

percentage_removed_decontam_sample # 0.133%
percentage_removed_decontam_total  # 0.156%

# Mean percent decontam by host genus
mean_decontam_by_genus <- decontam.removed.frame.sorted %>%
  mutate(decontam_percent = as.numeric(decontam_percent)) %>%
  group_by(host_genus) %>%
  summarise(decontam_percent = round(mean(decontam_percent, na.rm = TRUE), 2)) %>% 
  arrange(desc(decontam_percent))


as.data.frame(mean_decontam_by_genus)

# Show removed ASVs with genus name
top_decontam_removed <- names(sort(taxa_sums(decontam.removed), decreasing = TRUE)) # [1:50]

decontam.removed.top.ASVs <- data.frame(
  Abundance = taxa_sums(decontam.removed)[top_decontam_removed],
  Genus = tax_table(decontam.removed)[top_decontam_removed, "genus"])

decontam.removed.top.ASVs


# View the sorted dataframe
sink(file.path(out_dir, "00_data_pruned3_decontam.txt"))
print(decontam.removed.frame.sorted)
cat("Percent reads removed decontam total:", round(percentage_removed_decontam_total, 2), "%\n")
cat("Percent reads removed decontam sample:", round(percentage_removed_decontam_sample, 2), "%\n")
"Mean percent decontam by host genus"
as.data.frame(mean_decontam_by_genus)
"Show removed ASVs with genus name"
decontam.removed.top.ASVs
sink()

pdf(file.path(out_dir, "00_data_pruned3_decontam.pdf"), width=12, height=6)
decontam.neg
decontam.removed.neg.lowrank
decontam.removed.neg.bargraph 
decontam.pos
decontam.removed.pos.lowrank
decontam.removed.pos.bargraph
decontam.removed.rel.plot
decontam.percentremoved
dev.off()


### Check data.decontam
data.decontam
tail(tax_table(data.decontam))
table(tax_table(data.decontam)[, "phylum"], exclude = NULL)

## Validate sample names
data.decontam <- phyloseq_validate(data.decontam)


## 00 cleanup comparison  ------------------
# Overview of all filtering steps, will be exported later as csv file 

### Create filterframe to highlight all filtering steps
filterframe.df <- as(sample_data(data.decontam),"data.frame")
filterframe.df$A_data.comp <- sample_sums(data.comp.subset)
filterframe.df$B_data.bacteria <- sample_sums(data.fixed) # cyano removed
filterframe.df$C_data.prevfilter <- sample_sums(data.prevfilter) # prevfilter
filterframe.df$D_data.pruned <- sample_sums(data.pruned2) # spillover removed
filterframe.df$E_data.decontam <- sample_sums(data.decontam) # decontam removed

filterframe.df$B_Shannon <- estimate_richness(data.fixed, measures = "Shannon")$Shannon
filterframe.df$C_Shannon <- estimate_richness(data.prevfilter, measures = "Shannon")$Shannon
filterframe.df$D_Shannon <- estimate_richness(data.pruned2, measures = "Shannon")$Shannon
filterframe.df$E_Shannon <- estimate_richness(data.decontam, measures = "Shannon")$Shannon

filterframe.df$B_cyanoremoved <- sample_sums(cyano.removed.subset) 
filterframe.df$C_prevremoved <- sample_sums(prevfilter.removed)
filterframe.df$D_spillremoved <- sample_sums(spillover.removed)
filterframe.df$E_decontamremoved <- sample_sums(decontam.removed)

filterframe.df$B_cyanopercent <- ((filterframe.df$B_cyanoremoved) / (filterframe.df$A_data.comp)) * 100
filterframe.df$C_prevfilterpercent <- ((filterframe.df$C_prevremoved) / (filterframe.df$B_data.bacteria)) * 100
filterframe.df$D_prunedpercent <- ((filterframe.df$D_spillremoved) / (filterframe.df$C_data.prevfilter)) * 100
filterframe.df$E_decontampercent <- ((filterframe.df$E_decontamremoved) / (filterframe.df$D_data.pruned)) * 100

filterframe.df$filtersum <- (filterframe.df$A_data.comp - filterframe.df$E_data.decontam) /  filterframe.df$A_data.comp * 100

plot.cyano <- ggplot(filterframe.df, aes(x=B_Shannon, y=B_data.bacteria, size = B_cyanopercent , color=type))  + geom_point(alpha=0.7) + 
  ggtitle("Cyano / plants removed"  ) +  ylim(0, 75000)

plot.prevfilter <- ggplot(filterframe.df, aes(x=C_Shannon, y=C_data.prevfilter, size = C_prevfilterpercent , color=type))  + geom_point(alpha=0.7) + 
  ggtitle("Prevfilter (low stringent filtering)"  ) +   ylim(0, 75000)

plot.prune <- ggplot(filterframe.df, aes(x=D_Shannon, y=D_data.pruned, size = D_prunedpercent , color=type))  + geom_point(alpha=0.7) + 
  ggtitle("pos control removal"  ) +   ylim(0, 75000)

plot.decontam <- ggplot(filterframe.df, aes(x=E_Shannon, y=E_data.decontam, size = E_decontampercent , color=type))  + geom_point(alpha=0.7) + 
  ggtitle("decontam"  )  +   ylim(0, 75000)

# Sort by E_clean and create Rank
filterframe_long <- filterframe.df %>%
  rownames_to_column("Sample") %>%
  arrange(E_data.decontam) %>%
  mutate(Rank = row_number()) %>%
  pivot_longer(
    cols = c(A_data.comp, B_data.bacteria, C_data.prevfilter, D_data.pruned, E_data.decontam),
    names_to = "Stage",
    values_to = "ReadCount" )

# Library filter comparison
lib.comparison <- ggplot(data = filterframe_long, aes(x = Rank, y = ReadCount, color = subtype, shape = Stage)) +
  geom_point(size = 2) +
  geom_line(aes(group = Sample)) +
  geom_hline(yintercept = 2500, alpha = 0.5, linetype = 2) +
  #ylim(0, 50000) + xlim(0, 600) +
  coord_cartesian(ylim = c(0, 30000), xlim = c(0, 100)) + 
  ggtitle("Sample filter comparison") +
  labs(x = "Rank", y = "Read Count") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 20))


heatframe_long <- filterframe.df %>%
  rownames_to_column("Sample") %>%
  arrange(E_data.decontam) %>%
  mutate(
    Rank = row_number(),
    SampleLabel = paste0(Sample, " (", E_data.decontam, ")")
  ) %>%
  filter(type == "sample") %>%
  pivot_longer(
    cols = c(B_cyanopercent, C_prevfilterpercent, D_prunedpercent, E_decontampercent),
    names_to = "Variable",
    values_to = "Value"
  )

# Plot the heatmap
lib.heatmap <- ggplot(heatframe_long, aes(x = reorder(SampleLabel, Rank), y = Variable, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Sample filter comparison" ,  x = "Samples (sorted by sample sum)",
       y = "",
       fill = "% reads removed"  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  coord_flip(xlim = c(0, 60)) + geom_vline(xintercept = 10, linetype = 2) # mark last 10 samples


filterframe_sorted <- filterframe.df %>%
  arrange(filtersum) %>%
  select(E_data.decontam, B_cyanopercent, C_prevfilterpercent , D_prunedpercent,  E_decontampercent,  filtersum , host_genus) 
  
tail(filterframe_sorted, 50)

# View the sorted dataframe
sink(file.path(out_dir, "00_data_pruning_cleanup_comparison.txt"))
filterframe_sorted
sink()


# (optional) export sample names (filtersum >10%) to be pruned later
sample_names_to_prune <- filterframe.df %>%
  filter(type == "sample", filtersum > 10) %>%  
  rownames()   

sample_names_to_prune


pdf(file.path(out_dir, "00_data_pruning_cleanup_comparison.pdf"), width=12, height=6)
plot.cyano 
plot.prevfilter
plot.prune 
plot.decontam
lib.comparison
lib.heatmap 
dev.off()


### cleanup comparison diversity PCoA and sample sum

data.bacteria.rel <- transform_sample_counts(data.bacteria, function(x) x/sum(x))
data.bacteria.rel.PCoA <-  ordinate(data.bacteria.rel, method="PCoA", "bray")
data.bacteria.rel.PCoA.score <- cbind(
  as(sample_data(data.bacteria.rel), "data.frame"),
  as.data.frame(data.bacteria.rel.PCoA$vectors[, 1:3]) ) %>%
  mutate(logsum = log10(sample_sums(data.bacteria)))
data.bacteria.ordination <- ggplot(data.bacteria.rel.PCoA.score, aes(x = Axis.1, y = Axis.2, color = logsum, shape = type)) +
  geom_point(size = 6, alpha = 0.8) + scale_color_viridis_c(direction = -1) + theme_grid() + labs(title = "data.bacteria"  ) 
 
data.fixed.rel <- transform_sample_counts(data.fixed, function(x) x/sum(x))
data.fixed.rel.PCoA <-  ordinate(data.fixed.rel, method="PCoA", "bray")
data.fixed.rel.PCoA.score <- cbind(
  as(sample_data(data.fixed.rel), "data.frame"),
  as.data.frame(data.fixed.rel.PCoA$vectors[, 1:3]) ) %>%
  mutate(logsum = log10(sample_sums(data.fixed)))
data.fixed.ordination <- ggplot(data.fixed.rel.PCoA.score, aes(x = Axis.1, y = Axis.2, color = logsum, shape = type)) +
  geom_point(size = 6, alpha = 0.8) + scale_color_viridis_c(direction = -1) + theme_grid() + labs(title = "data.fixed"  ) 

data.prevfilter.rel <- transform_sample_counts(data.prevfilter, function(x) x/sum(x))
data.prevfilter.rel.PCoA <-  ordinate(data.prevfilter.rel, method="PCoA", "bray")
data.prevfilter.rel.PCoA.score <- cbind(
  as(sample_data(data.prevfilter.rel), "data.frame"),
  as.data.frame(data.prevfilter.rel.PCoA$vectors[, 1:3]) ) %>%
  mutate(logsum = log10(sample_sums(data.prevfilter)))
data.prevfilter.ordination <- ggplot(data.prevfilter.rel.PCoA.score, aes(x = Axis.1, y = Axis.2, color = logsum, shape = type)) +
  geom_point(size = 6, alpha = 0.8) + scale_color_viridis_c(direction = -1) + theme_grid() + labs(title = "data.prevfilter"  ) 

data.decontam.rel <- transform_sample_counts(data.decontam, function(x) x/sum(x))
data.decontam.rel.PCoA <-  ordinate(data.decontam.rel, method="PCoA", "bray")
data.decontam.rel.PCoA.score <- cbind(
  as(sample_data(data.decontam.rel), "data.frame"),
  as.data.frame(data.decontam.rel.PCoA$vectors[, 1:3]) ) %>%
  mutate(logsum = log10(sample_sums(data.decontam)))
data.decontam.ordination <- ggplot(data.decontam.rel.PCoA.score, aes(x = Axis.1, y = Axis.2, color = logsum, shape = type)) +
  geom_point(size = 6, alpha = 0.8) + scale_color_viridis_c(direction = -1) + theme_grid() + labs(title = "data.decontam"  ) 

data.high2000 = prune_samples(sample_sums(data.decontam)>=2000, data.decontam)
data.high2000.rel <- transform_sample_counts(data.high2000, function(x) x/sum(x))
data.high2000.rel.PCoA <-  ordinate(data.high2000.rel, method="PCoA", "bray")
data.high2000.rel.PCoA.score <- cbind(
  as(sample_data(data.high2000.rel), "data.frame"),
  as.data.frame(data.high2000.rel.PCoA$vectors[, 1:3]) ) %>%
  mutate(logsum = log10(sample_sums(data.high2000)))
data.high2000.ordination <- ggplot(data.high2000.rel.PCoA.score, aes(x = Axis.1, y = Axis.2, color = logsum, shape = type)) +
  geom_point(size = 6, alpha = 0.8) + scale_color_viridis_c(direction = -1) + theme_grid() + labs(title = "data.high2000"  ) 

data.high5000 = prune_samples(sample_sums(data.decontam)>=5000, data.decontam)
data.high5000.rel <- transform_sample_counts(data.high5000, function(x) x/sum(x))
data.high5000.rel.PCoA <-  ordinate(data.high5000.rel, method="PCoA", "bray")
data.high5000.rel.PCoA.score <- cbind(
  as(sample_data(data.high5000.rel), "data.frame"),
  as.data.frame(data.high5000.rel.PCoA$vectors[, 1:3]) ) %>%
  mutate(logsum = log10(sample_sums(data.high5000)))
data.high5000.ordination <- ggplot(data.high5000.rel.PCoA.score, aes(x = Axis.1, y = Axis.2, color = logsum, shape = type)) +
  geom_point(size = 6, alpha = 0.8) + scale_color_viridis_c(direction = -1) + theme_grid() + labs(title = "data.high5000"  ) 


# Alpha div

data.bacteria.rich  <- plot_richness(data.bacteria,x="host_subfamily", measures=c("Shannon","Observed")) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_point(size=4, aes(color=host_subfamily))  +geom_boxplot(aes(group = host_subfamily)) + ggtitle("data.bacteria")

data.fixed.rich  <- plot_richness(data.fixed,x="host_subfamily", measures=c("Shannon","Observed")) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_point(size=4, aes(color=host_subfamily))  +geom_boxplot(aes(group = host_subfamily)) + ggtitle("data.fixed")

data.prevfilter.rich  <- plot_richness(data.prevfilter,x="host_subfamily", measures=c("Shannon","Observed")) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_point(size=4, aes(color=host_subfamily))  +geom_boxplot(aes(group = host_subfamily)) + ggtitle("data.prevfilter")

data.decontam.rich  <- plot_richness(data.decontam,x="host_subfamily", measures=c("Shannon","Observed")) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_point(size=4, aes(color=host_subfamily))  +geom_boxplot(aes(group = host_subfamily)) + ggtitle("data.decontam")

data.high2000.rich  <- plot_richness(data.high2000,x="host_subfamily", measures=c("Shannon","Observed")) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_point(size=4, aes(color=host_subfamily))  +geom_boxplot(aes(group = host_subfamily)) +
  ggtitle("data.high2000") # + geom_label(aes(label = sampleID, color=host_tribe), size = 4) 

data.high5000.rich  <- plot_richness(data.high5000,x="host_subfamily", measures=c("Shannon","Observed")) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_point(size=4, aes(color=host_subfamily))  +geom_boxplot(aes(group = host_subfamily)) +
   ggtitle("data.high5000") # +geom_label(aes(label = sampleID, color=host_tribe), size = 4) 


pdf(file.path(out_dir, "00_data_pruning_cleanup_comparison_div.pdf"), width=12, height=6)
data.bacteria.ordination 
data.fixed.ordination 
data.prevfilter.ordination
data.decontam.ordination
data.high2000.ordination
data.high5000.ordination
data.bacteria.rich
data.fixed.rich
data.prevfilter.rich
data.decontam.rich
data.high2000.rich
data.high5000.rich
dev.off()


## 00 data.high  LT2000 cut-off ----------------------
# Set threshold and remove low throughput samples e.g. < LT2000

# Check sample abundance to find best cutoff LT2000 
sort(data.frame(sum = sample_sums(data.decontam)), decreasing = TRUE)

# Set cut-off LT2000
cutoff <- 2000

data.decontam.p <- tax_glom(data.decontam,taxrank="phylum") # speed up figure
data.decontam.p = prune_samples(sample_sums(data.decontam.p)<30000, data.decontam.p) # remove HT samples
clean.melt <- psmelt(data.decontam.p)

sample.sum.rank <-  ggplot(clean.melt, aes(x=reorder(Sample, Abundance), y=Abundance, fill=phylum)) +
  geom_bar(stat = "identity") +  geom_hline(yintercept=cutoff,linetype = 2) + #geom_vline(xintercept=33.5, linewidth=1) +
  coord_flip( )  + ggtitle( "Low throughput cut-off")

sample.sum.rank2 <-  ggplot(clean.melt, aes(x=reorder(Sample, Abundance), y=Abundance, fill=phylum)) +
  geom_bar(position="fill",  stat = "identity") +
  coord_flip()  + ggtitle( "Low throughput cut-off")


# sample sum vs shannon div vs Actinobacteria
head(data.decontam.df)
data.decontam.df$rich <- estimate_richness(data.decontam, measures=c("Observed", "Chao1", "Shannon", "Fisher"))
data.decontam.df$Actinobacteria <- sample_sums(subset_taxa(data.decontam.rel, phylum=="Actinobacteria" ))
data.decontam.df$Sphingomonas <- sample_sums(subset_taxa(data.decontam.rel, genus=="Sphingomonas" ))
data.decontam.df$Brevundimonas <- sample_sums(subset_taxa(data.decontam.rel, genus=="Brevundimonas" ))
data.decontam.df$Pseudomonas <- sample_sums(subset_taxa(data.decontam.rel, genus=="Pseudomonas" ))
data.decontam.df$Bacillus <- sample_sums(subset_taxa(data.decontam.rel, genus=="Bacillus" ))

# Sample sum vs Shannon diversity
plot.actino <- ggplot(data.decontam.df, aes(x=rich$Shannon, y=sump3, size = Actinobacteria, color=subtype))  + geom_point(alpha=0.7) + 
  geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) +  ggtitle(paste("LT:", cutoff)) + ylab("read counts")

plot.Brev <- ggplot(data.decontam.df, aes(x=rich$Shannon, y=sump3, size = Brevundimonas, color=subtype))  + geom_point(alpha=0.7) + 
  geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) +  ggtitle(paste("LT:", cutoff)) + ylab("read counts")

plot.Bacillus <- ggplot(data.decontam.df, aes(x=rich$Shannon, y=sump3, size = Bacillus, color=subtype))  + geom_point(alpha=0.7) + 
  geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) +  ggtitle(paste("LT:", cutoff)) + ylab("read counts")

sum.shannon.genus <- ggplot(data.decontam.df, aes(x=rich$Shannon, y=sump3, size = Sphingomonas, color=subtype, shape=type))  +
  geom_point(alpha=0.7)  +
  geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) +
  ggtitle(paste("LT:", cutoff)) + ylim(0, 10000)  +
  facet_wrap(~host_genus) + ylab("read counts")


# Remove LT2000 samples with defined cut-off e.g. 2000
data.high = prune_samples(sample_sums(data.decontam)>=cutoff, data.decontam)
data.low = prune_samples(sample_sums(data.decontam)<cutoff, data.decontam)
data.decontam # cleaned dataset
data.high # high throughput dataset
data.high.rel = transform_sample_counts(data.high, function(x) x/sum(x))

data.high.df <- as(sample_data(data.high),"data.frame")
data.high.df$samplesum <- sample_sums(data.high)
data.high.df$rich <- estimate_richness(data.high, measures=c("Observed", "Chao1", "Shannon", "Fisher"))
data.high.df$Sphingomonas <- sample_sums(subset_taxa(data.high.rel, genus=="Sphingomonas" ))
data.high.df$Bacillus <- sample_sums(subset_taxa(data.high.rel, genus=="Bacillus" ))

sumflop <- ggplot(data.decontam.df, aes(x=rich$Shannon, y=sump3, size = Sphingomonas, color=subtype))  + geom_point(alpha=0.7)  + ggtitle(
  "data.decontam"  ) + ylim(0, 75000) + xlim(0, 5) + geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) + ylab("read counts")
sumflophigh <- ggplot(data.high.df, aes(x=rich$Shannon, y=samplesum, size = Sphingomonas, color=subtype))  + geom_point(alpha=0.7) + ggtitle(
  "data.high (LT removal)"  ) + ylim(0, 75000) + xlim(0, 5) + geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) + ylab("read counts")


data.high.remainingsamples <- ggplot(data.high.df, aes(x=rich$Shannon, y=samplesum, color=host_genus))  +
  geom_point(alpha=0.7)  + ylim(0, 10000) + xlim(0, 5) +
  geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) +
  facet_wrap(~host_subfamily) + ggtitle("data.high") +
  geom_label(aes(label = sampleID), size = 4)

sampling.depth.per.genus <- ggplot(data.high.df, aes(x = fct_reorder(host_genus, -samplesum, .fun = median, na.rm = TRUE), y = samplesum)) +
  geom_boxplot() +
  geom_point(aes(color = host_subfamily), alpha = 0.5) +
  theme_line2()  + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +  labs(x="",y="sampling depth", title ="sampling depth per genus")


# Create a dataframe with sample sums and metadata
data.low.df <- data.frame(
  SampleSums = sample_sums(data.low),
  host_genus = sample_data(data.low)$host_genus,
  type = sample_data(data.low)$type
    )

# Sort the dataframe by SampleSums in decreasing order
sorted_data.low <- data.low.df[order(-data.low.df$SampleSums), ]

# Extract sample sums and add them to the phyloseq object sample data
sample_data(data.low)$samplesums <- sample_sums(data.low)

data.low.comp <- data.low %>%
  comp_barplot(tax_level = "genus", n_taxa = 30, merge_other = TRUE, sample_order = "bray",
               label = "samplesums") +
  facet_wrap(vars(host_genus), scales = "free") +
  coord_flip() +
  ggtitle(paste("data.low removed samples LT:", cutoff)) +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

LT.samples = prune_samples(sample_sums(data.high)<10000, data.high)
LT.samples # only LT samples
sample_data(LT.samples)$samplesums <- sample_sums(LT.samples)

data.high.low.comp <- LT.samples %>%
  #ps_filter(country=="Germany") %>%
  comp_barplot(tax_level = "genus", n_taxa = 30, merge_other = F, label = "samplesums", sample_order = "bray") +
  facet_wrap(vars(host_genus), scales = "free") +  
  coord_flip() + ggtitle(paste("data.high show sample reads between LT 10000 to LT",cutoff))

# Create a dataframe with sample sums and metadata
LT.samples.df <- data.frame(
  SampleSums = sample_sums(LT.samples),
  host_genus = sample_data(LT.samples)$host_genus,
  type = sample_data(LT.samples)$type )

# Sort the dataframe by SampleSums in decreasing order
sorted_LT.samples <- LT.samples.df[order(-LT.samples.df$SampleSums), ]

# samle sum before LT removal
genus_samplesum_beforeLT <- data.decontam %>%
  sample_data() %>%                # extract metadata
  data.frame() %>%                 # coerce to dataframe
  mutate(samplesum = sample_sums(data.decontam)) %>% 
  group_by(host_genus) %>%
  summarise(samplesum_beforeLT = round(mean(samplesum, na.rm = TRUE),0)) %>%
  arrange(desc(samplesum_beforeLT))

as.data.frame(genus_samplesum_beforeLT)

# Sample sum after LT removal
genus_samplesum_afterLT <- data.high %>%
  sample_data() %>%                # extract metadata
  data.frame() %>%                 # coerce to dataframe
  mutate(samplesum = sample_sums(data.high)) %>% 
  group_by(host_genus) %>%
  summarise(samplesum_afterLT = round(mean(samplesum, na.rm = TRUE),0)) %>%
  arrange(desc(samplesum_afterLT))

as.data.frame(genus_samplesum_afterLT)


# Show sample sums to identify best cut-off
sink(file.path(out_dir, "00_data_sample_cutoff_LT2000.txt"))
"Show sample sums to identify best cutoff"
print(sorted_LT.samples)
"cutoff"
cutoff
print(sorted_data.low)
"samplesum before LT2000"
as.data.frame(genus_samplesum_beforeLT)
"samplesum afterLT 2000"
as.data.frame(genus_samplesum_afterLT)
sink()

pdf(file.path(out_dir, "00_data_sample_cutoff_LT2000.pdf"), width=12, height=6)
sample.sum.rank
sample.sum.rank2
plot.actino
plot.Brev
plot.Bacillus
sum.shannon.genus
sumflop
sumflophigh
data.high.remainingsamples
sampling.depth.per.genus
data.low.comp 
data.high.low.comp
dev.off()

## sample.species (tax_glom genus)
## Multiple ASVs might represent the same bacterial genus, here they are collated 

data.high # ASV level including controls
sample.ASV <- subset_samples(data.high, type=="sample") #  ASV level only samples
sample.species <- tax_glom(sample.ASV,taxrank="genus") # tax_glom processing takes time
taxa_names(sample.species) <- tax_table(sample.species)[,"genus"]
sample.species  # Genus level only samples


## 00 (optional) sample.filter  -----------------------
# Optional low abundance filtering on genus level (only samples, no controls)
data.frame(sort(taxa_sums(sample.species), decreasing = T))

# Choose cutoff 
rel_ab_cutoff <- 0.035   # 0.035% results in 80 bacterial genera (good amount for phylo tree)
cutoff_min_total_abundance <- ((rel_ab_cutoff/100 ) * sum(sample_sums(sample.ASV))) # 0.035 percent cut-off
cutoff_min_total_abundance # value will be used to show cutoff in figure ,  equals ~2000 reads  
cutoff_min_percentage <- (cutoff_min_total_abundance*100)/sum(sample_sums(sample.ASV))
cutoff_min_percentage # equals 0.035 percentage

# Calculate min_prevalence e.g. 5% of samples
percent_of_samples<- 5 # set as 5% 
cutoff_min_prevalence <- (percent_of_samples*nsamples(sample.ASV))/100
cutoff_min_prevalence # equals about 8 samples --> set this number as min_prevalence
min_prev <- 8
# Calculated in % of samples
min_prev_cutoff <- (min_prev)/nsamples(sample.ASV)
min_prev_cutoff # value will be used to show cutoff in figure 

# Choose cutoff based on prevalence / abundance plot
prev_results3 <- calc_prevalence(sample.species, rank = "phylum")
prevdf3 <- subset(prev_results3$taxa_table, phylum %in% get_taxa_unique(sample.species, "phylum"))

sample.species.prevalence <- ggplot(prevdf3, aes(TotalAbundance, Prevalence / nsamples(sample.ASV),color=phylum)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = min_prev_cutoff, alpha = 0.5, linetype = 2) + # prev ca 5%
  geom_vline(xintercept = cutoff_min_total_abundance, alpha = 0.5, linetype = 2) + # abundance 2000 reads
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none") + 
  ggtitle(paste0("sample.species min_total_abundance cutoff ~ ", round(cutoff_min_total_abundance)))
  
# Filter on genus level >2000 reads = 80 taxa remaining
sample.g <- aggregate_taxa(sample.ASV, "genus")
data.frame(sort(taxa_sums(sample.g), decreasing = T))
sample.filter <- tax_filter(sample.ASV, tax_level = "genus", min_prevalence = min_prev, min_total_abundance = cutoff_min_total_abundance, min_sample_abundance = 1000 ) 
sample.filter.g <- aggregate_taxa(sample.filter, "genus")
data.frame(sort(taxa_sums(sample.filter.g), decreasing = T))

# percent of total reads remaining 
percent_retained <- (sum(otu_table(sample.filter)) / sum(otu_table(sample.ASV))) * 100
cat("Percent of reads retained:", round(percent_retained, 2), "%\n")

# Check what has been removed
genus.removed.taxa <- setdiff(taxa_names(sample.ASV), taxa_names(sample.filter))
genus.removed <- prune_taxa(genus.removed.taxa, sample.ASV)

# percentage of reads removed
percentage_removed_filter <- sum(sample_sums(genus.removed)) / sum(sample_sums(sample.ASV)) * 100
cat("Percent of reads removed:", round(percentage_removed_filter, 2), "%\n")

genus.removed.p <- tax_glom(genus.removed,taxrank="order")
sample.removed.melt <- psmelt(genus.removed.p)

sample.filter.lowrank <- ggplot(sample.removed.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=order)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 60)) + ggtitle("sample.filter")

sample.filter.bar <- ggplot(sample.removed.melt, aes(x=Sample, y=Abundance, fill=order)) + geom_bar(#position="fill",
  stat = "identity", linewidth = 5) + ggtitle("sample.filter") + theme(axis.text.x = element_blank()) + facet_wrap(~host_subfamily, scales="free_x")

# Compare data frames before after  
sample.filter.df <- as(sample_data(sample.filter),"data.frame")
sample.filter.df$sum_before <- sample_sums(sample.ASV)
sample.filter.df$sum_after <- sample_sums(sample.filter)
sample.filter.df$sumremoved <- sample_sums(genus.removed)
sample.filter.df$rich <- estimate_richness(sample.filter, measures=c("Shannon"))
sample.filter.df$percentfilter <- ((sample.filter.df$sumremoved) / (sample.filter.df$sum_before )) * 100

sample.filter.percent <- ggplot(sample.filter.df, aes(x=percentfilter , y=sum_after, shape=country, color=host_family, size = sumremoved)) +
  geom_point(alpha=0.7) + ylim(0, 50000) + xlim(0, 75) + xlab("percent removed filter") + ylab("read counts") + facet_wrap(~host_family)

# some samples might dropped below LT2000
sort(data.frame(sum = sample_sums(sample.filter)), decreasing = TRUE)

pdf(file.path(out_dir, "00_data_sample_filter_(optional).pdf"), width=12, height=6)
sample.species.prevalence
sample.filter.lowrank 
sample.filter.bar 
sample.filter.percent
dev.off()


sink(file.path(out_dir, "00_data_sample_filter_(optional).txt"))
"sample.filter"
sample.filter
"sample.filter.g"
sample.filter.g
"cutoff_min_percentage"
cutoff_min_percentage
"cutoff_min_total_abundance"
cutoff_min_total_abundance
cat("Percent of reads retained:", round(percent_retained, 2), "%\n")
cat("Percent of reads removed:", round(percentage_removed_filter, 2), "%\n")
sink()


### > cleanup pipeline / export sample data csv  -------------------
# Check filesize of all ps objects
# sapply(ls(), function(x) object.size(get(x))) %>% sort(decreasing = TRUE)
# Print size of ps objects
print(object.size(data.comp.subset), units = "auto")
print(object.size(sample.comp), units = "auto")
print(object.size(data.bacteria.rel), units = "auto")
print(object.size(data.bacteria), units = "auto")
print(object.size(data.fixed.rel), units = "auto")
print(object.size(sample.fixed), units = "auto")
print(object.size(prevfilter.removed ), units = "auto")
print(object.size(prevcheck.removed ), units = "auto")

# Clean pipeline delete large ps objects not needed anymore
rm(data.comp.subset)
rm(sample.comp)
rm(data.bacteria.rel)
rm(data.bacteria)
rm(data.fixed.rel)
rm(sample.fixed)
rm(prevfilter.removed)
rm(prevcheck.removed)
rm(data.prevfilter.rel)
rm(data.pruned2.rel)
rm(data.pruned2.pa)
gc()

# Export sample metadata as csv file
# sample_metadata <- as(sample_data(sample.ASV),"data.frame")
# ad column LT removed to filterframe.df
filterframe.df$LT_removed <- ifelse(filterframe.df$E_data.decontam < cutoff, "LT_remove", "HT_keep")
# select columns for metadata
colnames(filterframe.df)
sample_metadata <- filterframe.df[, c("chip", "sampleID", "collector", "host_order", "host_family", "host_subfamily", "host_tribe", "host_genus", "host_species", "BOLD_ID", "sex", "location", "sublocation", "country",
                                      "latitude" ,  "longitude"  , "date", "year", "kw", "elevation", "temp_week", "storage", "type", "PCR",
                                      "A_data.comp", "B_data.bacteria", "C_data.prevfilter", "D_data.pruned", "E_data.decontam", "LT_removed", "is.neg", "is.pos" , "B_cyanoremoved", "C_prevremoved", "D_spillremoved", "E_decontamremoved", "B_cyanopercent", "C_prevfilterpercent", "D_prunedpercent", "E_decontampercent", "filtersum", "Accession", "BioProject")]
#subset_metadata$TotalReads <- sample_sums(sample.ASV)
class(sample_metadata)
colnames(sample_metadata)
write.csv(sample_metadata, file.path(out_dir, "Suppl_table_sample_filter_metadata.csv"), row.names = T)


### > custom color palette -------
sample.species.df <- as(sample_data(sample.species),"data.frame")

# count group numbers for host palette
speciesCount = length(unique(sample.species.df$host_species))

# Set host palette for host colors
hostPalette = colorRampPalette(brewer.pal(12, "Set3")) # "Spectral" 

# Define specific fixed colors for host
# load with scale_color_manual(values = genus_colorfix) 
# load with scale_color_manual(values = hostPalette(speciesCount)) 

# Order samples manually by subfamily
table(sample_data(sample.species)$host_subfamily)
# order for figure
subfamily_ordered <- c("Satyrinae", "Dismorphiinae", "Pierinae", "Heliconiinae",  "Coliadinae", "Danainae", "Nymphalinae")
# order for color
subfamily_ordered_color <- c("Coliadinae" ,   "Dismorphiinae", "Heliconiinae" , "Nymphalinae"  , "Pierinae"  ,    "Satyrinae" , "Danainae" )

subfamily_names = sort(unique(sample.species.df$host_subfamily))
subfamily_count = length(subfamily_names)
subfamily_colorfix <- setNames(brewer.pal(n = 8, "Accent")[1:subfamily_count], subfamily_ordered_color) # max 8
#subfamily_colorfix <- setNames(colorRampPalette(brewer.pal(8, "Accent"))(subfamily_count), subfamily_ordered_color) # more than 8

tribe_names = sort(unique(sample.species.df$host_tribe))
tribe_count = length(tribe_names)
tribe_color = brewer.pal(n = 12, "Set3")[1:tribe_count]
tribe_color_sat = saturation(tribe_color, delta(+0.2))
tribe_colorfix <- setNames(tribe_color_sat, tribe_names)

genus_names = sort(unique(sample.species.df$host_genus))
genus_count = length(genus_names)
genus_color <- colorRampPalette(brewer.pal(n = 12, "Set3"))(genus_count)
genus_color_sat = saturation(genus_color, delta(+0.2))
genus_colorfix <- setNames(genus_color_sat, genus_names)

species_names = sort(unique(sample.species.df$host_species))
species_count = length(species_names)
species_colorfix <- setNames(colorRampPalette(brewer.pal(n = 12, "Set3"))(species_count), species_names)


country_names = sort(unique(sample.species.df$country))
country_count = length(country_names)
country_colorfix <- setNames(brewer.pal(n = 8, "Dark2")[1:country_count], country_names)


# Define specific fixed colors for microbial taxa
taxaPalette = colorRampPalette(brewer.pal(12, "Paired")) # main palette for core taxa
noncorePalette = colorRampPalette(brewer.pal(8, "Set1")) # secondary palette for non core taxa
orderPalette = colorRampPalette(brewer.pal(11, "Paired")) # main palette for order


## 01 sample comp all overview -------------------

# ASV level Transform to relative data
sample.ASV # Controls removed LT samples removed
sample.ASV.rel <- transform_sample_counts(sample.ASV, function(x) x/sum(x))
sample.filter # filtered from low abundant genera (but still ASV level)
sample.filter.rel <- transform_sample_counts(sample.filter, function(x) x/sum(x))

# Genus level Transform to relative data
sample.species # controls removed
sample.species.rel <- transform_sample_counts(sample.species, function(x) x/sum(x))

# Merge samples by subfamily
subfamily.merged = merge_samples(sample.species, "host_subfamily")
subfamily.merged.rel <- transform_sample_counts(subfamily.merged, function(x) x/sum(x))

# Aggregate taxa on family or order level for core analysis
sample.family <- aggregate_taxa(sample.species, "family")
sample.family.rel <- transform_sample_counts(sample.family, function(x) x/sum(x))
sample.order <- aggregate_taxa(sample.species, "order")
sample.order.rel <- transform_sample_counts(sample.order, function(x) x/sum(x))

# Subset country datasets
Peru.species <- subset_samples(sample.species, country=="Peru" )
Peru.species.rel <- subset_samples(sample.species.rel, country=="Peru" )
Germany.species <- subset_samples(sample.species, country=="Germany" )
Germany.species.rel <- subset_samples(sample.species.rel, country=="Germany" )

# Subset subfamilies
Satyrinae.rel <- subset_samples(sample.species.rel, host_subfamily=="Satyrinae")
Dismorphiinae.rel <- subset_samples(sample.species.rel, host_subfamily=="Dismorphiinae")
Pierinae.rel <- subset_samples(sample.species.rel, host_subfamily=="Pierinae")
Heliconiinae.rel <- subset_samples(sample.species.rel, host_subfamily=="Heliconiinae")
Coliadinae.rel <- subset_samples(sample.species.rel, host_subfamily=="Coliadinae")
Nymphalinae.rel <- subset_samples(sample.species.rel, host_subfamily=="Nymphalinae")

sample.order.comp.country <- sample.species %>%
  comp_barplot(tax_level = "order", n_taxa = 12, merge_other = T, label = "host_subfamily", sample_order = "bray") +
  facet_wrap(vars(country), scales = "free") +  
  coord_flip() + ggtitle( "sample.species")

# Comp barplot (merged) phylum
Sample.comp.merged.phylum <- sample.species%>%
  ps_select(host_genus,host_tribe,host_family,host_subfamily,type,country) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(type == "sample") %>%
  phyloseq::merge_samples(group = "host_subfamily") %>%
  comp_barplot(tax_level = "phylum", n_taxa = 5, merge_other = T, sample_order = subfamily_ordered , bar_width = 0.9) +
  labs(x = NULL, y = NULL) + ggtitle(""  ) +
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) #+ coord_flip()

# Comp barplot (merged) order
Sample.comp.merged.order <- sample.species%>%
  ps_select(host_genus,host_tribe,host_family,host_subfamily,type,country) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(type == "sample") %>%
  phyloseq::merge_samples(group = "host_subfamily") %>%
  comp_barplot(tax_level = "order", n_taxa = 10, merge_other = T, sample_order = subfamily_ordered , bar_width = 0.9) +
  labs(x = NULL, y = "rel. ab. [%]") + ggtitle(""  ) +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) 

# Comp barplot (merged) family
Sample.comp.merged.family <- sample.species%>%
  ps_select(host_genus,host_tribe,host_family,host_subfamily,type,country) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(type == "sample") %>%
  phyloseq::merge_samples(group = "host_subfamily") %>%
  comp_barplot(tax_level = "family", n_taxa = 12, merge_other = T, sample_order = subfamily_ordered, bar_width = 0.9) +
  labs(x = NULL, y = NULL) + ggtitle("Comp barplot merged family"  ) +
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) 

# Comp barplot (merged) genus
Sample.comp.merged.genus <- sample.species%>%
  ps_select(host_genus,host_tribe,host_family,host_subfamily,type,country) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(type == "sample") %>%
  phyloseq::merge_samples(group = "host_subfamily") %>%
  comp_barplot(tax_level = "genus", n_taxa = 12, merge_other = T, sample_order = subfamily_ordered, bar_width = 0.9) +
  labs(x = NULL, y = NULL) + ggtitle("Comp barplot merged genus"  ) +
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) )

pdf(file.path(out_dir, "01_sample_comp_all_overview.pdf"), width=12, height=6)
sample.order.comp.country
Sample.comp.merged.phylum
Sample.comp.merged.order
Sample.comp.merged.family
Sample.comp.merged.genus
dev.off()


## Overview about the main ps objects
sink(file.path(out_dir, "01_sample_comp_all_overview_median_sum.txt"))
"sample.ASV"
sample.ASV
"sample.species"
sample.species
"sample.filter"
sample.filter
"sample_sums(sample.species)"
sum(sample_sums(sample.species))
"median(sample_sums(sample.species))"
median(sample_sums(sample.species))
"mean(sample_sums(sample.species))"
mean(sample_sums(sample.species))
"host subfamily" 
table(sample_data(sample.species)$host_subfamily)
"host tribe" 
table(sample_data(sample.species)$host_tribe)
"host genus" 
table(sample_data(sample.species)$host_genus)
"host_genus number"
length(unique(sample_data(sample.species)$host_genus))
table(sample_data(sample.species)$host_species)
"host_species number"
length(unique(sample_data(sample.species)$host_species))
table(sample_data(sample.species)$country)
table(sample_data(sample.species)$location)
sink()




## 01 sample comp all rel abundance  -------------------

### Rel abundance all samples 
top_taxa_desc <- sort(rowMeans(otu_table(sample.species.rel)), decreasing = TRUE)[1:20]

relabundance <- ggplot(data.frame(genus = names(top_taxa_desc), mean_abundance = top_taxa_desc),
       aes(x = reorder(genus, -mean_abundance), y = mean_abundance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "rel. ab [%]") +
  theme_line2()


#### Top order + fam >1% mirrored  
#  merge by group select top order
sample.order.merged = merge_samples(sample.order, "host_subfamily")
data.frame(sort(taxa_sums(sample.order.merged), decreasing = F))
Top.order <- names(sort(taxa_sums(sample.order.merged), decreasing=T)[1:8])

# Subset Top order from sample.species dataset (to get genera within top order)
top_order_genera <- subset_taxa(sample.species.rel, order %in% Top.order)

percent_top_orders <- sum(sample_sums(top_order_genera))*100/nsamples(top_order_genera)
percent_top_orders # 87.84714 top 8 orders

# select top families (> 1% rel ab)
top_order_families <- tax_glom(top_order_genera, taxrank = "family")
taxa_names(top_order_families) <- tax_table(top_order_families)[,"family"]
data.frame(sort(taxa_sums(top_order_families)*100 / nsamples(top_order_families), decreasing = FALSE))
family_cutoff <- 0.01 * nsamples(top_order_families) # 0.01 is 1% rel ab
top_families <- taxa_names(top_order_families)[taxa_sums(top_order_families) > family_cutoff]

# select genera of top families
top.order.family <- prune_taxa(tax_table(top_order_genera)[, "family"] %in% top_families, top_order_genera)

percent_top_fam <- sum(sample_sums(top.order.family))*100/nsamples(top.order.family)
percent_top_fam # 86.72571 % top 8 order >1% family

# count for coloring
top_order_genera.melt <- psmelt(top.order.family)
orderCount = length(unique(top_order_genera.melt$order))
familyCount = length(unique(top_order_genera.melt$family))

# order by Top.order 
order_palette <- setNames(orderPalette(orderCount), Top.order)  # Assign colors in the abundance order
top_order_genera.melt$order <- factor(top_order_genera.melt$order, levels = Top.order) # order names by order

# normalize by country for rel ab
top_order_genera.melt.norm <- top_order_genera.melt %>%
  group_by(country) %>%
  mutate(Abundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

sum(top_order_genera.melt.norm$Abundance) # should equal group number

sample.top.tornado <- top_order_genera.melt.norm %>%
  group_by(country, order, family, genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(Abundance_mirrored = ifelse(country == "Peru", Abundance, -Abundance))

# Tornado top order / family min 1% rel abundance
order.fam.country.mirrored <- ggplot(sample.top.tornado , aes(family, Abundance_mirrored, fill= order)) +
  facet_grid(order~country, space = "free", scales = "free",switch = "y")+ # switch = "y" remove order names
  theme_grid()+
  geom_bar(stat="identity")+
  theme(strip.text.y.left = element_blank(),  # remove the order label
        axis.title.y = element_blank()) +    # remove "family" title
  scale_y_continuous(labels = function(x) scales::percent(abs(x))) + # show % as positive
  scale_fill_manual(values = order_palette) +
  coord_flip() + labs(y = "rel. ab. [%]")   
  
top_order_genera_agg <- top_order_genera.melt.norm %>%
  group_by(country, order, family, genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(Abundance_mirrored = ifelse(country == "Peru", Abundance, -Abundance))

# Tornado top order / family min 1% rel abundance / division by genus
order.fam.country.mirrored.genus <- ggplot(top_order_genera_agg, 
       aes(x = family, y = Abundance_mirrored, fill = order, group = genus)) +
  facet_grid(order ~ country, space = "free", scales = "free", switch = "y") +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3) +
  theme_grid() +
  theme(strip.text.y.left = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(labels = function(x) scales::percent(abs(x))) +
  coord_flip() +
  labs(y = "rel. ab. [%]") +
  scale_fill_manual(values = order_palette) 
  

# Total abundance of top order / families
top_order_total_stats <- top_order_genera.melt %>%
  group_by(Sample) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  summarise(
    mean_total_abundance = mean(total_abundance)*100,
    sd_total_abundance = sd(total_abundance)*100  )

top_order_total_stats # 86.7 
# crosscheck with previous calculated value
percent_top_fam       # 86.725

# taxa summary bacterial order (but here low abundant <1% families are missing)
top_order_family_stats <- top_order_genera.melt %>%
  group_by(order, Sample) %>%                      # keep Order grouping
  summarise(order_abundance = sum(Abundance), .groups = "drop") %>%
  group_by(order) %>%                              # then summarise per order
  summarise(
    mean_abundance = mean(order_abundance),
    sd_abundance = sd(order_abundance)
  )

top_order_family_stats # values lower than 

# taxa summary bacterial order for all samples
all_sample_order_stats <- as.data.frame(t(otu_table(sample.order.rel))) %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(-Sample, names_to = "order", values_to = "Abundance") %>%
  group_by(order) %>%
  summarise(mean_abundance = mean(Abundance)*100,
            sd_abundance   = sd(Abundance)*100,
            .groups = "drop") %>%
  slice_max(mean_abundance, n = 10)

all_sample_order_stats 

# taxa summary bacterial families for all samples
all_sample_family_stats <- as.data.frame(t(otu_table(sample.family.rel))) %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(-Sample, names_to = "family", values_to = "Abundance") %>%
  group_by(family) %>%
  summarise(mean_abundance = mean(Abundance)*100,
            sd_abundance   = sd(Abundance)*100,
            .groups = "drop") %>%
  slice_max(mean_abundance, n = 10)

all_sample_family_stats

sink(file.path(out_dir, "01_sample_comp_all_rel_abundance.txt"))
"sample.species"
table(tax_table(sample.species)[, "phylum"], exclude = NULL)
table(tax_table(sample.species)[, "order"], exclude = NULL)
"sample.species top 20" 
data.frame(sort(taxa_sums(sample.species), decreasing = T)[1:20])
"sample.species.rel"
data.frame(round(sort(taxa_sums(sample.species.rel)*100 / nsamples(sample.species.rel), decreasing = TRUE)[1:20],2 ))
"top 10 order"
as.data.frame(all_sample_order_stats)
"top 10 family"
as.data.frame(all_sample_family_stats)
sink()


pdf(file.path(out_dir, "01_sample_comp_all_rel_abundance.pdf"), width=12, height=6)
relabundance
order.fam.country.mirrored 
order.fam.country.mirrored.genus
dev.off()

# Cleanup pipeline
rm(top_order_genera.melt)
rm(top_order_genera.melt.norm)
rm(sample.top.tornado)

# 01 sample comp merged group barplot -------------

#### Group merged Top genus 
sample.species.merged = merge_samples(sample.species, "host_genus") # merge by host_genus or host_subfamily
sample.species.merged.rel <- transform_sample_counts(sample.species.merged, function(x) x/sum(x))
# Select top 10 genera
Top.genus <- names(sort(taxa_sums(sample.species.merged), decreasing=T)[1:12]) # top 10 taxa
# sample.species.top = prune_taxa(taxa_sums(sample.species.merged.rel)>0.20, sample.species.merged.rel)
sample.species.top <- subset_taxa(sample.species.merged.rel, taxa_names(sample.species.merged.rel)%in%Top.genus)
sample.species.merged.melt <- psmelt(sample.species.top)

# Calculate bacterial abundance (genus color order)
genus_abundance <- sample.species.merged.melt  %>%
  group_by(genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance))  # Sort by abundance (most abundant first)

# Color palette and color order
topgenusCount = length(unique(sample.species.merged.melt$genus))
topgenus_palette <- setNames(taxaPalette(topgenusCount), genus_abundance$genus)  # Assign colors in the abundance genus
sample.species.merged.melt$genus <- factor(sample.species.merged.melt$genus, levels = genus_abundance$genus) # genus names by order

sample.species.top.bar <- ggplot(sample.species.merged.melt,aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = TRUE), y=Abundance, fill = genus))+
  geom_bar(#position="fill",
    colour="black", linewidth=0.3, # remove black lines optional
    stat="identity") + 
  scale_fill_manual(values = topgenus_palette)+
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  theme(legend.title=element_text(size=11), legend.text=element_text(face="italic")) +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  labs(x="",y="rel. ab. [%]", fill = "top taxa")

### Top genera merge by group (host subfamily)
Top.genus <- names(sort(taxa_sums(sample.species.merged), decreasing=T)[1:50]) # select top taxa, only top 12 colored
abundant_group <- subset_taxa(subfamily.merged.rel, taxa_names(subfamily.merged.rel)%in%Top.genus)
group.melt <- psmelt(abundant_group)

# Color palette and color order for top 12 (other genera shown as NA)
group.melt$genus <- factor(group.melt$genus, levels = genus_abundance$genus) # genus names by order

sample.species.top.bar.group <- ggplot(group.melt,aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = TRUE), y=Abundance, fill = genus))+
  geom_bar( colour="black",
            stat="identity", linewidth=0.3)+
  scale_fill_manual(values = topgenus_palette)+
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1, face="italic") ) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  theme(legend.title=element_text(size=11), legend.text=element_text(size=10, face="italic")) + # legend text
  labs(x="",y="Relative Abundance [%]", fill = "top taxa") # scale_x_continuous(breaks=seq(0,6,1))

#### Group merged Top order 
sample.order.merged = merge_samples(sample.order, "host_subfamily")
sample.order.merged.rel <- transform_sample_counts(sample.order.merged, function(x) x/sum(x))
data.frame(sort(taxa_sums(sample.order.merged), decreasing = F))
Top.order <- names(sort(taxa_sums(sample.order.merged), decreasing=T)[1:8])
sample.order.top <- subset_taxa(sample.order.merged.rel, taxa_names(sample.order.merged.rel)%in%Top.order)
sample.order.top.melt <- psmelt(sample.order.top)
toporderCount = length(unique(sample.order.top.melt$order))

remaining_percent_order <- sum(sample_sums(sample.order.top))*100/nsamples(sample.order.top)
remaining_percent_order # 85.76097

# Set Sample as an ordered factor according to subfamily_ordered 
sample.order.top.melt$Sample <- factor(
  sample.order.top.melt$Sample,
  levels = subfamily_ordered)

# Calculate genus abundance for color order
order_abundance <- sample.order.top.melt  %>%
  group_by(order) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance))  # Sort by abundance (most abundant first)

# order bacterial taxa by Top.order same as order.fam.country.mirrored.genus
order_palette <- setNames(orderPalette(toporderCount), Top.order)  # Assign colors in the abundance order
sample.order.top.melt$order <- factor(sample.order.top.melt$order, levels = Top.order) # order names by order


sample.order.top.bar <- ggplot(sample.order.top.melt, aes(x = Sample, y = Abundance, fill = order)) +
  # ggplot(sample.order.top.melt,aes(x=fct_reorder(Sample, Abundance, .fun = sum, .desc = TRUE), y=Abundance, fill = order)) + # sort in descending order
  geom_bar(stat = "identity", colour = "black", linewidth = 0.1) + 
  scale_fill_manual(values = order_palette) +
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +
  labs(x = "", y = "Relative Abundance [%]")

# reverse sample order for coord_flip
sample.order.top.melt$Sample2 <- factor(sample.order.top.melt$Sample, levels = rev(c(subfamily_ordered)))

sample.order.top.bar.flip <- ggplot(sample.order.top.melt,aes(x = Sample2, y=Abundance, fill = order))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = order_palette)+
  theme_grid() + 
  scale_y_continuous(labels = scales::label_percent(scale = 100,  prefix = "", suffix = "" )) + 
  # guides(fill = guide_legend(reverse = TRUE)) + # inverse legend
  coord_flip() + labs(x="",y="rel. ab. [%]")


# Top 8 order abundance per host_subfamily
top_order_abundance <- sample.order.top.melt %>%
  group_by(Sample) %>%
  summarise(total_abundance = sum(Abundance)*100)

top_order_total_abundance <- top_order_abundance %>%
  summarise(
    mean_total_abundance = mean(total_abundance),
    sd_total_abundance = sd(total_abundance)   )

top_order_total_abundance # 85.8
remaining_percent_order   # 85.76

#### Group merged Top family 
sample.family.merged = merge_samples(sample.family, "host_subfamily")
sample.family.merged.rel <- transform_sample_counts(sample.family.merged, function(x) x/sum(x))
data.frame(sort(taxa_sums(sample.family.merged.rel)*100, decreasing = FALSE))

# select top 10 family
Top.family <- names(sort(taxa_sums(sample.family.merged), decreasing=T)[1:10]) 
sample.family.top <- subset_taxa(sample.family.merged.rel, taxa_names(sample.family.merged.rel)%in%Top.family)
sample.family.top.melt <- psmelt(sample.family.top)


# Calculate bacterial abundance for color order
family_abundance <- sample.family.top.melt  %>%
  group_by(family) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance))  # Sort by abundance (most abundant first)

topfamilyCount = length(unique(sample.family.top.melt$family))
family_palette <- setNames(taxaPalette(topfamilyCount), family_abundance$family)  # Assign colors in the abundance order
sample.family.top.melt$family <- factor(sample.family.top.melt$family, levels = family_abundance$family) # order names by order


sample.family.top.bar <- ggplot(sample.family.top.melt,aes(x=fct_reorder(Sample, Abundance, .fun = sum, .desc = T), y=Abundance, fill = family))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3)+
  scale_fill_manual(values = family_palette)+
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  #theme(legend.title=element_text(size=11), legend.text=element_text(size=10)) + # legend text
  labs(x="",y="Relative Abundance [%]") 

sample.family.top.bar.flip <- ggplot(sample.family.top.melt,aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = F), y=Abundance, fill = family))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = family_palette)+
  theme_grid() + #theme(axis.text.y = element_text(face="italic")) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  theme(legend.title=element_text(size=11)) + # legend.text=element_text(face="italic")
  #guides(fill = guide_legend(reverse = TRUE)) + # reverses bar fill
  coord_flip() + labs(x="",y="Relative Abundance [%]", fill="Top 12 families")

# Total abundance of top families per host_subfamily
top_family_abundance <- sample.family.top.melt %>%
  group_by(Sample) %>%
  summarise(total_abundance = sum(Abundance)*100)

top_family_total_abundance <- top_family_abundance %>%
  summarise(
    mean_total_abundance = mean(total_abundance),
    sd_total_abundance = sd(total_abundance)  )

pdf(file.path(out_dir, "01_sample_comp_merged_group.pdf"), width=12, height=6)
sample.species.top.bar
sample.species.top.bar.group
sample.order.top.bar
sample.order.top.bar.flip
sample.family.top.bar
sample.family.top.bar.flip
dev.off()

sink(file.path(out_dir, "01_sample_comp_merged_group.txt"))
"Top.order"
print(Top.order)
"top_order_total_abundance"
as.data.frame(top_order_total_abundance)
"top_order_abundance"
as.data.frame(top_order_abundance)
cat("\n")  # blank line
"Top.family"
Top.family 
"top_family_total_abundance"
as.data.frame(top_family_total_abundance)
"top_family_abundance"
as.data.frame(top_family_abundance)
sink()


## 01 sample core analysis -----------
rev_palette <- rev(brewer.pal(10, "RdYlBu"))
heat_core_palette <- rev_palette[-5] # Remove color 5

# ASV core
#ps1.rel <- microbiome::transform(sample.species, "compositional")
#Set different detection levels and prevalence
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- c(0.001, 0.003, 0.01, 0.05) # 0.1%, 0.5%, 1%, 5%

ASV.core <- plot_core(sample.ASV.rel, plot.type = "heatmap", 
                        colours = heat_core_palette, 
                        min.prevalence = 0.20, # 0.1 is 10%
                        prevalences = prevalences, 
                        detections = detections,
                        horizontal = F) +
  xlab("rel. ab. [%]") +theme_line2() +
  theme(axis.text.y= element_text(face="italic")) 
  #scale_fill_viridis_c(option = "inferno", name = "Prevalence") 

core.ASV <- core_members(sample.ASV.rel, detection = 0.01, prevalence = 20/100)

# Show genus names of ASVs
ASV.rel.tax_table <- tax_table(sample.ASV.rel)
core.ASV.genus <- as.character(ASV.rel.tax_table[core.ASV , "genus"])
core.ASV.as.genus <- data.frame(OTU = core.ASV , genus = core.ASV.genus)
core.ASV.as.genus

# Genus Core
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- c(0.002, 0.005, 0.01, 0.05) # 0.1%, 0.5%, 1%, 5%

genus.core <- plot_core(sample.species.rel, plot.type = "heatmap", 
                               colours = heat_core_palette, 
                               min.prevalence = 0.25, # 25%  
                               prevalences = prevalences, 
                               detections = detections,
                               horizontal = F               ) +
  xlab("rel. ab. [%]") +theme_line2() +
  theme(axis.text.y= element_text(face="italic")) 


core.data <- genus.core$data  
max_prev <- max(core.data$Prevalence, na.rm = TRUE)

# Change axis from 100% to max prevalence
genus.core.maxprev <- plot_core(sample.species.rel, plot.type = "heatmap", 
                        min.prevalence = 0.25, # 0.1 is 10% 
                        prevalences = prevalences, 
                        detections = detections,
                        horizontal = F               ) +
  xlab("rel. ab. [%]") +theme_line2() +
  scale_fill_gradientn(
    colours = heat_core_palette,
    limits = c(0, max_prev),
    oob = scales::squish,
    labels = scales::percent_format(accuracy = 1)  ) +
    theme(axis.text.y= element_text(face="italic")) 


# Reorder genus names based on prevalence at 1%
core_0.1 <- core.data %>% filter(DetectionThreshold == 0.01) #1% 0.01 
genus_order <- core_0.1 %>% arrange((Prevalence)) %>% pull(Taxa)

genus.core.order <- plot_core(sample.species.rel, plot.type = "heatmap", 
                         min.prevalence = 0.25, # 0.1 is 10%
                         prevalences = prevalences, 
                         detections = detections,
                         horizontal = F,
                         taxa.order = genus_order) +
  xlab("rel. ab. [%]") +theme_line2() +
  scale_fill_gradientn(
    colours = heat_core_palette,
    limits = c(0, max_prev),
    oob = scales::squish,
    labels = scales::percent_format(accuracy = 1) ) +
    theme(axis.text.y= element_text(face="italic")) 
    #theme(axis.text.x= element_text(size=8),axis.text.y= element_text(face="italic"),legend.text = element_text(size=8), legend.title= element_text(size=9))



### Core family analysis 
#Set different detection levels and prevalence
fam_prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
fam_detections <- c(0.002, 0.01, 0.05, 0.1) # 0.1%, 0.5%, 1%, 5% fam_detections <- c(0.001, 0.003, 0.01, 0.05)

family.core <- plot_core(sample.family.rel, plot.type = "heatmap", 
                         colours = heat_core_palette, 
                         min.prevalence = 0.3, # 0.1 is 10%
                         prevalences = fam_prevalences, 
                         detections = fam_detections,
                         horizontal = F) +
  xlab("rel. ab. [%]") +theme_line2() 
  

### Core order analysis 
#Set different detection levels and prevalence
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- c(0.002, 0.01, 0.05, 0.1) # 0.1%, 1%, 5%, 10%
order.core <- plot_core(sample.order.rel, plot.type = "heatmap", 
                         colours = heat_core_palette, 
                         min.prevalence = 0.25, # 0.1 is 10%
                         prevalences = prevalences, 
                         detections = detections,
                         horizontal = F) +
  xlab("rel. ab. [%]") +theme_line2() +
  theme(axis.text.x= element_text(size=8),legend.text = element_text(size=8), legend.title= element_text(size=9))


# Show Prevalence at 1% of top 20 taxa
prevalence(sample.species.rel, detection = 5/100, sort = TRUE)[1:20] # 5% high   ab; 10% prev
prevalence(sample.species.rel, detection = 1/100, sort = TRUE)[1:20] # 1% medium ab; 20% prev
prevalence(sample.species.rel, detection = 1/500, sort = TRUE)[1:20] # 0.2% low  ab; 40% prev

core.genus.high <- core_members(sample.species.rel, detection = 0.05, prevalence = 10/100)  # 11 taxa
core.genus      <- core_members(sample.species.rel, detection = 0.01, prevalence = 20/100)  # 12 taxa
core.genus.low  <- core_members(sample.species.rel, detection = 0.002, prevalence = 40/100) # 12 taxa

# Make a robust core by combining different definitions 
core_overlap <- list(high = core.genus.high, med = core.genus,
                    low = core.genus.low)
names(core_overlap) <- c("High abundance\nRA=5% Prev=10%", "Medium abundance\nRA=1% Prev=20%", "Low abundance\nRA=0.2% Prev=40%")

core_robust <- Reduce(intersect, core_overlap) # 9 taxa overlap from all three core definitions


#### Core Venn ------------
# Core genus per host_subfamily
Satyrinae.core <- core_members(Satyrinae.rel, detection = 0.01, prevalence = 20/100)
Dismorphiinae.core <- core_members(Dismorphiinae.rel, detection = 0.01, prevalence = 20/100)
Pierinae.core <- core_members(Pierinae.rel, detection = 0.01, prevalence = 20/100)
Heliconiinae.core <- core_members(Heliconiinae.rel, detection = 0.01, prevalence = 20/100)
Coliadinae.core <- core_members(Coliadinae.rel, detection = 0.01, prevalence = 20/100)
Nymphalinae.core <- core_members(Nymphalinae.rel, detection = 0.01, prevalence = 20/100)


core_genus_lists <- list(Satyrinae = Satyrinae.core, Dismorphiinae = Dismorphiinae.core,
                         Pierinae = Pierinae.core,  Heliconiinae = Heliconiinae.core,
                         Coliadinae = Coliadinae.core, Nymphalinae = Nymphalinae.core)

core_genus_lists2 <- list(Satyrinae = Satyrinae.core, Dismorphiinae = Dismorphiinae.core,
                          Pierinae = Pierinae.core,  Heliconiinae = Heliconiinae.core,
                          Coliadinae = Coliadinae.core)

core_genus_lists3 <- list(    Pierinae = Pierinae.core,  Heliconiinae = Heliconiinae.core,
                          Coliadinae = Coliadinae.core)

#### ggVennDiagram
# Six subfamilies
venn_core_plot <-  ggVennDiagram(core_genus_lists, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "none") +
  labs(title = "Core Venn")

core_genus_lists_overlap <- Reduce(intersect, core_genus_lists)

# Five subfamilies
venn_core_plot2 <-  ggVennDiagram(core_genus_lists2, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "none") +
  labs(title = "Core Venn")

core_genus_lists2_overlap <- Reduce(intersect, core_genus_lists2)

# Three subfamilies
venn_core_plot3 <-  ggVennDiagram(core_genus_lists3, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "none") +
  labs(title = "Core Venn")

core_genus_lists3_overlap <- Reduce(intersect, core_genus_lists3)

# Overlap of three core definitions
venn_core_overlap <-  ggVennDiagram(core_overlap, label = "count", label_alpha = 0, set_size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "none") +
  #labs(title = "Robust Core Definitions") +
  scale_x_continuous(expand = expansion(mult = 0.22)) +  # Make Venn smaller
  scale_y_continuous(expand = expansion(mult = 0.1))

# Modify line thickness in venn diagram
venn_core_overlap$layers[[2]]$aes_params$linewidth <- 0.2

# core_overlap between definitions
# high & medium = Entomomonas & Apibacter
# medium & low  = Acinetobacter
# only low      = Carnobacterium & Raoultella


#### Core boxplot subfamily-----

# Show percent of core genus per group
core.genus.rel <- prune_taxa(core.genus, sample.species.rel)
# same as: core.genus.rel2 <- core(sample.species.rel, detection = 0.01, prevalence = 20/100)
# same as: core.genus.rel3 <- subset_taxa(sample.species.rel, taxa_names(sample.species.rel)%in%core.genus)
# sample_sums(core.genus.rel) # core percentage per sample
tax_table(core.genus.rel)

# core genus ASV names (for phylo tree) 
core.genus.ASV.names <- taxa_names(sample.ASV)[tax_table(sample.ASV)[, "genus"] %in% core.genus]
#core.genus.ASVs <- prune_taxa(core.genus.ASV.names, sample.ASV)

# Sum relative abundance of core taxa per sample
core.genus.rel.df <- core.genus.rel %>%
  psmelt() %>%  
  group_by(Sample, host_genus, host_family, host_subfamily,host_tribe, country,elevation,kw,location, sublocation,sex,no,weight) %>%
  summarise(genus_core = sum(Abundance)) 

# Reorder host_genus based on mean Core_Abundance
core.genus.rel.df <- core.genus.rel.df %>%
  group_by(host_genus) %>%
  mutate(mediangenus_core = median(genus_core))  # calculate the mean core abundance

# Boxplot showing core relative abundance per host genus
core.genus.boxplot <- ggplot(core.genus.rel.df, aes(x = fct_reorder(host_genus, mediangenus_core), y = genus_core)) +
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 4, aes(color = host_subfamily, shape=country), alpha = 0.8) +
  theme_line2() + theme(axis.text.y = element_text(face="italic") ) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  scale_colour_manual(values=subfamily_colorfix) +
  coord_flip() + labs(x="",y="genus core [%]", color = "host subfamily") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )
  
# Boxplot showing core relative abundance per host subfamily
core.genus.boxplot.subfam <- ggplot(core.genus.rel.df, aes(x = fct_reorder(host_subfamily, mediangenus_core), y = genus_core)) +
  #geom_boxplot() +
  geom_violin(draw_quantiles = c(0.5), trim = T,  scale = "width") +
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 3, aes(color = host_subfamily, shape =country), alpha = 0.6) +
  theme_line2() + #theme(axis.text.y = element_text(face="italic") ) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  scale_colour_manual(values=subfamily_colorfix) +
  coord_flip() + labs(x="",y="genus core [%]", color = "host subfamily")


#### Core genus barplot (grouped by genus)
species.genus = merge_samples(sample.species, "host_genus") # merged by host_genus 
species.genus.rel <- transform_sample_counts(species.genus, function(x) x/sum(x))
species.genus.core <- subset_taxa(species.genus.rel, taxa_names(species.genus.rel)%in%core.genus)
species.genus.core.melt <- psmelt(species.genus.core)

# Calculate genus abundance for color order
core_genus_abundance <- species.genus.core.melt %>%
  group_by(genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance))  # Sort by abundance (most abundant first)

coregenusCount = length(unique(species.genus.core.melt$genus))
# core_palette definition
core_palette <- setNames(taxaPalette(coregenusCount), core_genus_abundance$genus)  # Assign colors in the abundance order
core_palette_alpha <- scales::alpha(core_palette, 0.8) # more pale
species.genus.core.melt$genus <- factor(species.genus.core.melt$genus, levels = core_genus_abundance$genus) # order genus names by abundance

core.genus.abundance <- ggplot(species.genus.core.melt,aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = F), y=Abundance, fill = genus))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = core_palette_alpha) +
  theme_grid() + theme(axis.text.y = element_text(face="italic")) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  coord_flip() + labs(x="",y="genus core [%]", fill="genus core") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))


### Core genus merged by host_subfamily  
subfamily.merged.core <- subset_taxa(subfamily.merged.rel, taxa_names(subfamily.merged.rel)%in%core.genus)
subfamily.merged.core.melt <- psmelt(subfamily.merged.core)

subfamily.merged.core.melt$genus <- factor(subfamily.merged.core.melt$genus, levels = core_genus_abundance$genus) # order genus names by abundance

subfamily.merged.core.melt$host_subfamily_ordered <- factor(subfamily.merged.core.melt$Sample, levels = rev(c(subfamily_ordered)))

core.genus.abundance.subfam <- ggplot(subfamily.merged.core.melt,aes(x = host_subfamily_ordered, y=Abundance, fill = genus))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3
    ) + 
  scale_fill_manual(values = core_palette_alpha) +
  theme_line2() + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""),limits = c(0, 1) ) + 
  coord_flip() + labs(x="",y="genus core [%]") 
  

core.genus.abundance.subfam2 <- ggplot(subfamily.merged.core.melt,aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = TRUE), y=Abundance, fill = genus))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = core_palette) +
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  labs(x="",y="genus core rel. abundance [%]")


# Total abundance of core genus per host_subfamily
genus_core_subfam_stats <- subfamily.merged.core.melt %>%
  group_by(Sample) %>%
  summarise(total_abundance = sum(Abundance))

genus_core_subfam_stats

subfam_stats <- genus_core_subfam_stats %>%
  summarise(
    mean_total_abundance = mean(total_abundance),
    sd_total_abundance = sd(total_abundance)
  )

#### Core sample  -----
core.melt<- psmelt(core.genus.rel)

# Calculate genus abundance for color order, or use core_genus_abundance$genus
#core_genus_order <- core.melt%>%
#  group_by(genus) %>%
#  summarise(total_abundance = sum(Abundance)) %>%
#  arrange(desc(total_abundance))  

sample.core.abundance <- ggplot(core.melt %>%
                                mutate(genus = factor(genus, levels = core_genus_abundance$genus))  # order genus names by abundance for legend 
                                ,aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = FALSE), y = Abundance, fill = genus)) +
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = core_palette) +
  theme_grid() + theme(axis.text.y = element_text(face="italic")) + 
  coord_flip() + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""),limits = c(0, 1)) +
  labs(x="",y="core [%]", fill="core") +
  facet_wrap(~country*host_subfamily, scales = "free") +
  scale_x_discrete(labels = setNames(core.melt$host_genus, core.melt$Sample)) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))

# Total abundance of core genus across all samples
genus_core_sample_stats <- core.melt %>%
  group_by(Sample) %>%
  summarise(total_abundance = sum(Abundance))  %>%
  summarise(
    mean_total_abundance = mean(total_abundance),
    sd_total_abundance = sd(total_abundance)
  )

#### Core family boxplot --------

# Show Prevalence at 1% of top 20 taxa
prevalence(sample.family.rel, detection = 1/100, sort = TRUE)[1:20]

# core family definition detection 1.0% (0.01)  prevalence 20% (20/100),  
core.family <- core_members(sample.family.rel, detection = 0.01, prevalence = 20/100) # 13 families

# core family on genus level 
#core.family.genus.names <- taxa_names(sample.species.rel)[tax_table(sample.species.rel)[, "family"] %in% core.family]
#core.family.genus <- prune_taxa(core.family.genus.names, sample.species.rel)
core.family.genus  <- subset_taxa(sample.species.rel, family %in% core.family) # same but simpler


# core family on ASV level
core.family.ASV.names <- taxa_names(sample.ASV)[tax_table(sample.ASV)[, "family"] %in% core.family]
core.family.ASVs <- prune_taxa(core.family.ASV.names, sample.ASV)

# Sum relative abundance of core taxa per sample
core.family.genus.df <- core.family.genus %>%
  psmelt() %>%  # convert phyloseq object to a dataframe
  group_by(Sample, host_genus, host_family, host_subfamily, host_tribe, country) %>%
  summarise(family_core = sum(Abundance)) 

# Reorder host_genus based on mean Core_Abundance
core.family.genus.df <- core.family.genus.df %>%
  group_by(host_genus) %>%
  mutate(meanfamily_core = mean(family_core))  # calculate the mean core abundance

# Boxplot showing core relative abundance per host genus
core.family.boxplot <- ggplot(core.family.genus.df, aes(x = fct_reorder(host_genus, meanfamily_core), y = family_core)) +
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 4, aes(color = host_subfamily), alpha = 0.8) +
  theme_line2() + theme(axis.text.y = element_text(face="italic") ) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  scale_colour_manual(values=subfamily_colorfix) +
  coord_flip() + labs(x="",y="family core [%]", color = "host subfamily")

core.family.boxplot.subfam <- ggplot(core.family.genus.df, aes(x = fct_reorder(host_subfamily, meanfamily_core), y = family_core)) +
  geom_violin( draw_quantiles = c(0.5),  trim = TRUE,
  #adjust = 1.5,         # Smoother curves (optional)
  scale = "width" ) +      # Keep widths consistent across groups 
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 4, aes(color = host_subfamily), alpha = 0.8) +
  theme_line2() +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  scale_colour_manual(values=subfamily_colorfix) +
  coord_flip() + labs(x="",y="family core [%]", color = "host subfamily")

# Core family genus list
data.frame(
  AvgRelAbundance = sort(taxa_sums(core.family.genus) / nsamples(core.family.genus) * 100, decreasing = TRUE))

core.family.genus.melt <- psmelt(core.family.genus)

# normalize across all samples (cum abundance to average rel abundance)
core.family.genus.avg <- core.family.genus.melt %>%
  group_by(genus, family, order, country, host_subfamily) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()

# Cleanup pieline
rm(core.family.genus.melt)

core.family.genus.relabundance_color <- ggplot(core.family.genus.avg, aes(x = mean_abundance, y = genus, fill = order)) +
  facet_grid(order*family~country, scales = "free_y", space = "free") +
  geom_bar(stat="identity")+
  #theme(strip.text.y.left = element_text(angle = 0))+
  #theme(axis.text.x = element_text(angle = 60, hjust = 1),axis.title.x = element_text(family = "sans", size = 15)) + 
  scale_fill_manual(values = order_palette) +
  scale_x_continuous(labels = scales::percent) +
  theme_grid() +  labs(x = "rel. ab.", y = "")  

core.family.genus.relabundance.fam <- ggplot(core.family.genus.avg, aes(x = mean_abundance, y = genus, fill = order)) +
  facet_grid(order*family~host_subfamily, scales = "free_y", space = "free") +
  geom_bar(stat="identity")+
  #theme(strip.text.y.left = element_text(angle = 0))+
  #theme(axis.text.x = element_text(angle = 60, hjust = 1),axis.title.x = element_text(family = "sans", size = 15)) + 
  scale_fill_manual(values = order_palette) +
  scale_x_continuous(labels = scales::percent) +
  theme_grid() +  labs(x = "rel. ab.", y = "") 

# core family per sample
sort(data.frame(sum = sample_sums(core.family.genus)), decreasing = T)

# Core family total rel abundance all samples
core_fam_relab <- mean(sample_sums(core.family.genus)) * 100 # percent of core family

# Core family total rel abundance with SD
family_core_sample_stats <- core.family.genus.df %>%
  group_by(Sample) %>%
  summarise(total_abundance = sum(family_core)) %>%
  summarise(
    mean_total_abundance = mean(total_abundance),
    sd_total_abundance = sd(total_abundance)   )


#### Core family barplot
### Core family merged by host_genus 
sample.family.merged.genus = merge_samples(sample.family, "host_genus")
sample.family.merged.genus.rel <- transform_sample_counts(sample.family.merged.genus, function(x) x/sum(x))
sample.family.merged.genus.core <- subset_taxa(sample.family.merged.genus.rel, taxa_names(sample.family.merged.genus.rel)%in%core.family)
sample.family.merged.genus.core.melt <- psmelt(sample.family.merged.genus.core)

# Calculate genus abundance for color order
fam_abundance <- sample.family.merged.genus.core.melt %>%
  group_by(family) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance))  # Sort by abundance (most abundant first)

corefamilyCount = length(unique(core.family))
fam_core_palette <- setNames(taxaPalette(corefamilyCount), fam_abundance$family)  # Assign colors in the abundance order
sample.family.merged.genus.core.melt$family <- factor(sample.family.merged.genus.core.melt$family, levels = fam_abundance$family) # order names by abundance 

core.family.abundance <- ggplot(sample.family.merged.genus.core.melt,aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = F), y=Abundance, fill = family))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = fam_core_palette) +
  theme_grid() + theme(axis.text.y = element_text(face="italic")) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  #scale_colour_manual(values=subfamily_color) +
  coord_flip() + labs(x="",y="rel. ab [%]", fill ="core families")


### Core family merged by host_subfamily 
sample.family.subfamily = merge_samples(sample.family, "host_subfamily")
sample.subfamily.rel <- transform_sample_counts(sample.family.subfamily, function(x) x/sum(x))
sample.subfamily.core <- subset_taxa(sample.subfamily.rel, taxa_names(sample.subfamily.rel)%in%core.family)
sample.subfamily.core.melt <- psmelt(sample.subfamily.core)

sample.subfamily.core.melt$family <- factor(sample.subfamily.core.melt$family, levels = fam_abundance$family) # order names by abundance 
# rev(fam_abundance$family)) to reverse order of stacked core families, include guides(fill = guide_legend(reverse = TRUE)) # reverse legend order

core.family.abundance.subfam <- ggplot(sample.subfamily.core.melt,aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = F), y=Abundance, fill = family))+
  geom_bar(#position="fill",F
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = fam_core_palette) +
  theme_grid() +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  coord_flip() + labs(x="",y="family core [%]") 


#### Core family group 'nonCore' in gray ----
# Melt dataset & classify families as core or non-core
sample.subfamily.melt <- psmelt(sample.subfamily.rel) %>%
  mutate(family_colored = ifelse(family %in% core.family, family, "nonCore")) %>%
  group_by(Sample) %>%
  mutate(core_abundance = sum(ifelse(family %in% core.family, Abundance, 0))) %>%  # Sum only core abundance
  ungroup()

# Set factor levels: Reverse for stacking, but keep "nonCore" last in legend
sample.subfamily.melt$family_colored <- factor(sample.subfamily.melt$family_colored, 
                                               levels = rev(c(fam_abundance$family, "nonCore")))

# Create stacked core family plot
core.family.abundance.others <- ggplot(sample.subfamily.melt, 
                                       aes(x = fct_reorder(Sample, core_abundance),  
                                           y = Abundance, 
                                           fill = family_colored)) +
  geom_bar(colour = "black", stat = "identity", linewidth = 0.1) + 
  scale_fill_manual(values = c(fam_core_palette, "nonCore" = "grey70")) +  # Directly define palette
  theme_line2() + 
  theme(axis.text.y = element_text(face="italic")) +  
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +  
  coord_flip() + 
  labs(x = "", y = "core family rel. abundance [%]", fill = "core families") +
  guides(fill = guide_legend(reverse = TRUE))  # Keep legend order correct



  
# Total abundance of top families per host_subfamily
family_core_subfam_stats <- sample.subfamily.core.melt %>%
    group_by(Sample) %>%
    summarise(total_abundance = sum(Abundance)) %>%
    summarise(
      mean_total_abundance = mean(total_abundance),
      sd_total_abundance = sd(total_abundance)
    )
  
family_core_subfam_stats_noNymph <- sample.subfamily.core.melt %>%
    filter(Sample != "Nymphalinae") %>%  # remove one subfamily
    group_by(Sample) %>%
    summarise(total_abundance = sum(Abundance)) %>%
    summarise(
      mean_total_abundance = mean(total_abundance),
      sd_total_abundance = sd(total_abundance)
    )
  

ggarrange.genus.core.subfamily <- ggarrange(core.genus.boxplot.subfam, core.genus.abundance.subfam, labels = c('a', 'b'),
                                  common.legend = F, legend = "right", ncol = 2,   nrow = 1, align = "h" )

ggarrange.fam.core.subfamily <- ggarrange(core.family.boxplot.subfam,  core.family.abundance.subfam, labels = c('a', 'b'),
                                common.legend = F, legend = "right", ncol = 2,   nrow = 1, align = "h" )

pdf(file.path(out_dir, "01_sample_core_analysis.pdf"), width=12, height=6)
venn_core_plot 
venn_core_plot2
venn_core_overlap
genus.core
core.genus.boxplot
core.genus.abundance
ggarrange.genus.core.subfamily
sample.core.abundance
family.core
core.family.boxplot
core.family.abundance
core.family.abundance.others
ggarrange.fam.core.subfamily
dev.off()

sink(file.path(out_dir, "01_sample_core_analysis.txt"))
"genus_core_sample_stats"
as.data.frame(genus_core_sample_stats)
"genus_core_subfam_stats" 
as.data.frame(genus_core_subfam_stats)
"family_core_sample_stats"
as.data.frame(family_core_sample_stats)
"family_core_subfam_stats"
as.data.frame(family_core_subfam_stats)
"family_core_subfam_stats_noNymph"
as.data.frame(family_core_subfam_stats_noNymph)
sink()



### > non-core genus group minor by 'other'-------------
# Get non-core taxa
noncore.samples <- subset_taxa(sample.species.rel, !taxa_names(sample.species.rel) %in% core.genus)
data.frame(sort(taxa_sums(noncore.samples), decreasing = TRUE)[1:30])

top_noncore_relab <- (taxa_sums(noncore.samples)*100 /  nsamples(noncore.samples)) %>%
  sort(decreasing = TRUE) %>%
  head(30) %>%
  as.data.frame()

top_noncore_relab # Select number of top non core taxa to show

noncore.samples.df <- psmelt(noncore.samples)

# Compute top taxa
top_taxa_samples <- noncore.samples.df %>%
  group_by(genus) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 16) %>%
  pull(genus)

# Label and ordering, remaining groups as 'others'
noncore.samples.df <- noncore.samples.df %>%
  mutate(label = if_else(genus %in% top_taxa_samples, as.character(genus), "Other"),
         label = fct_reorder(label, Abundance, .fun = sum, .desc = TRUE), # sort taxa by total abundance  
         label = fct_relevel(label, "Other", after = Inf)) # move "Other" to end

# noncore palette colors
noncore.samples_labels <- setdiff(unique(noncore.samples.df$label), "Other")
palette_colors <- noncorePalette(length(noncore.samples_labels))
noncore.samplescolors <- c(setNames(palette_colors, noncore.samples_labels), Other = "grey70")

single.sample.noncore.abundance <- ggplot(noncore.samples.df, aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = TRUE), y = Abundance, fill = label) ) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = noncore.samplescolors) +
  coord_flip() +  scale_y_reverse(labels = scales::label_percent(scale = 100), limits = c(1, 0)) +
  facet_wrap(~country * host_subfamily, scales = "free") +
  scale_x_discrete(labels = setNames(noncore.samples.df$host_genus,  noncore.samples.df$Sample)) +
  theme_grid() + theme(axis.text.y = element_text(face = "italic"),
                       legend.title = element_text(size = 10),
                       legend.text = element_text(size = 8),
                       legend.key.size = unit(0.5, "cm")) + labs(x = "", y = "non-core [%]", fill = "non-core")


# non-core taxa grouped by host_genus
noncore.genus <- subset_taxa(species.genus.rel, !taxa_names(species.genus.rel) %in% core.genus)
data.frame(sort(taxa_sums(noncore.genus), decreasing = TRUE)[1:30])

noncore.genus.df <- psmelt(noncore.genus)

# Compute top taxa
top_taxa_genus <- noncore.genus.df %>%
  group_by(genus) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 17) %>% # Color top taxa, rest as other
  pull(genus)

# Label and ordering, remaining groups as 'others'
noncore.genus.df <- noncore.genus.df %>%
  mutate(label = if_else(genus %in% top_taxa_genus, as.character(genus), "Other"),
         label = fct_reorder(label, Abundance, .fun = sum, .desc = TRUE), # sort taxa by total abundance  
         label = fct_relevel(label, "Other", after = Inf)) # move "Other" to end


noncore_labels <- setdiff(unique(noncore.genus.df$label), "Other")
#palette_colors <- taxaPalette(length(noncore_labels)+4)[5:(length(noncore_labels) + 4)] #  select from color 4 on
palette_colors <- noncorePalette(length(noncore_labels))
noncoregenuscolors <- c(setNames(palette_colors, noncore_labels), Other = "grey70")

# non core genus abundance by host_genus
noncore.genus.abundance <-  ggplot(noncore.genus.df,aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = TRUE), y=Abundance, fill = label))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = noncoregenuscolors) +
  theme_grid() + theme(axis.text.y = element_text(face="italic")) + 
  coord_flip() + scale_y_reverse(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""),limits = c(1, 0)) +
  labs(x="",y="non-core [%]", fill="non-core") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))


# Color based on total sample abundance (not host genus abundance)
# use 'noncore.samples_labels' list
noncoregenuscolors.sample <- c(setNames(noncorePalette(length(noncore.samples_labels)), noncore.samples_labels), Other = "grey70")

# New Label by 'noncore.samples_labels' list
noncore.genus.df <- noncore.genus.df %>%
  mutate(label2 = if_else(genus %in% noncore.samples_labels, as.character(genus), "Other"),
         label2 = fct_reorder(label2, Abundance, .fun = sum, .desc = TRUE), # sort taxa by total abundance  
         label2 = fct_relevel(label2, "Other", after = Inf)) # move "Other" to end

# Color based on total sample abundance (does not color biggest bars)
noncore.genus.abundance.samplecolor <- ggplot(noncore.genus.df,aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = TRUE), y=Abundance, fill = label2))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = noncoregenuscolors.sample) +
  theme_grid() + theme(axis.text.y = element_text(face="italic")) + 
  coord_flip() + scale_y_reverse(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""),limits = c(1, 0)) +
  labs(x="",y="non-core [%]", fill="non-core") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))


# non-core genus by subfamily 
noncore.subfam <- subset_taxa(subfamily.merged.rel, !taxa_names(subfamily.merged.rel) %in% core.genus)
noncore.subfam.df <- psmelt(noncore.subfam) 
data.frame(sort(taxa_sums(noncore.subfam), decreasing = TRUE)[1:30])

# Compute top taxa
top_taxa_subfam <- noncore.subfam.df  %>%
  group_by(genus) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 12) %>%
  pull(genus)

# Label and ordering, remaining groups as 'others'
noncore.subfam.df  <- noncore.subfam.df  %>%
  mutate(label = if_else(genus %in% top_taxa_subfam, as.character(genus), "Other"),
         label = fct_reorder(label, Abundance, .fun = sum, .desc = TRUE), # sort taxa by total abundance  
         label = fct_relevel(label, "Other", after = Inf)) # move "Other" to end

all_labels <- unique(noncore.subfam.df$label)
noncore_labels <- setdiff(all_labels, "Other")
palette_colors <- taxaPalette(length(noncore_labels))
noncorecolors <- c(setNames(palette_colors, noncore_labels), Other = "grey70")

noncore.genus.abundance.subfam <- ggplot(noncore.subfam.df,aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = TRUE), y=Abundance, fill = label))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = noncorecolors) +
  theme_grid() + #theme(axis.text.y = element_text(face="italic")) + 
  coord_flip() + scale_y_reverse(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +
  labs(x="",y="non core [%]",fill="non-core")

#### > non-core family group minor by 'other' -------------
noncore.genus.fam <- aggregate_taxa(noncore.genus, "family")
data.frame(sort(taxa_sums(noncore.genus.fam), decreasing = TRUE)[1:30])
noncore.genus.fam.df <- psmelt(noncore.genus.fam)

# Compute top taxa family level
top_taxa_fam <- noncore.genus.fam.df  %>%
  group_by(family) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 10) %>%
  pull(family)

# Label and ordering, remaining groups as 'others'
noncore.genus.fam.df  <- noncore.genus.fam.df  %>%
  mutate(label = if_else(family %in% top_taxa_fam, as.character(family), "Other"),
         label = fct_reorder(label, Abundance, .fun = sum, .desc = TRUE), # sort taxa by total abundance  
         label = fct_relevel(label, "Other", after = Inf)) # move "Other" to end

noncore_fam <- setdiff(unique(noncore.genus.fam.df$label), "Other")
fam_colors <- noncorePalette(length(noncore_fam))
noncorefamcolors <- c(setNames(fam_colors, noncore_fam), Other = "grey70")

# non core genera on family level
noncore.family.abundance <- ggplot(noncore.genus.fam.df,aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = TRUE), y=Abundance, fill = label))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = noncorefamcolors) +
  theme_grid() + theme(axis.text.y = element_text(face="italic")) + 
  coord_flip() + scale_y_reverse(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""),limits = c(1, 0)) +
  labs(x="",y="non-core [%]", fill="non-core") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))



# non core family grouped by subfamily
noncore.subfam.family <- aggregate_taxa(noncore.subfam, "family")
data.frame(sort(taxa_sums(noncore.subfam.family), decreasing = TRUE)[1:30])
noncore.subfam.family.df <- psmelt(noncore.subfam.family)

# Compute top taxa family level
top_subfam_fam <- noncore.subfam.family.df  %>%
  group_by(family) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 10) %>%
  pull(family)

# Label and ordering, remaining groups as 'others'
noncore.subfam.family.df  <- noncore.subfam.family.df  %>%
  mutate(label = if_else(family %in% top_subfam_fam, as.character(family), "Other"),
         label = fct_reorder(label, Abundance, .fun = sum, .desc = TRUE), # sort taxa by total abundance  
         label = fct_relevel(label, "Other", after = Inf)) # move "Other" to end


noncore_fam <- setdiff(unique(noncore.subfam.family.df$label), "Other")
fam_colors <- noncorePalette(length(noncore_fam))
noncorefamcolors <- c(setNames(fam_colors, noncore_fam), Other = "grey70")

noncore.family.subfamily <- ggplot(noncore.subfam.family.df,aes(x = fct_reorder(Sample, Abundance, .fun = sum, .desc = TRUE), y=Abundance, fill = label))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = noncorefamcolors) +
  theme_grid() + #theme(axis.text.y = element_text(face="italic")) + 
  coord_flip() + scale_y_reverse(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""),limits = c(1, 0)) +
  labs(x="",y="non-core [%]", fill="non-core taxa") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))

pdf(file.path(out_dir, "01_sample_core_noncore.pdf"), width=12, height=6)
single.sample.noncore.abundance 
noncore.genus.abundance
noncore.genus.abundance.samplecolor
noncore.genus.abundance.subfam
noncore.family.abundance
noncore.family.subfamily
dev.off()

rm(single.sample.noncore.abundance)
rm(noncore.samples.df)


## 02 sample div alpha Shannon --------------
 
ASV.rich  <- plot_richness(sample.ASV,x="host_subfamily", measures=c("Shannon","Observed","InvSimpson")) +
  geom_boxplot(aes(group = host_subfamily)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  ggtitle("sample.ASV") + scale_colour_manual(values=subfamily_colorfix) 

filter.rich  <- plot_richness(sample.filter,x="host_subfamily", measures=c("Shannon","Observed","InvSimpson")) +
  geom_boxplot(aes(group = host_subfamily)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  ggtitle("sample.filter") + scale_colour_manual(values=subfamily_colorfix) 

species.rich  <- plot_richness(sample.species,x="host_subfamily", measures=c("Shannon","Observed","InvSimpson")) +
  geom_boxplot(aes(group = host_subfamily)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  ggtitle("sample.species") + scale_colour_manual(values=subfamily_colorfix) 

family.rich  <- plot_richness(sample.family,x="host_subfamily", measures=c("Shannon","Observed","InvSimpson")) +
  geom_boxplot(aes(group = host_subfamily)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  ggtitle("sample.family") + scale_colour_manual(values=subfamily_colorfix) 

country.rich  <- plot_richness(sample.species,x="country", measures=c("Shannon","Observed","InvSimpson")) +
  geom_boxplot(aes(group = country)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  ggtitle("sample.species country") + scale_colour_manual(values=subfamily_colorfix) 

location.rich  <- plot_richness(sample.species,x="location", measures=c("Shannon","Observed","InvSimpson")) +
  geom_boxplot(aes(group = location)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  ggtitle("sample.species location") + scale_color_manual(values = subfamily_colorfix) 

sublocation.country.rich  <- plot_richness(sample.species,x="sublocation", measures=c("Shannon")) +
  geom_boxplot(aes(group = sublocation)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  ggtitle("sample.species country sublocation") + scale_color_manual(values = subfamily_colorfix) +
  facet_wrap(~ country, scales = "free_x")

genus.country.rich  <- plot_richness(sample.species,x="host_genus", measures=c("Shannon")) +
  geom_boxplot(aes(group = host_genus)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  ggtitle("sample.species country host") + scale_color_manual(values = subfamily_colorfix) +
  facet_wrap(~ country, scales = "free_x")

sample.rich.time  <- plot_richness(sample.species,x="kw", measures=c("Shannon","Observed","InvSimpson")) +
  theme_line2()+ geom_point(size=4, aes(color=host_subfamily)) + 
  stat_smooth(method="lm", color="black", linewidth=0.5,se=T) + 
  stat_regline_equation(label.y = 0.0, aes(label = after_stat(rr.label))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  ggtitle("sample.species time") + scale_color_manual(values = subfamily_colorfix) 

pdf(file.path(out_dir, "02_sample_div_alpha_Shannon.pdf"), width=12, height=6)
ASV.rich
filter.rich 
species.rich
family.rich
country.rich
location.rich
sublocation.country.rich
genus.country.rich 
sample.rich.time
dev.off()


#### alphaframe Shannon --------------
df.sample <- as(sample_data(sample.species), "data.frame") %>%
  rownames_to_column("Sample")

df.core <- data.frame(
  Sample = sample_names(core.genus.rel),
  genus_core = sample_sums(core.genus.rel))

df.alpha <- estimate_richness(sample.species, measures = c("Shannon", "Observed", "InvSimpson")) %>%
  rownames_to_column("Sample")

### New alphaframe simplify data frame for statistics 
alphaframe <- df.sample %>%
  left_join(df.core, by = "Sample") %>%
  left_join(df.alpha, by = "Sample")

alphaframe$hill_numbers <- exp(alphaframe$Shannon)
alphaframe$time <- alphaframe$kw
alphaframe$temperature <- alphaframe$temp_week
alphaframe$host_subfamily_ordered <- factor(alphaframe$host_subfamily, levels = subfamily_ordered)
#alphaframe$group <- alphaframe$host_subfamily # if there is no group defined


# adjust variables if necessary
rownames(alphaframe) <- alphaframe$Sample
rownames(alphaframe)
head(alphaframe)
class(alphaframe)
str(alphaframe) # variables can be converted as.numeric or as.factor
# use time point as.factor only for figure color
# use time point as integer for statistics
#alphaframe$host_genus <- as.factor(alphaframe$host_genus)
#alphaframe$kw <- as.numeric(alphaframe$kw)
#alphaframe$group <- factor(df.sample$group, levels = c("t0", "t1", "t2", "t3"))

# Quick correlations
pairs(alphaframe[, c("Shannon",  "Observed", "InvSimpson", "temperature",  "time")])

alphaframe.noAglais <- subset(alphaframe , host_genus != "Aglais")
alphaframe.Aglais <- subset(alphaframe , host_genus == "Aglais")

# alpha figures
alpha.shannon  <- ggplot(alphaframe, aes(x=host_subfamily_ordered, y=Shannon)) +
  geom_boxplot(aes(group = host_subfamily))+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_subfamily, shape = country), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="",y="Shannon Index", color = "host subfamily") +
  facet_wrap(~ host_family, scales = "free_x")
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

alpha.rich  <- ggplot(alphaframe, aes(x=host_subfamily_ordered, y=Observed)) +
  geom_boxplot(aes(group = host_subfamily))+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_subfamily, shape = country), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="",y="microbial richness (q=0)") 

alpha.rich.tribe  <- ggplot(alphaframe, aes(x=host_subfamily_ordered, y=Observed)) +
  geom_boxplot(aes(group = host_subfamily))+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_tribe, shape = country), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=tribe_colorfix) +
  labs(x="",y="microbial richness (q=0)") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

alpha.rich.tribe.sort  <- ggplot(alphaframe, aes(x = fct_reorder(host_tribe, Observed, .fun = median, na.rm = TRUE), y=Observed)) +
  geom_boxplot(aes(group = host_tribe))+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_subfamily, shape = country), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="",y="microbial richness (q=0)") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

lm.alpha.rich <- lm(Observed ~   host_tribe  , data = alphaframe) # + host_tribe host_subfamily+
summary(lm.alpha.rich) # R-squared:  0.60 
anova(lm.alpha.rich)

alpha.shannon2  <- ggplot(alphaframe, aes(x=Shannon, y=host_subfamily_ordered)) +
  #geom_violin(aes(group = host_subfamily), draw_quantiles = c(0.5), trim = F,  scale = "width")+
  geom_boxplot(aes(group = host_subfamily))+
  geom_point(position = position_jitter(h = 0.1, w = 0), size=4, aes(color=host_subfamily, shape = country), alpha = 0.8) +
  theme_line2()+ scale_y_discrete(limits = rev) + # reverse order
  scale_colour_manual(values=subfamily_colorfix) +
  labs(y="",x="Shannon Index", color ="host subfamily") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

alpha.shannon.tribe.sort <- ggplot(alphaframe, aes(x = fct_reorder(host_tribe, Shannon, .fun = median, na.rm = TRUE), y=Shannon)) +
  geom_boxplot(aes(group = host_tribe))+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_subfamily, shape = country), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="",y="Shannon Index") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

alpha.shannon.country  <- ggplot(alphaframe, aes(x=country, y=Shannon)) +
  #geom_violin(aes(group = country), draw_quantiles = c(0.5), trim = T,  scale = "area")+ 
  geom_boxplot(aes(group = country), outlier.shape = NA)+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_subfamily, shape=country), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="",y="Shannon Index") + # , subtitle="country"
  theme(plot.subtitle = element_text(hjust = 0.5))

alpha.subfamily.country  <- ggplot(alphaframe, aes(x=country, y=Shannon)) +
  #geom_violin(aes(group = country), draw_quantiles = c(0.5), trim = T,  scale = "width")+ # scale = "width" for equal width
  geom_boxplot(aes(group = country), outlier.shape = NA)+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_subfamily, shape=country), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="",y="Shannon Index", color="host subfamily") +
  facet_wrap(~ host_subfamily, scales = "free_x") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

alpha.subfamily.country.tribe <- ggplot(alphaframe, aes(x=country, y=Shannon)) +
  geom_violin(aes(group = country), draw_quantiles = c(0.5), trim = T,  scale = "width")+ # scale = "width" for equal width
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_tribe), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=tribe_colorfix) +
  labs(x="",y="Shannon Index", color = "host subfamily") +
  facet_wrap(~ host_subfamily, scales = "free_x")

alpha.location  <- ggplot(alphaframe, aes(x=location, y=Shannon)) +
  geom_boxplot(aes(group = location))+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_subfamily, shape=country), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="",y="Shannon Index", color ="host subfamily") +
  facet_wrap(~ host_subfamily, scales = "free_x") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

alpha.location.tribe  <- ggplot(alphaframe, aes(x=location, y=Shannon)) +
  geom_boxplot(aes(group = location))+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_tribe), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=tribe_colorfix) +
  labs(x="",y="Shannon Index") +
  facet_wrap(~ host_subfamily, scales = "free_x")


pdf(file.path(out_dir, "02_sample_div_alphaframe_groups.pdf"), width=12, height=6)
alpha.shannon
alpha.rich
alpha.rich.tribe
alpha.rich.tribe.sort
alpha.shannon2
alpha.shannon.tribe.sort
alpha.shannon.country
alpha.subfamily.country 
alpha.subfamily.country.tribe
alpha.location
alpha.location.tribe 
dev.off()


#### alphaframe stats lm Shannon -------------

# Statistic notes: Quick plot data groups that will be tested
boxplot(alphaframe$Shannon ~ alphaframe$group) # chr
plot(alphaframe$Shannon ~ alphaframe$kw) # integer

# Quick test for 'group'
summary(aov(alphaframe$Shannon ~ alphaframe$group))
TukeyHSD(aov(alphaframe$Shannon ~ alphaframe$group)) 
anova(lm(Shannon ~ group, data=alphaframe))
oneway.test(alphaframe$Shannon ~ alphaframe$group)
kruskal.test(alphaframe$Shannon ~ alphaframe$group)
pairwise.wilcox.test(alphaframe$Shannon, alphaframe$group, p.adjust="BH")

# Check normal distribution of data
#hist(alphaframe$Shannon)

# Bartlett test of homogeneity of variances (too sensitive)
bartlett.test(Shannon ~ group, data=alphaframe)
# Fligner-Killeen’s test of homogeneity of variances (less sensitive)
fligner.test(Shannon ~ group, data=alphaframe)
# Levene's Test for homogeneity of variances (needs 'library(car)')
leveneTest(Shannon ~ group, data=alphaframe)

# optional Welch’s ANOVA (robust to unequal variance)
oneway.test(Shannon ~ group, data = alphaframe, var.equal = FALSE)

# Alternative kruskal Wallis test with pairwise Wilcox test
kruskal.test(Shannon ~ host_subfamily, data = alphaframe)
pairwise.wilcox.test(alphaframe$Shannon, alphaframe$host_subfamily, p.adjust="BH")

# Two sample Wilcox test
wilcox.test(Shannon ~  country, data = alphaframe)
t.test(Shannon ~  country, data = alphaframe)

# Illustrate 2 way test design for two way anova
#boxplot(Shannon ~ group1 * group2, data=alphaframe, frame = FALSE, 
#        col = c("#00AFBB", "#E7B800"), ylab="Shannon")

# In case of violation of normality, transform data, check distribution
hist(alphaframe$Shannon)
hist(sqrt(alphaframe$Shannon))
hist(log(alphaframe$Shannon))
bartlett.test(Shannon ~ group, data=alphaframe)
bartlett.test(sqrt(Shannon) ~ group, data=alphaframe)
bartlett.test(log(Shannon) ~ group, data=alphaframe)
fligner.test(Shannon ~ group, data=alphaframe)
fligner.test(sqrt(Shannon) ~ group, data=alphaframe)
fligner.test(log(Shannon) ~ group, data=alphaframe)
# if necessary transform data in dataframe
#alphaframe <- mutate(alphaframe, sqrtshannon = sqrt(Shannon))


## For Shannon diversity use lm/aov or lmer, but not glm or glmer.
boxplot(alphaframe$Shannon ~ alphaframe$host_subfamily) 
lm.shannon <- lm(Shannon ~ host_subfamily + location , data = alphaframe)
summary(lm.shannon) 
summary(lm.shannon)$adj.r.squared # R-squared:   0.323 
# Shannon diversity was tested by a linear model (‘lm’) and the ‘Anova’ function applied to the fitted model.
anova(lm.shannon)
# location not significant when host_subfamily is in the model!!
Anova(lm.shannon, type=2) # independent order
Anova(lm.shannon, type=2, white.adjust = TRUE) # if unequal variances across groups
vif(lm.shannon) # variance inflation factor should be <3  
plot(lm.shannon,2)
hist(lm.shannon$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon)) # Run Shapiro-Wilk test on residuals 

# Univariate tests: Shannon

# Quick model comparison test 
predictors <- c("host_family", "host_subfamily", "host_tribe",
                "host_genus", "country", "location", "sublocation",
                "temperature" , "year", "kw")

# Linear model fit for Shannon on predictors, extract F-values, p-values, R² and AIC
lm.results.summary <- map_dfr(predictors, function(p) {
  mod <- lm(reformulate(p, response = "Shannon"), data = alphaframe)
  s <- summary(mod)
  f <- s$fstatistic
  tibble(
    predictor = p,
    F_value   = unname(f[1]),
    p_value   = pf(f[1], f[2], f[3], lower.tail = FALSE),
    R2_adj    = s$adj.r.squared,
    AIC       = AIC(mod) ) }) %>%
  arrange(AIC)

# Show results as list 
as.data.frame(lm.results.summary)
# Shannon best AIC: genus > subfamily > tribe > location  > sublocation > country
# Main effect of host taxonomy followed by location

# Multivariate tests: Shannon, combine taxonomy with location
lm.comb1 <- lm(Shannon ~ host_subfamily + location, data = alphaframe)
lm.comb2 <- lm(Shannon ~ host_subfamily + location + temperature, data = alphaframe)
lm.comb3 <- lm(Shannon ~ host_tribe + location, data = alphaframe)
lm.comb4 <- lm(Shannon ~ host_tribe + location + temperature , data = alphaframe)
lm.comb5 <- lm(Shannon ~ host_genus + location  , data = alphaframe)
lm.comb6 <- lm(Shannon ~ host_genus + location + temperature , data = alphaframe)
lm.comb7 <- lm(Shannon ~ host_subfamily + sublocation , data = alphaframe)
lm.comb8 <- lm(Shannon ~ host_tribe + sublocation , data = alphaframe)
lm.comb9 <- lm(Shannon ~ host_tribe + country , data = alphaframe)
lm.comb10 <- lm(Shannon ~ host_subfamily + country , data = alphaframe)

AIC(lm.comb1,lm.comb2,lm.comb3,lm.comb4,lm.comb5, lm.comb6, lm.comb7, lm.comb8, lm.comb9, lm.comb10 )
# subfamily + location better than 3 predictors

# For univariate tests with host genus
# Robustness Check against over fitting  (filter genus with ≥2 samples)
genus_to_keep <- alphaframe %>%
  count(host_genus) %>%
  filter(n >= 2) %>%
  pull(host_genus)

# Filter the dataframe
alphaframe.robust <- alphaframe %>%
  filter(host_genus %in% genus_to_keep)

boxplot(alphaframe.robust$Shannon ~ alphaframe.robust$host_genus)
lm.shannon.genus.robust <- lm(Shannon ~ host_genus  , data = alphaframe.robust)
summary(lm.shannon.genus.robust) # R-squared:  0.3848
Anova(lm.shannon.genus.robust, type=2)

lm.shannon.robust.combine <- lm(Shannon ~ host_genus + location  , data = alphaframe.robust)
summary(lm.shannon.robust.combine) # R-squared:  0.397 
Anova(lm.shannon.robust.combine, type=2)
# check colinearity
#table(alphaframe.robust$host_genus, alphaframe.robust$location)
vif(lm.shannon.robust.combine, type = "predictor") # variance inflation factor should be <3 
alias(lm.shannon.robust.combine) 
# Aliased coefficients in the model: location + host_genus !! Cannot be combined!


# Shannon best final models: host subfamily alone
boxplot(alphaframe$Shannon ~ alphaframe$host_subfamily)
lm.shannon.subfam <- lm(Shannon ~ host_subfamily , data = alphaframe)
summary(lm.shannon.subfam) 
summary(lm.shannon.subfam)$adj.r.squared 
Anova(lm.shannon.subfam, type=2)
plot(lm.shannon.subfam,2)
plot(lm.shannon.subfam, which = 1)
hist(lm.shannon.subfam$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon.subfam)) # Run Shapiro-Wilk test on residuals 

# Shannon best final models: host genus alone
summary(lm.shannon.genus.robust) 
summary(lm.shannon.genus.robust)$adj.r.squared 
Anova(lm.shannon.genus.robust, type=2)
plot(lm.shannon.genus.robust,2)
plot(lm.shannon.genus.robust, which = 1)
hist(lm.shannon.genus.robust$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon.genus.robust)) # Run Shapiro-Wilk test on residuals 

# Shannon best final models: host_subfamily  + location 
lm.shannon.subfam.location <- lm(Shannon ~  host_subfamily  + location     , data=alphaframe)
summary(lm.shannon.subfam.location)
summary(lm.shannon.subfam.location)$adj.r.squared 
Anova(lm.shannon.subfam.location, type=2)
# location not significant when host_subfamily is in the model!!
vif(lm.shannon.subfam.location) # variance inflation factor should be <3 okay
alias(lm.shannon.subfam.location) # all good
table(alphaframe$host_subfamily, alphaframe$location) # check colinearity
plot(lm.shannon.subfam.location,2)
plot(lm.shannon.subfam.location, which = 1)
hist(lm.shannon.subfam.location$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon.subfam.location)) # Run Shapiro-Wilk test on residuals 

# variance partitioning for lm model on Type II ANOVA
testme <- lm.shannon.subfam.location
anovaII_out <- Anova(testme, type = 2)
ss_res <- sum(residuals(testme)^2)
# Partial R² per term
partial_R2 <- anovaII_out$`Sum Sq` / (anovaII_out$`Sum Sq` + ss_res)
names(partial_R2) <- rownames(anovaII_out)
partial_R2

# plot variance explained barplot
barplot(partial_R2[-length(partial_R2)], # exclude residuals
        las=2, col="skyblue",
        ylab="Proportion of variance explained",
        main="Variance partitioning")

# Partial R2 or eta squared (library effectsize)
eta_squared(anovaII_out, partial = TRUE)
# same as partial_R2


lm.shannon.tribe.country <- lm(Shannon ~  host_tribe  + country , data=alphaframe)
summary(lm.shannon.tribe.country)
summary(lm.shannon.tribe.country)$adj.r.squared 
Anova(lm.shannon.tribe.country, type=2)
# Minor effect of country when including tribe!
vif(lm.shannon.tribe.country) # variance inflation factor should be <3
alias(lm.shannon.tribe.country) # all good
table(alphaframe$host_tribe, alphaframe$country) # check colinearity

# variance partitioning for lm model on Type II ANOVA
testme <- lm.shannon.tribe.country
anovaII_out <- Anova(testme, type = 2)
ss_res <- sum(residuals(testme)^2)
# Partial R² per term
partial_R2 <- anovaII_out$`Sum Sq` / (anovaII_out$`Sum Sq` + ss_res)
names(partial_R2) <- rownames(anovaII_out)
partial_R2

# plot variance explained barplot
barplot(partial_R2[-length(partial_R2)], # exclude residuals
        las=2, col="skyblue",
        ylab="Proportion of variance explained",
        main="Variance partitioning")

# Partial R2 or eta squared (library effectsize)
eta_squared(anovaII_out, partial = TRUE)


# Test individual subfamilies Shannon country
lm.shannon.coliadinae <- lm(Shannon ~ country, data = subset(alphaframe, host_subfamily == "Coliadinae"))
summary(lm.shannon.coliadinae)
anova(lm.shannon.coliadinae)
plot(lm.shannon.coliadinae,2)
hist(lm.shannon.coliadinae$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon.coliadinae)) # Run Shapiro-Wilk test on residuals

# Two sample tests
wilcox.test(Shannon ~  country, data = subset(alphaframe, host_subfamily == "Coliadinae"))
t.test(Shannon ~  country, data = subset(alphaframe, host_subfamily == "Coliadinae"))


lm.shannon.dismorphinae <- lm(Shannon ~ country, data = subset(alphaframe, host_subfamily == "Dismorphiinae"))
summary(lm.shannon.dismorphinae)
anova(lm.shannon.dismorphinae)
plot(lm.shannon.dismorphinae,2)
hist(lm.shannon.dismorphinae$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon.dismorphinae)) # Run Shapiro-Wilk test on residuals

wilcox.test(Shannon ~ country, data = subset(alphaframe, host_subfamily == "Dismorphiinae"))  
t.test(Shannon ~ country, data = subset(alphaframe, host_subfamily == "Dismorphiinae"))  


lm.shannon.pierinae <- lm(Shannon ~ country, data = subset(alphaframe, host_subfamily == "Pierinae"))
summary(lm.shannon.pierinae)
anova(lm.shannon.pierinae)
plot(lm.shannon.pierinae,2)
hist(lm.shannon.pierinae$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon.pierinae)) # Run Shapiro-Wilk test on residuals 

wilcox.test(Shannon ~ country, data = subset(alphaframe, host_subfamily == "Pierinae"))
t.test(Shannon ~ country, data = subset(alphaframe, host_subfamily == "Pierinae"))


lm.shannon.heliconiinae <- lm(Shannon ~ location, data = subset(alphaframe, host_subfamily == "Heliconiinae"))
summary(lm.shannon.heliconiinae)
anova(lm.shannon.heliconiinae)
plot(lm.shannon.heliconiinae,2)
hist(lm.shannon.heliconiinae$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon.heliconiinae)) # Run Shapiro-Wilk test on residuals


sink(file.path(out_dir, "02_sample_div_alphaframe_groups_stats.txt"))
as.data.frame(lm.results.summary)
cat("\n")  # blank line
"Shannon host subfamily"
summary(lm.shannon.subfam)$adj.r.squared 
Anova(lm.shannon.subfam, type=2)
cat("\n")  # blank line
"Shannon genus robust"
summary(lm.shannon.genus.robust)$adj.r.squared 
Anova(lm.shannon.genus.robust, type=2)
cat("\n")  # blank line
"Shannon host subfamily + location"
summary(lm.shannon.subfam.location)$adj.r.squared 
Anova(lm.shannon.subfam.location, type=2)
cat("\n")  # blank line
"Shannon host tribe + country"
summary(lm.shannon.tribe.country)$adj.r.squared # 0.3230627
Anova(lm.shannon.tribe.country, type=2)
cat("\n")  # blank line
"Host taxa"
summary(lm.shannon.coliadinae)
summary(lm.shannon.dismorphinae)
summary(lm.shannon.pierinae)
summary(lm.shannon.heliconiinae)
sink()



#### alpha split country Peru / Germany---------------
alphaframe.Peru <- subset(alphaframe, country == "Peru")


# high correlation temp and elevation
lm.temp.el <- lm(temperature ~ elevation , data = alphaframe.Peru)
summary(lm.temp.el)$adj.r.squared # adjusted R-squared: 
Anova(lm.temp.el, type=2)

temp.elevation <- ggplot(alphaframe.Peru, aes(x = temperature , y = elevation)) +
  theme_line2() + 
  geom_point(shape = 17, size = 4, alpha = 0.8, aes(color=host_subfamily)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(label.y = 10, aes(label = after_stat(rr.label)),size=3) + #stat_cor(label.y = 200, size=3, p.accuracy = 0.00001) +
  scale_color_manual(values=subfamily_colorfix) +
  labs(x="temperature [°C]",y="Elevation [m.a.s.l.]", fill="host subfamily") 

# Peru alpha elevation
alpha.elevation.Peru  <- plot_richness(Peru.species,x="elevation", measures=c("Shannon","Observed", "InvSimpson", "Fisher")) +
  geom_point(size=4 , alpha=0.6, shape = 16 )  +
  stat_smooth(method="lm", color="black", linewidth=0.5,se=T) +
  stat_regline_equation(aes(label = after_stat(rr.label)), size=3) + stat_cor(label.y.npc = 0.90,size=3, p.accuracy = 0.00001) +
  scale_color_manual(values=subfamily_colorfix) +
  labs(x="elevation [m]",y="alpha div measures",  subtitle="Peru samples only" ) +
  theme_line2()

sample.shannon.elevation.Peru <-
  ggplot(alphaframe.Peru, aes(x = elevation, y = Shannon)) +
  theme_line2() + 
  geom_point(shape = 17, size = 4, alpha = 0.8, aes(color=host_subfamily)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(label.y = 4.0, aes(label = after_stat(rr.label)),size=3) + stat_cor(label.y = 3.7, size=3, p.accuracy = 0.00001) +
  scale_color_manual(values=subfamily_colorfix) +
  labs(x="elevation [m]",y="Shannon Index", fill="host subfamily",  subtitle="Peru samples" ) +
  coord_cartesian(ylim = c(0, 4), expand = TRUE,  default = FALSE,   clip = "on" ) 
 
sample.rich.elevation.Peru <-
  ggplot(alphaframe.Peru, aes(x = elevation, y = Observed)) +
  theme_line2() + 
  geom_point(shape = 17, size = 4, alpha = 0.8, aes(color=host_subfamily)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(label.y = 20, aes(label = after_stat(rr.label)),size=3) + stat_cor(label.y = 10, size=3, p.accuracy = 0.00001) +
  scale_color_manual(values=subfamily_colorfix) +
  labs(x="Elevation [m]",y="microbial richness (q=0)", fill="host subfamily") 
  coord_cartesian(ylim = c(0, 4), expand = TRUE,  default = FALSE,   clip = "on" ) 

sample.shannon.temp.Peru <-
    ggplot(alphaframe.Peru, aes(x = temperature, y = Shannon)) +
    theme_line2() + 
    geom_point(shape = 17, size = 4, alpha = 0.8, aes(color=host_subfamily)) +
    geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = TRUE) +
    stat_regline_equation(label.y = 4.0, aes(label = after_stat(rr.label)),size=3) + stat_cor(label.y = 3.7, size=3, p.accuracy = 0.00001) +
    scale_color_manual(values=subfamily_colorfix) +
    labs(x="temperature [°C]",y="Shannon Index", fill="host subfamily") +
    coord_cartesian(ylim = c(0, 4), expand = TRUE,  default = FALSE,   clip = "on" ) 
  
 
# Quick model comparison Shannon Peru 
peru_predictors <- c("host_family", "host_subfamily", "host_tribe",
                     "host_genus", "location", "sublocation",
                     "elevation",  "temperature" , "year", 
                     "host_subfamily + temperature" , "host_subfamily + location"  )

lm.results.peru <- map_dfr(peru_predictors, function(p) {
  mod <- lm(reformulate(p, response = "Shannon"), data = alphaframe.Peru)
  s <- summary(mod)
  f <- s$fstatistic
  tibble(
    predictor = p,
    F_value   = unname(f[1]),
    p_value   = pf(f[1], f[2], f[3], lower.tail = FALSE),
    R2_adj    = s$adj.r.squared,
    AIC       = AIC(mod) ) }) %>%
  arrange(AIC)

# Show model comparison Peru as list 
as.data.frame(lm.results.peru)
# best univariate predictor:  subfamily > genus > tribe > temperature > elevation > location
# best multivariate predictor:  subfamily + temperature



# Best model Shannon Peru
lm.peru <- lm(Shannon ~ host_subfamily +  temperature, data=alphaframe.Peru)
summary(lm.peru)
Anova(lm.peru, type=2) # independent order
plot(lm.peru,2)
hist(lm.peru$residuals) # extract the residuals
shapiro.test(residuals(lm.peru)) # Run Shapiro-Wilk test on residuals 
vif(lm.peru) # variance inflation factor should be <3
alias(lm.peru) # all good

# variance partitioning for lm model on Type II ANOVA
testme <- lm.peru
anovaII_out <- Anova(testme, type = 2)
ss_res <- sum(residuals(testme)^2)
# Partial R² per term
partial_R2 <- anovaII_out$`Sum Sq` / (anovaII_out$`Sum Sq` + ss_res)
names(partial_R2) <- rownames(anovaII_out)
partial_R2



#### alpha split country Germany 
alphaframe.Germany <- subset(alphaframe, country == "Germany")

# Modify alphaframe to remove Aglais
alphaframe.Germany.noAglais <- subset(alphaframe.Germany , host_genus != "Aglais")
alphaframe.Germany.Aglais <- subset(alphaframe.Germany , host_genus == "Aglais")
alphaframe.Germany.rbind <- rbind(alphaframe.Germany.noAglais, alphaframe.Germany.Aglais)

# correlation temp and calendar week
lm.temp.kw <- lm(temperature ~ poly(kw, 2)  , data = alphaframe.Germany)
summary(lm.temp.kw)$adj.r.squared # adjusted R-squared: 
Anova(lm.temp.kw, type=2)

kw.temp <- ggplot(alphaframe.Germany, aes(x = kw, y = temperature)) +
  theme_line2() + 
  geom_point(size = 4, alpha = 0.8, aes(color=host_subfamily)) +
  #geom_smooth(method = "auto", color = "black", linewidth = 0.5, se = TRUE) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(label.y = 2.0,formula = y ~ poly(x, 2), aes(label = after_stat(rr.label)),size=3) +# stat_cor(label.y = 4,size=3, p.accuracy = 0.001) +
  scale_color_manual(values = subfamily_colorfix) +
  labs(x="calendar week",y="temperature [°C]",
       color = "host subfamily") 
  

Shannon.temp.Germany <- ggplot(alphaframe.Germany.noAglais, aes(x = temperature, y = Shannon)) +
  theme_line2() + 
  geom_point(size = 4, alpha = 0.8, aes(color=host_subfamily)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(label.y = 4.0, aes(label = after_stat(rr.label)),size=3) + stat_cor(label.y = 3.7 ,size=3, p.accuracy = 0.00001) +
  #stat_regline_equation(label.y = 2.0,formula = y ~ poly(x, 2), aes(label = after_stat(rr.label)),size=3) +# stat_cor(label.y = 4,size=3, p.accuracy = 0.001) +
  scale_color_manual(values = subfamily_colorfix) +
  labs(x="temp",y="Shannon",
       color = "host subfamily") 


shannon.kw.Germany.noAglais <-
  ggplot(alphaframe.Germany.noAglais, aes(x = kw, y = Shannon)) +
  theme_line2() + 
  geom_point(size = 4, alpha = 0.8, aes(color=host_subfamily)) +
  #stat_smooth(data = alphaframe.Germany.noAglais, method="auto", color="black", linewidth=0.5,se=T) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(label.y = 4.0, aes(label = after_stat(rr.label)),size=3) + stat_cor(label.y = 3.7 ,size=3, p.accuracy = 0.00001) +
  geom_point(data = alphaframe.Germany.Aglais, aes(x = kw, y = Shannon, shape = "Aglais", fill=host_subfamily), color="black", stroke = 1, size = 4) +
  scale_shape_manual(values = c("Aglais" = 21)) +
  scale_color_manual(values = subfamily_colorfix) +
  scale_fill_manual(values = subfamily_colorfix, guide = "none") + # guide hide second legend
  labs(x="calendar week",y="Shannon Index",
           color = "host subfamily", shape = "excluded outgroup", subtitle="Germany samples" ) +
  coord_cartesian(ylim = c(0, 4), expand = TRUE,  default = FALSE,   clip = "on" ) +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

shannon.kw.Germany.noAglais2 <- ggplot(alphaframe.Germany.noAglais, aes(x = kw, y = Shannon)) +
  theme_line2() + 
  geom_point(size = 4, alpha = 0.8, aes(color = host_subfamily)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(label.y = 4.0, aes(label = after_stat(rr.label)), size = 3) + 
  stat_cor(label.y = 3.7, size = 3, p.accuracy = 0.00001) +
  geom_point(data = alphaframe.Germany.Aglais, 
             aes(x = kw, y = Shannon, fill = host_subfamily), 
             shape = 21, color = "black", stroke = 1, size = 4) +
  scale_color_manual(values = subfamily_colorfix) +
  scale_fill_manual(values = subfamily_colorfix) +
  labs(x = "Calendar Week", y = "Shannon Index", 
       color = "host subfamily", fill = "Aglais (outgroup)") +
  coord_cartesian(ylim = c(0, 4), expand = TRUE, default = FALSE, clip = "on")


# Shannon Germany model fit
ger_predictors <- c("host_family", "host_subfamily", "host_tribe",
                     "host_genus", "location", "sublocation",
                     "kw",  "temperature" )

# Quick check model fit for Shannon on predictors, extract F-values, p-values, R² and AIC
lm.results.germany <- map_dfr(ger_predictors, function(p) {
  mod <- lm(reformulate(p, response = "Shannon"), data = alphaframe.Germany.noAglais)
  s <- summary(mod)
  f <- s$fstatistic
  tibble(
    predictor = p,
    F_value   = unname(f[1]),
    p_value   = pf(f[1], f[2], f[3], lower.tail = FALSE),
    R2_adj    = s$adj.r.squared,
    AIC       = AIC(mod) ) }) %>%
  arrange(AIC)

# Show results as list 
as.data.frame(lm.results.germany)
# best predictor: kw > host_subfamily > host_tribe

lm.shannon.kw.subfam <- lm(Shannon ~ kw + host_subfamily, data=alphaframe.Germany.noAglais) 
summary(lm.shannon.kw.subfam)
Anova(lm.shannon.kw.subfam, type=2)
# Calendar week explains Shannon diversity independently of host identity
vif(lm.shannon.kw.subfam) # variance inflation factor should be <3

# variance partitioning for lm model on Type II ANOVA
testme <- lm.shannon.kw.subfam
anovaII_out <- Anova(testme, type = 2)
ss_res <- sum(residuals(testme)^2)
# Partial R² per term
partial_R2 <- anovaII_out$`Sum Sq` / (anovaII_out$`Sum Sq` + ss_res)
names(partial_R2) <- rownames(anovaII_out)
partial_R2

# Partial R2 or eta squared (library effectsize)
eta_squared(anovaII_out, partial = TRUE)
# Calendar week more important than host!

lm.shannon.kw.genus <- lm(Shannon ~ kw +host_genus , data=alphaframe.Germany.noAglais) 
summary(lm.shannon.kw.genus)
Anova(lm.shannon.kw.genus, type=2)
# No evidence that the effect of kw depends on host_genus turnover 

AIC(lm.shannon.kw.subfam, lm.shannon.kw.genus)
# subfamily better than genus!


pdf(file.path(out_dir, "02_sample_div_alphaframe_linear_split_country.pdf"), width=12, height=6)
temp.elevation
alpha.elevation.Peru 
sample.shannon.elevation.Peru
sample.rich.elevation.Peru
sample.shannon.temp.Peru 
Shannon.temp.Germany
shannon.kw.Germany.noAglais
shannon.kw.Germany.noAglais2
dev.off()


sink(file.path(out_dir, "02_sample_div_alphaframe_linear_split_country_stats.txt"))
"Peru model fit"
as.data.frame(lm.results.peru)
cat("\n")  # blank line
"Peru subfam temp"
summary(lm.peru)$adj.r.squared 
Anova(lm.peru, type=2) 
cat("\n")  # blank line
"Germany model fit"
as.data.frame(lm.results.germany)
cat("\n")  # blank line
"Germany subfam kw"
summary(lm.shannon.kw.subfam)$adj.r.squared
Anova(lm.shannon.kw.subfam, type=2)
sink()




#### alpha shannon temp / core  ----

temp.shannon.noAglais <- ggplot(alphaframe.noAglais, aes(x = temperature, y = Shannon)) +
  theme_line2() + 
  geom_point(size=4, alpha = 0.8, aes(color=host_subfamily)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(formula = y ~ x, label.y = 3, aes(label = after_stat(rr.label)),size=3) +
  stat_cor(method = "spearman", label.y = 0.6, size = 3, p.accuracy = 0.001) + # Spearman
  labs(x="temperature [°C]",y="Shannon", title="", color="host subfamily" ) +
  #facet_wrap(~country) +
  scale_colour_manual(values=subfamily_colorfix)


temp.shannon.noAglais2  <- ggplot(alphaframe.noAglais, aes(x = temperature, y = Shannon)) +
  theme_line2() + 
  geom_point(size = 4, alpha = 0.8, aes(color=host_subfamily, shape=country)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(label.y = 4.0, aes(label = after_stat(rr.label)),size=3) + stat_cor(label.y = 3.7 ,size=3, p.accuracy = 0.00001) +
  geom_point(data = alphaframe.Aglais, aes(x = temperature, y = Shannon, shape = "Aglais", fill=host_subfamily), stroke = 1, size = 4) +
  scale_shape_manual(
    values = c("Germany" = 19, "Peru" = 17, "Aglais" = 21),
    breaks = c("Germany", "Peru"),  # This controls the first shape legend
    guide = guide_legend(title = "Country", order = 1) ) +
  scale_fill_manual(
    values = subfamily_colorfix,
    guide = guide_legend(
      override.aes = list(shape = 21, size = 4, color = "black"),
      title = "Aglais (outgroup)",
      order = 2 ) ) +
  scale_color_manual(values = subfamily_colorfix) +
  labs(x="temperature [°C]",y="Shannon Index",
       color = "host subfamily") +
  coord_cartesian(ylim = c(0, 4), expand = TRUE,  default = FALSE,   clip = "on" ) +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )


# Shannon model fit full dataset noAglais 
predictors_all <- c("host_family", "host_subfamily", "host_tribe",
                    "host_genus", "location", "sublocation",
                    "elevation",  "temperature" ,  "kw",
                    "host_subfamily + temperature" , "host_subfamily + location"  )

# Quick check model fit for Shannon on predictors, extract F-values, p-values, R² and AIC
lm.results.all <- map_dfr(predictors_all, function(p) {
  mod <- lm(reformulate(p, response = "Shannon"), data = alphaframe.noAglais)
  s <- summary(mod)
  f <- s$fstatistic
  tibble(
    predictor = p,
    F_value   = unname(f[1]),
    p_value   = pf(f[1], f[2], f[3], lower.tail = FALSE),
    R2_adj    = s$adj.r.squared,
    AIC       = AIC(mod) ) }) %>%
  arrange(AIC)

# Show results as list 
as.data.frame(lm.results.all)
# final model subfamily * temp


## Shannon complete dataset from both countries host_subfamily * temperature
lm.temp.subfam <- lm(Shannon ~  host_subfamily * temperature , data=alphaframe.noAglais) 
summary(lm.temp.subfam)
summary(lm.temp.subfam)$adj.r.squared 
Anova(lm.temp.subfam, type=2) # independent order
vif(lm.temp.subfam, type = "predictor") # type = predictor for interaction models, 1 should be fine
alias(lm.temp.subfam) # only due to interaction effect, host_subfamily + temperature is okay
plot(lm.temp.subfam,2)
hist(lm.temp.subfam$residuals) # extract the residuals
shapiro.test(residuals(lm.temp.subfam)) # Run Shapiro-Wilk test on residuals 

# Partial R2 variance partitioning for lm model on Type II ANOVA
testme <- lm.temp.subfam
anovaII_out <- Anova(testme, type = 2)
ss_res <- sum(residuals(testme)^2)
# Partial R² per term
partial_R2 <- anovaII_out$`Sum Sq` / (anovaII_out$`Sum Sq` + ss_res)
names(partial_R2) <- rownames(anovaII_out)
partial_R2

# Partial R2 or eta squared (library effectsize)
eta_squared(anovaII_out, partial = TRUE)

## Host turnover at different temperatures
ggplot(alphaframe.noAglais, aes(x = host_subfamily, y = temperature)) +
  theme_line2() + geom_point(size = 4, alpha = 0.8, aes(color=host_genus)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), color = "black", linewidth = 0.5, se = TRUE)


lm.temp.host.turnover <- lm(temperature ~  host_subfamily , data=alphaframe.noAglais) 
hist(alphaframe.noAglais$temperature)
summary(lm.temp.host.turnover)
summary(lm.temp.host.turnover)$adj.r.squared
Anova(lm.temp.host.turnover, type=2) # independent order
plot(lm.temp.host.turnover,2)
hist(lm.temp.host.turnover$residuals) # extract the residuals
shapiro.test(residuals(lm.temp.host.turnover)) # Run Shapiro-Wilk test on residuals 

kruskal.test(temperature ~ host_subfamily, data = alphaframe.noAglais)


# Observed richness
temp.q0 <- ggplot(alphaframe, aes(x = temperature, y = Observed)) +
  theme_line2() + 
  geom_point(size=4, alpha = 0.8, aes(color=host_subfamily, shape=country)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = T) +
  stat_regline_equation(formula = y ~ x, label.y = 140, aes(label = after_stat(rr.label)),size=3) +
  stat_cor(label.y = 130,size=3, p.accuracy = 0.0001) +
  #stat_cor(method = "spearman", label.y = 1, size = 3, p.accuracy = 0.001) + # Spearman
  labs(x="temperature [°C]",y="microbial richness (q=0)", title="", color="host subfamily" ) +
  #facet_wrap(~country) +
  scale_colour_manual(values=subfamily_colorfix) +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

#plot(alphaframe$Observed ~ alphaframe$temperature) 
lm.rich.temp <- lm(Observed ~     host_subfamily  * temperature   , data=alphaframe) 
summary(lm.rich.temp)
summary(lm.rich.temp)$adj.r.squared
Anova(lm.rich.temp, type=2) # independent order
vif(lm.rich.temp, type = "predictor") # variance inflation factor should be <3
plot(lm.rich.temp,2)
hist(lm.rich.temp$residuals) # extract the residuals
shapiro.test(residuals(lm.rich.temp)) # Run Shapiro-Wilk test on residuals 

# variance partitioning for lm model on Type II ANOVA
testme <- lm.rich.temp
anovaII_out <- Anova(testme, type = 2)
ss_res <- sum(residuals(testme)^2)
# Partial R² per term
partial_R2 <- anovaII_out$`Sum Sq` / (anovaII_out$`Sum Sq` + ss_res)
names(partial_R2) <- rownames(anovaII_out)
partial_R2

# eta squared
eta_squared(anovaII_out, partial = TRUE)


# Interaction effect on temp and subfamily
temp.q0_subfamily <- ggplot(alphaframe, aes(x = temperature, y = Observed, color=host_subfamily)) +
  theme_line2() + 
  geom_point(size=3, alpha = 0.8, aes(shape=country)) +
  geom_smooth(method = "lm", linewidth = 1, se = F) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = T) +
  labs(x="temperature [°C]",y="microbial richness (q=0)", title="", color="host subfamily" ) +
  scale_colour_manual(values=subfamily_colorfix) +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )


temp.q1 <- ggplot(alphaframe, aes(x = temperature, y = hill_numbers )) +
  theme_line2() + 
  geom_point(size = 4, alpha = 0.8, aes(color=host_subfamily)) +
  geom_smooth(method = "auto", color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(formula = y ~ x, label.y = 3, aes(label = after_stat(rr.label)),size=3) +
  #annotate("text", x = 0.0, y = 0.2, label = label_text, size = 3, hjust = 0, parse = TRUE) +
  stat_cor(method = "spearman", label.y = 0.6, size = 3, p.accuracy = 0.001) + # Spearman
  labs(x="temperature [°C]",y="exp(Shannon)(q=1)", title="", color="host subfamily" ) +
  scale_colour_manual(values=subfamily_colorfix)

temp.q2 <- ggplot(alphaframe, aes(x = temperature, y = InvSimpson )) +
  theme_line2() + 
  geom_point(size = 4, alpha = 0.8, aes(color=host_subfamily)) +
  geom_smooth(method = "auto", color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(formula = y ~ poly(x, 2), label.y = 5, aes(label = after_stat(rr.label)),size=3) +
  #annotate("text", x = 0.0, y = 0.2, label = label_text, size = 3, hjust = 0, parse = TRUE) +
  stat_cor(method = "spearman", label.y = 1, size = 3, p.accuracy = 0.001) + # Spearman
  labs(x="temperature [°C]",y="InvSimpson (q=2)", title="", color="host subfamily" ) +
  scale_colour_manual(values=subfamily_colorfix)


#### alpha shannon core 
core.shannon.auto <- ggplot(alphaframe, aes(x = genus_core, y = Shannon )) +
  theme_line2() + 
  geom_point(size = 4, alpha = 0.8, aes(color=host_subfamily)) +
  geom_smooth(method = "auto", color = "black", linewidth = 0.5, se = TRUE) +
  #stat_regline_equation(label.y = 0.4, aes(label = after_stat(rr.label)),size=3) + # only for lm
  #stat_cor(label.y = 0.6,size=3, p.accuracy = 0.001) + # pearson assumes a linear relationship
  stat_cor(method = "spearman", label.y = 0.6, size = 3, p.accuracy = 0.001) + # Spearman rank correlation Non-parametric not necessarily linear
  labs(x="genus core [%]",y="Shannon Index", title="", color="host subfamily" ) +
  scale_colour_manual(values=subfamily_colorfix)

lm.quad <- lm(Shannon ~ poly(genus_core, 2) , data=alphaframe) # quadratic fit 
summary(lm.quad)
summary(lm.quad)$adj.r.squared
anova(lm.quad)
plot(lm.quad,2)
hist(lm.quad$residuals) # extract the residuals
shapiro.test(residuals(lm.quad))  

fstat <- summary(lm.quad)$fstatistic
f_val <- round(fstat[1], 2)
p_val <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
p_round <- ifelse(p_val < 0.0001, "< 0.0001", paste0("= ", round(p_val, 4)))
label_text <- paste0(
  "F == ", f_val, "~','~italic(p) ", 
  ifelse(p_val < 0.0000000001, "< 0.0000000001", paste0("== ", p_val_rounded)))

core.shannon <- ggplot(alphaframe, aes(x = genus_core, y = Shannon )) +
  theme_line2() + 
  geom_point(size = 4, alpha = 0.8, aes(color=host_subfamily, shape=country)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2),  color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(formula = y ~ poly(x, 2), label.y = 0.4, aes(label = after_stat(rr.label)),size=3) +
  annotate("text", x = 0.0, y = 0.2, label = label_text, size = 3, hjust = 0, parse = TRUE) +
  scale_x_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  labs(x="genus core [%]",y="Shannon Index", title="", color="host subfamily" ) +
  scale_colour_manual(values=subfamily_colorfix)


core.shannon.temp <- ggplot(alphaframe, aes(x = genus_core, y = Shannon )) +
  theme_line2() + 
  geom_point(alpha = 0.8, aes(color=host_subfamily, size=temperature)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2),  color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(formula = y ~ poly(x, 2), label.y = 0.4, aes(label = after_stat(rr.label)),size=3) +
  annotate("text", x = 0.0, y = 0.2, label = label_text, size = 3, hjust = 0, parse = TRUE) +
  labs(x="genus core [%]",y="Shannon Index", title="", color="host subfamily" ) +
  scale_colour_manual(values=subfamily_colorfix)


pdf(file.path(out_dir, "02_sample_div_alphaframe_temp_core.pdf"), width=12, height=6)
temp.shannon.noAglais
temp.shannon.noAglais2
temp.q0 
temp.q0_subfamily 
temp.q1
temp.q2
core.shannon.auto 
core.shannon 
core.shannon.temp 
dev.off()


sink(file.path(out_dir, "02_sample_div_alphaframe_temp_core_stats.txt"))
"Model fit noAglais"
as.data.frame(lm.results.all)
cat("\n")  # blank line
"Shannon Temperature subfamily"
summary(lm.temp.subfam)$adj.r.squared 
Anova(lm.temp.subfam, type=2) 
cat("\n")  # blank line
"Richness Temperature subfamily"
summary(lm.rich.temp)$adj.r.squared
Anova(lm.rich.temp, type=2) 
cat("\n")  # blank line
"Shannon Core "
summary(lm.quad)$adj.r.squared 
anova(lm.quad)
sink()


## 02 sample div beta (NMDS) -------------------

# ASV-level
# sample.ASV       #sample.ASV.rel    #LT removed (samples only) ASV level
# sample.filter    #sample.filter.rel 

# genus-level
# sample.species   #sample.species.rel #LT removed (samples only) genus level

# family-level
# sample.family    #sample.family.rel


sample.nmds <- ordinate(sample.species.rel, method="NMDS",distance = "bray", k=3, trymax=200)
sample.nmds.ASV <- ordinate(sample.ASV.rel, method="NMDS",distance = "bray", k=3, trymax=200)
# sample.nmds.family <- ordinate(sample.family.rel, method="NMDS",distance = "bray", k=3, trymax=200)

# Check goodness of fit with Shepards diagram
#stressplot(sample.nmds)
sample.nmds.stress <- round(sample.nmds$stress, 3)
sample.nmds.stress
# dimension : 2 Stress:     0.282 
# dimensions: 3 Stress:     0.214 # better 3 dimensions
# NMDS stress level 0.05	Excellent, 0.05–0.1 Good, 0.1–0.2 Fair, 0.2–0.3 Weak, > 0.3 Very poor

# NMDS
nmds.ASV <- plot_ordination(sample.ASV.rel,sample.nmds.ASV, color="host_subfamily", shape="country")+
  geom_point(size=4) + theme_grid() + stat_ellipse(aes(group = country))  +  ggtitle("sample.nmds.ASV"  )
nmds.species <- plot_ordination(sample.species.rel,sample.nmds, color="host_subfamily", shape="country")+
  geom_point(size=4) + theme_grid() + stat_ellipse(aes(group = factor(country))) +  ggtitle("sample.nmds.species" ) +
  annotate("text", x = Inf, y = -Inf, label = paste("Stress =", sample.nmds.stress), 
           hjust = 1.1, vjust = -0.5, size = 4)

# NMDS parameter
#sample_data(sample.species.rel)$weight <- as.numeric(as.character(sample_data(sample.species.rel)$weight))
sample_data(sample.species.rel)$year <- as.factor(sample_data(sample.species.rel)$year)

nmds.subfamily <- plot_ordination(sample.species.rel,sample.nmds, color="host_tribe", shape="year")+
  geom_point(size=4) + theme_grid() + facet_wrap(~host_subfamily)

#### beta div plots (NMDS) binder 
# combine NMDS sites with data frame to include custom alpha levels in figure
scores(sample.nmds)$sites
rownames(scores(sample.nmds)$sites) == row.names(sample_data(sample.species.rel)) # check if true before cbind
nmds.binder <- cbind(scores(sample.nmds)$sites,sample_data(sample.species.rel))

nmds.binder$samplesum <- sample_sums(sample.species)
nmds.binder$logsum <- log10(sample_sums(sample.species))
nmds.binder$rich <- estimate_richness(sample.species, measures=c("Observed", "InvSimpson", "Shannon"))
nmds.binder$Actinobacteria <- sample_sums(subset_taxa(sample.species.rel, phylum=="Actinobacteria" ))
nmds.binder$entero <- sample_sums(subset_taxa(sample.species.rel, order=="Enterobacterales" ))
nmds.binder$familycore <- sample_sums(core.family.genus)
nmds.binder$genuscore <- sample_sums(core.genus.rel)
nmds.binder$LT7000 <- ifelse(sample_sums(sample.species) > 7000, 'HT', 'LT7000')
nmds.binder$LT5000 <- ifelse(sample_sums(sample.species) > 5000, 'HT', 'LT5000')
nmds.binder$LT3000 <- ifelse(sample_sums(sample.species) > 3000, 'HT', 'LT3000')

#str(nmds.binder)

NMDS.binder.shannon.logsum <- ggplot(nmds.binder, aes(x=NMDS1,y=NMDS2,shape=LT5000, color=rich$Shannon, size=logsum))+
  stat_ellipse(linewidth=1, alpha=0.8, aes(group = factor(country))) +
  geom_point( alpha=0.8)+
  theme(plot.title =element_text(size=10, face='bold')) +
  scale_colour_viridis_c(direction = 1) +
  theme_grid() #+ labs(title="sample.species") +labs(color="Shannon")
  facet_wrap(~country)

NMDS.binder.wrap.shannon <-
  ggplot(nmds.binder , aes(x=NMDS1,y=NMDS2,color=rich$Shannon, shape=LT5000, size=Actinobacteria))+
  #stat_ellipse(linewidth=1, alpha=0.5 ) +
  geom_point( alpha=0.8)+
  theme(plot.title =element_text(size=10, face='bold')) +
  scale_colour_viridis_c(direction = 1) + facet_wrap(~country) +
  theme_grid() #labs(title="all time points") +labs(color="sampling\ntime point") +

NMDS.binder.genuscore <- ggplot(nmds.binder, aes(x=NMDS1,y=NMDS2, color=genuscore))+
  stat_ellipse(linewidth=1, alpha=0.8, aes(group = factor(host_tribe))) +
  geom_point( alpha=0.8)+
  theme(plot.title =element_text(size=10, face='bold')) +
  scale_colour_viridis_c(direction = -1) +
  theme_grid() + facet_wrap(~country) + geom_label(aes(label = sampleID), size = 3)


NMDS.binder.tribe <- ggplot(nmds.binder, aes(x=NMDS1,y=NMDS2))+
  stat_ellipse(level = 0.9, linewidth = 1, aes(group = host_tribe, color = host_tribe)) +
  geom_point(size = 5, aes(shape = country, color = host_tribe), alpha = 0.8) +
  facet_wrap(~country) +
  scale_colour_manual(values = tribe_colorfix) +
  theme_grid()

NMDS.binder.tribe.wrap <- ggplot(nmds.binder, aes(x=NMDS1,y=NMDS2))+
  stat_ellipse(level = 0.9, linewidth = 1, aes(group = host_tribe, color = host_tribe)) +
  geom_point(size = 5, aes(shape = country, color = host_tribe), alpha = 0.8) +
  facet_wrap(~host_subfamily) +
  scale_colour_manual(values = tribe_colorfix) +
  theme_grid()

NMDS.binder.genus.wrap <- ggplot(nmds.binder, aes(x=NMDS1,y=NMDS2))+
  stat_ellipse(level = 0.9, linewidth = 1, aes(group = host_genus, color = host_tribe)) +
  geom_point(size = 5, aes(shape = country, color = host_tribe), alpha = 0.8) +
  facet_wrap(~host_genus) +
  scale_colour_manual(values = tribe_colorfix) +
  theme_grid()


pdf(file.path(out_dir, "02_sample_div_beta_NMDS_binder.pdf"), width=8, height=8)
nmds.ASV
nmds.species
nmds.subfamily
NMDS.binder.shannon.logsum
NMDS.binder.wrap.shannon
NMDS.binder.genuscore
NMDS.binder.tribe
NMDS.binder.tribe.wrap
NMDS.binder.genus.wrap
dev.off()



# 02 sample div beta (PCoA)  ------------------
sample.PCoA <- ordinate(sample.species.rel, method="PCoA",distance = "bray")
sample.PCoA.ASV <- ordinate(sample.ASV.rel, method="PCoA",distance = "bray")
# sample.PCoA.family <- ordinate(sample.family.rel, method="PCoA",distance = "bray")

# Check variance explained PCoA
explained <- sample.PCoA$values$Relative_eig
round(explained[1:5] * 100, 2) 

# PCoA
PCoA.ASV <- plot_ordination(sample.ASV.rel,sample.PCoA.ASV, color="host_subfamily", shape = "country")+
  geom_point(size=4)+theme_grid() + stat_ellipse(aes(group = country)) +    ggtitle("sample.PCoA.ASV"  )

PCoA.species <- plot_ordination(sample.species.rel,sample.PCoA, color="host_subfamily", shape = "country")+
  geom_point(size=4)+theme_grid() + stat_ellipse(aes(group = country))  +  ggtitle("sample.PCoA.species"  )

# PCoA parameter
sample.PCoA.genus <- plot_ordination(sample.species.rel,sample.PCoA, color="host_genus", shape = "country")+
  geom_point(size=4)+theme_grid() + stat_ellipse(aes(group = country)) 
sample.PCoA.subfamily <- plot_ordination(sample.species.rel,sample.PCoA, color="host_subfamily", shape = "country")+
  geom_point(size=4)+theme_grid() + stat_ellipse(aes(group = host_subfamily))
sample.PCoA.location <- plot_ordination(sample.species.rel,sample.PCoA, color="host_tribe", shape = "country")+
  geom_point(size=4)+theme_grid() + facet_wrap(~location)
sample.PCoA.species <- plot_ordination(sample.species.rel,sample.PCoA, shape = "country")+
  geom_point(size=4, aes(color = host_species))+theme_grid() + stat_ellipse(aes(group = country)) +
  scale_color_manual(values = species_colorfix) 


# simple plot_ordination PCoA plots
pdf(file.path(out_dir, "02_sample_div_beta_PCoA.pdf"), width=12, height=6)
PCoA.ASV
PCoA.species
sample.PCoA.genus
sample.PCoA.subfamily
sample.PCoA.location 
sample.PCoA.species
dev.off()


#### beta div plots (PCoA) binder
# Export pcoa scores 
pcoa_scores_3 <- as.data.frame(sample.PCoA$vectors[, 1:3]) %>%
  rownames_to_column("Sample")

all(df.sample$Sample %in% df.core$Sample) # check if true before left_join
all(df.core$Sample %in% df.sample$Sample) # check if true before left_join

# Combine with sample data and core data into pcoa.binder # left_join better than cbind
pcoa.binder <- df.sample %>% 
  left_join(pcoa_scores_3, by = "Sample") %>%
  left_join(df.core, by = "Sample")

# Set rownames to Sample, keep Sample column for left_join
rownames(pcoa.binder) <- pcoa.binder$Sample

axes_var <- sample.PCoA$values$Relative_eig * 100 # Export axis % variation 
xlab_pcoa1 <- paste0("PCoA1 (", round(axes_var[1], 1), "%)")
ylab_pcoa2 <- paste0("PCoA2 (", round(axes_var[2], 1), "%)")
ylab_pcoa3 <- paste0("PCoA3 (", round(axes_var[3], 1), "%)")

PCoA.host.overview  <- ggplot(pcoa.binder, aes(Axis.1, Axis.2)) +
  stat_ellipse(level = 0.9, linewidth = 0.7, alpha = 0.5, aes(group = host_subfamily, color = host_subfamily)) + #linetype = 5,
  geom_point(size = 4, aes(shape = country, color = host_subfamily), alpha = 0.8) +
  facet_wrap(~host_family) +
  scale_colour_manual(values = subfamily_colorfix) +
  labs(x = xlab_pcoa1, y = ylab_pcoa2, color = "host subfamily") +
  theme_grid() +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

PCoA.axis13  <- ggplot(pcoa.binder, aes(Axis.1, Axis.3)) +
  stat_ellipse(level = 0.9, linewidth = 0.7, alpha = 0.5, aes(group = host_subfamily, color = host_subfamily)) +
  geom_point(size = 5, aes(shape = country, color = host_subfamily), alpha = 0.8) +
  facet_wrap(~host_family) +
  scale_colour_manual(values = subfamily_colorfix) +
  labs(x = xlab_pcoa1, y = ylab_pcoa3) +
  theme_grid()

PCoA.axis23  <- ggplot(pcoa.binder, aes(Axis.2, Axis.3)) +
  stat_ellipse(level = 0.9, linewidth = 0.7, alpha = 0.5, aes(group = host_subfamily, color = host_subfamily)) +
  geom_point(size = 5, aes(shape = country, color = host_subfamily), alpha = 0.8) +
  facet_wrap(~host_family) +
  scale_colour_manual(values = subfamily_colorfix) +
  labs(x = ylab_pcoa2, y = ylab_pcoa3) +
  theme_grid()

PCoA.host.overview.tribe  <- ggplot(pcoa.binder, aes(Axis.1, Axis.2)) +
  stat_ellipse(type = "t",level = 0.9, linewidth = 0.7, alpha = 0.5, aes(group = host_tribe, color = host_tribe)) +
  geom_point(size = 4, aes(shape = country, color = host_tribe), alpha = 0.8) +
  facet_wrap(~host_family) +
  scale_colour_manual(values = tribe_colorfix) +
  labs(x = xlab_pcoa1, y = ylab_pcoa2) +
  theme_grid()

PCoA.subfamily.wrap  <- ggplot(pcoa.binder, aes(Axis.1, Axis.2)) +
  stat_ellipse(level = 0.9, linewidth = 0.7, alpha = 0.5, aes(group = host_subfamily, color = host_tribe)) +
  geom_point(size = 5, aes(shape = country, color = host_tribe), alpha = 0.8) +
  facet_wrap(~host_subfamily) +
  scale_colour_manual(values = tribe_colorfix) +
  labs(x = xlab_pcoa1, y = ylab_pcoa2) +
  theme_grid() #+ geom_label(aes(label = sampleID), size = 4)

PCoA.tribe.wrap <- ggplot(pcoa.binder, aes(Axis.1, Axis.2)) +
  stat_ellipse(level = 0.9, linewidth = 0.7, alpha = 0.5, aes(color = host_tribe)) +
  geom_point(size=4, alpha = 0.8, aes(shape = country, color = host_tribe)) + 
  facet_wrap(~host_tribe) +
  scale_color_manual(values = tribe_colorfix) +
  labs(x = xlab_pcoa1, y = ylab_pcoa2) + theme_grid() #+ geom_label(aes(label = sampleID), size = 4) 


# Group unique genera together (<2 specimen)
table(pcoa.binder$host_genus) # count the number of samples per host_genus
low_rep_genera <- names(which(table(pcoa.binder$host_genus) == 1))

pcoa.binder <- pcoa.binder %>%
  mutate(host_genus_grouped = ifelse(host_genus %in% low_rep_genera, "unique genera", as.character(host_genus)))

genus_level_order <- sort(unique(pcoa.binder$host_genus_grouped))
genus_level_order <- c(setdiff(genus_level_order, "unique genera"), "unique genera")
pcoa.binder$host_genus_grouped <- factor(pcoa.binder$host_genus_grouped, levels = genus_level_order)

# genus wrap with grouped genera
PCoA.genus.wrap <- ggplot(pcoa.binder, aes(Axis.1, Axis.2)) +
  stat_ellipse(level = 0.9,  alpha = 0.5, aes(group = host_subfamily, color = host_subfamily)) +
  geom_point(size=4, alpha = 0.8, aes(shape = country, color = host_subfamily)) + 
  facet_wrap(~host_genus_grouped) +
  scale_color_manual(values = subfamily_colorfix) +
  labs(x = xlab_pcoa1, y = ylab_pcoa2) + theme_grid() #+ geom_label(aes(label = sampleID), size = 4) 
  
PCoA.genus.wrap.tribe <- ggplot(pcoa.binder, aes(Axis.1, Axis.2), color = host_tribe)+
  stat_ellipse(linewidth = 0.7, alpha = 0.5, aes(group = host_genus, color = host_tribe))  + 
  geom_point(size=4, alpha = 0.8, aes(shape = country, color = host_tribe)) +
  facet_wrap(~host_genus_grouped) +
  scale_color_manual(values = tribe_colorfix) +
  labs(x = xlab_pcoa1, y = ylab_pcoa2) +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2)) +
  theme_grid()

PCoA.host.country.wrap  <- ggplot(pcoa.binder, aes(Axis.1, Axis.2) )+
  stat_ellipse(level = 0.9, linewidth = 0.5, aes(group = host_subfamily, color = host_subfamily)) +
  geom_point(size = 5, aes(shape = location, color = host_subfamily), alpha = 0.8) +
  facet_wrap(~country) +
  scale_colour_manual(values = subfamily_colorfix) +
  labs(x = xlab_pcoa1, y = ylab_pcoa2, color = "host subfamily") +
  theme_grid() +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )


PCoA.host.country  <- ggplot(pcoa.binder, aes(Axis.1, Axis.2)) +
  stat_ellipse(level = 0.9, linewidth = 0.5, aes(group = country)) +
  geom_point(size = 5, aes(shape = country, color = host_subfamily), alpha = 0.8) +
  #facet_wrap(~country) +
  scale_colour_manual(values = subfamily_colorfix) +
  labs(x = xlab_pcoa1, y = ylab_pcoa2, color = "host subfamily") +
  theme_grid() +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )


pcoa.binder.core <-
  ggplot(pcoa.binder, aes(x=Axis.1,y=Axis.2, shape=country))+
  stat_ellipse(linewidth=1, alpha=0.2, aes(group = factor(country)) ) +
  geom_point(alpha=0.8, aes(size=genus_core, color=host_subfamily))+
  theme_grid() +
  scale_colour_manual(values = subfamily_colorfix) +
  scale_size_continuous(range = c(2, 10))
  

pdf(file.path(out_dir, "02_sample_div_beta_PCoA_binder.pdf"), width=12, height=6)
PCoA.host.overview
PCoA.axis13
PCoA.axis23
PCoA.host.overview.tribe
PCoA.subfamily.wrap
PCoA.tribe.wrap 
PCoA.genus.wrap
PCoA.genus.wrap.tribe
PCoA.host.country.wrap 
PCoA.host.country 
pcoa.binder.core 
dev.off()


#### betaframe betadisper ------------
sample.vegdist <- vegan::vegdist(t(otu_table(sample.species.rel)), index="bray")
#sample.distance <- phyloseq::distance(otu_table(sample.species.rel), method = "bray") # identical outcome

### New betaframe 
betaframe <- pcoa.binder
all(rownames(betaframe) == names(betadisper(sample.vegdist, group = betaframe$host_subfamily)$distances)) # check if true

# Beta dispersion for groups of interest
beta_family      <- betadisper(sample.vegdist, group = betaframe$host_family)          
beta_subfamily   <- betadisper(sample.vegdist, group = betaframe$host_subfamily)
beta_tribe       <- betadisper(sample.vegdist, group = betaframe$host_tribe)
beta_genus       <- betadisper(sample.vegdist, group = betaframe$host_genus)
beta_country     <- betadisper(sample.vegdist, group = betaframe$country)
beta_location    <- betadisper(sample.vegdist, group = betaframe$location)
beta_sublocation <- betadisper(sample.vegdist, group = betaframe$sublocation)

setdiff(betaframe$Sample, names(beta_family$distances)) # optional: check if samples are missing

# Put distances in betaframe for plotting and linear model
betaframe$beta_family    <- beta_family$distances[match(betaframe$Sample, names(beta_family$distances))]
betaframe$beta_subfamily <- beta_subfamily$distances[match(betaframe$Sample, names(beta_subfamily$distances))]
betaframe$beta_tribe     <- beta_tribe$distances[match(betaframe$Sample, names(beta_tribe$distances))]
betaframe$beta_genus     <- beta_genus$distances[match(betaframe$Sample, names(beta_genus$distances))]
betaframe$beta_country   <- beta_country$distances[match(betaframe$Sample, names(beta_country$distances))]
betaframe$beta_location  <- beta_location$distances[match(betaframe$Sample, names(beta_location$distances))]
betaframe$beta_sublocation <- beta_sublocation$distances[match(betaframe$Sample, names(beta_sublocation$distances))]
betaframe$host_subfamily_ordered <- factor(betaframe$host_subfamily, levels = subfamily_ordered)
betaframe$temperature <- betaframe$temp_week

# Simple base R plot
plot(betadisper(sample.vegdist, group = betaframe$host_subfamily))
plot(betadisper(sample.vegdist, group = betaframe$location))
plot(betadisper(sample.vegdist, group = betaframe$country))


# Betadisper family level
betadisp.family <- ggplot(betaframe, aes(x=host_family, y=beta_family))+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  scale_colour_manual(values=subfamily_colorfix) +  theme_line2() +
  #geom_smooth(method="lm", color="black", linewidth=0.5) +
  labs(y="distance to centroid", x='') 

permutest(beta_family, permutations = 999) 
plot(beta_family)

# Betadisper subfamily level
betadisp.subfamily.country <- ggplot(betaframe, aes(x=host_subfamily_ordered, y=beta_subfamily))+
  geom_violin( trim = T, adjust = 0.8, scale = "width" ) + # draw_quantiles = c(0.5),
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_tribe, shape=country),alpha = 0.8) +
  scale_colour_manual(values=tribe_colorfix) +  
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y="distance to centroid", x='', color ="host tribe") + #subtitle="distance to centroid", 
  facet_wrap(~ country, scales = "free_x") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2))

betadisp.subfamily <- ggplot(betaframe, aes(x=host_subfamily_ordered, y=beta_subfamily))+
  #geom_boxplot() +
  geom_violin( trim = T, adjust = 0.8, scale = "width" ) + # draw_quantiles = c(0.5),
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_tribe, shape=country),alpha = 0.8) +
  scale_colour_manual(values=tribe_colorfix) +  
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y="distance to centroid", x='', color ="host tribe") + #subtitle="beta dispersion", 
  facet_wrap(~ host_family, scales = "free_x") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2))

permutest(beta_subfamily, permutations = 999) 
plot(beta_subfamily)
# Difference in beta dispersion indicates variation in within-group community heterogeneity


# sort by subfamily
betadisp.subfamily.sorted <- ggplot(betaframe, aes(x = fct_reorder(host_subfamily, beta_subfamily, .fun = median, na.rm = TRUE), y = beta_subfamily)) +
  #geom_boxplot() +
  #geom_violin(#draw_quantiles = c(0.5), trim = F,  scale = "width"
  #           ) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 4, aes(color = host_tribe, shape = country), alpha = 0.8) +
  scale_colour_manual(values = tribe_colorfix) +  
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(subtitle = "beta dispersion",  y = "distance to centroid",  x = '',  color = "host tribe") +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2) )

betadisp.subfamily.ordered <- ggplot(betaframe, aes(x=beta_subfamily, y=host_subfamily_ordered))+
  #geom_boxplot(aes(group = host_subfamily))+
  geom_violin( trim = T, adjust = 0.6 , scale = "width" ) + # draw_quantiles = c(0.5),
  geom_point(position = position_jitter(h = 0.1, w = 0), size=4, aes(color = host_tribe, shape=country),alpha = 0.8) +
  scale_colour_manual(values=tribe_colorfix) +  
  theme_line2()+ scale_y_discrete(limits = rev) + # reverse order
  labs(subtitle="", y="", x='distance to centroid', color ="host tribe") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2)) 

# Betadisper tribe level
betadisp.tribe <- ggplot(betaframe, aes(x=host_tribe, y=beta_tribe))+
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_tribe),alpha = 0.8) +
  scale_colour_manual(values=tribe_colorfix) +  theme_line2() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y="distance to centroid", x='')

permutest(beta_tribe, permutations = 999) 
plot(beta_tribe)

betadisp.tribe.sorted <- ggplot(betaframe, aes(x=fct_reorder(host_tribe, beta_tribe, .fun = median, na.rm = TRUE), y=beta_tribe))+
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  scale_colour_manual(values=subfamily_colorfix) +  theme_line2() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y="distance to centroid", x='')

# Betadisper country
betadisp.country <- ggplot(betaframe, aes(x=country, y=beta_country))+
  geom_boxplot(aes(group = country))+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  scale_colour_manual(values=subfamily_colorfix) +  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  #geom_smooth(method="lm", color="black", linewidth=0.5) +
  labs(y="beta distance", x='')

permutest(beta_country, permutations = 999) 
plot(beta_country)
# No variation in within-group community heterogeneity by country

# Betadisper location
betadisp.location <- ggplot(betaframe , aes(x=location, y=beta_location))+
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_tribe),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  facet_wrap(~ country, scales = "free_x") +
  scale_colour_manual(values=tribe_colorfix)

permutest(beta_location, permutations = 999) 
plot(beta_location)
# No variation in within-group community heterogeneity by location 


### Beta Dispersion genus 
# Simplify dataset / remove host_genus with only 1 count (as beta dispersion would be zero anyway)
genus_to_keep <- df.sample %>%
  group_by(host_genus) %>%
  filter(n() >= 2) %>%   # keep genera with >= 2 samples
  pull(host_genus)

# Subset the phyloseq object to keep only the selected host genera
genus.robust <- subset_samples(sample.species.rel, host_genus %in% genus_to_keep)
genus.robust.df  <- data.frame(sample_data(genus.robust))

# Tests if groups have different beta dispersion
genus.vegdist <- vegan::vegdist(t(otu_table(genus.robust)), index="bray")
beta_genus_robust <- betadisper(genus.vegdist, group = genus.robust.df$host_genus)
permutest(beta_genus_robust, permutations = 9999)
plot(beta_genus_robust)

# Put distances in data frame for plotting and linear model
genus.robust.df$beta <- beta_genus_robust$distances

betadisp.genus   <- ggplot(genus.robust.df  , aes(x=fct_reorder(host_genus, beta, .fun = median, na.rm = TRUE), y=beta))+
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_tribe, shape=country),alpha = 0.8) +
  scale_colour_manual(values=tribe_colorfix) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(subtitle="Within-Group Variation per Host Genus", y="distance to centroid", x='', color ="host tribe") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )


hist(genus.robust.df$beta)
beta.lm.select <- lm(beta ~ host_genus, data=genus.robust.df)
summary(beta.lm.select)
summary(beta.lm.select)$adj.r.squared 
Anova(beta.lm.select)
plot(beta.lm.select,2)
hist(beta.lm.select$residuals) # extract the residuals
shapiro.test(residuals(beta.lm.select))

### betadisper sublocation 
# Simplify dataset / remove sublocations with 3 or less samples 
sublocation_to_keep <- df.sample %>%
  group_by(sublocation) %>%
  filter(n() > 3) %>%   # keep sublocation with >3 samples
  pull(sublocation)

# Subset the phyloseq object to keep only the selected host sublocation
sublocation.robust <- subset_samples(sample.species.rel, sublocation %in% sublocation_to_keep)
sublocation.robust.df  <- data.frame(sample_data(sublocation.robust))

# Tests if groups have different beta dispersion
sublocation.vegdist <- vegan::vegdist(t(otu_table(sublocation.robust)), index="bray")
beta_sublocation_robust <- betadisper(sublocation.vegdist, group = sublocation.robust.df$sublocation)
permutest(beta_sublocation_robust, permutations = 9999)

# Put distances in data frame for plotting and linear model
sublocation.robust.df$beta <- beta_sublocation_robust$distances

betadisp.sublocation   <- ggplot(sublocation.robust.df, aes(x=fct_reorder(sublocation, beta, .fun = median, na.rm = TRUE), y=beta))+
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_tribe, shape=country),alpha = 0.8) +
  scale_colour_manual(values=tribe_colorfix) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(subtitle="Within-Group Variation per Sublocation", y="distance to centroid", x='', color ="host tribe") +
  facet_wrap(~ country, scales = "free_x") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )


pdf(file.path(out_dir, "02_sample_div_beta_plot_betadisper.pdf"), width=12, height=6)
betadisp.family
betadisp.subfamily.country
betadisp.subfamily 
betadisp.subfamily.sorted
betadisp.subfamily.ordered
betadisp.tribe
betadisp.tribe.sorted 
betadisp.genus
betadisp.country
betadisp.location 
betadisp.sublocation 
dev.off()


#### betaframe adonis2 --------------------------------
# PERMANOVA
# by margin: Tests each term after accounting for all other terms (order does not matter)
# by terms: Adds terms sequentially from first to last (order matters)

# Quick check colinearity of variables with contingency table
table(betaframe$year, betaframe$collector)
table(betaframe$year, betaframe$storage)
table(betaframe$host_genus, betaframe$country)
table(betaframe$host_genus, betaframe$location)

set.seed(781)
# test adonis individual taxomomic ranks
adonis.fam  <- adonis2(phyloseq::distance(sample.species.rel, method = "bray")~   host_family , data=betaframe,  permutations = 9999, by = "margin")
adonis.fam  # 3.4%
#              Df SumOfSqs      R2     F Pr(>F)    
# host_family   1    2.239 0.03598 6.159  1e-04 ***

adonis.subfam  <- adonis2(phyloseq::distance(sample.species.rel, method = "bray")~   host_subfamily , data=betaframe,  permutations = 9999, by = "margin")
adonis.subfam  # 17%
#              Df SumOfSqs      R2     F Pr(>F)    
# host_subfamily   6   10.784 0.17334 5.5917  1e-04 ***

adonis.tribe  <-  adonis2(phyloseq::distance(sample.species.rel, method = "bray")~   host_tribe , data=betaframe,  permutations = 9999, by = "margin")
adonis.tribe # 24%
#              Df SumOfSqs      R2     F Pr(>F)    
# host_tribe  11   15.131 0.2432 4.5283  1e-04 ***

adonis.genus  <- adonis2(phyloseq::distance(sample.species.rel, method = "bray")~   host_genus , data=betaframe,  permutations = 999, by = "margin")
adonis.genus  # 36%
#              Df SumOfSqs      R2     F Pr(>F)    
# host_genus  23   22.297 0.3584 3.473  0.001 ***

adonis.species  <- adonis2(phyloseq::distance(sample.species.rel, method = "bray")~   host_species , data=betaframe,  permutations = 999, by = "margin")
adonis.species # 43%
#              Df SumOfSqs      R2     F Pr(>F)    
# host_species  34   26.502 0.42598 2.8811  0.001 ***


adonis.storage  <- adonis2(phyloseq::distance(sample.species.rel, method = "bray")~   host_tribe + storage + location  , data=betaframe,  permutations = 999, by = "margin")
adonis.storage
# host_tribe  11   12.570 0.20205 3.9454  0.001 ***
# storage      1    0.370 0.00595 1.2791  0.184    
# location     2    0.618 0.00994 1.0673  0.340 
# storage, location not significant when tribe is in the model

adonis.country  <- adonis2(phyloseq::distance(sample.species.rel, method = "bray")~ country , data=betaframe,  permutations = 9999, by = "margin")
adonis.country
# country    1    3.592 0.05773 10.11  1e-04 *** # alone 6%

adonis.location  <- adonis2(phyloseq::distance(sample.species.rel, method = "bray")~ location  , data=betaframe,  permutations = 9999, by = "margin")
adonis.location
# location   3    4.809 0.07729 4.5511  1e-04 *** # alone 8% 

adonis.sublocation  <- adonis2(phyloseq::distance(sample.species.rel, method = "bray")~ sublocation  , data=betaframe,  permutations = 999, by = "margin")
adonis.sublocation
# location   3    4.813 0.07737 4.5561  0.001 *** # alone 8% 


# Interaction full model, test host taxa and environmental predictors
full.model.adonis <- adonis2(phyloseq::distance(sample.species.rel, method = "bray") ~ host_subfamily + temp_week + sublocation   , 
                        data = betaframe, permutations = 999)
full.model.adonis

# adonis only Peru dataset
Peru.df <- as(sample_data(Peru.species.rel),"data.frame")
adonis.Peru.only <-  adonis2(phyloseq::distance(Peru.species.rel, method = "bray") ~   year +  elevation  +  location + host_genus  , data = Peru.df , permutations = 999, by = "margin")
adonis.Peru.only 
# year        1   0.3007 0.00964 1.1317  0.318    
# elevation   1   0.1999 0.00641 0.7525  0.741    
# location    1   0.1917 0.00615 0.7216  0.789    
# host_genus 16   9.6830 0.31041 2.2780  0.001 ***
# No influence of sampling year, elevation or location!

# adonis only Germany dataset
Germany.df <- as(sample_data(Germany.species.rel),"data.frame")
adonis.Germany.only <-  adonis2(phyloseq::distance(Germany.species.rel, method = "bray") ~   year +  kw  +  location + host_genus  , data = Germany.df , permutations = 999, by = "margin")
adonis.Germany.only
# year        2   0.6991 0.02550 1.2405  0.166    
# kw          1   0.4384 0.01599 1.5558  0.075 .  
# location    1   0.1485 0.00542 0.5271  0.957    
# host_genus  6   4.4572 0.16258 2.6364  0.001 ***
# No influence of sampling year, kw or location!

#### Robustness Check against over fitting  (filter species with ≥2 samples)
host_species_counts <- betaframe %>%
  group_by(host_species) %>%
  summarise(n_samples = n())

species_to_keep <- host_species_counts %>%
  filter(n_samples >= 2) %>%
  pull(host_species)

species.robust <- subset_samples(sample.species.rel, host_species %in% species_to_keep)
species.robust.df <- data.frame(sample_data(species.robust))

set.seed(123)
adonis.species.robust  <- adonis2(phyloseq::distance(species.robust, method = "bray")~ host_species  , data=species.robust.df,  permutations = 999)
adonis.species.robust
# Model     26   23.328 0.39512 3.3164  1e-04 ***  


set.seed(456)
adonis.species.robust.full  <- adonis2(phyloseq::distance(species.robust, method = "bray")~ host_species +  location + temp_week +  elevation   , data=species.robust.df,  permutations = 999)
adonis.species.robust.full #    
# species                 Model     26   23.303 0.39482 3.3122  0.001 ***
# species +location       Model     28   23.753 0.40246 3.127  0.001 ***
# species+ location +temp Model     29   24.159 0.40932 3.0825  0.001 ***
# spc+ loc +temp + elev   Model     30   24.428 0.41389 3.0129  0.001 ***
# Adding environmental factors to species model increases R2 only minimal 


# Robustness Check against over fitting  (filter genus with <2 samples)
set.seed(321)
adonis.genus.robust  <- adonis2(phyloseq::distance(genus.robust, method = "bray")~ host_genus , data=genus.robust.df,  permutations = 9999)
adonis.genus.robust
# Model     17   19.786 0.33141 4.1696  1e-04 ***
# robust genus 33% 

adonis.genus.robust.location  <- adonis2(phyloseq::distance(genus.robust, method = "bray")~  host_genus + location , data=genus.robust.df,  permutations = 999, by = "margin")
adonis.genus.robust.location
# host_genus  16   15.660 0.26229 3.5172  0.001 ***
# location     2    0.680 0.01140 1.2226  0.133   

adonis.robust.storage  <- adonis2(phyloseq::distance(genus.robust, method = "bray")~  host_tribe +  storage  + location , data=genus.robust.df,  permutations = 999, by = "margin")
adonis.robust.storage
# location   3    4.807 0.08051 4.5824  0.001 ***

# Quick predictor correlation matrix
cor(genus.robust.df [, c("temp_week", "kw", "elevation")]) # quick check predictors
cor_matrix <- cor(genus.robust.df[, c("temp_week", "kw", "elevation")])
r2_matrix <- cor_matrix^2
r2_matrix # Within the full dataset no major correlation of temp, elevation and calendar week 

# Check for colinearity using VIF
genus.robust.df$dummy <- seq_len(nrow(genus.robust.df)) 
lm_dummy <- lm(dummy ~  host_tribe +  location + elevation +  temp_week + kw    , data = genus.robust.df )
#summary(lm_dummy)
vif(lm_dummy)   # If variance inflation factor VIF > 5–10 indicate strong collinearity, consider removing variables
vif(lm_dummy, type = "predictor") # GVIF^(1/2Df) adjusted GVIF for categorical variables with multiple levels
alias(lm_dummy)
# working combinations without colinearity:
# species + temp, kw,
# genus   + temp, kw, elev
# tribe   + temp, kw, elev, location
# subfam  + temp, kw, elev, location

#                  GVIF Df GVIF^(1/(2*Df))
#host_tribe  74.231157 10        1.240307
#location   102.654525  3        2.163863
#elevation    8.472244  1        2.910712
#temp_week    9.849418  1        3.138378 max VIF 3.1 is still acceptable
#kw           2.028766  1        1.424347


# Use tribe + location + elevation +  temp_week + kw 
# Use all predictors in multivariate model
adonis.robust.tribe  <- adonis2(phyloseq::distance(genus.robust, method = "bray")~ host_tribe + location  +  elevation + kw  + temp_week , data=genus.robust.df,  permutations = 999, by = "margin")
adonis.robust.tribe
# host_tribe  10   10.996 0.18418 3.8660  0.001 ***
# location     3    1.173 0.01966 1.3752  0.036 *  
# elevation    1    0.296 0.00496 1.0414  0.390    
# kw           1    0.383 0.00642 1.3466  0.144    
# temp_week    1    0.405 0.00678 1.4223  0.099 .  
adonis.genus.robust.tribe.full  <- adonis2(phyloseq::distance(genus.robust, method = "bray")~ host_tribe + location +  elevation + kw + temp_week , data=genus.robust.df,  permutations = 999)


#### Marginal R2 / Variance Partitioning  ----
r2_df_clean <- adonis.robust.tribe %>%            
  as.data.frame() %>%                             # Extract marginal (partial) R2 values
  rownames_to_column("term") %>%                  
  filter(!term %in% c("Residual", "Total")) %>%   
  select(term, R2, p = `Pr(>F)`)                  # select and rename columns

# optional: rename factors 
r2_df_clean$term <- factor(r2_df_clean$term,
                           levels = c("host_tribe", "location", "kw", "elevation", "temp_week"),
                           labels = c("host tribe", "location", "calendar week", "elevation", "temperature"))

# Sort by explained variance
r2_df_clean <- r2_df_clean %>%
  arrange(desc(R2)) %>%  # sort by explained variance
  mutate(
    term = factor(term, levels = term),  # preserve sorted order for plotting
    signif = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***","**","*","")) )


# Show effect size as marginal (partial) R2
marginalR2.plot <- ggplot(r2_df_clean, aes(x = term, y = R2)) +
  geom_col(fill = "grey70") +
  geom_text(aes(label = signif), 
            vjust = -0.5, # place above the bar
            size = 6) +    # adjust size
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.1))) + # add space above bars
  labs(x = NULL, y = "partial R² (marginal)" ) +
  theme_line2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  )


#### Variance partitioning beta diversity
community.matrix <- as(otu_table(genus.robust), "matrix")
if (taxa_are_rows(genus.robust)) {
  community.matrix <- t(community.matrix) }

# Hellinger transformation
comm_hel <- decostand(community.matrix, method = "hellinger")
stopifnot(identical(rownames(comm_hel), rownames(genus.robust.df)))
genus.robust.df <- genus.robust.df[rownames(comm_hel), , drop = FALSE]

# Use only predictors without colinearity (checked above vif(lm_dummy))
X_host     <- genus.robust.df[, "host_tribe", drop = FALSE]
#X_location <- genus.robust.df[, "location", drop = FALSE]
#X_env  <- genus.robust.df[, c("temp_week"), drop = FALSE]
#X_env  <- genus.robust.df[, c("temp_week", "kw"), drop = FALSE]
#X_env  <- genus.robust.df[, c("temp_week", "kw", "elevation"), drop = FALSE]
X_env  <- genus.robust.df[, c("temp_week", "kw",  "location"), drop = FALSE]
#X_env  <- genus.robust.df[, c("temp_week", "kw",  "location" , "elevation"), drop = FALSE]

# Simple overview: Lists unique values of each predictor
list(
  host = sapply(X_host, function(x) length(unique(x))),
  env  = sapply(X_env, function(x) length(unique(x))) )

# Variance partitioning using varpart() function to decompos shared vs unique fractions visually
vp1 <- varpart(comm_hel, X_host, X_env)
vp1 # automatic warns if colinearity is detected
#plot(vp1) # gives a Venn diagram-style plot
plot(vp1,
     digits = 1,
     bg = c("skyblue", "salmon"),
     Xnames = c("Host", "Environment"))

# subfamily 10% both 6%  env 5% # (temp, kw, location)  = okay
# tribe     14% both 6%  env 4% # (temp, kw, location)  = okay 
# genus     16% both 10  env 1% # (temp, kw, location)  = collinearity detected with genus and location!!
# genus     19% both 7%  env 1% # (temp, kw, elevation) = okay, genus without location
# species   21% both 7%  env 1% # (temp, kw, elevation) = okay 

# using only temperature as env, shows minor influence of temp
# subfamily 12% both 4% temp 0%
# tribe     17% both 3% temp 1%
# genus     22% both 4% temp 0%
# species   24% both 4% temp 0%

# possible to check variance inflation, but less reliable than adjusted GVIF^(1/(2*Df))
vif.cca(rda(comm_hel ~ ., data = cbind(X_host, X_env)))

# Testing varpart fractions with rda (redundancy analysis) e.g. host controlling for environment
rda_host <- rda(comm_hel ~ host_tribe + Condition(kw + temp_week + location), data = genus.robust.df)
anova(rda_host, permutations = 999)
RsquareAdj(rda_host)

# Test environment controlling for host
rda_env <- rda(comm_hel ~ kw + temp_week + location +  Condition(host_tribe), data = genus.robust.df)
anova(rda_env, permutations = 999)
RsquareAdj(rda_env)

sink(file.path(out_dir, "02_sample_div_beta_stats_adonis_varpart.txt"))
adonis.fam
adonis.subfam
adonis.tribe
adonis.genus.robust
adonis.species.robust
adonis.storage
adonis.country
adonis.location
adonis.sublocation
adonis.robust.tribe
"Show effect size as marginal R2"
r2_df_clean
"X_host"
colnames(X_host)
"X_env"
colnames(X_env)
vp1
sink()


## 03 taxa abundance core / endo / family stripchart --------------

#### Core abundance subfamily # manually rebuild stripchart
core_stripchart.rebuild <- ggplot(core.melt, aes(x = factor(host_subfamily, levels = rev(subfamily_ordered)), y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +  # boxplot alpha = 0.8
  geom_jitter(aes(color = host_subfamily, shape = country),width = 0.2, size = 3, alpha = 0.8) +  # points on top
  facet_wrap(~genus) + # , scales = "free"
  coord_flip() +
  theme_grid() + theme(axis.text.x = element_text(size = 8)) +
  scale_colour_manual(values = subfamily_colorfix) + # alternative color country
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +
  labs(x = "", y = "genus core [%]", color = "host subfamily", fill = "host subfamily") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

# Stats kruskal.test core genera abundance per subfamily 
core_stats_subfam  <- run_kruskal(core.melt, "host_subfamily")
core_stats_subfam # BH correction for adj p

# posthoc pairwise tests for core genera difference among subfamilies
with(subset(core.melt, genus == "Apibacter"), 
     pairwise.wilcox.test(Abundance, host_subfamily, p.adjust.method = "BH", exact = FALSE))

with(subset(core.melt, genus == "Bartonella"), 
     pairwise.wilcox.test(Abundance, host_subfamily, p.adjust.method = "BH", exact = FALSE))

with(subset(core.melt, genus == "Orbus"), 
     pairwise.wilcox.test(Abundance, host_subfamily, p.adjust.method = "BH", exact = FALSE))

# Run kruskal.test function per country 
core_stats_country <- run_kruskal(core.melt, "country")
core_stats_country

sig_genera <- core_stats_country %>% filter(p_adj < 0.05)
sig_genera

# Make label for figure
core_country_label <- core_stats_country %>%
  mutate(label = paste0("p = ", p_adj, " ", signif))

# sort core genus names by their country ratio
sorted_taxa <- core.melt %>%
  group_by(genus, country) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>% 
  pivot_wider(names_from = country, values_from = total, values_fill = 0) %>%
  mutate(Ratio = Peru / (Germany + Peru)) %>%
  arrange(Ratio) %>%
  pull(genus)

core_stripchart.country.stats <- ggplot(core.melt, aes(x = country, y = Abundance)) +
  geom_jitter(aes(color = country),width = 0.2, size = 3, alpha = 0.6) +  # points on top
  geom_violin() +
  stat_summary(fun = mean,  geom = "crossbar", width = 0.2,             
               color = "black",  linewidth = 0.5  ) +
  facet_wrap(~factor(genus, levels = sorted_taxa)) +  # use sorted_taxa for facet order #facet_wrap(~genus)
  theme_grid() + theme(axis.text.x = element_text(size = 8)) +
  scale_colour_manual(values = country_colorfix) +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +
  labs(x = "", y = "rel ab [%]", color = "country") +
  geom_text(data = core_country_label,
    aes(x = 1.5, y = 0.95, label = label ), inherit.aes = FALSE, size = 3  ) 

core_stripchart.rebuild.tribe <- ggplot(core.melt, aes(x = genus, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +  # boxplot alpha = 0.8
  geom_jitter(aes(color = host_subfamily),width = 0.2, size = 3, alpha = 0.8) +  # points on top
  facet_wrap(~host_tribe) + # , scales = "free"
  #facet_grid(~host_subfamily~genus, scales = "free", space = "free")
  coord_flip() +
  theme_grid() + theme(axis.text.x = element_text(size = 8)) +
  scale_colour_manual(values = subfamily_colorfix) +
  #scale_fill_manual(values = subfamily_colorfix) +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +
  labs(x = "", y = "genus core [%]", color = "host subfamily", fill = "host subfamily")


#### ASV core abundance-------------------

sample.ASV.core <- subset_taxa(sample.ASV.rel, genus %in% core.genus)
data.frame(sort(taxa_sums(sample.ASV.core),decreasing = TRUE))
#sample.ASV.core.top <- prune_taxa(taxa_sums(sample.ASV.core) > 0.04, sample.ASV.core) # cumulative abundance 4% (not clear)
#sample.ASV.core.top <- prune_taxa(taxa_sums(sample.ASV.core) / nsamples(sample.ASV.core) >= 0.01,  sample.ASV.core) # ≥ 1% relative abundance in entire dataset
#sample.ASV.core.top <- prune_taxa(taxa_sums(sample.ASV.core) / nsamples(sample.ASV.core) >= 0.00024,  sample.ASV.core) # same as cum 4% or 0.024% relative abundance in entire dataset
sample.ASV.core.top <- prune_taxa(apply(otu_table(sample.ASV.core), 1, max) >= 0.02, sample.ASV.core) # ≥ 2% relative abundance in at least one sample # better!

sample.ASV.core.top.melt <- psmelt(sample.ASV.core.top)

ASV.bubble.country <- ggplot(sample.ASV.core.top.melt,aes(sampleID,fct_reorder(OTU, genus, .desc=TRUE))) +
  geom_point(aes (size = Abundance, color = factor(genus, levels = core_genus_abundance$genus)), alpha = 0.8) +
  theme_line2() +  
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  facet_grid(~country, scales = "free", space = "free") +
  scale_color_manual(values = core_palette, name = "genus") +
  scale_size(name = "Rel. Abundance [%]") +
  guides(color = guide_legend(override.aes = list(size = 5))) + # bigger dots in legend 
  xlab("sample ID") + ylab("") 

ASV.bubble.core.country <- ggplot(sample.ASV.core.top.melt, aes(sampleID,  fct_reorder(OTU, genus, .desc = TRUE) )) +
  geom_point(aes (size = Abundance, color = factor(genus, levels = core_genus_abundance$genus)), alpha = 0.8) +
  theme_line2() +
  facet_wrap(~genus + country, scales = "free", ncol = 4) +
  scale_color_manual(values = core_palette, name = "genus") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank() ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +  # bigger dots in legend
  labs(x="",y="", fill="genus core", size="rel. ab. [%]") 

pdf(file.path(out_dir, "03_sample_taxa_abundance_core_ASV.pdf"), width=12, height=6)
core_stripchart.rebuild
core_stripchart.country.stats
core_stripchart.rebuild.tribe
ASV.bubble.country
ASV.bubble.core.country
dev.off()

# Cleanup pipeline 
rm(sample.ASV.core.top.melt)

sink(file.path(out_dir, "03_sample_taxa_abundance_core_ASV.txt"))
"Kruskal.test core genus per subfamily"
data.frame(core_stats_subfam)
"Kruskal.test core genus per country"
data.frame(core_stats_country)
sink()


#### Endosymbiont abundance  ----------------

endosymbiont.rel <- subset_taxa(sample.species.rel , genus == "Wolbachia" |
                                  #genus=="Spiroplasma" |
                                  #genus=="Spiroplasma-like" |
                                  #genus=="Arsenophonus" |
                                  #genus=="Flavobacterium"|
                                  #genus=="Pectobacterium"|
                                  #genus=="Mesoplasma"  |
                                  #genus=="Rickettsia" |
                                  genus=="Kinetoplastibacterium"  |
                                  genus=="Sodalis"  |
                                  #genus=="Rickettsiaceae spc"|
                                  family=="Holosporaceae"|
                                  family=="Anaplasmataceae"|
                                  #family=="Rickettsiaceae"|
                                  #order == "Rickettsiales"|
                                  phylum=="Tenericutes")

data.frame(sort(taxa_sums(endosymbiont.rel), decreasing = FALSE))
data.frame(sort(sample_sums(endosymbiont.rel), decreasing = FALSE))

# Remove low abundances for figure, set everything below 1% to zero
otu_table(endosymbiont.rel)[otu_table(endosymbiont.rel )<0.01 ]<-0


endosymbiont.rel.melt <- psmelt(endosymbiont.rel)
# alternative: filter low abundant taxa
# endosymbiont.rel.filtered <- endosymbiont.rel.melt %>%   filter(Abundance >= 0.01)


endosymbiont.rel_stripchart <- ggstripchart(endosymbiont.rel.melt, "host_genus", "Abundance", 
                                            facet.by = "genus", color = "host_subfamily", 
                                                   ) +
                                            geom_boxplot(aes(fill = host_subfamily, color = host_subfamily), 
                                            alpha = 0.6,  outlier.shape = NA) +
                                            theme_grid() +  coord_flip()+
                                            scale_colour_manual(values=subfamily_colorfix) + scale_fill_manual(values=subfamily_colorfix)
  

endosymbiont.rel_stripchart1 <- ggstripchart(endosymbiont.rel.melt, "host_subfamily", "Abundance", 
                                            facet.by = "genus", color = "host_subfamily",
                                            add = "boxplot") +
                                            theme_grid() +  coord_flip()+
                                            scale_colour_manual(values=subfamily_colorfix)


sym_abundance <- endosymbiont.rel.melt %>%
  group_by(genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance))  # Sort by abundance (most abundant first)

symbiontCount = length(unique(endosymbiont.rel.melt$genus))
symb_palette <- setNames(taxaPalette(symbiontCount), sym_abundance$genus)  
                                          
                                             
endo.bar2 <- ggplot(endosymbiont.rel.melt,aes(x=Sample, y=Abundance, fill = genus))+
  geom_bar( stat="identity")+
  scale_fill_manual(values = symb_palette)+ 
  theme_line2() + #theme(axis.text.x = element_text(angle = 60, hjust = 1, face="italic") ) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + # remove sample names 
  labs(x="",y="rel. ab. [%]", fill="endosymbiont") +
  #facet_wrap(~host_genus, scale = "free_y" ) + coord_flip() +
  facet_wrap(~host_genus+country, scale = "free_y" ) + coord_flip() +
  #facet_grid(~host_genus, scales = "free", space = "free")
  guides(fill = guide_legend(reverse = TRUE))  # reverse legend order
  

# Average abundance per genus
avg_ab_genus <- endosymbiont.rel.melt %>%
  group_by(host_genus, genus, country) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(host_genus) %>%
  mutate(total_abundance = sum(mean_abundance)) %>%
  ungroup()


endo.bar3 <- ggplot(avg_ab_genus ,aes(x= mean_abundance, y=fct_rev(fct_reorder(host_genus, total_abundance)), fill = genus)) + #fct_reorder(genus, Abundance)
  geom_bar(#position="fill",
    stat="identity") +
  scale_fill_manual(values=symb_palette) +
  theme_line2() + theme(axis.text.y = element_text(face = 'italic') ) +
  scale_x_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +
  labs(x="rel. ab. [%]",y="", fill="endosymbiont") +
  guides(fill = guide_legend(reverse = TRUE)) # reverse legend order

pdf(file.path(out_dir, "03_sample_taxa_abundance_endosymbiont.pdf"), width=12, height=6)
endosymbiont.rel_stripchart
endosymbiont.rel_stripchart1
endo.bar2
endo.bar3
dev.off()


## 05 ggtree ASV alignment all taxa -----------------
# use reduced dataset sample.filter for tree

# percent of total reads remaining 
percent_retained <- (sum(otu_table(sample.filter)) / sum(otu_table(sample.species))) * 100
percent_retained # still represents 98.0337% of total dataset

# Check if core taxa remain in sample.filter
sample.filter.rel.g <- aggregate_taxa(sample.filter.rel, "genus")
data.frame(sort(taxa_sums(sample.filter.rel.g)))

filter.core.genus <- core_members(sample.filter.rel.g , detection = 0.01, prevalence = 20/100)
filter.core.genus.high <- core_members(sample.filter.rel.g, detection = 0.05, prevalence = 10/100) 
filter.core.genus.low <- core_members(sample.filter.rel.g, detection = 0.002, prevalence = 40/100) 

# sample.filter check individual taxa
table(tax_table(sample.filter)[, "phylum"], exclude = NULL)
check_taxa <- subset_taxa(sample.filter, genus == "Pseudomonas")
sum(taxa_sums(check_taxa))
sort(data.frame(sum = taxa_sums(check_taxa)), decreasing = FALSE)
sort(data.frame(sum = sample_sums(check_taxa)), decreasing = FALSE)

# select most abundant ASVs for each genus (manual version)
asv_data <- data.frame(
  ASV = taxa_names(sample.filter), 
  reads = taxa_sums(sample.filter),   
  genus = tax_table(sample.filter)[, "genus"] )

selected_asvs <- asv_data %>%
  group_by(genus) %>%
  slice_max(order_by = reads, n = 1, with_ties = FALSE) %>% 
  ungroup()

abundant_ASVs <- selected_asvs$ASV

# Alternative: Collapse OTU table by genus and select most abundant ASV from ps object
sample.filter.genus <- tax_glom(sample.filter, taxrank = "genus", NArm=F) # glom by genus but keeps ASV name info
data.frame(sort(taxa_sums(sample.filter.genus)))
asv_names <- rownames(otu_table(sample.filter.genus))

# compare ASV name list both methods should give same outcome
setequal(asv_names, abundant_ASVs)

# Import sequences and rename to ASV1 etc.
asv_sequences <- readDNAStringSet("asvs.merge.fa")
names(asv_sequences) <- gsub(";size=.*", "", names(asv_sequences))
head(names(asv_sequences))

# Match most abundant ASVs per genus to sequences
matched_all_seqs <- asv_sequences[names(asv_sequences) %in% asv_names]
sample.filter.seq <- merge_phyloseq(sample.filter.genus, matched_all_seqs)
data.frame(sort(taxa_sums(sample.filter.genus)))
data.frame(sort(taxa_sums(sample.filter.seq)))

# Export ASV sequences to make tree
writeXStringSet(matched_all_seqs , file = "tree/ASV_export.fasta")
# Use to make IQ tree or RAxML tree with SINA 1.2.12 (https://www.arb-silva.de/aligner/)
# Parameter: SSU, min identity with query 0.90, Number of neighbours 1
# tree building RAxML with GTR model and Gamma likelihood
# rename tree file in nwk and reimport into R (arb-silva.nwk)

# Find outgroup for root tree
check_taxa <- subset_taxa(sample.filter.seq , phylum == "Bacteroidetes") # Fusobacteria Deinococcus-Thermus Verrucomicrobia
data.frame(genus = as.character(tax_table(check_taxa)[names(sort(taxa_sums(check_taxa), decreasing = T)), "genus"]),
  Taxa_Sum = sort(taxa_sums(check_taxa), decreasing = T))
OG <- "ASV99" # define outgroup ASV name # ASV3819 # ASV937

# Load IQ tree with references
# Combine with reference sequences, run IQ tree, import nwk tree here
IQ_tree <- read.tree("tree/IQtree.nwk")
IQ_tree.root <- root(IQ_tree, outgroup = OG, resolve.root = TRUE) 
sample.filter.IQtree <- merge_phyloseq(sample.filter.seq, IQ_tree.root )
IQ.tree <-  ggtree(sample.filter.IQtree, layout="circular") + geom_tiplab(size = 3, align=TRUE, aes(label=genus, color=phylum)) + ggtitle("IQ tree")
#IQ.tree.ASV <-  ggtree(sample.filter.IQtree, layout="circular") + geom_tiplab(size = 3, align=TRUE, aes(color=phylum)) 

# Load SINA tree with references
# Combine with SILVA reference sequences, run RAxML tree, import nwk tree here
SINA_tree <- read.tree("tree/arb-silva.nwk")
SINA_tree.root <- root(SINA_tree, outgroup = OG, resolve.root = TRUE) 
sample.filter.SINA_tree <- merge_phyloseq(sample.filter.seq, SINA_tree.root )
SINA.tree <-  ggtree(sample.filter.SINA_tree, layout="circular") + geom_tiplab(size = 3, align=TRUE, aes(label=genus, color=phylum)) + ggtitle("SINA RAxML tree")

# Choose best tree clean up tip labels
tree.root  <- drop.tip(SINA_tree.root, setdiff(SINA_tree.root$tip.label, asv_names))
sample.filter.tree <- merge_phyloseq(sample.filter.seq, tree.root )

### ggtree ad core info 
tip_data <- data.frame(label = tree.root$tip.label)
tip_data$core <- ifelse(tip_data$label %in% core.genus.ASV.names, "Core", "Non-Core")
core_tips <- tip_data[tip_data$core == "Core", ]

# Ad phylum genus names for core tips
taxa_data <- as.data.frame(tax_table(sample.filter.tree))
core_tips$genus_name <- taxa_data$genus[match(core_tips$label, rownames(taxa_data))]
core_tips$phylum_name <- taxa_data$phylum[match(core_tips$label, rownames(taxa_data))]
core_tips

tree.core <- ggtree(tree.root, layout="circular") + 
  geom_tiplab(size = 3, align = TRUE, aes(label = ifelse(label %in% core_tips$label, 
                                 core_tips$genus_name[match(label, core_tips$label)], ""))) # Only show core genera
  
# Select bee associated taxa 
select_bee_genus_list <- c("Gilliamella", "Snodgrassella", "Lactobacillus", "Bifidobacterium", "Frischella", "Bartonella", "Apibacter", "Apilactobacillus", "Commensalibacter")
bee_taxa <- taxa_data[rownames(taxa_data) %in% select_bee_genus_list | taxa_data[,"genus"] %in% select_bee_genus_list , ]

bee_tips <- data.frame(genus = bee_taxa[,"genus"],  # Extract bee core genus names
                       label = rownames(bee_taxa) ) # Extract ASV/tip labels


tree.bee <- ggtree(sample.filter.tree, layout="circular") +
  geom_tree() + geom_tiplab(size = 3, align = TRUE, aes(label = ifelse(label %in% bee_tips$label, genus, ""))) +  # Only show core genera
  theme(legend.position = "right") 

tree.bee.asteriks <- ggtree(sample.filter.tree, layout="circular") +
  geom_tree() + geom_tiplab(size = 3, align = TRUE, aes(label = ifelse(label %in% bee_tips$label, "*", "")), color = "red") +  # Show * for matching tips
  theme(legend.position = "right")


pdf(file.path(out_dir, "05_ggtree_core.pdf"), width=12, height=6)
IQ.tree
SINA.tree
tree.core
tree.bee
tree.bee.asteriks
dev.off()

# Cleanup pipeline
rm(IQ.tree)
rm(SINA.tree)
rm(tree.bee)
rm(tree.bee.asteriks)
rm(asv_sequences)
rm(matched_all_seqs)
rm(sample.filter.seq)

rm(sample.filter.IQtree)
rm(sample.filter.SINA_tree)


####  ggtree round tree  ---------------
taxonomy_df <- as.data.frame(tax_table(sample.filter.tree))
tree_tips <- phy_tree(sample.filter.tree)$tip.label
all(tree_tips %in% rownames(taxonomy_df)) # Check if all okay
taxonomy_df <- taxonomy_df[tree_tips, , drop = FALSE]

# Adding label, total, relative and log abundance
taxonomy_df$label <- rownames(taxonomy_df)
taxonomy_df$dummy <- seq_len(nrow(taxonomy_df))
otu_data <- as.matrix(otu_table(sample.filter.tree))
total_abundance <- rowSums(otu_data)
taxonomy_df$total_ab <- total_abundance

# ad rel ab (rough calculation without sample normalization)
taxonomy_df$rel_ab   <- total_abundance / sum(total_abundance) * 100

# ad rel ab normalized per sample (to account for different sampling depth)
otu_rel <- as.matrix(otu_table(transform_sample_counts(sample.filter.tree, function(x) x / sum(x)) ))
total_abundance_rel <- rowSums(otu_rel)
taxonomy_df$rel_ab2 <- total_abundance_rel / sum(total_abundance_rel) * 100

# ad log ab
taxonomy_df$log_ab <- log10(total_abundance)

head(taxonomy_df)



sink(file.path(out_dir, "05_ggtree_layers.txt"))
taxonomy_df
sink()

# phyla color by order for ggtree
phylum_order <- names(sort(tapply(taxonomy_df$total_ab, taxonomy_df$phylum, sum), decreasing = TRUE))
phylumGradient <- colorRampPalette(brewer.pal(12, "Paired")[2:10])(length(phylum_order)) # gradient with selected bright colors
color_phyla <- setNames(phylumGradient, phylum_order) # scale_fill_manual(values = color_phyla)

# prevalence sample.filter.tree
prev_cutoff <- 20 # Prevalence cutoff 20 reads minimum
prevalence_df <- data.frame(
  ASV = rownames(otu_table(sample.filter.tree)),
  prevalence = apply(otu_table(sample.filter.tree), 1, function(x) sum(x > prev_cutoff)) )
# ad relative prevalence in percent
prevalence_df$rel_prev <- prevalence_df$prevalence / ncol(otu_table(sample.filter.tree))
head(prevalence_df)


Satyrinae <- subset_samples(sample.filter.tree, host_subfamily=="Satyrinae")
Dismorphiinae <- subset_samples(sample.filter.tree, host_subfamily=="Dismorphiinae")
Pierinae <- subset_samples(sample.filter.tree, host_subfamily=="Pierinae")
Heliconiinae <- subset_samples(sample.filter.tree, host_subfamily=="Heliconiinae")
Coliadinae <- subset_samples(sample.filter.tree, host_subfamily=="Coliadinae")
Nymphalinae <- subset_samples(sample.filter.tree, host_subfamily=="Nymphalinae")

Dismorph_Nymph <- merge_phyloseq(Dismorphiinae, Nymphalinae)

Satyrinae_prevalence_df <- data.frame(
  ASV = rownames(otu_table(Satyrinae)),
  prevalence = apply(otu_table(Satyrinae), 1, function(x) sum(x > prev_cutoff)) )
Satyrinae_prevalence_df$rel_prev <- Satyrinae_prevalence_df$prevalence / ncol(otu_table(Satyrinae))

Pierinae_prevalence_df <- data.frame(
  ASV = rownames(otu_table(Pierinae)),
  prevalence = apply(otu_table(Pierinae), 1, function(x) sum(x > prev_cutoff)) )
Pierinae_prevalence_df$rel_prev <- Pierinae_prevalence_df$prevalence / ncol(otu_table(Pierinae))

Heliconiinae_prevalence_df <- data.frame(
  ASV = rownames(otu_table(Heliconiinae )),
  prevalence = apply(otu_table(Heliconiinae ), 1, function(x) sum(x > prev_cutoff)) )
Heliconiinae_prevalence_df$rel_prev <- Heliconiinae_prevalence_df$prevalence / ncol(otu_table(Heliconiinae))

Coliadinae_prevalence_df <- data.frame(
  ASV = rownames(otu_table(Coliadinae)),
  prevalence = apply(otu_table(Coliadinae), 1, function(x) sum(x > prev_cutoff)) )
Coliadinae_prevalence_df$rel_prev <- Coliadinae_prevalence_df$prevalence / ncol(otu_table(Coliadinae))

Nymphalinae_prevalence_df <- data.frame(
  ASV = rownames(otu_table(Nymphalinae)),
  prevalence = apply(otu_table(Nymphalinae), 1, function(x) sum(x > prev_cutoff)) )
Nymphalinae_prevalence_df$rel_prev <- Nymphalinae_prevalence_df$prevalence / ncol(otu_table(Nymphalinae))

Dismorphiinae_prevalence_df <- data.frame(
  ASV = rownames(otu_table(Dismorphiinae)),
  prevalence = apply(otu_table(Dismorphiinae), 1, function(x) sum(x > prev_cutoff)) )
Dismorphiinae_prevalence_df$rel_prev <- Dismorphiinae_prevalence_df$prevalence / ncol(otu_table(Dismorphiinae))

combined_rel_prev.df <- data.frame(
  Satyrinae = Satyrinae_prevalence_df$rel_prev,
  #Dismorphiinae = Dismorphiinae_prevalence_df$rel_prev,
  Pierinae = Pierinae_prevalence_df$rel_prev,
  Heliconiinae = Heliconiinae_prevalence_df$rel_prev,
  Coliadinae = Coliadinae_prevalence_df$rel_prev,
  Nymphalinae = Nymphalinae_prevalence_df$rel_prev )

rownames(combined_rel_prev.df) <- rownames(taxonomy_df)


# Core Tree fan including layers
fan.core <- ggtree(tree.root, layout="fan", size=0.2, open.angle=30) +
  geom_tiplab(aes(label = ""),  align = TRUE, offset = 0.01, linetype = "dotted")

fan.prev_all <- gheatmap(fan.core, combined_rel_prev.df, 
                            width = 0.3,
                            offset = 0.03,
                            colnames_position = "bottom",
                            colnames_angle =90, font.size = 2.8,
                            hjust = 1,
                            colnames_offset_y = 0) +
  scale_fill_viridis_c(option = "inferno", name = "Prevalence") 
  
fan.core.bar <- fan.prev_all + new_scale_fill() +
  geom_fruit(
    data=taxonomy_df, 
    geom=geom_bar, 
    mapping=aes(y=label, x=rel_ab, fill = phylum),
    pwidth=0.5, offset = 0.44, # 0.2
    orientation="y", stat="identity",
    axis.params = list(axis = "x", #title = "rel", #title.size = 3,
      text.size = 2.8, hjust = 0.5, vjust = 1, line.color = "grey"  ),
    grid.params = list(size = 0.2, color = "grey")) + scale_fill_manual(values = color_phyla)
 

fan.core.label <- fan.core.bar  + geom_tiplab(size = 3, align = T, offset = 0.19, #nudge_y = 0.1, # offset 0.13
            aes(label = ifelse(label %in% core_tips$label, core_tips$genus_name[match(label, core_tips$label)], "")), linetype = NA)   # Only show core genera

fan.core.point <- fan.core.label + geom_point(aes(color = taxonomy_df$phylum, x = x + 0.01),
                                              data = ~ .x[.x$isTip, ],  
                                              size = 2,
                                              show.legend = FALSE) +
                                              scale_color_manual(values = color_phyla)
  
fan.core.bee <- fan.core.point  + 
    geom_tiplab(size = 6, align = TRUE, offset = 0.001, # nudge_y = -0.3, # offset = 0.14
                aes(label = ifelse(label %in% bee_tips$label, "\u2022", "")), # * or bold asteriks ✱ \u2022
                linetype = NA, # NA dotted
                family = "mono")   # Only show * for matching tips  
  
fan.core.annotation <- fan.core.bee + annotate("text",  x = 0,  y = 0, 
  label = "rel. ab [%]", size = 2.8, angle = 0, hjust = -4,  vjust = 2.5) +
  theme(legend.key.size = unit(0.6, "lines"),  # smaller legend (usually 1 line)
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.5), # place legend closer to plot
        legend.justification = c(0, 0.5),     # Right-center of the legend box
        plot.margin = margin(0, 40, 0, 0) )   # increase right margin for legend
        # figure border highlight for diagnostics
        # theme(panel.border = element_rect(color = "blue", fill = NA, linewidth = 1)) +
        # theme(plot.background = element_rect(color = "red", fill = NA, linewidth = 1)) 
  
  
#### ggtree square tree  -----
#square.core <- ggtree(tree.root, size=0.2)
head(taxonomy_df)
all(tree.root$tip.label %in% rownames(taxonomy_df)) #check if labels are okay
all(tree.root$nodes %in% rownames(taxonomy_df)) #check if labels are okay

# Ad a row label before merge tree with tax data
square_taxa_data <- tax_table(sample.filter.tree) %>% 
  as.data.frame() %>%
  rename_with(~ paste0("tree_", .x)) %>%  # prefix to avoid conflicts
  rownames_to_column("label")             # add tip labels for merge

# Ad taxdata to tree and ad tax info as tip label
square.core <- ggtree(tree.root, size = 0.2) %<+% square_taxa_data +
  geom_tiplab(aes(label = tree_genus, color = tree_phylum),
    align = TRUE, size = 2.5, hjust = 1, offset = 0.08, #family='mono'#, linesize = .7, offset = 0.035 0.005
    linetype = NA, # NA  "dotted"  linetype = "dotted"
    show.legend = FALSE, lineend = "round", lineheight = 0.9  ) +
  scale_color_manual(values = color_phyla)

square.prev_all <- gheatmap(square.core, combined_rel_prev.df, 
                            width = 0.2,
                            offset = 0.084, # 0.035, 0.06
                            colnames_position = "top",
                            colnames_angle =60, font.size = 3,
                            hjust = 0,
                            colnames_offset_y = 0.5) +
  scale_fill_viridis_c(option = "inferno", name = "Prevalence") +
  theme(  plot.margin = margin(t = 45, r = 10, b = 10, l = 10) # Add space at the top
    ) + coord_cartesian(clip = "off")


square.core.bar <- square.prev_all + new_scale_fill() +
  geom_fruit(data=taxonomy_df, geom=geom_bar, 
    mapping=aes(y=label, x=rel_ab, fill=phylum), 
    pwidth=0.5, # 0.5 50% space for geom bar
    offset = 0.45, # 0.25 0.18 0.31
    orientation="y", stat="identity",
    axis.params = list(axis = "x", title = NULL, text.size = 3, hjust = 0.5, vjust = 1.5, line.color = "black" ),
    grid.params = list(size = 0.2, color = "black") ) +  # coord_cartesian(clip = "off") + 
  scale_fill_manual(values = color_phyla) +
  labs(x = "Rel. Abundance [%]") +
  theme(axis.title.x = element_text(size = 10, vjust = -1, hjust = 0.85 )) 


square.label <- square.core.bar  + geom_tiplab(size = 2.5, align = T, offset = 0.19, #nudge_y = 0.3 , # 0.13 0.085
                                   aes(label = ifelse(label %in% core_tips$label, core_tips$genus_name[match(label, core_tips$label)], "")), linetype = NA)  # Only show core genera

square.bee.core <- square.label + geom_tiplab(size = 5, align = TRUE, offset = 0.08, #nudge_y = -0.4,
              aes(label = ifelse(label %in% bee_tips$label, "\u2022", "")),
              linetype = NA, family = "mono")   # Only show * for matching tips

square.point <- square.bee.core + geom_point(aes(color = taxonomy_df$phylum, x = x + 0.005),
                                              data = ~ .x[.x$isTip, ],  
                                              size = 2, show.legend = FALSE) +
                                              scale_color_manual(values = color_phyla)

square.annotation <-  square.point +   geom_treescale(fontsize=3, linesize=0.1, x=0, y=-2) +
  theme(legend.position.inside=c(1.3, 0.5),
        #legend.background=element_rect(fill=NA), #legend.title=element_blank(),
        #legend.text=element_text(size=12), #legend.spacing.y = unit(0.2, "cm"),
  ) 

# convert inch to cm -> 16 / 2.54 
pdf(file.path(out_dir, "05_ggtree_layers.pdf"), width= 20 / 2.54 , height= 16 / 2.54)
fan.core.label
fan.core.annotation
square.annotation 
dev.off()


## 06 Country Color Chord Alluvial  --------------

# Network Chord diagramm country vs core taxa
countrymat <- otu_table(merge_samples(core.genus.rel,group="country"))
plotweb(data.frame(t(otu_table(countrymat))))
#chordDiagram(countrymat) # works inverse taxa country
matcountry <- t(as.matrix(countrymat))
#chordDiagram(matcountry) # does not work

# Make long format for chord diagram and alluvial plot 
df_long <- as.data.frame(as.table(as.matrix(matcountry))) 
colnames(df_long) <- c("Taxa", "Sample", "Abundance")
#chordDiagram(df_long) # country taxa works

# ad customn color
countryColors <- brewer.pal(n = min(12, nrow(countrymat)), name = "Dark2") # Colors for country Dark2
names(countryColors) <- c(rownames(countrymat))
sectorColorsCountry <- c(core_palette, countryColors)

# Network connection dprime
dfun(t(matcountry)) # degree distribution, how connected nodes are (number of partners)
# dprime (d′) Normalized Specialization Index. It ranges from 0 (no specialization) to 1 (perfect specialist).

H2fun(t(matcountry), H2_integer = F) 
# 𝐻′₂ Standardized Specialization Index (generalization to the entire interaction web)


### 06 Country Alluvial plot 
# Sort dataframe like in plotweb Compute the ratio of abundance per taxon across countries 
taxa_abundance <- as.data.frame(matcountry) %>%
  mutate(Ratio = Peru / (Germany + Peru)) %>%
  arrange(Ratio)

sorted_taxa <- rownames(taxa_abundance)

df_sorted <- df_long %>% 
mutate(Taxa = factor(Taxa, levels = sorted_taxa),  # Sort Taxa based on sorted_taxa
       Sample = factor(Sample, levels = sort(levels(Sample))) ) %>% # Alphabetical order for countries
       arrange(Sample, Taxa)  # Arrange by Sample first, then Taxa


# Sorted alluvial plot
Allu_sorted <- ggplot(df_sorted, aes(axis1 = Taxa, axis2 = Sample, y = Abundance)) +
  geom_alluvium(aes(fill = Taxa), width = 0.2, knot.pos = 0.3) +  # Flow lines
  geom_stratum(fill = "gray80", color = "black") +  # Nodes
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) + # Labels
  scale_x_discrete(limits = c("Taxa", "Sample")) +  # Set axis names
  theme_heat() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),panel.border = element_blank())  
  
AlluTaxaCountry <- ggplot(df_sorted, aes(axis1 = Taxa, axis2 = Sample, y = Abundance)) +
  geom_alluvium(aes(fill = Taxa), width = 0.3, knot.pos = 0.3, alpha = 0.6) +  
  geom_stratum(fill = NA, color = "black", width = 0.3)  +  #geom_stratum(fill = "gray80", color = "black") +  
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +  
  scale_fill_manual(values=core_palette) +
  scale_x_discrete(limits = c("Taxa", "Sample"), expand = c(0.2, 0.2)) +
  theme_heat() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),panel.border = element_blank(),
        legend.position = "none",
        plot.margin = margin(2, 1, 2, 2)) 
 
AlluTaxaCountryreverse <- ggplot(df_sorted, aes(axis1 = Sample, axis2 = Taxa, y = Abundance, fill = Taxa)) +
  geom_alluvium(aes(fill = Taxa), width = 0.3, knot.pos = 0.3, alpha = 0.6) +  
  geom_stratum(color = "black", width = 0.3, alpha = 0.6)  +  
  geom_text(stat = "stratum", aes(label = after_stat(stratum), angle = after_stat(ifelse(after_stat(x) < 1.5, 90, 0))  # rotate only axis1 labels
    ),  hjust = 0.5,   size = 4 ) +
  scale_fill_manual(values=core_palette) + 
  scale_x_discrete(limits = c("Taxa", "Sample"), expand = c(0.2, 0.2), guide = guide_axis(n.dodge = 1)) +
  theme_heat() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2) )

AlluTaxaCountry2 <- ggplot(df_sorted, aes(axis1 = Taxa, axis2 = Sample, y = Abundance)) +
  geom_alluvium(aes(fill = Taxa), width = 0.3, knot.pos = 0.3, alpha = 0.6) +  
  geom_stratum(fill = NA, color = "black", width = 0.3)  +    
  geom_stratum(
    data = subset(df_sorted, !is.na(Taxa)),
    aes(fill = Taxa),
    width = 0.3,
    color = "black" ) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
  scale_fill_manual(values=core_palette) +
  scale_x_discrete(limits = c("Taxa", "Sample")) +
  theme_heat() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),panel.border = element_blank(), legend.position = "none") 


# For dual plot normalize for rel abundance
df_sorted_norm <- df_sorted %>%
  group_by(Sample, Taxa) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  complete(Sample, Taxa, fill = list(RelAbundance = 0)) %>%  # Ensure all taxa appear in all groups
  filter(RelAbundance > 0)  # Remove any zero values to clean the visualization

AlluTaxaCountry_dual <- ggplot(df_sorted_norm, aes(x = Sample, stratum = Taxa, alluvium = Taxa, y = RelAbundance, fill = Taxa, label = Taxa)) +
  geom_flow(alpha = 0.5) +  geom_stratum(alpha = 0.8, color = "black") +  geom_text(stat = "stratum", size = 3) + 
  theme_line2() + scale_fill_manual(values=core_palette) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),panel.border = element_blank(), legend.position = "none") +
  labs(x = "", y = "", fill = "Taxa")

AlluCountry <- ggplot(df_sorted, aes(axis1 = Taxa, axis2 = Sample, y = Abundance)) +
  geom_alluvium(aes(fill = Sample), width = 0.2, knot.pos = 0.3) +  # Color flows by Sample
  geom_stratum(fill = "gray80", color = "black", width = 0.2)  +  # Color strata by Sample
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
  scale_fill_manual(values = sectorColorsCountry) +  # Define colors based on Sample
  scale_x_discrete(limits = c("Taxa", "Sample")) +
  theme_heat() +  
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank()) 
  #guides(fill = "none")  # Remove legend if not needed

pdf(file.path(out_dir, "06_CountryAlluvial_color.pdf"), width=8, height=8)
Allu_sorted
AlluTaxaCountry
AlluTaxaCountryreverse
AlluTaxaCountry2
AlluTaxaCountry_dual
AlluCountry
dev.off()



## 07 Final Figures  -----
empty.plot <- ggplot() + theme_void()

betadisp.subfamily.theme <- betadisp.subfamily.ordered  + theme(
  legend.key.size = unit(0.9, "lines") )

PCoA.genus.wrap.tribe_noleg <- PCoA.genus.wrap.tribe + theme(legend.position = "none")

core.genus.boxplot.subfam.no.text <- core.genus.boxplot.subfam + theme(axis.text.y = element_blank()) 

# arrange with ggarrange
fig1.arrange <- ggarrange(empty.plot , alpha.shannon2,  PCoA.host.country ,    labels = c('a', 'b', 'c'),
                          common.legend = T, legend = "right", ncol = 3,   nrow = 1, align = "h", widths = c(1, 1, 1) )

fig1.arrangeb <- ggarrange(empty.plot , alpha.shannon2,  alpha.shannon.tribe.sort ,    labels = c('a', 'b', 'c'),
                          common.legend = T, legend = "right", ncol = 3,   nrow = 1, align = "h", widths = c(1, 1, 1) )


fig2.betadisp.tribe.wrap <- ggarrange( betadisp.subfamily.theme, PCoA.host.overview.tribe ,   labels = c('a', 'b'),
                                       common.legend = T, legend = "right", ncol = 2,   nrow = 1, align = "h", widths = c(1.2 ,2) )

# arrange with patchwork
fig3.gen.core.abundance <- (order.core + theme(legend.position = "right"))+ 
  (sample.order.top.bar.flip + theme(legend.position = "right") ) +
  (order.fam.country.mirrored.genus + theme(legend.position = "none") ) +
  plot_layout(ncol = 3, guides = "collect", widths = c(0.5,0.8, 2)) &  # , guides = "collect"
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14, hjust = 0),
        plot.tag.position = c(0.01, 0.98),
        plot.margin = margin(3, 3, 3, 3),
        legend.key.size = unit(0.9, "lines") )      # key box size

# Fig 5 net core
fig5.net.core.comp <- 
  (AlluTaxaCountryreverse) +
  (core.genus.abundance.subfam + theme(legend.position = "none") ) +
  (core.genus.boxplot.subfam.no.text + theme(legend.position = "none") ) +
  (core.shannon+ theme(legend.position = "none") ) +
  plot_layout(ncol = 4, widths = c(0.9, 0.7, 0.7, 1.4),nrow=1, height = c(1,1)) &  
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14, hjust = 0),
        plot.tag.position = c(0.01, 0.98),
        plot.margin = margin(4, 4, 4, 4))  # styling tags




# supplement S1 alpha.country patchwork
s1.plot_alpha.country <- (alpha.shannon.tribe.sort + theme(legend.position = "none"))+ 
  (alpha.shannon.country  + theme(legend.position = "none") ) +
  (alpha.subfamily.country  + theme(legend.position = "right") ) +
  plot_layout(ncol = 3, widths = c(2, 1,2)) &  # , guides = "collect"
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14, hjust = 0),
        plot.tag.position = c(0.01, 0.98),
        plot.margin = margin(3, 3, 3, 3),
        legend.spacing.x = unit(0.3, 'cm'),        # horizontal spacing
        legend.spacing.y = unit(0.2, 'cm') )        # vertical spacing


# supplement S2 alpha elevation patchwork
s2.plot_alpha.elevation <- (
  (sample.shannon.elevation.Peru  + theme(legend.position = "none"))+ 
    (shannon.kw.Germany.noAglais + theme(legend.position = "none") ) 
)/( 
  (temp.shannon.noAglais2 + theme(legend.position = "right") ) +
    (temp.q0 + theme(legend.position = "none"))
)+
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14, hjust = 0),
        plot.margin = margin(3, 3, 3, 3),
        legend.spacing.x = unit(0.3, 'cm'),        # horizontal spacing
        legend.spacing.y = unit(0.2, 'cm')   )      # vertical spacing


# supplement S3 PCoA genus tribe
s3.betadisp.tribe.wrap2 <- ggarrange(betadisp.genus, PCoA.genus.wrap.tribe_noleg, marginalR2.plot,   labels = c('a', 'b', 'c'),
                                     common.legend = F, ncol = 3,   nrow = 1, align = "h", widths = c(1,1.2,0.4) ) # , legend = "right"



# Venn core definition 3 area subplots A/B|C
pA <- genus.core.maxprev + theme(legend.position = "right")
pB <- venn_core_overlap + theme(legend.position = "none")

pD <- noncore.genus.abundance + 
  theme(
    axis.title.y = element_blank(),    # remove y-axis title
    axis.text.y = element_blank(),     # remove y-axis tick labels
    axis.ticks.y = element_blank()     # remove y-axis ticks
  )

left_col <- (pA / pB) + plot_layout(heights = c(1.3, 1))  # customize if needed

right_col <- (core.genus.abundance | pD) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# supplement S4 patch_venn_core_noncore_plot
s4.patch_venn_core_noncore_plot <- (left_col | right_col) +
  plot_layout(widths = c(0.7,1.2)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 14, hjust = 0),
    #plot.tag.position = c(0.02, 0.98),
    #plot.margin = margin(3, 1, 3, 1),
    legend.key.size = unit(0.8, "lines"),      # key box size legend
  )



# supplement S7 Endosymbiont
s8.endosymbiont.arrange <- ggarrange(endo.bar3, endo.bar2,    labels = c('a', 'b'),
          common.legend = T, legend = "right", ncol = 2,   nrow = 1, align = "h", widths = c(1,2) )


# supplement fig
s9.stripchart.arrange <- ggarrange( core.genus.boxplot , core_stripchart.rebuild, labels = c('a', 'b'),
                                    common.legend = T, legend = "none", ncol = 2,   nrow = 1, align = "h", widths = c(1, 2) )


# Supplement family core
fam.core.abundance2 <- ggarrange(family.core,  core.family.abundance, core.family.boxplot, labels = c('a', 'b', 'c'),
                                 common.legend = F, legend = "right", ncol = 3,   nrow = 1, align = "h" )


### ggsave Figures -------------

# figure dimensions
# 1/3 page  20 x 7 cm
# 2/3 page  20 x 20 cm
# full page 20 x 25 cm
# 1 column  10 x 12 cm

ggsave(file.path(out_dir, "Fig1.png"), plot = fig1.arrange, width = 30, height = 9, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "Fig1.pdf"), plot = fig1.arrange, width = 30, height = 9, units = "cm", dpi = 300)

ggsave(file.path(out_dir, "Fig2.png"), plot = fig2.betadisp.tribe.wrap , width = 30, height = 9, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "Fig2.pdf"), plot = fig2.betadisp.tribe.wrap , width = 30, height = 9, units = "cm", dpi = 300)

ggsave(file.path(out_dir, "Fig3.png"), plot = fig3.gen.core.abundance  , width = 30, height = 10, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "Fig3.pdf"), plot = fig3.gen.core.abundance  , width = 30, height = 10, units = "cm", dpi = 300)

ggsave(file.path(out_dir, "Fig4.png"), plot = fan.core.annotation, width = 20, height = 16, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "Fig4.pdf"), plot = fan.core.annotation, width = 20, height = 16, units = "cm", dpi = 300)

ggsave(file.path(out_dir, "Fig5.png"), plot = fig5.net.core.comp, width = 30, height = 10, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "Fig5.pdf"), plot = fig5.net.core.comp, width = 30, height = 10, units = "cm", dpi = 300)


### ggsave Supplemental Figures -----

ggsave(file.path(out_dir, "FigS1_alpha.country.png"), plot = s1.plot_alpha.country , width = 24, height = 14, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "FigS1_alpha.country.pdf"), plot = s1.plot_alpha.country , width = 24, height = 14, units = "cm", dpi = 300)

ggsave(file.path(out_dir, "FigS2_alphaelevation.png"), plot = s2.plot_alpha.elevation , width = 24, height = 20, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "FigS2_alphaelevation.pdf"), plot = s2.plot_alpha.elevation , width = 24, height = 20, units = "cm", dpi = 300)

ggsave(file.path(out_dir, "FigS3_beta_genus.png"), plot = s3.betadisp.tribe.wrap2 , width = 30, height = 13, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "FigS3_beta_genus.pdf"), plot = s3.betadisp.tribe.wrap2 , width = 30, height = 13, units = "cm", dpi = 300)

ggsave(file.path(out_dir, "FigS4_corevenn_alternative.png"), plot = s4.patch_venn_core_noncore_plot , width = 34, height = 22, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "FigS4_corevenn_alternative.pdf"), plot = s4.patch_venn_core_noncore_plot , width = 34, height = 22, units = "cm", dpi = 300)

ggsave(file.path(out_dir, "FigS5_core_square.png"), plot = square.annotation , width = 25, height = 25, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "FigS5_core_square.pdf"), plot = square.annotation , width = 25, height = 25, units = "cm", dpi = 300)

ggsave(file.path(out_dir, "FigS6_country_core.png"), plot = core_stripchart.country.stats, width = 24, height = 15, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "FigS6_country_core.pdf"), plot = core_stripchart.country.stats, width = 24, height = 15, units = "cm", dpi = 300)

ggsave(file.path(out_dir, "FigS7_bubble.country.png"), plot = ASV.bubble.core.country , width = 25, height = 30, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "FigS7_bubble.country.pdf"), plot = ASV.bubble.core.country , width = 25, height = 30, units = "cm", dpi = 300)

ggsave(file.path(out_dir, "FigS8_endo.png"), plot = s8.endosymbiont.arrange , width = 25, height = 12, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "FigS8_endo.pdf"), plot = s8.endosymbiont.arrange , width = 25, height = 12, units = "cm", dpi = 300)

ggsave(file.path(out_dir, "FigS9_corestrip.png"), plot = s9.stripchart.arrange , width = 24, height = 12, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "FigS9_corestrip.pdf"), plot = s9.stripchart.arrange , width = 24, height = 12, units = "cm", dpi = 300)

ggsave(file.path(out_dir, "FigSx_core-fam_genus.png"), plot = fam.core.abundance2, width = 35, height = 15, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "FigSx_core-fam_genus.pdf"), plot = fam.core.abundance2, width = 35, height = 15, units = "cm", dpi = 300)


# Cleanup pipeline to optimize .RData file size
rm(ASV.bubble.core.country)
rm(ASV.bubble.country)
rm(patch_venn_core_noncore_plot)
rm(Pierinae)
rm(Coliadinae)
rm(Heliconiinae)
rm(Satyrinae)
rm(Nymphalinae)
rm(Dismorphiinae)
rm(Dismorph_Nymph)

# optional: Run extended workflow for extra figures
if (file.exists("R_16S_AW_extension.R")) {
  message("Running extended workflow...")
  source("R_16S_AW_extension.R", echo = TRUE, local = TRUE) }

message("Pipeline finished at ", Sys.time())
