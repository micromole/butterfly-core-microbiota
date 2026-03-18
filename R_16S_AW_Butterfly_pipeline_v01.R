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
library(decontam); packageVersion("decontam") 
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
setwd("../16S_AW_Butterfly21_Peru_pipeline")
#setwd("../LRZ Sync+Share/data/16S_AW_Butterfly21_Peru_pipeline")
getwd()

# Main steps of the analysis pipeline:
# data.comp       # Raw project sample data file
# data.bacteria   # Cyano reads / plant reads / unresolved taxa removed
# data.fixed      # Low quality samples removed (low PCR / high cyano reads)
# data.prevfilter # Prevalence cutoff / low stringent filtering / rare taxa removed <0.01% 
# data.pruned     # positive control / spike-in taxa removed (mock community)
# data.decontam   # decontam applied on cleared controls 
# data.high       # LT2000 Low throughput samples removed 
# sample.ASV      # controls removed / only samples on ASV level 
# sample.species  # controls removed / only samples on genus / final dataset for analysis
# sample.filter   # optional: filter low abundant genera to simplify phylo tree

###  Load custom themes and functions 
source('./R_16S_AW_functions.R')

# Create output folder for plots
if(!dir.exists("plots_peru")) {
  dir.create("plots_peru")
}

sink("sessionInfo.txt")
sessionInfo()
sink()

## 00 data.comp  -------
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

# Cleanup
# rm(data.ps)
# rm(data.otu)

## Export ps data.comp
# data.comp.df <- as(sample_data(data.comp),"data.frame") 
# write.csv(data.comp.df, "data/data.comp.metadata.csv", row.names = TRUE)
# data.comp.otu <- as.data.frame(otu_table(data.comp))
# write.csv(data.comp.otu, "data/data.comp.otu.csv", row.names = TRUE)
# data.comp.tax <- as.data.frame(tax_table(data.comp))
# write.csv(data.comp.tax, "data/data.comp.tax.csv", row.names = TRUE)

# Reimport ps data.comp
re_otu <- read.csv("data/data.comp.otu.csv", row.names = 1)
otu_tab <- otu_table(as.matrix(re_otu), taxa_are_rows = T)
re_tax <- read.csv("data/data.comp.tax.csv", row.names = 1)
tax_tab <- tax_table(as.matrix(re_tax))
re_meta <- read.csv("data/data.comp.metadata.csv", row.names = 1)
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

# Cleanup
rm(re_otu)
rm(otu_tab)

## 00 data.bacteria / Cyano removal ------------------------
sample.comp <- subset_samples(data.comp, type=="sample")

# Histogram of sample read counts, Make a data frame for the read counts of each sample
data.comp.df <- data.frame(sample_data(data.comp))
data.comp.df$LibrarySize <- sample_sums(data.comp)
data.comp.df <- data.comp.df[order(data.comp.df$LibrarySize), ]
data.comp.df$Rank <- seq_len(nrow(data.comp.df))

first.view.histogram <- ggplot(data.comp.df, aes(x = LibrarySize)) +
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  ggtitle("Distribution of sample read counts") +
  xlab("Read counts") +
  scale_x_continuous(breaks = seq(0, 100000, 5000)) +
  theme_line2() 

first.view.lib <- ggplot(data.comp.df, aes(x = Rank, y = LibrarySize, color = type)) +
  geom_point() + theme_grid() 

sampling.depth <- ggplot(data.comp.df, aes(x = chip, y = LibrarySize)) +
  geom_boxplot() +
  geom_point(aes(color = type), alpha = 0.5) +
  theme_line2() +
  ggtitle("sampling depth seq chips")

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


data.bacteria.p <- tax_filter(data.bacteria, min_prevalence = 2, min_total_abundance = 50, min_sample_abundance = 10)
second.view.comp <- data.bacteria.p %>%
  #ps_filter(type == "sample") %>%
  comp_barplot(tax_level = "phylum", n_taxa = 15, sample_order = "asis",
               label = NULL ) +
  facet_wrap(vars(country), scales = "free") +  
  coord_flip() + ggtitle("data.bacteria","Cyano / plants removed") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

# Check what has been removed
bacteria.taxa <- names(sort(taxa_sums(data.bacteria),TRUE))
data.comp.taxa <- names(sort(taxa_sums(data.comp),TRUE))
cyano.taxa <- data.comp.taxa[!(data.comp.taxa %in% bacteria.taxa)]
cyano.removed = prune_taxa(cyano.taxa, data.comp)

# Show top removed taxa with genus names
# Try to identify d:Bacteria_spc_spc_spc_spc with Blast search
cyano.taxa.abundance <- data.frame(sort(taxa_sums(cyano.removed), decreasing = TRUE)[1:50])
cyano.taxa.genus  <- data.frame(
  genus = as.character(tax_table(cyano.removed)[names(sort(taxa_sums(cyano.removed), decreasing = TRUE)[1:50]), "genus"]),
  Taxa_Sum = sort(taxa_sums(cyano.removed), decreasing = TRUE)[1:50])
cyano.taxa.genus

# Make dataframe with percent removed Cyano reads per sample
cyano.removed.frame <- data.frame(
  cyano = sample_sums(cyano.removed),
  samplesum = sample_sums(data.comp),
  host_genus = sample_data(cyano.removed)$host_genus,
  PCR = sample_data(cyano.removed)$PCR,
  percent_removed = (sample_sums(cyano.removed) / sample_sums(data.comp)) * 100)
# Sort the dataframe in decreasing order
cyano.removed.frame.sorted <- cyano.removed.frame[order(cyano.removed.frame$cyano), ]

# print tail cyano removed 
tail(cyano.removed.frame.sorted, 50)


# Cyano removed per sample
cyano.removed.sample <- subset_samples(cyano.removed, type=="sample")

# Calculate percentage of cyano reads removed
percentage_removed_cyano_sample <- (sum(sample_sums(cyano.removed.sample)) / 
                                sum(sample_sums(sample.comp ))) * 100

percentage_removed_cyano_total <- (sum(sample_sums(cyano.removed)) / 
                               sum(sample_sums(data.comp))) * 100

percentage_removed_cyano_sample # 0.599%
percentage_removed_cyano_total # 0.628%

# Mean percentage of cyano removed by genus
mean_cyanoremoved_by_genus <- cyano.removed.frame.sorted %>%
  mutate(percent_removed = as.numeric(percent_removed)) %>%
  group_by(host_genus) %>%
  summarise(mean_cyano = mean(percent_removed, na.rm = TRUE)) %>%
  arrange(desc(mean_cyano))

as.data.frame(mean_cyanoremoved_by_genus)

# Show removed taxa
# simplify to higher tax rank to speed up figures (and safe space)
cyano.removed.p <- tax_glom(cyano.removed,taxrank="phylum")
cyano.removed.melt <- psmelt(cyano.removed.p)
cyano.removed.lowrank <- ggplot(cyano.removed.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 60)) + ggtitle("Cyano / plants removed") 

cyano.removed.plot <- ggplot(cyano.removed.melt, aes(x=Sample, y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity") + ggtitle("Cyano / plants removed") + theme(axis.text.x = element_blank()) + facet_wrap(~host_genus, scales="free_x", nrow=1)

# Compare data frames before after  
data.bacteria.df <- as(sample_data(data.bacteria),"data.frame")
data.bacteria.df$sumsafter <- sample_sums(data.bacteria)
data.bacteria.df$sumsbefore <- sample_sums(data.comp)
data.bacteria.df$cyanoremoved <- sample_sums(cyano.removed)
data.bacteria.df$percentcyano <-  ((data.bacteria.df$cyanoremoved) / (data.bacteria.df$sumsbefore)) * 100
data.bacteria.df$rich <- estimate_richness(data.bacteria, measures=c("Observed", "Chao1", "Shannon", "Fisher"))
str(data.bacteria.df)
data.bacteria.df$weight <- as.numeric(as.character(data.bacteria.df$weight))

cyano.percentremoved <- ggplot(data.bacteria.df  , aes(x=percentcyano , y=sumsbefore, shape=study, color=host_family, size = percentcyano)) +
    geom_point(alpha=0.7) + ylim(0, 50000) + xlim(0, 50) + xlab("percent removed") + ylab("read counts") + facet_wrap(~host_family) 

cyano.percent.boxplot <- ggplot(data.bacteria.df, aes(x = host_family, y = percentcyano, fill = host_family)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  
  geom_jitter(width = 0.2, shape = 21, size = 2,  alpha = 0.7) +  
  labs(x = "", y = "percent reads removed", title = "Cyano / plants removed") +
  scale_fill_brewer(palette = "Set2") 

cyano.boxplot <- ggplot(data.bacteria.df, aes(x = host_family, y = cyanoremoved, fill = host_family)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  
  geom_jitter(width = 0.2, shape = 21, size = 2,  alpha = 0.7) +  
  labs(x = "", y = "reads removed", title = "Cyano / plants removed") +
  scale_fill_brewer(palette = "Set2")

data.bacteria.div <- ggplot(data.bacteria.df, aes(x=rich$Shannon, y=sumsafter, size = cyanoremoved, color=host_family, shape=country))  + geom_point(alpha=0.7)  + ggtitle(
  "data.bacteria"  )  + geom_hline(yintercept = 2000, alpha = 0.5, linetype = 2)  + ylim(0, 75000) + xlim(0, 6)


pdf("plots_peru/00_data_comp_first_view_cyano.pdf", width=12, height=6)
first.view.histogram
first.view.lib
sampling.depth
first.view.comp
second.view.comp
cyano.removed.lowrank
cyano.removed.plot
cyano.percentremoved
cyano.percent.boxplot
cyano.boxplot
data.bacteria.div
dev.off()


options(max.print = 4000) 
sink("plots_peru/00_data_comp_first_view_cyano.txt")
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
"percentage_removed_cyano_total [%]"
percentage_removed_cyano_total
"percentage_removed_cyano_sample [%]"
percentage_removed_cyano_sample 
"mean_cyanoremoved_by_genus"
as.data.frame(mean_cyanoremoved_by_genus)
sink()


# Inspect sample community composition for outlier 
bacteria.comp.plot <- data.bacteria.p %>% 
  ps_filter(host_genus == "Aglais") %>% 
  comp_barplot(tax_level = "genus", n_taxa = 30,merge_other = F, sample_order = "bray",
               label = "SAMPLE" ) +
  facet_wrap(vars(host_genus), scales = "free") +  
  coord_flip() + ggtitle("data.bacteria") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))


### > Single sample check / Taxon Check----------------
# inspect unidentified ASVs in a single sample for manual BLAST and classification
single.sample.check <- prune_samples(sample_names(data.comp) == "AW.B1.G067_S67", data.comp)  
sample_sums(single.sample.check) # list total sample sum
# Show genus names and taxa sums of top 30 ASVs in a sample
data.frame(
  genus = as.character(tax_table(single.sample.check)[names(sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:30]), "genus"]),
  Taxa_Sum = sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:30])

# Check samples with a specific ASV
# e.g. ASV341,d:Bacteria,p:uncultured_bacteria 
tag_ASV <- c("ASV13") 
singleASV <- prune_taxa(tag_ASV, data.comp)
sort(data.frame(sum = sample_sums(singleASV)), decreasing = FALSE)
singleASV %>%   tax_top(n = 5, rank = "genus")

### > Taxon check 
# e.g. p:Cyanobacteria, p:Abditibacteriota p:Thermotogae p:Gemmatimonadetes p:Verrucomicrobia 
#check_taxa <- subset_taxa(data.bacteria, genus == "p:Proteobacteria_spc_spc_spc")
check_taxa <- subset_taxa(data.bacteria, genus == "o:Enterobacterales_spc_spc")
check_taxa <- subset_taxa(data.bacteria, genus == "f:Enterobacteriaceae_spc")
check_taxa <- subset_taxa(data.bacteria, phylum == "p:Deinococcus-Thermus")

# Samples with taxon
sort(data.frame(sum = sample_sums(check_taxa),
  host_genus = sample_data(check_taxa)$host_genus), decreasing = FALSE)

# ASV abundance
sort(data.frame(sum = taxa_sums(check_taxa)), decreasing = FALSE)

## 00 data.fixed / PCR depth / delete troublemaker -----------
data.bacteria

# PCR quality high med low
PCR.depth <- data.bacteria.p %>%
  #ps_filter(type == "sample") %>%
  comp_barplot(tax_level = "genus", n_taxa = 30, sample_order = "bray",
               label = "host_genus" ) +
  facet_wrap(vars(PCR), scales = "free") +  
  coord_flip() + ggtitle("data.bacteria","PCR depth") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

# Select low amplification PCR samples
data.bacteria.low = subset_samples(data.bacteria, PCR=="low")
data.frame(sort(sample_sums(data.bacteria.low), decreasing = F))

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
sample.bacteria.low = subset_samples(data.bacteria.low, type=="sample")
# Exclude Aglais to keep it in the dataset
sample.bacteria.low.select = subset_samples(sample.bacteria.low, host_genus!="Aglais")
PCR.low.names <- sample_names(sample.bacteria.low.select)

# Remove samples with more than 500 cyano reads 
cyano.removed.sample = subset_samples(cyano.removed, type=="sample")
high_cyano_reads  <- sample_names(cyano.removed.sample)[sample_sums(cyano.removed.sample) > 500]
high_cyano.reads.subset <- subset_samples(data.comp, rownames(sample_data(data.comp)) %in% high_cyano_reads)
high_cyano.reads.subset.p <- tax_filter(high_cyano.reads.subset , min_prevalence = 3, min_total_abundance = 10, min_sample_abundance = 10)

# show samples with > 500 cyano reads 
high.cyano.samples.comp <- 
  high_cyano.reads.subset.p %>%
  comp_barplot(tax_level = "genus", n_taxa = 30,merge_other = F, sample_order = "bray",
               label = "SAMPLE" ) +
  facet_wrap(vars(host_genus), scales = "free") +  
  coord_flip() + ggtitle("data.bacteria", "samples with > 500 cyano reads") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

# Remove all samples with more then 500 cyano reads
# Combine with low PCR samples (partly overlapping)
only_cyano_names <- setdiff(high_cyano_reads, PCR.low.names) # only cyano list
tag_troublemaker <- union(high_cyano_reads, PCR.low.names) # combines with low PCR

#### Remove troublemaker (high cyano and low PCR) from dataset
only_troublemaker <- subset_samples(data.bacteria, rownames(sample_data(data.bacteria)) %in% tag_troublemaker)

troublemaker.comp <- 
  only_troublemaker %>%
  comp_barplot(tax_level = "genus", n_taxa = 30,merge_other = F, sample_order = "bray",
               label = "SAMPLE" ) +
  facet_wrap(vars(host_genus), scales = "free") +  
  coord_flip() + ggtitle("only_troublemaker") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

# Coordinate to highlight troublemaker
PCR.all.low.names <- sample_names(data.bacteria.low)
all_troublemaker <- union(high_cyano_reads, PCR.all.low.names) 

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
 
data.comp.subset <- subset_samples(data.comp, !(rownames(sample_data(data.comp)) %in% tag_troublemaker))
data.bacteria.subset <- subset_samples(data.bacteria, !(rownames(sample_data(data.bacteria)) %in% tag_troublemaker))
cyano.removed.subset <- subset_samples(cyano.removed, !(rownames(sample_data(cyano.removed)) %in% tag_troublemaker))


# Many ASVs might contain zero counts (since hyper diverse samples are now removed)
data.frame(sort(taxa_sums(data.bacteria.subset), decreasing = F)[1:1000])
# Remove zero counts from dataset (singletons are removed later)
data.condensed = prune_taxa(taxa_sums(data.bacteria.subset )>0, data.bacteria.subset)
data.frame(sort(taxa_sums(data.condensed), decreasing = F)[1:1000])

### Make taxa labels nice for plots 

# removes 'd:' 'p:' 'o:' in taxa names
data.condensed <- replace_tax_prefixes(data.condensed)

### Check the names
tail(tax_table(data.condensed))
rank_names(data.condensed)


# Validate data with tax_fix, as replace_tax_prefix makes shorter names
data.condensed <- phyloseq_validate(data.condensed)
# tax_fix_interactive(data.condensed) # adjust conditions with tax interactive
data.fixed <- tax_fix(data.condensed,
                      min_length = 4,
                      unknowns = c("NA"),
                      sep = " ", anon_unique = TRUE,)

phyloseq_validate(data.fixed,
                  remove_undetected = TRUE,
                  min_tax_length = 4,
                  verbose = TRUE )

rm(data.condensed)

pdf("plots_peru/00_data_delete_troublemaker.pdf", width=12, height=6)
PCR.depth
PCR.low.samplesums 
high.cyano.samples.comp
troublemaker.comp
troublemaker.ordinate.plot
dev.off()



# Pre-Check which samples have unique and low abundant ASVs, optionally remove samples from which more than xx% would be removed
data.prevcheck <- tax_filter(data.fixed , min_prevalence = 3, min_total_abundance = 50, min_sample_abundance = 5 )

prevcheck.removed <- prune_taxa(
  setdiff(names(taxa_sums(data.fixed)), names(taxa_sums(data.prevcheck))),
  data.fixed)

prevcheck.removed.frame <- data.frame(
  prevfilter = sample_sums(prevcheck.removed),
  samplesum = sample_sums(data.fixed),
  host_genus = sample_data(prevcheck.removed)$host_genus,
  prevfilter_percent = (sample_sums(prevcheck.removed) / sample_sums(data.fixed)) * 100 )
# Sort the dataframe in decreasing order
prevcheck.removed.frame.sorted <- prevcheck.removed.frame[order(prevcheck.removed.frame$prevfilter_percent), ]
  
print(prevcheck.removed.frame.sorted)
# Pre check to identify hyperdiverse samples (from which more then xx% ASVs would be removed)


## 00 data.prevfilter / Prevalence Filtering ---------------------
data.fixed
sample.fixed <- subset_samples(data.fixed, type=="sample")

# Low stringent filtering to remove spurious phyla e.g. Acidobacteria, Tenericutes, Verrucomicrobia, candidate division WPS−1 etc.
# Show prevalence of phyla 
prev_results <- calc_prevalence(data.fixed, rank = "phylum")
prevdf1 <- subset(prev_results$taxa_table, phylum %in% get_taxa_unique(data.fixed, "phylum"))

### Prevalence Plot 
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
cutoff_min_total_abundance # equals 600 reads  

cutoff_percent <- cutoff_min_total_abundance*100/(sum(data.sums))  
cutoff_percent # 0.01% cutoff


# Insepct data: Only the first 1000 ASVs have more than 100 reads, 
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


### Low stringent filtering ASV
# For Peru set min_prevalence = 3, min_total_abundance = 60 0.001 percent (60 reads)
data.frame(sort(taxa_sums(data.fixed), decreasing = F)[1:1000])
data.filter.ASV <- tax_filter(data.fixed,  min_prevalence = 3, min_total_abundance = 60, min_sample_abundance = 10 )
data.frame(sort(taxa_sums(data.filter.ASV), decreasing = F)[1:500])

ASV_percentage_removed <- 100 * (1 - sum(sample_sums(data.filter.ASV)) / sum(sample_sums(data.fixed)))
ASV_percentage_removed # 0.8% 

# Low stringent filtering genus
data.filter.g <- aggregate_taxa(data.filter.ASV, "genus")
data.frame(sort(taxa_sums(data.filter.g))) # remove genus below  0.01 percent (600 reads)
data.filter.genus <- tax_filter(data.filter.ASV, tax_level = "genus", min_prevalence = 3, min_total_abundance = 600, min_sample_abundance = 10 )
data.filter.g <- aggregate_taxa(data.filter.genus, "genus")
data.frame(sort(taxa_sums(data.filter.g)))

genus_percentage_removed <- 100 * (1 - sum(sample_sums(data.filter.genus)) / sum(sample_sums(data.filter.ASV)))
genus_percentage_removed  # 0.39

data.filter.f <- aggregate_taxa(data.filter.genus, "family")
data.frame(sort(taxa_sums(data.filter.f))) 

data.filter.o <- aggregate_taxa(data.filter.genus, "order")
data.frame(sort(taxa_sums(data.filter.o))) # optional remove unresolved order as well

data.filter.p <- aggregate_taxa(data.filter.genus, "phylum")
data.frame(sort(taxa_sums(data.filter.p))) 

data.prevfilter <- data.filter.genus
#sort(data.frame(sum = sample_sums(data.prevfilter)), decreasing = TRUE)

data.filter.prune		= prune_taxa(taxa_sums(data.prevfilter)<500, data.prevfilter) # speed up figure
ASV_sums_df_postfilter <- data.frame(ASV_reads = taxa_sums(data.filter.prune))

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

# Prevfilter removed list 
prevfilter.removed.frame <- data.frame(
  prevfilter = sample_sums(prevfilter.removed),
  samplesum = sample_sums(data.fixed),
  host_genus = sample_data(prevfilter.removed)$host_genus,
  PCR = sample_data(prevfilter.removed)$PCR,
  prevfilter_percent = (sample_sums(prevfilter.removed) / sample_sums(data.fixed)) * 100)
# Sort the dataframe in decreasing order
prevfilter.removed.frame.sorted <- prevfilter.removed.frame[order(prevfilter.removed.frame$prevfilter_percent), ]

print(prevfilter.removed.frame.sorted)
#total_percentage_prevfilter_removed <-  (sum(prevfilter.removed.frame.sorted$prevfilter) / sum(prevfilter.removed.frame.sorted$samplesum)) * 100

# Check prevfilter Removed Aglais
prevfilter.removed.frame.sorted %>%
  filter(host_genus == "Aglais") %>%
  summarise(mean_prevfilter = mean(as.numeric(prevfilter_percent)))

# Genus list of removed reads
mean_prevfilter_by_genus <- prevfilter.removed.frame.sorted %>%
  mutate(prevfilter_percent = as.numeric(prevfilter_percent)) %>%
  group_by(host_genus) %>%
  summarise(mean_prevfilter = mean(prevfilter_percent, na.rm = TRUE)) %>%
  arrange(desc(mean_prevfilter))

as.data.frame(mean_prevfilter_by_genus)

# Check what had been removed for individual sample
single.sample.check <- prune_samples(sample_names(prevfilter.removed) == "AW.B1.G032_S33", prevfilter.removed)
sample_sums(single.sample.check) # list total sample sum
data.frame(
  genus = as.character(tax_table(single.sample.check)[names(sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:30]), "genus"]),
  Taxa_Sum = sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:30])

# Show low stringent filter removed taxa
prevfilter.removed.p <- tax_glom(prevfilter.removed,taxrank="phylum") # simplify to phylum rank to speed up figures
prevfilter.removed.melt <- psmelt(prevfilter.removed.p)

prevfilter.lowrank <- ggplot(prevfilter.removed.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 50)) + ggtitle("removed taxa low stringent filtering")

prevfilter.removed.plot <- ggplot(prevfilter.removed.melt, aes(x=Sample, y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity") + ggtitle("removed taxa low stringent filtering") + theme(axis.text.x = element_blank()) + facet_wrap(~host_genus, scales="free_x", nrow=1)


# Compare data frames before after  
data.prevfilter.df <- as(sample_data(data.prevfilter),"data.frame")
data.prevfilter.df$sumsfilter <- sample_sums(data.prevfilter)
data.prevfilter.df$sumsfixed <- sample_sums(data.fixed)
data.prevfilter.df$sumsremoved <- sample_sums(prevfilter.removed)
data.prevfilter.df$percentremoved <- ((data.prevfilter.df$sumsremoved) / (data.prevfilter.df$sumsfixed)) * 100

prevfilter.percentremoved <- ggplot(data.prevfilter.df , aes(x=percentremoved, y=sumsfilter, shape=study, color=host_family, size = percentremoved)) +
  geom_point(alpha=0.7) + facet_wrap(~host_family) +
 ylim(0, 80000) + xlim(0, 50) + xlab("percent removed") + ylab("read counts") 

# Show filtered ASVs and genus names
prevfilter.removed.sums <- data.frame(
  Genus = tax_table(prevfilter.removed)[names(sort(taxa_sums(prevfilter.removed), decreasing = TRUE)[1:50]), "genus"],
  Abundance = sort(taxa_sums(prevfilter.removed), decreasing = TRUE)[1:50] )

prevfilter.removed.sums

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
data.minorphyla
data.minorphyla.p <- tax_glom(data.minorphyla,taxrank="phylum") # compress ps object
tail(tax_table(data.minorphyla.p))
data.minorphyla.melt <- psmelt(data.minorphyla.p)


minorphyla.bar <- ggplot(data.minorphyla.melt, aes(x = subtype, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", position = "stack") + ggtitle("samples with rare phyla") 
  theme_line2() 

minorphyla <- ggplot(data.minorphyla.melt, aes(x=Sample, y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity", linewidth = 5) + ggtitle("samples with rare phyla") + theme(axis.text.x = element_blank()) + facet_wrap(~host_subfamily, scales="free_x", nrow=1)

minorphyla.lowrank <- ggplot(data.minorphyla.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 30)) + ggtitle("samples with rare phyla")

pdf("plots_peru/00_data_filter_prevalence.pdf", width=12, height=6)
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
minorphyla.bar 
minorphyla
minorphyla.lowrank
dev.off()

options(max.print=4000)
sink("plots_peru/00_data_filter_prevalence.txt")
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
"ASVs filtered percent"
ASV_percentage_removed
"genus filtered percent"
genus_percentage_removed
"percentage removed sum sample"
percentage_removed_sample
"percentage removed sum total"
percentage_removed_total
"mean_prevfilter_by_genus"
as.data.frame(mean_prevfilter_by_genus)
"top filtered ASVs"
prevfilter.removed.sums
sink()


#### > Check genus or ASV on plate frame ------------------
# check spillover of genera from positive control e.g. Listeria, Limosilactobacillus
check.genus <- subset_taxa(data.prevfilter, genus=="Listeria" ) 
 #tag_ASV <- c("ASV62") 
 #check.genus <- prune_taxa(tag_ASV, data.prevfilter)
sort(data.frame(sum = sample_sums(check.genus)), decreasing = FALSE)
reads_check.genus <- sample_sums(check.genus)
plate.frame.check.genus <- cbind(as.data.frame(sample_data(check.genus)), reads = reads_check.genus)
# plate frame removed of "low stringent filtering"
plate.frame.check.genus %>%
  separate(well, into = c("Row", "Column"), sep = "(?<=\\D)(?=\\d)", convert = TRUE) %>%
  arrange(Row, Column) %>%
  mutate(Row = fct_relevel(Row, "H", "G", "F", "E", "D", "C", "B", "A")) %>%
  ggplot(aes(x = Column, y = Row, fill = reads, label = reads)) +
  geom_tile() +
  geom_text(color = "black", size = 2) +  # Add text annotations
  labs(title = "Plate frame 'check single genus for spillover removal'",
       x = "Column",
       y = "Row") +
  facet_wrap(~ plate) +  # Create a facet for each plate
  #scale_fill_continuous(trans = "log10") + 
  scale_fill_viridis_c(direction = -1) +  
  theme_minimal()


# 00 data.pruned1  ----------------
# remove positive control taxa from dataset 
data.prevfilter
sample.prevfilter <- subset_samples(data.prevfilter, type=="sample")

# Identify Mock taxa spillover on ASV level
data.mocks <- subset_samples(data.prevfilter, type == "positive") # select mocks
sample_names(data.mocks)
sample_sums(data.mocks)

# Select top 10 ASVs from data.mocks
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


# (optional) check individual samples that only mock community gets removed
single.sample.check <- prune_samples(sample_names(spillover.removed1) == "AW.B2.P80_S234", spillover.removed1)
sample_sums(single.sample.check) # list total sample sum
# Show genus names and taxa sums of top 30 ASVs in a sample
data.frame(
  genus = as.character(tax_table(single.sample.check)[names(sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:10]), "genus"]),
  Taxa_Sum = sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:10])

# Should only remove mock community from positive controls
spillover.melt <- psmelt(spillover.removed1)
spillover.bar  <- ggplot(spillover.melt, aes(x=Sample, y=Abundance, fill=genus)) + geom_bar(#position="fill",
  stat = "identity", linewidth = 5) + ggtitle("pos control removal") + theme(axis.text.x = element_blank()) + facet_wrap(~subtype, scales="free_x", nrow=1)

spillover.lowrank <- ggplot(spillover.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=genus)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 30))

# Show spillover on plate frame
reads_spillover <- sample_sums(spillover.removed1)
# Add sample abundance information to the original data
plate.frame.spillover <- cbind(as.data.frame(sample_data(spillover.removed1)), reads = reads_spillover)
# plate frame removed of "low stringent filtering"
spillover.plate <- plate.frame.spillover %>%
  separate(well, into = c("Row", "Column"), sep = "(?<=\\D)(?=\\d)", convert = TRUE) %>%
  arrange(Row, Column) %>%
  mutate(Row = fct_relevel(Row, "H", "G", "F", "E", "D", "C", "B", "A")) %>%
  ggplot(aes(x = Column, y = Row, fill = reads, label = reads)) +
  geom_tile() +
  geom_text(color = "black", size = 2) +  # Add text annotations
  labs(title = "Plate frame 'spillover mock community'",
       x = "Column",
       y = "Row") +
  facet_wrap(~ plate) +  # Create a facet for each plate
  #scale_fill_continuous(trans = "log10") + 
  scale_fill_viridis_c(direction = -1) +  
  theme_minimal() 

# Remove spillover of mock taxa from dataset (top mock ASVs)
data.pruned1 <- prune_taxa(!(taxa_names(data.prevfilter) %in% mock.select), data.prevfilter)

data.prevfilter # pre-filter
data.pruned1 # mock ASVs removed


# Percentage of spillover removed per sample
spillover.removed1.frame <- data.frame(
  spillover = sample_sums(spillover.removed1),
  samplesum = sample_sums(data.prevfilter),
  host_genus = sample_data(spillover.removed1)$host_genus,
  spillover_percent = (sample_sums(spillover.removed1) / sample_sums(data.prevfilter)) * 100)

# Sort the dataframe in decreasing order
spillover.removed1.frame.sorted <- spillover.removed1.frame[order(spillover.removed1.frame$spillover_percent), ]
print(spillover.removed1.frame.sorted)


# Calculate percentage of spillover removed
percentage_spillover_sample <- (sum(sample_sums(spillover.removed1.sample)) / 
                                  sum(sample_sums(sample.prevfilter ))) * 100

percentage_spillover_total <- (sum(sample_sums(spillover.removed1)) / 
                                 sum(sample_sums(data.prevfilter))) * 100

percentage_spillover_sample # 0.14
percentage_spillover_total # 1.68

# Collect mock genera names for comp barplot
mock_genera <- spillover.removed1%>% tax_top(n = 20, rank = "genus")
mock_genera

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


pdf("plots_peru/00_data_pruned1_pos_control_removal.pdf", width=12, height=6)
spillover.bar
spillover.lowrank
spillover.plate
spillover.prefilter.comp
spillover.postfilter.comp
dev.off()


sink("plots_peru/00_data_pruned1_pos_control_removal.txt")
"mock.top.genus"
mock.top.genus
"spillover.removed1.frame.sorted"
print(spillover.removed1.frame.sorted)
"data.prevfilter "
data.prevfilter 
"data.pruned1"
data.pruned1 
"percentage_spillover_sample"
percentage_spillover_sample
"percentage_spillover_total"
percentage_spillover_total
sink()

## 00 data.pruned2  -----------------
## pos control removal fine tuning (remove remaining ASVs from mock controls)
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
# Remove also unresolved taxa: Actinobacteria spc, Firmicutes spc, Proteobacteria spc


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
percentage_spillover2_total # total percent removed

check.spillover2 <- subset_taxa(spillover.removed2, genus=="Pseudomonas" )
sort(data.frame(sum = sample_sums(check.spillover2 )), decreasing = FALSE)

# (optional) check individual samples if only mock gets removed, if not put ASV in keep_ASV2
single.sample.check <- prune_samples(sample_names(spillover.removed2) == "AW.B1.G029_S29", spillover.removed2) # AW.B2.P80_S234 AW.B2.P82_S237
# single.sample.check <- subset_samples(spillover.removed2, study=="Peru") # look at entire projects  
sample_sums(single.sample.check) # list total sample sum
# Show genus names and taxa sums of top 30 ASVs in a sample
data.frame(
  genus = as.character(tax_table(single.sample.check)[names(sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:10]), "genus"]),
  Taxa_Sum = sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:10])

spillover2.melt <- psmelt(spillover.removed2)
spillover.bar2 <- ggplot(spillover2.melt, aes(x=Sample, y=Abundance, fill=genus)) + geom_bar(#position="fill",
  stat = "identity") + ggtitle("pos control removal finetuning") + theme(axis.text.x = element_blank()) + facet_wrap(~subtype, scales="free_x", nrow=1)

spillover.lowrank2 <- ggplot(spillover2.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=genus)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 30))

# Show spillover2 on plate frame
reads_spillover2 <- sample_sums(spillover.removed2)
# Add sample abundance information to the original data
plate.frame.spillover2 <- cbind(as.data.frame(sample_data(spillover.removed2)), reads = reads_spillover2)
# plate frame removed of "low stringent filtering"
spillover.plate2 <- plate.frame.spillover2 %>%
  separate(well, into = c("Row", "Column"), sep = "(?<=\\D)(?=\\d)", convert = TRUE) %>%
  arrange(Row, Column) %>%
  mutate(Row = fct_relevel(Row, "H", "G", "F", "E", "D", "C", "B", "A")) %>%
  ggplot(aes(x = Column, y = Row, fill = reads, label = reads)) +
  geom_tile() +
  geom_text(color = "black", size = 2) +  # Add text annotations
  labs(title = "Plate frame 'pos control removal finetuning'",
       x = "Column",
       y = "Row") +
  facet_wrap(~ plate) +  # Create a facet for each plate
  #scale_fill_continuous(trans = "log10") + 
  scale_fill_viridis_c(direction = -1) +  
  theme_minimal() 

# Remove spillover of mock taxa from dataset
data.pruned2 <- prune_taxa(!(taxa_names(data.pruned1) %in% mock.taxa.fine), data.pruned1)

# merge both spillover.removed
spillover.removed <- merge_phyloseq(spillover.removed1, spillover.removed2)
sort(data.frame(sum = sample_sums(spillover.removed)), decreasing = FALSE)

data.prevfilter # pre-filter
data.pruned1 # pos control removal
data.pruned2 # pos control removal fine tuning

# Collect mock genera2 names for comp barplot
mock_genera2 <- spillover.removed2%>% tax_top(n = 20, rank = "genus")
mock_genera2

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

pdf("plots_peru/00_data_pruned2_pos_control_removal_finetuning.pdf", width=12, height=6)
mock.fine.select
spillover.bar2
spillover.lowrank2
spillover.plate2
spillover.prefilter.comp2
spillover.postfilter.comp2
dev.off()
# Cleanup pieline
rm(spillover.bar2)
rm(spillover.lowrank2)
rm(spillover2.melt)

## 00 data.decontam  -------------------
data.pruned2
sample.pruned2 <- subset_samples(data.pruned2, type=="sample")

#### decontam.neg 
# identify ASVs as contamination using negative controls 
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

# Top ASVs and genus names in decontam.removed.neg
data.frame(
  genus = as.character(tax_table(decontam.removed.neg)[names(sort(taxa_sums(decontam.removed.neg), decreasing = TRUE)[1:10]), "genus"]),
  Taxa_Sum = sort(taxa_sums(decontam.removed.neg), decreasing = TRUE)[1:10])
# manually put in "keep_decontam" if no contamination

#                           genus Taxa_Sum
# ASV62               Pseudomonas     3191 pattern on plateframe
# ASV241        Clostridiales spc     1464 pattern on plateframe
# ASV335        Clostridiales spc      807 pattern on plateframe
# ASV737                   Saezia      408 pattern on plateframe
# ASV452        Clostridiales spc      316 pattern on plateframe
# ASV49               Spiroplasma      275 pattern on plateframe
# ASV1991      Proteobacteria spc      207 
# ASV3180      Proteobacteria spc      202
# ASV732                 Bacillus      168
# ASV994            Acinetobacter      166


# check individual ASVs in dataset if contamination or not
tag_ASV <- c("ASV241") 
singleASV <- prune_taxa(tag_ASV, data.pruned2)
sort(data.frame(sum = sample_sums(singleASV)), decreasing = FALSE)

data.frame(sort(sample_sums(decontam.removed.neg), decreasing = T)[1:10])
sum(taxa_sums(decontam.removed.neg))

single.sample.check <- prune_samples(sample_names(decontam.removed.neg) == "AW.B1.G032_S33", decontam.removed.neg)
data.frame(
  genus = as.character(tax_table(single.sample.check)[names(sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:10]), "genus"]),
  Taxa_Sum = sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:10])


# Show to be removed taxa
decontam.removed.neg.melt <- psmelt(decontam.removed.neg)

decontam.removed.neg.lowrank <- ggplot(decontam.removed.neg.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 60)) + ggtitle("decontam neg")

decontam.removed.neg.bargraph <- ggplot(decontam.removed.neg.melt, aes(x=Sample, y=Abundance, fill=order)) + geom_bar(#position="fill",
  stat = "identity", linewidth = 5) + ggtitle("decontam neg") + theme(axis.text.x = element_blank()) + facet_wrap(~group, scales="free_x", nrow=1)

# Show decontam.neg on plate frame
reads_decontam.neg <- sample_sums(decontam.removed.neg)
# Add sample abundance information to the original data
plate.frame.decontam.neg <- cbind(as.data.frame(sample_data(decontam.removed.neg)), reads = reads_decontam.neg)
plate.decontam.neg <- plate.frame.decontam.neg %>%
  separate(well, into = c("Row", "Column"), sep = "(?<=\\D)(?=\\d)", convert = TRUE) %>%
  arrange(Row, Column) %>%
  mutate(Row = fct_relevel(Row, "H", "G", "F", "E", "D", "C", "B", "A")) %>%
  ggplot(aes(x = Column, y = Row, fill = reads, label = reads)) +
  geom_tile() +
  geom_text(color = "black", size = 2) +  # Add text annotations
  labs(title = "Plate frame 'decontam.neg mock community'",
       x = "Column",
       y = "Row") +
  facet_wrap(~ plate) +  # Create a facet for each plate
  #scale_fill_continuous(trans = "log10") + 
  scale_fill_viridis_c(direction = -1) +  
  theme_minimal() 

#### decontam pos 
# identify ASVs as contamination using positive controls
sample_data(data.pruned2)$is.pos <- sample_data(data.pruned2)$type == "positive" # 
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
                        "ASV121", "ASV194", "ASV3076", "ASV12115","ASV8279","ASV2235" ) # ASV13250 "ASV575"  
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
  genus = as.character(tax_table(single.sample.check)[names(sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:20]), "genus"]),
  Taxa_Sum = sort(taxa_sums(single.sample.check), decreasing = TRUE)[1:20])



# Show to be removed taxa
decontam.removed.pos.melt <- psmelt(decontam.removed.pos)

decontam.removed.pos.lowrank <- ggplot(decontam.removed.pos.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=family)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 60)) + ggtitle("decontam pos")

decontam.removed.pos.bargraph <- ggplot(decontam.removed.pos.melt, aes(x=Sample, y=Abundance, fill=order)) + geom_bar(#position="fill",
  stat = "identity", linewidth = 5) + ggtitle("decontam pos") + theme(axis.text.x = element_blank()) + facet_wrap(~group, scales="free_x", nrow=1)

# Show decontam.pos on plate frame
reads_decontam.pos <- sample_sums(decontam.removed.pos)
# Add sample abundance information to the original data
plate.frame.decontam.pos <- cbind(as.data.frame(sample_data(decontam.removed.pos)), reads = reads_decontam.pos)
plate.decontam.pos <- plate.frame.decontam.pos %>%
  separate(well, into = c("Row", "Column"), sep = "(?<=\\D)(?=\\d)", convert = TRUE) %>%
  arrange(Row, Column) %>%
  mutate(Row = fct_relevel(Row, "H", "G", "F", "E", "D", "C", "B", "A")) %>%
  ggplot(aes(x = Column, y = Row, fill = reads, label = reads)) +
  geom_tile() +
  geom_text(color = "black", size = 2) +  # Add text annotations
  labs(title = "Plate frame 'decontam.neg mock community'",
       x = "Column",
       y = "Row") +
  facet_wrap(~ plate) +  # Create a facet for each plate
  #scale_fill_continuous(trans = "log10") + 
  scale_fill_viridis_c(direction = -1) +  
  theme_minimal() 

decontam_me_combined <- union(decontam_me_neg, decontam_me_pos)

decontam.removed <- prune_taxa(decontam_me_combined, data.pruned2)

# Remove decontam taxa from dataset
data.decontam <- prune_taxa(!(taxa_names(data.pruned2) %in% decontam_me_combined), data.pruned2)

data.pruned2
data.decontam # decontam

# Relative amount removed by decontam
data.pruned2.rel <- transform_sample_counts(data.pruned2, function(x) x / sum(x))
decontam.removed.rel <- subset_taxa(data.pruned2.rel, taxa_names(data.pruned2.rel) %in% taxa_names(decontam.removed))
decontam.removed.rel <- tax_glom(decontam.removed.rel,taxrank="order") # compress
decontam.removed.rel.melt <- psmelt(decontam.removed.rel)

decontam.removed.rel.plot <- ggplot(decontam.removed.rel.melt, aes(x=Sample, y=Abundance, fill=order)) + geom_bar(#position="fill",
  stat = "identity", linewidth = 5) + ggtitle("percent decontam removed") + theme(axis.text.x = element_blank()) + facet_wrap(~host_genus, scales="free_x", nrow=1)


# Make dataframe with percent removed decontam per sample
decontam.removed.frame <- data.frame(
  decontam = sample_sums(decontam.removed),
  samplesum = sample_sums(data.pruned2),
  host_genus = sample_data(decontam.removed)$host_genus,
  study = sample_data(decontam.removed)$study,
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

decontam.percentremoved <- ggplot(data.decontam.df, aes(x=percentdecontam, y=sump3, shape=study, color=host_family, size = percentdecontam)) +
  geom_point(alpha=0.7) + ylim(0, 50000) + xlim(0, 75) + xlab("percent removed") + ylab("read counts") + facet_wrap(~host_family)


## Check what has been decontamed
decontam.removed.sample <- subset_samples(decontam.removed, type=="sample")

# Divide each count by the total count and multiply by 100 to get the percentage
percentage_removed_decontam_sample <- (sum(sample_sums(decontam.removed.sample)) / 
                                         sum(sample_sums(sample.pruned2))) * 100
percentage_removed_decontam_total <- (sum(sample_sums(decontam.removed)) / 
                                     sum(sample_sums(data.pruned2))) * 100
percentage_removed_decontam_sample # 0.13
percentage_removed_decontam_total # 0.15


# Mean decontam by host genus
mean_decontam_by_genus <- decontam.removed.frame.sorted %>%
  mutate(decontam_percent = as.numeric(decontam_percent)) %>%
  group_by(host_genus) %>%
  summarise(mean_decontam = mean(decontam_percent, na.rm = TRUE)) %>%
  arrange(desc(mean_decontam))

as.data.frame(mean_decontam_by_genus)


# Show filtered ASVs with genus name
top_decontam_removed <- names(sort(taxa_sums(decontam.removed), decreasing = TRUE)) # [1:50]

decontam.removed.sums <- data.frame(
  Abundance = taxa_sums(decontam.removed)[top_decontam_removed],
  Genus = tax_table(decontam.removed)[top_decontam_removed, "genus"])

decontam.removed.sums


# View the sorted dataframe
sink("plots_peru/00_data_pruned3_decontam.txt")
print(decontam.removed.frame.sorted)
"percentage_removed_decontam_sample"
percentage_removed_decontam_sample
"percentage_removed_decontam_total"
percentage_removed_decontam_total
"mean_decontam_by_genus"
as.data.frame(mean_decontam_by_genus)
"decontam.removed.sums"
decontam.removed.sums
sink()



pdf("plots_peru/00_data_pruned3_decontam.pdf", width=12, height=6)
decontam.neg
decontam.removed.neg.lowrank
decontam.removed.neg.bargraph 
decontam.pos
decontam.removed.pos.lowrank
decontam.removed.pos.bargraph
decontam.removed.rel.plot
decontam.percentremoved
dev.off()



# perform additional clean-up steps if necessary data.pruned4 etc.
# Rename final data set as data.clean
data.clean <- data.decontam

### Check data.clean
data.clean
tail(tax_table(data.clean))
table(tax_table(data.clean)[, "phylum"], exclude = NULL)

mean_samplesum_by_genus <- data.clean %>%
  sample_data() %>%                # extract metadata
  data.frame() %>%                 # coerce to dataframe
  mutate(samplesum = sample_sums(data.clean)) %>% 
  group_by(host_genus) %>%
  summarise(mean_samplesum = round(mean(samplesum, na.rm = TRUE),0)) %>%
  arrange(desc(mean_samplesum))

as.data.frame(mean_samplesum_by_genus)

## Validate sample names
data.clean <- phyloseq_validate(data.clean)


### cleanup comparison  ------------------

### filterframe to highlight all filtering steps
filterframe.df <- as(sample_data(data.clean),"data.frame")
filterframe.df$A_data.comp <- sample_sums(data.comp.subset)
filterframe.df$B_data.bacteria <- sample_sums(data.fixed) # cyano removed
filterframe.df$C_data.prevfilter <- sample_sums(data.prevfilter) # prevfilter
filterframe.df$D_data.pruned <- sample_sums(data.pruned2) # spillover removed
filterframe.df$E_data.decontam <- sample_sums(data.clean) # decontam removed

filterframe.df$Arich <- estimate_richness(data.comp.subset, measures=c("Shannon"))
filterframe.df$Brich <- estimate_richness(data.fixed, measures=c("Shannon"))
filterframe.df$Crich <- estimate_richness(data.prevfilter, measures=c("Shannon"))
filterframe.df$Drich <- estimate_richness(data.pruned2, measures=c("Shannon"))
filterframe.df$Erich <- estimate_richness(data.clean, measures=c("Shannon"))

filterframe.df$cyanoremoved <- sample_sums(cyano.removed.subset) 
filterframe.df$prevremoved <- sample_sums(prevfilter.removed)
filterframe.df$spillremoved <- sample_sums(spillover.removed)
filterframe.df$decontamremoved <- sample_sums(decontam.removed)

filterframe.df$B_cyanopercent <- ((filterframe.df$cyanoremoved) / (filterframe.df$A_data.comp)) * 100
filterframe.df$C_prevfilterpercent <- ((filterframe.df$prevremoved) / (filterframe.df$B_data.bacteria)) * 100
filterframe.df$D_prunedpercent <- ((filterframe.df$spillremoved) / (filterframe.df$C_data.prevfilter)) * 100
filterframe.df$E_decontampercent <- ((filterframe.df$decontamremoved) / (filterframe.df$D_data.pruned)) * 100

filterframe.df$filtersum <- (filterframe.df$A_data.comp - filterframe.df$E_data.decontam) /  filterframe.df$A_data.comp * 100


plot.cyano <- ggplot(filterframe.df, aes(x=Brich$Shannon, y=B_data.bacteria, size = B_cyanopercent , color=type))  + geom_point(alpha=0.7) + 
  ggtitle("Cyano"  ) +  ylim(0, 75000)

plot.prevfilter <- ggplot(filterframe.df, aes(x=Crich$Shannon, y=C_data.prevfilter, size = C_prevfilterpercent , color=type))  + geom_point(alpha=0.7) + 
  ggtitle("Prevfilter"  ) +   ylim(0, 75000)

plot.prune <- ggplot(filterframe.df, aes(x=Drich$Shannon, y=D_data.pruned, size = D_prunedpercent , color=type))  + geom_point(alpha=0.7) + 
  ggtitle("Prune"  ) +   ylim(0, 75000)

plot.decontam <- ggplot(filterframe.df, aes(x=Erich$Shannon, y=E_data.decontam, size = E_decontampercent , color=type))  + geom_point(alpha=0.7) + 
  ggtitle("decontam"  )  +   ylim(0, 75000)


# Sort by E_clean and create Rank
sortframe.df <- filterframe.df %>%
  arrange(E_data.decontam) %>%
  mutate(Rank = row_number())

# Melt the dataframe for easier plotting with ggplot2
filterframe_long <- sortframe.df %>%
  rownames_to_column(var = "Sample") %>%
  select(Sample, Rank, A_data.comp, B_data.bacteria, C_data.prevfilter, D_data.pruned, E_data.decontam, type, subtype) %>%
  pivot_longer(cols = c(A_data.comp, B_data.bacteria, C_data.prevfilter, D_data.pruned, E_data.decontam
                         ), names_to = "Stage", values_to = "ReadCount")

# Library filter comparison
lib.comparison <- ggplot(data = filterframe_long, aes(x = Rank, y = ReadCount, color = subtype, shape = Stage)) +
  geom_point(size = 2) +
  geom_line(aes(group = Sample)) +
  geom_hline(yintercept = 2500, alpha = 0.5, linetype = 2) +
  #ylim(0, 50000) + xlim(0, 600) +
  coord_cartesian(ylim = c(0, 30000), xlim = c(0, 100)) + 
  ggtitle("Sample processing comparison") +
  labs(x = "Rank", y = "Read Count") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 20))

# Select only samples
sampleframe.df <- sortframe.df %>% filter(type == 'sample')

# Create long format data for the heatmap
heatframe_long <- sampleframe.df %>%
  rownames_to_column(var = "Sample") %>%
  select(Sample, E_data.decontam, B_cyanopercent, C_prevfilterpercent, D_prunedpercent, E_decontampercent,  Rank) %>%
  pivot_longer(cols = c(B_cyanopercent, C_prevfilterpercent, D_prunedpercent, E_decontampercent), names_to = "Variable", values_to = "Value")


# Create a combined Sample label with Sample name and E_clean value
heatframe_long$SampleLabel <- paste(heatframe_long$Sample, "(", heatframe_long$F_data.decontam, ")", sep = "")

# Find best cut-off value for geom_vline (here around the last 10 samples) 
sampleframe.df$E_data.decontam # equals ca 2500 reads 

# Plot the heatmap
lib.heatmap <- ggplot(heatframe_long, aes(x = reorder(SampleLabel, Rank), y = Variable, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Sample processing comparison" ,  x = "Samples (sorted by sample sum)",
       y = "Filtering steps",
       fill = "% removed"  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  coord_flip(xlim = c(0, 80)) + geom_vline(xintercept = 10, alpha = 0.5, linetype = 2)


sorted_samples <- filterframe.df %>%
  arrange(filtersum) %>%
  select(E_data.decontam, B_cyanopercent, C_prevfilterpercent , D_prunedpercent,  E_decontampercent,  filtersum , host_genus) 
  
tail(sorted_samples, 50)

# View the sorted dataframe
sink("plots_peru/00_data_pruning_cleanup_comparison.txt")
sorted_samples
sink()


# (optional) export sample names (filtersum >10%) to be pruned later
exportframe.df <- rownames_to_column(filterframe.df, var = "Sample")
sample_names_to_prune <- exportframe.df %>%
  filter(filtersum > 10) %>%
  pull(Sample)


# Combine ordination with filtersum
# Remove genera with only one count
data.clean.count <- table(sample_data(data.clean)$host_genus)
data.clean.mod <- subset_samples(data.clean, host_genus %in% names(data.clean.count[data.clean.count > 1]))
data.clean.rel <- transform_sample_counts(data.clean.mod, function(x) x/sum(x))
data.clean.PCoA <-  ordinate(data.clean.rel, method="PCoA", "bray")

data.clean.ordination <- plot_ordination(data.clean.rel, data.clean.PCoA, color="subtype")+
  geom_point(size=6)+theme_grid() + geom_label(aes(label = sample_sums(data.clean.mod)), size = 3) +
  labs(title="sample_sums(data.clean)")

data.clean.ordination.genuswrap <- plot_ordination(data.clean.rel, data.clean.PCoA, color="host_genus")+
  geom_point(size=4)+theme_grid() + #stat_ellipse(aes(group = host_genus))  +
  facet_wrap(~host_genus) +
  geom_label(aes(label = sampleID), size = 2)
  
scoreframe.df <- as(sample_data(data.clean.mod),"data.frame")
pcoa_scores <- data.clean.PCoA$vectors
pcoa_scores_df <- as.data.frame(pcoa_scores)
pcoa_scores_subset <- pcoa_scores_df[, 1:3]
pcoa_scores_subset$sampleID <- scoreframe.df$sampleID 
pcoa_scores_subset$E_clean <- sample_sums(data.clean.mod)
pcoa_scores_subset$logsum <- log10(sample_sums(data.clean.mod))
pcoa_scores_subset$type <- scoreframe.df$type
pcoa_scores_subset$host_family <- scoreframe.df$host_family
pcoa_scores_subset$host_genus <- scoreframe.df$host_genus
pcoa_scores_subset$LT2000 <- ifelse(sample_sums(data.clean.mod) > 2000, 'HT', 'LT2000')


data.clean.ordination.samplesum <-
  ggplot(pcoa_scores_subset, aes(x=Axis.1,y=Axis.2,color=logsum,shape=LT2000))+
  geom_point(size=6, alpha=0.8)+
  theme_grid()+ 
  theme(plot.title =element_text(size=10, face='bold')) +
  scale_colour_viridis_c(direction = -1) +
  geom_label(aes(label = sampleID), size = 2) +
  facet_wrap(~host_genus)



pdf("plots_peru/00_data_pruning_cleanup_comparison.pdf", width=12, height=6)
plot.cyano 
plot.prevfilter
plot.prune 
plot.decontam
lib.comparison
lib.heatmap 
data.clean.ordination
data.clean.ordination.genuswrap 
data.clean.ordination.samplesum
dev.off()

### cleanup comparison diversity

sample.comp.rel <- transform_sample_counts(sample.comp, function(x) x/sum(x))
sample.comp.rel.PCoA <-  ordinate(sample.comp.rel, method="PCoA", "bray")
sample.comp.ordination <- plot_ordination(sample.comp.rel, sample.comp.rel.PCoA, color="host_subfamily", shape="country")+
  geom_point(size=6, alpha=0.8)+theme_grid() + stat_ellipse(aes(group = country)) + labs(title = "data.comp"  )  

sample.bacteria <- subset_samples(data.bacteria, type=="sample")
sample.bacteria.rel <- transform_sample_counts(sample.bacteria, function(x) x/sum(x))
sample.bacteria.rel.PCoA <-  ordinate(sample.bacteria.rel, method="PCoA", "bray")
sample.bacteria.ordination <- plot_ordination(sample.bacteria.rel, sample.bacteria.rel.PCoA, color="host_subfamily", shape="country")+
  geom_point(size=6, alpha=0.8)+theme_grid() + stat_ellipse(aes(group = country)) + labs(title = "data.bacteria"  )  


sample.fixed.rel <- transform_sample_counts(sample.fixed, function(x) x/sum(x))
sample.fixed.rel.PCoA <-  ordinate(sample.fixed.rel, method="PCoA", "bray")
sample.fixed.ordination <- plot_ordination(sample.fixed.rel, sample.fixed.rel.PCoA, color="host_subfamily", shape="country")+
  geom_point(size=6, alpha=0.8)+theme_grid() + stat_ellipse(aes(group = country)) + labs(title = "data.fixed"  )   


sample.prevfilter.rel <- transform_sample_counts(sample.prevfilter, function(x) x/sum(x))
sample.prevfilter.rel.PCoA <-  ordinate(sample.prevfilter.rel, method="PCoA", "bray")
sample.prevfilter.ordination <- plot_ordination(sample.prevfilter.rel, sample.prevfilter.rel.PCoA, color="host_subfamily", shape="country")+
  geom_point(size=6, alpha=0.8)+theme_grid() + stat_ellipse(aes(group = country)) + labs(title = "data.prevfilter"  )  

sample.clean <- subset_samples(data.clean, type=="sample")
sample.clean.rel <- transform_sample_counts(sample.clean, function(x) x/sum(x))
sample.clean.rel.PCoA <-  ordinate(sample.clean.rel, method="PCoA", "bray")
sample.clean.ordination <- plot_ordination(sample.clean.rel, sample.clean.rel.PCoA, color="host_subfamily", shape="country")+
  geom_point(size=6, alpha=0.8)+theme_grid() + stat_ellipse(aes(group = country)) + labs(title = "data.clean"  )  

data.high2000 = prune_samples(sample_sums(data.clean)>=2000, data.clean)
sample.high2000 <- subset_samples(data.high2000, type=="sample")
sample.high2000.rel <- transform_sample_counts(sample.high2000, function(x) x/sum(x))
sample.high2000.rel.PCoA <-  ordinate(sample.high2000.rel, method="PCoA", "bray")
sample.high2000.ordination <- plot_ordination(sample.high2000.rel, sample.high2000.rel.PCoA, color="host_subfamily", shape="country")+
  geom_point(size=6, alpha=0.8)+theme_grid() + stat_ellipse(aes(group = country)) + labs(title = "data.high LT2000"  )  

data.high5000 = prune_samples(sample_sums(data.clean)>=5000, data.clean)
sample.high5000 <- subset_samples(data.high5000, type=="sample")
sample.high5000.rel <- transform_sample_counts(sample.high5000, function(x) x/sum(x))
sample.high5000.rel.PCoA <-  ordinate(sample.high5000.rel, method="PCoA", "bray")
sample.high5000.ordination <- plot_ordination(sample.high5000.rel, sample.high5000.rel.PCoA, color="host_subfamily", shape="country")+
  geom_point(size=6, alpha=0.8)+theme_grid() + stat_ellipse(aes(group = country)) + labs(title = "data.high LT5000"  ) 

data.comp.rich  <- plot_richness(data.comp,x="host_subfamily", measures=c("Shannon","Observed")) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_point(size=4, aes(color=host_tribe))  +geom_boxplot(aes(group = host_subfamily)) + ggtitle("data.comp")


data.bacteria.rich  <- plot_richness(data.bacteria,x="host_subfamily", measures=c("Shannon","Observed")) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_point(size=4, aes(color=host_tribe))  +geom_boxplot(aes(group = host_subfamily)) + ggtitle("data.bacteria")


data.fixed.rich  <- plot_richness(data.fixed,x="host_subfamily", measures=c("Shannon","Observed")) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_point(size=4, aes(color=host_tribe))  +geom_boxplot(aes(group = host_subfamily)) + ggtitle("data.fixed")


data.prevfilter.rich  <- plot_richness(data.prevfilter,x="host_subfamily", measures=c("Shannon","Observed")) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_point(size=4, aes(color=host_tribe))  +geom_boxplot(aes(group = host_subfamily)) + ggtitle("data.prevfilter")

data.clean.rich  <- plot_richness(data.clean,x="host_subfamily", measures=c("Shannon","Observed")) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_point(size=4, aes(color=host_tribe))  +geom_boxplot(aes(group = host_subfamily)) + ggtitle("data.clean")

sample.high2000.rich  <- plot_richness(data.high2000,x="host_subfamily", measures=c("Shannon","Observed")) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_point(size=4, aes(color=host_tribe))  +geom_boxplot(aes(group = host_subfamily)) +
  geom_label(aes(label = sampleID, color=host_tribe), size = 4) + ggtitle("data.high2000")

sample.high5000.rich  <- plot_richness(data.high5000,x="host_subfamily", measures=c("Shannon","Observed")) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_point(size=4, aes(color=host_tribe))  +geom_boxplot(aes(group = host_subfamily)) +
  geom_label(aes(label = sampleID, color=host_tribe), size = 4) + ggtitle("data.high5000")



pdf("plots_peru/00_data_pruning_cleanup_comparison_div.pdf", width=12, height=6)
sample.comp.ordination
sample.fixed.ordination 
sample.prevfilter.ordination
sample.clean.ordination
sample.high2000.ordination
sample.high5000.ordination
data.comp.rich
data.bacteria.rich
data.fixed.rich
data.prevfilter.rich
data.clean.rich
sample.high2000.rich
sample.high5000.rich
dev.off()




## 00 data.high cut-off LT500 (Low Throughput) ----------------------

# Check sample abundance to find best cutoff LT500 
sort(data.frame(sum = sample_sums(data.clean)), decreasing = TRUE)

## (optional) Label samples with low throughput e.g. LT500
data.clean <- label_low_throughput(data.clean , 2000) # optimal Peru 2000-2500 
sample_names(data.clean)

# Set cut-off LT500
cutoff <- 2000

data.clean.p <- tax_glom(data.clean,taxrank="phylum") # speed up figure
clean.melt <- psmelt(data.clean.p)

sample.sum.rank <-  ggplot(clean.melt, aes(x=reorder(Sample, Abundance), y=Abundance, fill=phylum)) +
  geom_bar(stat = "identity") +  geom_hline(yintercept=cutoff,linetype = 2) + #geom_vline(xintercept=33.5, linewidth=1) +
  coord_flip( ylim = c(0, 30000),xlim = c(0, 100))  + ggtitle( "Low throughput cut-off")

sample.sum.rank2 <-  ggplot(clean.melt, aes(x=reorder(Sample, Abundance), y=Abundance, fill=phylum)) +
  geom_bar(position="fill",  stat = "identity") +
  coord_flip(xlim = c(0, 100))  + ggtitle( "Low throughput cut-off")

# Tag Top and Flop ASVs based on relative abundance data
data.clean.rel = transform_sample_counts(data.clean, function(x) x/sum(x))
#top.ASV.names = names(sort(taxa_sums(data.clean.rel), decreasing=T)[1:50]) # Select Top ASVs
#data.top.ASV <- subset_taxa(data.clean.rel, taxa_names(data.clean.rel)%in%top.ASV.names)

# sample sum vs shannon div vs Actinobacteria
data.clean.df <- as(sample_data(data.clean),"data.frame")
data.clean.df$pruned2 <- sample_sums(data.pruned2)
data.clean.df$samplesum <- sample_sums(data.clean)
data.clean.df$rich <- estimate_richness(data.clean, measures=c("Observed", "Chao1", "Shannon", "Fisher"))
data.clean.df$actino <- sample_sums(subset_taxa(data.clean.rel, phylum=="Actinobacteria" ))
data.clean.df$Sphingomonas <- sample_sums(subset_taxa(data.clean.rel, genus=="Sphingomonas" ))
data.clean.df$Brevundimonas <- sample_sums(subset_taxa(data.clean.rel, genus=="Brevundimonas" ))
data.clean.df$Pseudomonas <- sample_sums(subset_taxa(data.clean.rel, genus=="Pseudomonas" ))
data.clean.df$Bacillus <- sample_sums(subset_taxa(data.clean.rel, genus=="Bacillus" ))

# Sample sum vs Shannon diversity
plot.actino <- ggplot(data.clean.df, aes(x=rich$Shannon, y=samplesum, size = actino, color=subtype))  + geom_point(alpha=0.7) + 
  geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) +  ggtitle(paste("LT:", cutoff)) 

plot.Sphingo <- ggplot(data.clean.df, aes(x=rich$Shannon, y=samplesum, size = Sphingomonas, color=subtype))  + geom_point(alpha=0.7) + 
  geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) +  ggtitle(paste("LT:", cutoff)) 

plot.Brev <- ggplot(data.clean.df, aes(x=rich$Shannon, y=samplesum, size = Brevundimonas, color=subtype))  + geom_point(alpha=0.7) + 
  geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) +  ggtitle(paste("LT:", cutoff)) 

plot.Bacillus <- ggplot(data.clean.df, aes(x=rich$Shannon, y=samplesum, size = Bacillus, color=subtype))  + geom_point(alpha=0.7) + 
  geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) +  ggtitle(paste("LT:", cutoff)) 

sum.shannon.genus <- ggplot(data.clean.df, aes(x=rich$Shannon, y=samplesum, size = Sphingomonas, color=host_genus, shape=study))  +
  geom_point(alpha=0.7)  +
  geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) +
  ggtitle(paste("LT:", cutoff)) + ylim(0, 10000)  +
  facet_wrap(~host_genus)


# Remove LT500 samples with defined cut-off e.g. 2000
data.high = prune_samples(sample_sums(data.clean)>=cutoff, data.clean)
data.low = prune_samples(sample_sums(data.clean)<cutoff, data.clean)
data.clean # cleaned dataset
data.high # high throughput dataset
sample_names(data.clean)
sample_names(data.high) # check removal of LT samples
data.high.rel = transform_sample_counts(data.high, function(x) x/sum(x))

data.high.df <- as(sample_data(data.high),"data.frame")
data.high.df$samplesum <- sample_sums(data.high)
data.high.df$rich <- estimate_richness(data.high, measures=c("Observed", "Chao1", "Shannon", "Fisher"))
data.high.df$Sphingomonas <- sample_sums(subset_taxa(data.high.rel, genus=="Sphingomonas" ))
data.high.df$Bacillus <- sample_sums(subset_taxa(data.high.rel, genus=="Bacillus" ))

sumflop <- ggplot(data.clean.df, aes(x=rich$Shannon, y=samplesum, size = Sphingomonas, color=subtype))  + geom_point(alpha=0.7)  + ggtitle(
  "data.clean"  ) + ylim(0, 75000) + xlim(0, 5) + geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2)
sumflophigh <- ggplot(data.high.df, aes(x=rich$Shannon, y=samplesum, size = Sphingomonas, color=subtype))  + geom_point(alpha=0.7) + ggtitle(
  "data.high (LT removal)"  ) + ylim(0, 75000) + xlim(0, 5) + geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2)

Shann_rich <- ggplot(data.clean.df, aes(x=rich$Shannon, y=rich$Observed))  + geom_point(aes(size = samplesum,  color=subtype), alpha = 0.7)  + ggtitle(
  "data.clean"  ) + xlim(0, 5) + stat_smooth(method="lm", color="black", se=T) +  stat_regline_equation(label.y = 300.0, aes(label = after_stat(rr.label))) 

Shann_rich_high <- ggplot(data.high.df, aes(x=rich$Shannon, y=rich$Observed))  + geom_point(aes(size = samplesum,  color=subtype), alpha = 0.7)  + ggtitle(
  "data.high (LT removal)"  ) +xlim(0, 5) + stat_smooth(method="lm", color="black", se=T) +  stat_regline_equation(label.y = 300.0, aes(label = after_stat(rr.label))) 

data.high.remainingsamples <- ggplot(data.high.df, aes(x=rich$Shannon, y=samplesum, color=host_genus))  +
  geom_point(alpha=0.7)  + ylim(0, 10000) + xlim(0, 5) +
  geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) +
  facet_wrap(~host_genus) + ggtitle("data.high") +
  geom_label(aes(label = sampleID), size = 4)

## Ordination data.low
data.low.rel <- transform_sample_counts(data.low, function(x) x/sum(x))
data.low.PCoA <-  ordinate(data.low.rel, method="PCoA", "bray")


data.low.PCoA.plot <- plot_ordination(data.low.rel, data.low.PCoA, color="host_subfamily", shape="study")+
  geom_point(size=6)+theme_grid()+ ggtitle(paste("data.low LT:", cutoff))# +geom_label(label=sample_names(data.low.rel))

# Create a dataframe with sample sums and metadata
data.low.df <- data.frame(
  SampleSums = sample_sums(data.low),
  host_genus = sample_data(data.low)$host_genus,
  study = sample_data(data.low)$study,
  subtype = sample_data(data.low)$subtype
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

data.high.LT <- label_low_throughput(data.high , 5000)
sample_names(data.high.LT)

data.high.LT.genus.comp <- data.high.LT %>%
  ps_filter(type=="sample") %>%
  comp_barplot(tax_level = "genus", n_taxa = 30, merge_other = T, label = "SAMPLE", sample_order = "bray") +
  facet_wrap(vars(country), scales = "free") +
  coord_flip() + ggtitle( "data.high show LT5000")


LT.samples = prune_samples(sample_sums(data.high)<10000, data.high)
LT.samples # only low samples
sample_data(LT.samples)$samplesums <- sample_sums(LT.samples)

data.high.low.comp <- LT.samples %>%
  #ps_filter(country=="Germany") %>%
  comp_barplot(tax_level = "genus", n_taxa = 30, merge_other = F, label = "samplesums", sample_order = "bray") +
  facet_wrap(vars(host_genus), scales = "free") +  
  coord_flip() + ggtitle( "data.high.LT2500 to LT10000")

# Create a dataframe with sample sums and metadata
LT.samples.df <- data.frame(
  SampleSums = sample_sums(LT.samples),
  host_genus = sample_data(LT.samples)$host_genus,
  study = sample_data(LT.samples)$study,
  subtype = sample_data(LT.samples)$subtype
)

# Sort the dataframe by SampleSums in decreasing order
sorted_LT.samples <- LT.samples.df[order(-LT.samples.df$SampleSums), ]


# Show sample sums to identify best cut-off
sink("plots_peru/00_data_sample_cutoff_LT500.txt")
"Show sample sums to identify best cutoff"
"data.clean"
data.frame(sort(sample_sums(data.clean), decreasing = TRUE))
print(sorted_LT.samples)
"cutoff"
cutoff
print(sorted_data.low)
sink()



pdf("plots_peru/00_data_sample_cutoff_LT500.pdf", width=12, height=6)
sample.sum.rank
sample.sum.rank2
plot.actino
plot.Sphingo 
plot.Brev
plot.Bacillus
sum.shannon.genus
sumflop
sumflophigh
Shann_rich 
Shann_rich_high 
data.high.remainingsamples
data.low.PCoA.plot
data.low.comp 
data.high.LT.genus.comp
data.high.low.comp
dev.off()


# (optional) rarefaction as alternative for normalizing sampling depth differences
#set.seed(111) # keep result reproductive
#ps.rarefied = rarefy_even_depth(data.high, rngseed=1, sample.size=2000, replace=F)

# Use data.high for follow up analysis as data.ASV
data.ASV <- data.high

## data.species (tax_glom genus)
## Multiple ASVs might represent the same species, here they are collated (takes time)
data.species <- tax_glom(data.ASV,taxrank="genus")
taxa_names(data.species) <- tax_table(data.species)[,"genus"]


data.ASV # ASV level
sample.ASV <- subset_samples(data.ASV, type=="sample")
data.species # species level
sample.species <- subset_samples(data.species, type=="sample")

mean_data.species_samplesum <- data.species %>%
  sample_data() %>%                # extract metadata
  data.frame() %>%                 # coerce to dataframe
  mutate(samplesum = sample_sums(data.species)) %>% 
  group_by(host_genus) %>%
  summarise(mean_samplesum = round(mean(samplesum, na.rm = TRUE),0)) %>%
  arrange(desc(mean_samplesum))

as.data.frame(mean_data.species_samplesum)

#### > Check which taxa have been removed with samples
# Identify removed samples
all.removed.names <- setdiff(sample_names(data.comp), sample_names(data.ASV))
all.removed <- prune_samples(all.removed.names, data.comp)
samples.removed <- subset_samples(all.removed, type=="sample") # sample
samples.removed.species <- tax_glom(samples.removed,taxrank="genus") # genus level
taxa_names(samples.removed.species) <- tax_table(samples.removed.species)[,"genus"]
samples.removed.high = prune_samples(sample_sums(samples.removed.species)>=500, samples.removed.species) #remove LT1000 samples
samples.removed.rel <- transform_sample_counts(samples.removed.high, function(x) x/sum(x))
samples.removed.rel  %>% tax_top(n = 20, rank = "genus")

# Inspect genus of interest
removed_taxa <- subset_taxa(samples.removed.rel, genus == "g:Pseudomonas")
removed_taxa <- subset_taxa(samples.removed.rel, genus == "g:Wolbachia")
removed_taxa <- subset_taxa(samples.removed.rel, genus == "g:Entomomonas")
removed_taxa <- subset_taxa(samples.removed.rel, genus == "g:Apibacter")
removed_taxa <- subset_taxa(samples.removed.rel, genus == "g:Bacillus")

# Show rel. ab. of  taxa in removed samples
sort(data.frame(sum = sample_sums(removed_taxa)*100,
                host_genus = sample_data(removed_taxa)$host_genus), decreasing = FALSE)





## 00 (optional) sample.filter genus -----------------------
# optional low abundance filtering on genus level
data.frame(sort(taxa_sums(sample.ASV), decreasing = F)[1:100])

# Choose cutoff 
rel_ab_cutoff <- 0.035   # 0.035%
cutoff_min_total_abundance <- ((rel_ab_cutoff/100 ) * sum(sample_sums(sample.ASV))) # 0.04 percent cut-off
cutoff_min_total_abundance # equals ~2000 reads  
cutoff_min_percentage <- (cutoff_min_total_abundance*100)/sum(sample_sums(sample.ASV))
cutoff_min_percentage # equals 0.035 percentage

# Calculate min_prevalence e.g. 5% of samples
cutoff_min_prevalence <- (5*nsamples(sample.ASV))/100
cutoff_min_prevalence # equals about 8 samples --> set this number as min_prevalence
# 'min_prevalence = 17' calculated in % of samples
percent_of_samples <- (8*100)/nsamples(sample.ASV)
percent_of_samples


# Choose good cutoff based on prevalence / abundance plot
prev_results3 <- calc_prevalence(sample.species, rank = "phylum")
prevdf3 <- subset(prev_results3$taxa_table, phylum %in% get_taxa_unique(sample.species, "phylum"))

sample.species.prevalence <- ggplot(prevdf3, aes(TotalAbundance, Prevalence / nsamples(sample.ASV),color=phylum)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + # prev 0.05%
  geom_vline(xintercept = 2000.0, alpha = 0.5, linetype = 2) + # abundance 2000 reads
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none") + ggtitle(
    "sample.species min_total_abundance = 2000"  ) 

# Filter on genus level >2000 reads = 80 taxa remaining
sample.g <- aggregate_taxa(sample.ASV, "genus")
data.frame(sort(taxa_sums(sample.g), decreasing = T)[1:149])
sample.filter <- tax_filter(sample.ASV, tax_level = "genus", min_prevalence = 8, min_total_abundance = 2000, min_sample_abundance = 1000 ) 
sample.filter.g <- aggregate_taxa(sample.filter, "genus")
data.frame(sort(taxa_sums(sample.filter.g)))

# percent of total reads remaining 
percent_retained <- (sum(otu_table(sample.filter)) / sum(otu_table(sample.ASV))) * 100
percent_retained # 98.14692 of total reads remaining

#sample.ASV.filter <- tax_filter(sample.ASV, tax_level = "genus", min_prevalence = 9, min_total_abundance = 2000, min_sample_abundance = 1000 ) # 1800
#sample.filter <- tax_glom(sample.ASV.filter,taxrank="genus")
#taxa_names(sample.filter) <- tax_table(sample.filter)[,"genus"]


# Check what has been removed
genus.removed.taxa <- setdiff(taxa_names(sample.ASV), taxa_names(sample.filter))
genus.removed <- prune_taxa(genus.removed.taxa, sample.ASV)

# percentage of reads removed
percentage_species_removed <- sum(sample_sums(genus.removed)) / sum(sample_sums(sample.ASV)) * 100
percentage_species_removed # 1.85% removed

genus.removed.f <- tax_glom(genus.removed,taxrank="order")
sample.removed.melt <- psmelt(genus.removed.f)
sample.filter.lowrank <- ggplot(sample.removed.melt, aes(x=reorder(Sample, -Abundance), y=Abundance, fill=order)) + geom_bar(#position="fill",
  stat = "identity") +
  coord_flip(xlim = c(0, 60)) + ggtitle("sample.filter")

sample.filter.bar <- ggplot(sample.removed.melt, aes(x=Sample, y=Abundance, fill=phylum)) + geom_bar(#position="fill",
  stat = "identity", linewidth = 5) + ggtitle("sample.filter") + theme(axis.text.x = element_blank()) + facet_wrap(~host_genus, scales="free_x", nrow=1)

# Compare data frames before after  
sample.filter.df <- as(sample_data(sample.filter),"data.frame")
sample.filter.df$sum_before <- sample_sums(sample.ASV)
sample.filter.df$sum_after <- sample_sums(sample.filter)
sample.filter.df$sumremoved <- sample_sums(genus.removed)
sample.filter.df$rich <- estimate_richness(sample.filter, measures=c("Shannon"))
sample.filter.df$percentfilter <- ((sample.filter.df$sumremoved) / (sample.filter.df$sum_before )) * 100


sample.filter.percent <- ggplot(sample.filter.df, aes(x=percentfilter , y=sum_after, shape=country, color=host_family, size = percentfilter)) +
  geom_point(alpha=0.7) + ylim(0, 50000) + xlim(0, 75) + xlab("percent removed filter") + ylab("read counts") + facet_wrap(~host_family)


# some samples might dropped below LT2000
sort(data.frame(sum = sample_sums(sample.filter)), decreasing = TRUE)

pdf("plots_peru/00_data_sample_filter_(optional).pdf", width=12, height=6)
sample.species.prevalence
sample.filter.lowrank 
sample.filter.bar 
sample.filter.percent
dev.off()


###  cleanup pipeline / export sample data csv  -------------------
# Check filesize of all ps objects
# sapply(ls(), function(x) object.size(get(x))) %>% sort(decreasing = TRUE)
# Print size of ps objects
print(object.size(data.comp), units = "auto")
print(object.size(data.comp.subset), units = "auto")
print(object.size(sample.comp), units = "auto")
print(object.size(sample.comp.rel), units = "auto")
print(object.size(sample.bacteria.rel), units = "auto")
print(object.size(troublemaker.comp), units = "auto")
print(object.size(data.bacteria), units = "auto")
print(object.size(data.bacteria.subset), units = "auto")
print(object.size(sample.bacteria), units = "auto")
print(object.size(sample.fixed.rel), units = "auto")


# Clean pipeline delete large ps objects not needed anymore
rm(data.comp)
rm(data.comp.subset)
rm(sample.comp)
rm(sample.comp.rel)
rm(sample.bacteria.rel)
rm(troublemaker.comp)
rm(data.bacteria)
rm(data.bacteria.subset)
rm(sample.bacteria)
rm(sample.fixed.rel)
gc()

# export sample metadata as csv file
#sample_metadata <- as(sample_data(sample.ASV),"data.frame")
# ad column LT removed to filterframe.df
filterframe.df$LT_removed <- ifelse(filterframe.df$E_data.decontam < cutoff, "LT_remove", "HT_keep")
# select columns for metadata
sample_metadata <- filterframe.df[, c("chip", "sampleID", "collector", "host_order", "host_family", "host_subfamily", "host_tribe", "host_genus", "host_species", "BOLD_ID", "sex", "location", "sublocation", "country",
                                      "date", "year", "kw", "elevation", "temp_week", "storage", "type", "PCR",
                                      "A_data.comp", "B_data.bacteria", "C_data.prevfilter", "D_data.pruned", "E_data.decontam", "LT_removed", "cyanoremoved", "prevremoved", "spillremoved", "decontamremoved", "B_cyanopercent", "C_prevfilterpercent", "D_prunedpercent", "E_decontampercent", "filtersum", "Accession", "BioProject")]
#subset_metadata$TotalReads <- sample_sums(sample.ASV)
class(sample_metadata)
write.csv(sample_metadata, "plots_peru/Suppl_table_sample_filter_metadata.csv", row.names = T)



## 01 sample comp overview  -------------------

# ASV level
data.ASV # LT samples removed
sample.ASV %>% tax_top(n = 20, rank = "genus")
sample.filter # filtered from low abundant genera

# genus level
data.species # LT samples removed
sample.species # controls removed


#Transform to relative data
# ASV level data sets
data.ASV.rel <- transform_sample_counts(data.ASV, function(x) x/sum(x))
sample.ASV.rel <- transform_sample_counts(sample.ASV, function(x) x/sum(x))
sample.filter.rel <- transform_sample_counts(sample.filter, function(x) x/sum(x))

# Genus level data sets
data.species.rel <- transform_sample_counts(data.species, function(x) x/sum(x))
sample.species.rel <- transform_sample_counts(sample.species, function(x) x/sum(x))

# Aggregate on family level for core analysis
sample.family <- aggregate_taxa(sample.species, "family")
sample.family.rel <- transform_sample_counts(sample.family, function(x) x/sum(x))

sample.order <- aggregate_taxa(sample.species, "order")
sample.order.rel <- transform_sample_counts(sample.order, function(x) x/sum(x))

# Subset countries
Peru.species <- subset_samples(sample.species, country=="Peru" )
Peru.species.rel <- subset_samples(sample.species.rel, country=="Peru" )

Germany.species <- subset_samples(sample.species, country=="Germany" )
Germany.species.rel <- subset_samples(sample.species.rel, country=="Germany" )

# subset subfamilies
Satyrinae.rel <- subset_samples(sample.species.rel, host_subfamily=="Satyrinae")
Dismorphiinae.rel <- subset_samples(sample.species.rel, host_subfamily=="Dismorphiinae")
Pierinae.rel <- subset_samples(sample.species.rel, host_subfamily=="Pierinae")
Heliconiinae.rel <- subset_samples(sample.species.rel, host_subfamily=="Heliconiinae")
Coliadinae.rel <- subset_samples(sample.species.rel, host_subfamily=="Coliadinae")
Nymphalinae.rel <- subset_samples(sample.species.rel, host_subfamily=="Nymphalinae")

## sample.ASV comp overview
sink("plots_peru/01_All_sample_comp_overview_median_sum.txt")
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
table(sample_data(sample.species)$year)
table(sample_data(sample.species)$sublocation)
"mean_samplesum data.clean before LT500"
as.data.frame(mean_samplesum_by_genus)
"mean_data.species_samplesum after LT500"
as.data.frame(mean_data.species_samplesum)
sink()


# Taxonomy overview samples
sample.species.rel.melt <- psmelt(sample.species.rel)

# taxa summary bacterial order
all_sample_order_stats <- sample.species.rel.melt %>%
  group_by(order) %>%
  summarise(
    mean_abundance = mean(tapply(Abundance, Sample, sum)),
    sd_abundance   = sd(tapply(Abundance, Sample, sum)),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance)) %>% 
  slice_head(n = 10)             

# taxa summary bacterial families
all_sample_family_stats <- sample.species.rel.melt %>%
  group_by(family) %>%
  summarise(
    mean_abundance = mean(tapply(Abundance, Sample, sum)),
    sd_abundance   = sd(tapply(Abundance, Sample, sum)),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 10)

sink("plots_peru/01_All_sample_comp_overview_taxa.txt")
"phyla sample.species"
table(tax_table(sample.species)[, "phylum"], exclude = NULL)
table(tax_table(sample.species)[, "order"], exclude = NULL)
"sample.species" 
data.frame(sort(taxa_sums(sample.species), decreasing = T)[1:50])
"sample.species.rel"
data.frame(round(sort(taxa_sums(sample.species.rel)*100 / nsamples(sample.species.rel), decreasing = TRUE)[1:20],2 ))
"top 10 order"
as.data.frame(all_sample_order_stats)
"top 10 family"
as.data.frame(all_sample_family_stats)
sink()

rm(sample.species.rel.melt)






#### custom color palette -------

# count group numbers for host palette
sample.species.melt <- psmelt(sample.species)
speciesCount = length(unique(sample.species.melt$host_species))

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

subfamily_names = sort(unique(sample.species.melt$host_subfamily))
subfamily_count = length(subfamily_names)
subfamily_colorfix <- setNames(brewer.pal(n = 8, "Accent")[1:subfamily_count], subfamily_ordered_color) # max 8
#subfamily_colorfix <- setNames(colorRampPalette(brewer.pal(8, "Accent"))(subfamily_count), subfamily_ordered_color) # more than 8

tribe_names = sort(unique(sample.species.melt$host_tribe))
tribe_count = length(tribe_names)
tribe_color = brewer.pal(n = 12, "Set3")[1:tribe_count]
tribe_color_sat = saturation(tribe_color, delta(+0.2))
tribe_colorfix <- setNames(tribe_color_sat, tribe_names)

genus_names = sort(unique(sample.species.melt$host_genus))
genus_count = length(genus_names)
genus_color <- colorRampPalette(brewer.pal(n = 12, "Set3"))(genus_count)
genus_color_sat = saturation(genus_color, delta(+0.2))
genus_colorfix <- setNames(genus_color_sat, genus_names)

species_names = sort(unique(sample.species.melt$host_species))
species_count = length(species_names)
species_colorfix <- setNames(colorRampPalette(brewer.pal(n = 12, "Set3"))(species_count), species_names)


country_names = sort(unique(sample.species.melt$country))
country_count = length(country_names)
country_colorfix <- setNames(brewer.pal(n = 8, "Dark2")[1:country_count], country_names)


# Define specific fixed colors for microbial taxa
taxaPalette = colorRampPalette(brewer.pal(12, "Paired")) # main palette for core taxa
noncorePalette = colorRampPalette(brewer.pal(8, "Set1")) # secondary palette for non core taxa
orderPalette = colorRampPalette(brewer.pal(11, "Paired")) # main palette for order


## 01 Cum Rel Abundance  -------------------

### Rel abundance all samples 
top_taxa_desc <- sort(rowMeans(otu_table(sample.species.rel)), decreasing = TRUE)[1:20]

relabundance <- ggplot(data.frame(genus = names(top_taxa_desc), mean_abundance = top_taxa_desc),
       aes(x = reorder(genus, -mean_abundance), y = mean_abundance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "rel. ab [%]") +
  theme_line2()


#### Top order + fam >1% mirrored  -----------
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
percent_top_fam # 86.72571 % top 8 order >1% fam

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
  mutate(Abundance_mirrored = ifelse(country == "Peru", Abundance, -Abundance))

# Tornado top order / family min 1% rel abundance
order.fam.country.mirrored <- ggplot(sample.top.tornado , aes(family, Abundance_mirrored, fill= order)) +
  facet_grid(order~country, space = "free", scales = "free",switch = "y")+ # switch = "y" remove order names
  theme_grid()+
  geom_bar(stat="identity")+
  #theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
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
    mean_total_abundance = mean(total_abundance),
    sd_total_abundance = sd(total_abundance)  )

top_order_total_stats

top_order_family_stats <- top_order_genera.melt %>%
  group_by(order, Sample) %>%                      # keep Order grouping
  summarise(order_abundance = sum(Abundance), .groups = "drop") %>%
  group_by(order) %>%                              # then summarise per order
  summarise(
    mean_abundance = mean(order_abundance),
    sd_abundance = sd(order_abundance)
  )

top_order_family_stats

pdf("plots_peru/01_Cum_Rel_Abundance.pdf", width=8, height=6)
relabundance
order.fam.country.mirrored 
order.fam.country.mirrored.genus
dev.off()



### 01 Group merged comp barplot -------------

# Comp barplot (merged) phylum
Top.comp.merged.phylum <- sample.species%>%
  ps_select(host_genus,host_tribe,host_family,host_subfamily,type,country) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(type == "sample") %>%
  phyloseq::merge_samples(group = "host_subfamily") %>%
  comp_barplot(tax_level = "phylum", n_taxa = 5, merge_other = T, sample_order = subfamily_ordered , bar_width = 0.9) +
  labs(x = NULL, y = NULL) + ggtitle(""  ) +
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) #+ coord_flip()


# Comp barplot (merged) order
Top.comp.merged.order <- sample.species%>%
  ps_select(host_genus,host_tribe,host_family,host_subfamily,type,country) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(type == "sample") %>%
  phyloseq::merge_samples(group = "host_subfamily") %>%
  comp_barplot(tax_level = "order", n_taxa = 10, merge_other = T, sample_order = subfamily_ordered , bar_width = 0.9) +
  labs(x = NULL, y = "rel. ab. [%]") + ggtitle(""  ) +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) 


# Comp barplot (merged) family
Top.comp.merged.family <- sample.species%>%
  ps_select(host_genus,host_tribe,host_family,host_subfamily,type,country) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(type == "sample") %>%
  phyloseq::merge_samples(group = "host_subfamily") %>%
  comp_barplot(tax_level = "family", n_taxa = 14, merge_other = T, sample_order = subfamily_ordered, bar_width = 0.9) +
  labs(x = NULL, y = NULL) + ggtitle("Comp barplot merged family"  ) +
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) 

# Comp barplot (merged) genus
Top.comp.merged.genus <- sample.species%>%
  ps_select(host_genus,host_tribe,host_family,host_subfamily,type,country) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(type == "sample") %>%
  phyloseq::merge_samples(group = "host_subfamily") %>%
  comp_barplot(tax_level = "genus", n_taxa = 12, merge_other = T, sample_order = subfamily_ordered, bar_width = 0.9) +
  labs(x = NULL, y = NULL) + ggtitle("Comp barplot merged genus"  ) +
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) )


#### Group merged Top genus ------
sample.species.merged = merge_samples(sample.species, "host_genus") # merge by host_genus or host_subfamily
sample.species.merged.rel <- transform_sample_counts(sample.species.merged, function(x) x/sum(x))
# Select top 10 genera
Top.genus <- names(sort(taxa_sums(sample.species.merged), decreasing=T)[1:12]) # top 10 taxa
# sample.species.top = prune_taxa(taxa_sums(sample.species.merged.rel)>0.20, sample.species.merged.rel)
sample.species.top <- subset_taxa(sample.species.merged.rel, taxa_names(sample.species.merged.rel)%in%Top.genus)
sample.species.top.melt <- psmelt(sample.species.top)

# Calculate total abundance per sample (sample sorting)
sample.species.top.melt <- sample.species.top.melt %>%
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))

# Calculate genus abundance (genus color order)
genus_abundance <- sample.species.top.melt  %>%
  group_by(genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance))  # Sort by abundance (most abundant first)

# Color palette and color order
topgenusCount = length(unique(sample.species.top.melt$genus))
topgenus_palette <- setNames(taxaPalette(topgenusCount), genus_abundance$genus)  # Assign colors in the abundance genus
sample.species.top.melt$genus <- factor(sample.species.top.melt$genus, levels = genus_abundance$genus) # genus names by order

sample.species.top.bar <- ggplot(sample.species.top.melt,aes(x = fct_reorder(Sample, -total_abundance), y=Abundance, fill = genus))+
  geom_bar(#position="fill",
    colour="black", linewidth=0.3, # remove black lines optional
    stat="identity") + 
  scale_fill_manual(values = topgenus_palette)+
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  theme(legend.title=element_text(size=11), legend.text=element_text(face="italic")) +# x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  labs(x="",y="rel. ab. [%]", fill = "top taxa")


#### Top 20 genera dotplot
Top.genus <- names(sort(taxa_sums(sample.species.rel ), decreasing=T)[1:20]) # top 20 taxa
# sample.species.top = prune_taxa(taxa_sums(sample.species.merged.rel)>0.20, sample.species.merged.rel)
sample.species.top <- subset_taxa(sample.species.rel, taxa_names(sample.species.rel)%in%Top.genus)
sample.species.top.melt <- psmelt(sample.species.top)

### Dotplot (Top 20 genera) 
dotplot.top.genus <- ggplot(sample.species.top.melt,aes(x=host_subfamily, y=Abundance)) +
  geom_point(position = position_jitter(w = 0.1, h = 0),aes(color=host_subfamily),  size=3, alpha=1) +theme_line2() +
  theme(strip.text = element_text(colour = "black", size=10, face="italic"), legend.text = element_text(face="italic")) +
  stat_smooth(method="lm", color="black", linewidth=0.5,se=T) +
  facet_wrap(~genus) + 
  #scale_x_continuous(breaks=seq(0,6,1)) + scale_y_continuous(breaks=seq(0,1.0,0.25), labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +
  stat_regline_equation(label.y = 0.9, aes(label = after_stat(rr.label)),size =3) + stat_cor(label.y = 0.75, size =3, p.accuracy = 0.001) +
  scale_color_manual(values = subfamily_colorfix) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),axis.title.x = element_text(family = "sans", size = 15)) +
  labs(x="",y="Relative Abundance [%]")


### Top genera merge by group (host subfamily)
Top.genus <- names(sort(taxa_sums(sample.species.merged), decreasing=T)[1:12]) # top 10 taxa
subfamily.merged = merge_samples(sample.species, "host_subfamily")
subfamily.merged.rel <- transform_sample_counts(subfamily.merged, function(x) x/sum(x))
abundant_group <- subset_taxa(subfamily.merged.rel, taxa_names(subfamily.merged.rel)%in%Top.genus)
group.melt <- psmelt(abundant_group)

# Calculate total abundance per sample (sample sorting)
group.melt <- group.melt %>%
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))

# Color palette and color order
group.melt$genus <- factor(group.melt$genus, levels = genus_abundance$genus) # genus names by order

sample.species.top.bar.group <- ggplot(group.melt,aes(x=fct_reorder(Sample, -total_abundance), y=Abundance, fill = genus))+
  geom_bar( colour="black",
     stat="identity", linewidth=0.3)+
  scale_fill_manual(values = topgenus_palette)+
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1, face="italic") ) + # x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  theme(legend.title=element_text(size=11), legend.text=element_text(size=10, face="italic")) + # legend text
  labs(x="",y="Relative Abundance [%]", fill = "top taxa") # scale_x_continuous(breaks=seq(0,6,1))


#### Group merged Top order -----

# Top order merge by group
sample.order.merged = merge_samples(sample.order, "host_subfamily")
sample.order.merged.rel <- transform_sample_counts(sample.order.merged, function(x) x/sum(x))
data.frame(sort(taxa_sums(sample.order.merged), decreasing = F))
Top.order <- names(sort(taxa_sums(sample.order.merged), decreasing=T)[1:8])
sample.order.top <- subset_taxa(sample.order.merged.rel, taxa_names(sample.order.merged.rel)%in%Top.order)
sample.order.top.melt <- psmelt(sample.order.top)
toporderCount = length(unique(sample.order.top.melt$order))

remaining_percent <- sum(sample_sums(sample.order.top))*100/nsamples(sample.order.top)

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


sample.order.top.bar <- ggplot(sample.order.top.melt,
       aes(x = Sample, y = Abundance, fill = order)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.1) + 
  scale_fill_manual(values = order_palette) +
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +
  #coord_flip() +
  labs(x = "", y = "Relative Abundance [%]")


# Calculate total abundance per sample
sample.order.top.melt <- sample.order.top.melt %>%
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))

# sort samples in descending order
sample.order.top.bar.ab <- ggplot(sample.order.top.melt,aes(x=fct_reorder(Sample, -total_abundance), y=Abundance, fill = order))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3)+
  scale_fill_manual(values = order_palette)+ #coord_flip() +
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + # x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  theme(legend.title=element_text(size=11)) + # legend text
  labs(x="",y="Relative Abundance [%]") # scale_x_continuous(breaks=seq(0,6,1))

# reverse sample order for coord_flip
sample.order.top.melt$Sample2 <- factor(sample.order.top.melt$Sample, levels = rev(c(subfamily_ordered)))

sample.order.top.bar.flip <- ggplot(sample.order.top.melt,aes(x = Sample2, y=Abundance, fill = order))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = order_palette)+
  theme_grid() + 
  scale_y_continuous(labels = scales::label_percent(scale = 100,)) + # y-axis 100% prefix = "", suffix = ""
  #guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() + labs(x="",y="rel. ab. [%]")


# Total abundance of top 8 order per host_subfamily
top_order_abundance <- sample.order.top.melt %>%
  group_by(Sample) %>%
  summarise(total_abundance = sum(Abundance))

order_stats <- top_order_abundance %>%
  summarise(
    mean_total_abundance = mean(total_abundance),
    sd_total_abundance = sd(total_abundance)
  )


#### Group merged Top family ----
sample.family.merged = merge_samples(sample.family, "host_subfamily")
sample.family.merged.rel <- transform_sample_counts(sample.family.merged, function(x) x/sum(x))
data.frame(sort(taxa_sums(sample.family.merged.rel)*100, decreasing = FALSE))

# select top 10 family
Top.family <- names(sort(taxa_sums(sample.family.merged), decreasing=T)[1:18]) #top 15 no Bacillaceae
sample.family.top <- subset_taxa(sample.family.merged.rel, taxa_names(sample.family.merged.rel)%in%Top.family)
sample.family.top.melt <- psmelt(sample.family.top)
topfamilyCount = length(unique(sample.family.top.melt$family))

# Calculate total abundance per sample
sample.family.top.melt <- sample.family.top.melt %>%
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))

# Calculate genus abundance for color order
family_abundance <- sample.family.top.melt  %>%
  group_by(family) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance))  # Sort by abundance (most abundant first)

family_palette <- setNames(taxaPalette(topfamilyCount), family_abundance$family)  # Assign colors in the abundance order
sample.family.top.melt$family <- factor(sample.family.top.melt$family, levels = family_abundance$family) # order names by order

sample.family.top.bar.flip <- ggplot(sample.family.top.melt,aes(x = fct_reorder(Sample, total_abundance), y=Abundance, fill = family))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = family_palette)+
  theme_grid() + #theme(axis.text.y = element_text(face="italic")) + # x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  theme(legend.title=element_text(size=11)) + # legend.text=element_text(face="italic")
  #guides(fill = guide_legend(reverse = TRUE)) + # reverses bar fill
  coord_flip() + labs(x="",y="Relative Abundance [%]", fill="Top 12 families")

sample.family.top.bar <- ggplot(sample.family.top.melt,aes(x=fct_reorder(Sample, -total_abundance), y=Abundance, fill = family))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3)+
  scale_fill_manual(values = family_palette)+
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + # x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  #theme(legend.title=element_text(size=11), legend.text=element_text(size=10)) + # legend text
  labs(x="",y="Relative Abundance [%]") 

# Total abundance of top families per host_subfamily
top_family_abundance <- sample.family.top.melt %>%
  group_by(Sample) %>%
  summarise(total_abundance = sum(Abundance))

top_family_stats <- top_family_abundance %>%
  summarise(
    mean_total_abundance = mean(total_abundance),
    sd_total_abundance = sd(total_abundance)
  )



pdf("plots_peru/01_Group_merged_comp_barplot.pdf", width=8, height=6)
Top.comp.merged.phylum
Top.comp.merged.order
Top.comp.merged.family
Top.comp.merged.genus
sample.species.top.bar
sample.species.top.bar.group
sample.order.top.bar.flip
sample.order.top.bar
sample.family.top.bar.flip
sample.family.top.bar
dev.off()

sink("plots_peru/01_Group_merged_comp_barplot.txt")
"top_order_abundance"
as.data.frame(top_order_abundance)
"top_family_abundance"
as.data.frame(top_family_abundance)
"summary_stats order"
as.data.frame(order_stats)
"summary_stats family"
as.data.frame(top_family_stats)
sink()



## 01 Core analysis -----------

rev_palette <- rev(brewer.pal(10, "RdYlBu"))
# Remove color 5
heat_core_palette <- rev_palette[-5]

# ASV core
#ps1.rel <- microbiome::transform(sample.species, "compositional")
#Set different detection levels and prevalence
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- c(0.001, 0.003, 0.01, 0.05) # 0.1%, 0.5%, 1%, 5%

ASV.core <- plot_core(sample.ASV.rel, plot.type = "heatmap", 
                        colours = heat_core_palette, #rev(brewer.pal(10, "RdYlBu")), #RdBu
                        min.prevalence = 0.20, # 0.1 is 10%
                        prevalences = prevalences, 
                        detections = detections,
                        horizontal = F) +
  xlab("rel. ab. [%]") +theme_line2() +
  theme(axis.text.x= element_text(size=8),axis.text.y= element_text(face="italic"),legend.text = element_text(size=8), legend.title= element_text(size=9)) 
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
                               colours = heat_core_palette, # RdYlBu RdBu Spectral
                               min.prevalence = 0.25, # 25%  
                               prevalences = prevalences, 
                               detections = detections,
                               horizontal = F               ) +
  xlab("rel. ab. [%]") +theme_line2() +
  theme(axis.text.x= element_text(size=8),axis.text.y= element_text(face="italic"),legend.text = element_text(size=8), legend.title= element_text(size=9)) 
  #scale_fill_viridis_c(option = "inferno")


core.data <- genus.core$data  
max_prev <- max(core.data$Prevalence, na.rm = TRUE)

genus.core.maxprev <- plot_core(sample.species.rel, plot.type = "heatmap", 
                        min.prevalence = 0.25, # 0.1 is 10% 0.35
                        prevalences = prevalences, 
                        detections = detections,
                        horizontal = F               ) +
  xlab("rel. ab. [%]") +theme_line2() +
  scale_fill_gradientn(
    colours = heat_core_palette,
    limits = c(0, max_prev),
    oob = scales::squish,
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme(axis.text.x= element_text(size=8),axis.text.y= element_text(face="italic"),legend.text = element_text(size=8), legend.title= element_text(size=9))

# Reorder genus names based on prevalence 
core_0.1 <- core.data %>% filter(DetectionThreshold == 0.01) #1% 0.01 
genus_order <- core_0.1 %>% arrange((Prevalence)) %>% pull(Taxa)

genus.core.order <- plot_core(sample.species.rel, plot.type = "heatmap", 
                         #colours = rev(brewer.pal(10, "RdYlBu")), # RdYlBu RdBu Spectral YlGnBu
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
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme(axis.text.x= element_text(size=8),axis.text.y= element_text(face="italic"),legend.text = element_text(size=8), legend.title= element_text(size=9))



### Core family analysis 
#Set different detection levels and prevalence
fam_prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
fam_detections <- c(0.002, 0.01, 0.05, 0.1) # 0.1%, 0.5%, 1%, 5% fam_detections <- c(0.001, 0.003, 0.01, 0.05)

family.core <- plot_core(sample.family.rel, plot.type = "heatmap", 
                         colours = heat_core_palette, #rev(brewer.pal(10, "RdYlBu")),
                         min.prevalence = 0.3, # 0.1 is 10%
                         prevalences = fam_prevalences, 
                         detections = fam_detections,
                         horizontal = F) +
  xlab("rel. ab. [%]") +theme_line2() #+
  #theme(axis.text.x= element_text(size=8),axis.text.y= element_text(face="italic"),legend.text = element_text(size=8), legend.title= element_text(size=9))


### Core order analysis 
#Set different detection levels and prevalence
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- c(0.002, 0.01, 0.05, 0.1) # 0.1%, 1%, 5%, 10%
order.core <- plot_core(sample.order.rel, plot.type = "heatmap", 
                         colours = heat_core_palette, #rev(brewer.pal(10, "RdYlBu")),
                         min.prevalence = 0.25, # 0.1 is 10%
                         prevalences = prevalences, 
                         detections = detections,
                         horizontal = F) +
  xlab("rel. ab. [%]") +theme_line2() +
  theme(axis.text.x= element_text(size=8),legend.text = element_text(size=8), legend.title= element_text(size=9))


# core genus definition detection 1.0% (0.01)  prevalence  20% (20/100)
check_taxa <- subset_taxa(sample.species.rel, genus == "Wolbachia")
check_taxa <- subset_taxa(sample.species.rel, genus == "Entomomonas")
check_taxa <- subset_taxa(sample.species.rel, genus == "Apibacter") 

sort(data.frame(sum = sample_sums(check_taxa)*100,
                host_genus = sample_data(check_taxa)$host_genus), decreasing = FALSE)

# Show Prevalence at 1% of top 20 taxa
prevalence(sample.species.rel, detection = 5/100, sort = TRUE)[1:20] # 5% high ab; 10% prev
prevalence(sample.species.rel, detection = 1/100, sort = TRUE)[1:20] # 1% medium ab; 20% prev
prevalence(sample.species.rel, detection = 1/500, sort = TRUE)[1:20] # 0.2% low ab; 40% prev

core.genus.high <- core_members(sample.species.rel, detection = 0.05, prevalence = 10/100) # 11 taxa
core.genus <- core_members(sample.species.rel, detection = 0.01, prevalence = 20/100) # 12 taxa
core.genus.low <- core_members(sample.species.rel, detection = 0.002, prevalence = 40/100) # 12 taxa

core.genus.low2 <- core_members(sample.species.rel, detection = 0.005, prevalence = 25/100) # 11 taxa

# Make a robust core by combining different definitions 
core_overlap <- list(high = core.genus.high, med = core.genus,
                    low = core.genus.low)
names(core_overlap) <- c("High abundance\nRA=5% Prev=10%", "Medium abundance\nRA=1% Prev=20%", "Low abundance\nRA=0.2% Prev=40%")

core_robust <- Reduce(intersect, core_overlap) # 9 taxa from all core definitons

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
venn_core_plot <-  ggVennDiagram(core_genus_lists, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "none") +
  labs(title = "Core Venn")

core_genus_lists_overlap <- Reduce(intersect, core_genus_lists)

venn_core_plot2 <-  ggVennDiagram(core_genus_lists2, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "none") +
  labs(title = "Core Venn")

core_genus_lists2_overlap <- Reduce(intersect, core_genus_lists2)

venn_core_plot3 <-  ggVennDiagram(core_genus_lists3, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "none") +
  labs(title = "Core Venn")

core_genus_lists3_overlap <- Reduce(intersect, core_genus_lists3)

venn_core_overlap <-  ggVennDiagram(core_overlap, label = "count", label_alpha = 0, set_size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "none") +
  #labs(title = "Robust Core Definitions") +
  scale_x_continuous(expand = expansion(mult = 0.22)) +  # Make Venn smaller
  scale_y_continuous(expand = expansion(mult = 0.1))

# Modify line thickness in venn diagram
venn_core_overlap$layers[[2]]$aes_params$linewidth <- 0.2

# core_overlap
# high & medium  = Entomomonas & Apibacter
# medium & low =  Acinetobacter
# low = Carnobacterium Raoultella




#### Core boxplot subfam-----

# Show percent of core genus per group
core.genus.rel <- prune_taxa(core.genus, sample.species.rel)
# same as: core.genus.rel2 <- core(sample.species.rel, detection = 0.01, prevalence = 20/100)
# same as: core.genus.rel3 <- subset_taxa(sample.species.rel, taxa_names(sample.species.rel)%in%core.genus)
tax_table(core.genus.rel)

#core.abundance.per.sample <- sample_sums(core(sample.species.rel, detection = .01, prevalence = .2))

# core genus on ASV level 
core.genus.ASV.names <- taxa_names(sample.ASV)[tax_table(sample.ASV)[, "genus"] %in% core.genus]
core.genus.ASVs <- prune_taxa(core.genus.ASV.names, sample.ASV)


# Sum relative abundance of core taxa per sample
core.genus.rel.df <- core.genus.rel %>%
  psmelt() %>%  # convert phyloseq object to a dataframe
  group_by(Sample, host_genus, host_family, host_subfamily,host_tribe, country,elevation,kw,location, sublocation,sex,no,weight) %>%
  summarise(genus_core = sum(Abundance)) 

core.genus.rel.df$weight <- as.numeric(as.character(core.genus.rel.df$weight))


# Reorder host_genus based on mean Core_Abundance
core.genus.rel.df <- core.genus.rel.df %>%
  group_by(host_genus) %>%
  mutate(mediangenus_core = median(genus_core))  # calculate the mean core abundance

# Boxplot showing core relative abundance per host genus
core.genus.boxplot <- ggplot(core.genus.rel.df, aes(x = fct_reorder(host_genus, mediangenus_core), y = genus_core)) +
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 4, aes(color = host_subfamily, shape=country), alpha = 0.8) +
  theme_line2() + theme(axis.text.y = element_text(face="italic") ) + # x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
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
  theme_line2() + #theme(axis.text.y = element_text(face="italic") ) + # x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  scale_colour_manual(values=subfamily_colorfix) +
  coord_flip() + labs(x="",y="genus core [%]", color = "host subfamily")




#### Core genus barplot (grouped by genus)
### Core genus merged by host_genus 
species.genus = merge_samples(sample.species, "host_genus")
species.genus.rel <- transform_sample_counts(species.genus, function(x) x/sum(x))
species.genus.core <- subset_taxa(species.genus.rel, taxa_names(species.genus.rel)%in%core.genus)
species.genus.core.melt <- psmelt(species.genus.core)

# Calculate total abundance per sample
species.genus.core.melt <- species.genus.core.melt %>%
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))

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

core.genus.abundance <- ggplot(species.genus.core.melt,aes(x = fct_reorder(Sample, total_abundance), y=Abundance, fill = genus))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = core_palette_alpha) +
  theme_grid() + theme(axis.text.y = element_text(face="italic")) + # x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  coord_flip() + labs(x="",y="genus core [%]", fill="genus core") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))




### Core genus merged by host_subfamily  
species.subfamily = merge_samples(sample.species, "host_subfamily")
species.subfamily.rel <- transform_sample_counts(species.subfamily, function(x) x/sum(x))
species.subfamily.core <- subset_taxa(species.subfamily.rel, taxa_names(species.subfamily.rel)%in%core.genus)
species.subfamily.core.melt <- psmelt(species.subfamily.core)

# Calculate total abundance per sample
species.subfamily.core.melt <- species.subfamily.core.melt %>%
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))

#species.subfamily.core.melt$genus <- factor(species.subfamily.core.melt$genus, levels = genus_ordered) # order genus names by abundance
species.subfamily.core.melt$genus <- factor(species.subfamily.core.melt$genus, levels = core_genus_abundance$genus) # order genus names by abundance

species.subfamily.core.melt$host_subfamily_ordered <- factor(species.subfamily.core.melt$Sample, levels = rev(c(subfamily_ordered)))

core.genus.abundance.subfam <- ggplot(species.subfamily.core.melt,aes(x = host_subfamily_ordered, y=Abundance, fill = genus))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3
    ) + 
  scale_fill_manual(values = core_palette_alpha) +
  theme_line2() + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""),limits = c(0, 1) ) + # y-axis 100%
  coord_flip() + labs(x="",y="genus core [%]") 
  

core.genus.abundance.subfam2 <- ggplot(species.subfamily.core.melt,aes(x = fct_reorder(Sample, -total_abundance), y=Abundance, fill = genus))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = core_palette) +
  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  labs(x="",y="genus core rel. abundance [%]")


# Total abundance of core genus per host_subfamily
genus_core_subfam_stats <- species.subfamily.core.melt %>%
  group_by(Sample) %>%
  summarise(total_abundance = sum(Abundance))

genus_core_subfam_stats

subfam_stats <- genus_core_subfam_stats %>%
  summarise(
    mean_total_abundance = mean(total_abundance),
    sd_total_abundance = sd(total_abundance)
  )






#### Core sample  -----
core.genus.rel.melt <- psmelt(core.genus.rel)

# Calculate total abundance per sample
core.genus.rel.melt <- core.genus.rel.melt  %>% 
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))

# Calculate genus abundance for color order
core_genus_order <- core.genus.rel.melt %>%
  group_by(genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance))  # Sort by abundance (most abundant first)

# use global core_palette definition
core.genus.rel.melt$genus <- factor(core.genus.rel.melt$genus, levels = core_genus_order$genus) # order genus names by abundance

sample.core.abundance <- ggplot(core.genus.rel.melt,aes(x = fct_reorder(Sample, total_abundance), y=Abundance, fill = genus))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = core_palette) +
  theme_grid() + theme(axis.text.y = element_text(face="italic")) + 
  coord_flip() + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""),limits = c(0, 1)) +
  labs(x="",y="core [%]", fill="core") +
  facet_wrap(~country*host_subfamily, scales = "free") +
  scale_x_discrete(labels = setNames(core.genus.rel.melt$host_genus, core.genus.rel.melt$Sample)) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))

# Total abundance of core genus across all samples
genus_core_sample_stats <- core.genus.rel.melt %>%
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
  theme_line2() + theme(axis.text.y = element_text(face="italic") ) + # x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  scale_colour_manual(values=subfamily_colorfix) +
  coord_flip() + labs(x="",y="family core [%]", color = "host subfamily")

core.family.boxplot.subfam <- ggplot(core.family.genus.df, aes(x = fct_reorder(host_subfamily, meanfamily_core), y = family_core)) +
  geom_violin( draw_quantiles = c(0.5),  trim = TRUE,
  #adjust = 1.5,         # Smoother curves (optional)
  scale = "width" ) +      # Keep widths consistent across groups 
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 4, aes(color = host_subfamily), alpha = 0.8) +
  theme_line2() +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
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

# Calculate total abundance per sample
sample.family.merged.genus.core.melt <- sample.family.merged.genus.core.melt %>%
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))

# Calculate genus abundance for color order
fam_abundance <- sample.family.merged.genus.core.melt %>%
  group_by(family) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance))  # Sort by abundance (most abundant first)

corefamilyCount = length(unique(sample.family.merged.genus.core.melt$family))
fam_core_palette <- setNames(taxaPalette(corefamilyCount), fam_abundance$family)  # Assign colors in the abundance order
sample.family.merged.genus.core.melt$family <- factor(sample.family.merged.genus.core.melt$family, levels = fam_abundance$family) # order genus names by 

core.family.abundance <- ggplot(sample.family.merged.genus.core.melt,aes(x = fct_reorder(Sample, total_abundance), y=Abundance, fill = family))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = fam_core_palette) +
  theme_grid() + theme(axis.text.y = element_text(face="italic")) + # x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  #scale_colour_manual(values=subfamily_color) +
  coord_flip() + labs(x="",y="rel. ab [%]", fill ="core families")


### Core family merged by host_subfamily 
sample.family.subfamily = merge_samples(sample.family, "host_subfamily")
sample.subfamily.rel <- transform_sample_counts(sample.family.subfamily, function(x) x/sum(x))
sample.subfamily.core <- subset_taxa(sample.subfamily.rel, taxa_names(sample.subfamily.rel)%in%core.family)
sample.subfamily.core.melt <- psmelt(sample.subfamily.core)


# Calculate total abundance per sample
sample.subfamily.core.melt <- sample.subfamily.core.melt %>%
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))

sample.subfamily.core.melt$family <- factor(sample.subfamily.core.melt$family, levels = fam_abundance$family) # order genus names by 
# rev(fam_abundance$family)) to reverse order of stacked core families, include guides(fill = guide_legend(reverse = TRUE)) # reverse legend order

core.family.abundance.subfam <- ggplot(sample.subfamily.core.melt,aes(x = fct_reorder(Sample, total_abundance), y=Abundance, fill = family))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = fam_core_palette) +
  theme_grid() +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  coord_flip() + labs(x="",y="family core [%]") 
  

#### Core family group 'others' in gray ----
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

pdf("plots_peru/01_sample_core_analysis.pdf", width=8, height=6)
venn_core_plot 
venn_core_plot2
venn_core_overlap
genus.core
family.core
core.genus.boxplot
core.genus.boxplot.subfam 
core.genus.abundance
core.genus.abundance.subfam 
sample.core.abundance
core.family.boxplot
core.family.boxplot.subfam
core.family.abundance
core.family.abundance.subfam 
core.family.abundance.others
ggarrange.genus.core.subfamily
ggarrange.fam.core.subfamily
dev.off()

sink("plots_peru/01_sample_core_analysis.txt")
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



### non-core genus group minor by 'other'-------------
# Get non-core taxa
noncore.samples <- subset_taxa(sample.species.rel, !taxa_names(sample.species.rel) %in% core.genus)
data.frame(sort(taxa_sums(noncore.samples), decreasing = TRUE)[1:30])

# Quick inspect taxa of noncore.samples
# Remove samples with low non core amount (below 40%)
noncore.filtered <- prune_samples(sample_sums(otu_table(noncore.samples)) >= 0.4,  noncore.samples)
data.frame(sort(taxa_sums(noncore.filtered), decreasing = TRUE)[1:30])

noncore.filtered.comp <- noncore.filtered %>%
  #ps_filter(type == "sample") %>%
  comp_barplot(tax_level = "genus", n_taxa = 18, merge_other = F, label = "host_genus", sample_order = "bray") +
  facet_wrap(vars(country), scales = "free") +  
  coord_flip() + ggtitle( "noncore.filtered")

top_noncore_relab <- (taxa_sums(noncore.samples)*100 /  nsamples(noncore.samples)) %>%
  sort(decreasing = TRUE) %>%
  head(30) %>%
  as.data.frame()

top_noncore_relab # Select number of top non core taxa to show

noncore.samples.df <- psmelt(noncore.samples)

# Compute top taxa
top_taxa_samples <- noncore.samples.df %>%
  group_by(Taxon = genus) %>%   
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  mutate(label = if_else(row_number() <= 16, as.character(Taxon), "Other")) # Show top 20, rest as other

noncore.samples.df <- noncore.samples.df %>%
  left_join(top_taxa_samples %>% select(Taxon, label), by = c("genus" = "Taxon"))

noncore.samples.df <- noncore.samples.df %>%
  mutate(label = fct_reorder(label, Abundance, .fun = sum, .desc = TRUE)) %>%  # sort taxa by total abundance
  mutate(label = fct_relevel(label, "Other", after = Inf))  # move "Other" to end

noncore.samples.df.sort <- noncore.samples.df  %>% # Sort samples by abundance
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))

# noncore palette colors
noncore.samples_labels <- setdiff(unique(noncore.samples.df$label), "Other")
palette_colors <- noncorePalette(length(noncore.samples_labels))
noncore.samplescolors <- c(setNames(palette_colors, noncore.samples_labels), Other = "grey70")

single.sample.noncore.abundance <- ggplot(noncore.samples.df.sort,aes(x = fct_reorder(Sample, -total_abundance), y=Abundance, fill = label))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = noncore.samplescolors) +
  theme_grid() + theme(axis.text.y = element_text(face="italic")) + 
  coord_flip() + scale_y_reverse(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""),limits = c(1, 0)) +
  labs(x="",y="non-core [%]", fill="non-core") +
  #facet_wrap(~host_subfamily*country, scales = "free") +
  facet_wrap(~country*host_subfamily, scales = "free") +
  scale_x_discrete(labels = setNames(noncore.samples.df.sort$host_genus, noncore.samples.df.sort$Sample)) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))


# non-core taxa grouped by host_genus
noncore.genus <- subset_taxa(species.genus.rel, !taxa_names(species.genus.rel) %in% core.genus)
data.frame(sort(taxa_sums(noncore.genus), decreasing = TRUE)[1:30])

noncore.genus.df <- psmelt(noncore.genus)

# Compute top taxa
top_taxa <- noncore.genus.df %>%
  group_by(Taxon = genus) %>%   
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  mutate(label = if_else(row_number() <= 17, as.character(Taxon), "Other")) # Show top 20, rest as other

noncore.genus.df <- noncore.genus.df %>%
  left_join(top_taxa %>% select(Taxon, label), by = c("genus" = "Taxon"))


noncore.genus.df <- noncore.genus.df %>%
  mutate(label = fct_reorder(label, Abundance, .fun = sum, .desc = TRUE)) %>%  # sort taxa by total abundance
  mutate(label = fct_relevel(label, "Other", after = Inf))  # move "Other" to end


noncore.genus.df.sort <- noncore.genus.df  %>% # Sort samples by abundance
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))


noncore_labels <- setdiff(unique(noncore.genus.df$label), "Other")
#palette_colors <- taxaPalette(length(noncore_labels)+4)[5:(length(noncore_labels) + 4)] #  select from color 4 on
palette_colors <- noncorePalette(length(noncore_labels))
noncoregenuscolors <- c(setNames(palette_colors, noncore_labels), Other = "grey70")

# non core genus abundance by host_genus
noncore.genus.abundance <- ggplot(noncore.genus.df.sort,aes(x = fct_reorder(Sample, -total_abundance), y=Abundance, fill = label))+
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
noncore.sample.taxa.df <- noncore.genus.df %>%
  mutate(label = if_else(genus %in% noncore.samples_labels, genus, "Other"))

noncore.sample.taxa.df <- noncore.sample.taxa.df %>%
  mutate(label = fct_reorder(label, Abundance, .fun = sum, .desc = TRUE)) %>%  # sort taxa by total abundance
  mutate(label = fct_relevel(label, "Other", after = Inf))  # move "Other" to end

noncore.sample.taxa.df.sort <- noncore.sample.taxa.df   %>% # Sort samples by abundance
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))

noncoregenuscolors.sample <- c(setNames(noncorePalette(length(noncore.samples_labels)), noncore.samples_labels), Other = "grey70")

# Color based on total sample abundance (does not color biggest bars)
noncore.genus.abundance.samplecolor <- ggplot(noncore.sample.taxa.df.sort,aes(x = fct_reorder(Sample, -total_abundance), y=Abundance, fill = label))+
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
noncore.taxa.fam <- subset_taxa(species.subfamily.rel, !taxa_names(species.subfamily.rel) %in% core.genus)
noncore.taxa.fam.df <- psmelt(noncore.taxa.fam) 
data.frame(sort(taxa_sums(noncore.taxa.fam), decreasing = TRUE)[1:30])

top_taxa <- noncore.taxa.fam.df %>%
  group_by(Taxon = genus) %>%   
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  mutate(label = if_else(row_number() <= 12, as.character(Taxon), "Other")) # Show top 20, rest as other

noncore.taxa.fam.df <- noncore.taxa.fam.df %>%
  left_join(top_taxa %>% select(Taxon, label), by = c("genus" = "Taxon"))

noncore.taxa.fam.df <- noncore.taxa.fam.df %>%
  mutate(label = fct_reorder(label, Abundance, .fun = sum, .desc = TRUE)) %>%  # sort taxa by total abundance
  mutate(label = fct_relevel(label, "Other", after = Inf))  # move "Other" to end

noncore.taxa.fam.df.sort <- noncore.taxa.fam.df  %>% # Sort samples by abundance
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))

all_labels <- unique(noncore.taxa.fam.df$label)
noncore_labels <- setdiff(all_labels, "Other")
palette_colors <- taxaPalette(length(noncore_labels))
noncorecolors <- c(setNames(palette_colors, noncore_labels), Other = "grey70")

noncore.genus.abundance.subfam <- ggplot(noncore.taxa.fam.df.sort,aes(x = fct_reorder(Sample, -total_abundance), y=Abundance, fill = label))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = noncorecolors) +
  theme_grid() + #theme(axis.text.y = element_text(face="italic")) + # x-axis angle
  coord_flip() + scale_y_reverse(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +
  labs(x="",y="non core [%]",fill="non-core")

#### non-core grouped by family -------------

noncore.taxa.family <- aggregate_taxa(noncore.genus, "family")
data.frame(sort(taxa_sums(noncore.taxa.family), decreasing = TRUE)[1:30])
noncore.taxa.family.df <- psmelt(noncore.taxa.family)

top_family <- noncore.taxa.family.df %>%
  group_by(Taxon = family) %>%   
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  mutate(label = if_else(row_number() <= 10, as.character(Taxon), "Other")) # Show top 20, rest as other

noncore.taxa.family.df <- noncore.taxa.family.df %>%
  left_join(top_family %>% select(Taxon, label), by = c("family" = "Taxon"))

noncore.taxa.family.df <- noncore.taxa.family.df %>%
  mutate(label = fct_reorder(label, Abundance, .fun = sum, .desc = TRUE)) %>%  # sort taxa by total abundance
  mutate(label = fct_relevel(label, "Other", after = Inf))  # move "Other" to end

noncore.taxa.family.df.sort <- noncore.taxa.family.df  %>% # Sort samples by abundance
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))


noncore_fam <- setdiff(unique(noncore.taxa.family.df$label), "Other")
fam_colors <- noncorePalette(length(noncore_fam))
noncorefamcolors <- c(setNames(fam_colors, noncore_fam), Other = "grey70")

# non core genera on family level
noncore.family.abundance <- ggplot(noncore.taxa.family.df.sort,aes(x = fct_reorder(Sample, -total_abundance), y=Abundance, fill = label))+
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
noncore.subfamily <- subset_taxa(subfamily.merged.rel, !taxa_names(subfamily.merged.rel) %in% core.genus)
noncore.subfamily.family <- aggregate_taxa(noncore.subfamily, "family")
data.frame(sort(taxa_sums(noncore.subfamily.family), decreasing = TRUE)[1:30])

noncore.subfamily.family.df <- psmelt(noncore.subfamily.family)

top_fam.sub <- noncore.subfamily.family.df %>%
  group_by(Taxon = family) %>%   
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  mutate(label = if_else(row_number() <= 10, as.character(Taxon), "Other")) # Show top 20, rest as other

noncore.subfamily.family.df <- noncore.subfamily.family.df %>%
  left_join(top_fam.sub %>% select(Taxon, label), by = c("family" = "Taxon"))

noncore.subfamily.family.df <- noncore.subfamily.family.df %>%
  mutate(label = fct_reorder(label, Abundance, .fun = sum, .desc = TRUE)) %>%  # sort taxa by total abundance
  mutate(label = fct_relevel(label, "Other", after = Inf))  # move "Other" to end

noncore.subfamily.family.df.sort <- noncore.subfamily.family.df  %>% # Sort samples by abundance
  group_by(Sample) %>%
  mutate(total_abundance = sum(Abundance))


noncore_fam <- setdiff(unique(noncore.subfamily.family.df$label), "Other")
fam_colors <- noncorePalette(length(noncore_fam))
noncorefamcolors <- c(setNames(fam_colors, noncore_fam), Other = "grey70")

noncore.family.subfamily <- ggplot(noncore.subfamily.family.df.sort,aes(x = fct_reorder(Sample, -total_abundance), y=Abundance, fill = label))+
  geom_bar(#position="fill",
    colour="black", stat="identity", linewidth=0.3) + 
  scale_fill_manual(values = noncorefamcolors) +
  theme_grid() + #theme(axis.text.y = element_text(face="italic")) + 
  coord_flip() + scale_y_reverse(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""),limits = c(1, 0)) +
  labs(x="",y="non-core [%]", fill="non-core taxa") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))


pdf("plots_peru/01_sample_core_noncore.pdf", width=8, height=6)
single.sample.noncore.abundance 
noncore.genus.abundance
noncore.genus.abundance.samplecolor
noncore.genus.abundance.subfam
noncore.family.abundance
noncore.family.subfamily
dev.off()

rm(single.sample.noncore.abundance)
rm(noncore.samples.df.sort)
rm(noncore.samples.df)




## 02 sample div alpha Shannon --------------

#### alpha div plots  
# categorical data for box plot
ASV.rich  <- plot_richness(sample.ASV,x="host_subfamily", measures=c("Shannon","Observed","Simpson")) +
  geom_boxplot(aes(group = host_subfamily)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  ggtitle("sample.ASV") + scale_colour_manual(values=subfamily_colorfix) 

species.rich  <- plot_richness(sample.species,x="host_subfamily", measures=c("Shannon","Observed","Simpson")) +
  geom_boxplot(aes(group = host_subfamily)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  ggtitle("sample.species") + scale_colour_manual(values=subfamily_colorfix) 

filter.rich  <- plot_richness(sample.filter,x="host_subfamily", measures=c("Shannon","Observed","Simpson")) +
  geom_boxplot(aes(group = host_subfamily)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  ggtitle("sample.filter") + scale_colour_manual(values=subfamily_colorfix) 

family.rich  <- plot_richness(sample.family,x="host_subfamily", measures=c("Shannon","Observed","Simpson")) +
  geom_boxplot(aes(group = host_subfamily)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) +
  ggtitle("sample.family") + scale_colour_manual(values=subfamily_colorfix) 

country.rich  <- plot_richness(sample.species,x="country", measures=c("Shannon","Observed","Simpson")) +
  geom_boxplot(aes(group = country)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  ggtitle("sample.species country") + scale_colour_manual(values=subfamily_colorfix) 

location.rich  <- plot_richness(sample.species,x="location", measures=c("Shannon","Observed","Simpson")) +
  geom_boxplot(aes(group = location)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  ggtitle("sample.species location") + scale_color_manual(values = subfamily_colorfix) 

sublocation.country.rich  <- plot_richness(sample.species,x="sublocation", measures=c("Shannon")) +
  geom_boxplot(aes(group = sublocation)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  ggtitle("sample.species sublocation") + scale_color_manual(values = subfamily_colorfix) +
  facet_wrap(~ country, scales = "free_x")

genus.country.rich  <- plot_richness(sample.species,x="host_genus", measures=c("Shannon")) +
  geom_boxplot(aes(group = host_genus)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  ggtitle("sample.species country") + scale_color_manual(values = subfamily_colorfix) +
  facet_wrap(~ country, scales = "free_x")

# continuous data for linear regression
sample.rich.time  <- plot_richness(sample.species,x="kw", measures=c("Shannon","Observed","Simpson")) +
  theme_line2()+ geom_point(size=4, aes(color=host_subfamily)) + 
  stat_smooth(method="lm", color="black", linewidth=0.5,se=T) + 
  stat_regline_equation(label.y = 0.0, aes(label = after_stat(rr.label))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  ggtitle("sample.species") + scale_color_manual(values = subfamily_colorfix) 



pdf("plots_peru/02_sample_div_alpha_Shannon.pdf", width=8, height=6)
ASV.rich
species.rich
filter.rich 
family.rich
country.rich
location.rich
sublocation.country.rich
genus.country.rich 
sample.rich.time
dev.off()



#### alphaframe Shannon --------------
#df.sample <- as(sample_data(sample.species),"data.frame")
df.sample <- as(sample_data(sample.species), "data.frame") %>%
  rownames_to_column("Sample")

df.sample$richness <- estimate_richness(sample.species, split=TRUE, measures=c("Shannon","Observed","InvSimpson"))

df.core <- core.genus.rel %>%
  psmelt() %>%  # convert phyloseq object to a dataframe
  group_by(Sample) %>%
  summarise(genus_core = sum(Abundance))

### New alphaframe simplify data frame for statistics 
alphaframe <- left_join(df.sample, df.core, by = "Sample")
alphaframe$Shannon <- alphaframe$richness$Shannon
alphaframe$Observed <- alphaframe$richness$Observed
alphaframe$InvSimpson<- alphaframe$richness$InvSimpson
alphaframe$hill_numbers <- exp(alphaframe$richness$Shannon)
alphaframe$time <- alphaframe$kw
alphaframe$temperature <- alphaframe$temp_week
alphaframe$host_subfamily_ordered <- factor(alphaframe$host_subfamily, levels = subfamily_ordered)
#alphaframe$group <- alphaframe$host_subfamily # if there is no group 


# adjust variables if necessary
str(alphaframe) # variables can be converted as.numeric or as.factor
# use time point as.factor only for figure color
# use time point as integer for statistics
#alphaframe$host_genus <- as.factor(alphaframe$host_genus)
#alphaframe$kw <- as.numeric(alphaframe$kw)
#alphaframe$group <- factor(df.sample$group, levels = c("t0", "t1", "t2", "t3"))

alphaframe.noAglais <- subset(alphaframe , host_genus != "Aglais")
alphaframe.Aglais <- subset(alphaframe , host_genus == "Aglais")

alphaframe <- alphaframe %>%
  group_by(host_tribe) %>%
  mutate(median_richness = median(Observed))  # calculate the median richness

alphaframe <- alphaframe %>%
  group_by(host_tribe) %>%
  mutate(median_Shannon = median(Shannon))  # calculate the median richness


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

alpha.rich.tribe2  <- ggplot(alphaframe, aes(x = fct_reorder(host_tribe, median_richness), y=Observed)) +
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
  geom_boxplot(aes(group = host_subfamily))+
  geom_point(position = position_jitter(h = 0.1, w = 0), size=4, aes(color=host_subfamily, shape = country), alpha = 0.8) +
  theme_line2()+ scale_y_discrete(limits = rev) + # reverse order
  scale_colour_manual(values=subfamily_colorfix) +
  labs(y="",x="Shannon Index", color ="host subfamily") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

alpha.shannon3  <- ggplot(alphaframe, aes(x=Shannon, y=host_subfamily_ordered)) +
  geom_violin(aes(group = host_subfamily), draw_quantiles = c(0.5), trim = F,  scale = "width")+
  geom_point(position = position_jitter(h = 0.1, w = 0), size=4, aes(color=host_subfamily, shape = country), alpha = 0.8) +
  theme_line2()+ scale_y_discrete(limits = rev) + # reverse order
  scale_colour_manual(values=subfamily_colorfix) +
  labs(y="",x="Shannon Index", color ="host subfamily") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )


alpha.shannon.tribe2 <- ggplot(alphaframe, aes(x = fct_reorder(host_tribe, median_Shannon), y=Shannon)) +
  geom_boxplot(aes(group = host_tribe))+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_subfamily, shape = country), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="",y="Shannon Index") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

alpha.rich.tribe.country <- ggplot(alphaframe, aes(x = fct_reorder(host_tribe, median_richness), y=Observed)) +
  geom_boxplot(aes(group = country))+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_subfamily, shape = country), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="",y="microbial richness (q=0)") +
  facet_wrap(~ country, scales = "free_x") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )


alpha.shannon.country  <- ggplot(alphaframe, aes(x=country, y=Shannon)) +
  geom_boxplot(aes(group = country), outlier.shape = NA)+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_subfamily, shape=country), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="",y="Shannon Index") + # , subtitle="country"
  theme(plot.subtitle = element_text(hjust = 0.5))

alpha.shannon.country2  <- ggplot(alphaframe, aes(x=country, y=Shannon)) +
  geom_violin(aes(group = country), draw_quantiles = c(0.5), trim = T,  scale = "area")+ 
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_subfamily, shape=country), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="",y="Shannon Index", color="host subfamily") + 
  theme(plot.subtitle = element_text(hjust = 0.5))

alpha.subfamily.country  <- ggplot(alphaframe, aes(x=country, y=Shannon)) +
  geom_boxplot(aes(group = country), outlier.shape = NA)+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_subfamily, shape=country), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="",y="Shannon Index", color="host subfamily") +
  facet_wrap(~ host_subfamily, scales = "free_x") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

alpha.subfamily.country2 <- ggplot(alphaframe, aes(x=country, y=Shannon)) +
  geom_violin(aes(group = country), draw_quantiles = c(0.5), trim = T,  scale = "width")+ # scale = "width" for equal width
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color=host_subfamily), alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="",y="Shannon Index", color = "host subfamily") +
  facet_wrap(~ host_subfamily, scales = "free_x")

alpha.subfamily.country3 <- ggplot(alphaframe, aes(x=country, y=Shannon)) +
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

violin.country.subfamilies <- ggplot(alphaframe, aes(x = country, y = Shannon, fill = country)) +
  geom_violin(trim = F,  scale = "width", alpha = 0.6) + # draw_quantiles = c(0.5),
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.5, aes(fill=country)) +
  #stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black", fatten = 2) +
  scale_fill_manual(values=country_colorfix) +
  theme_line2() +
  geom_point(position = position_jitter(w = 0.01, h = 0), size=3, alpha = 0.2,  aes(fill=country)) +
  facet_wrap(~ host_subfamily, scales = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x=NULL ,y="Shannon Index") 


pdf("plots_peru/02_sample_div_alpha_Shannon_01.pdf", width=8, height=6)
alpha.shannon
alpha.rich
alpha.rich.tribe
alpha.rich.tribe2
alpha.shannon2
alpha.shannon3 
alpha.shannon.tribe2
alpha.shannon.country
alpha.shannon.country2 
alpha.subfamily.country 
alpha.subfamily.country2 
alpha.subfamily.country3
alpha.location
alpha.location.tribe 
violin.country.subfamilies
dev.off()

#### alpha Shannon stats lm lmer -------------

# Show data groups that will be tested
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
Anova(lm.shannon, type=2) # independent order
Anova(lm.shannon, type=2, white.adjust = TRUE) # if unequal variances across groups
# location not significant when host_subfamily is in the model!!
vif(lm.shannon) # variance inflation factor should be <3  
plot(lm.shannon,2)
hist(lm.shannon$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon)) # Run Shapiro-Wilk test on residuals 

# Univariate tests Shannon
lm.shannon.family <- lm(Shannon ~ host_family, data=alphaframe)
summary(lm.shannon.family)$adj.r.squared # R-squared:  0.03161
Anova(lm.shannon.family, type=2)

lm.shannon.subfam <- lm(Shannon ~ host_subfamily , data = alphaframe)
summary(lm.shannon.subfam)$adj.r.squared # R-squared:  0.3088
Anova(lm.shannon.subfam, type=2)

lm.shannon.tribe <- lm(Shannon ~ host_tribe , data = alphaframe)
summary(lm.shannon.tribe)$adj.r.squared # R-squared:  0.3178
Anova(lm.shannon.tribe, type=2)

lm.shannon.genus <- lm(Shannon ~ host_genus  , data = alphaframe)
summary(lm.shannon.genus)$adj.r.squared # R-squared:  0.3916
Anova(lm.shannon.genus, type=2)

lm.shannon.country <- lm(Shannon ~ country, data = alphaframe)
summary(lm.shannon.country)$adj.r.squared # R-squared:  0.1491 
Anova(lm.shannon.country, type=2)

lm.shannon.location <- lm(Shannon ~ location, data = alphaframe)
summary(lm.shannon.location)$adj.r.squared # Adjusted R-squared:  0.191
Anova(lm.shannon.location, type=2)

lm.shannon.sublocation <- lm(Shannon ~ sublocation, data = alphaframe)
summary(lm.shannon.sublocation)$adj.r.squared # Adjusted R-squared:  0.2059 
Anova(lm.shannon.sublocation, type=2)

lm.shannon.temp <- lm(Shannon ~ temperature, data = alphaframe)
summary(lm.shannon.temp)$adj.r.squared # Adjusted R-squared:  0.05995  
Anova(lm.shannon.temp, type=2)

lm.shannon.year <- lm(Shannon ~ year , data = alphaframe)
summary(lm.shannon.year)$adj.r.squared # Adjusted R-squared:  0.01247  
Anova(lm.shannon.year, type=2)


# Compare model fit. Lower AIC values or higher adjusted R-squared values indicate a better fit.
AIC(lm.shannon.family,lm.shannon.subfam , lm.shannon.tribe, lm.shannon.genus, lm.shannon.country, lm.shannon.location, lm.shannon.sublocation, lm.shannon.temp, lm.shannon.year)
# Shannon best AIC genus, subfamily, tribe
# But subfamily better F value
# location > country > sublocation

# Quick model comparison
predictors <- c("host_family", "host_subfamily", "host_tribe",
                "host_genus", "country", "location", "sublocation",
                "temperature" , "year", "kw", "storage")

# Create a named list of models automatically
lm.list <- lapply(predictors, function(p) {
  formula <- as.formula(paste("Shannon ~", p))
  lm(formula, data = alphaframe) })

names(lm.list) <- predictors  # optional: set names for easy reference

# Extract F-values, p-values, and R²
lm.results <- lapply(lm.list, function(mod) {
  s <- summary(mod)
  f <- s$fstatistic
  F_value <- f[1]
  p_value <- pf(f[1], f[2], f[3], lower.tail = FALSE)
  R2 <- s$r.squared
  AIC_val <- AIC(mod)  # add AIC here
  c(F_value = F_value, p_value = p_value, R2 = R2, AIC = AIC_val) })

# Combine into a data frame
lm.results.summary <- do.call(rbind, lm.results)
lm.results.summary <- as.data.frame(lm.results.summary)
lm.results.summary$predictor <- rownames(lm.results.summary)
rownames(lm.results.summary) <- NULL

# Optional: sort by F-value (strongest predictors first)
lm.results.summary <- lm.results.summary[order(-lm.results.summary$R2), ]

lm.results.summary
# main effect of host taxonomy, but also location

# Multivariate tests, combine taxonomy with location
lm.shannon.combine1 <- lm(Shannon ~ host_subfamily + location, data = alphaframe)
lm.shannon.combine2 <- lm(Shannon ~ host_subfamily + location + temperature, data = alphaframe)
lm.shannon.combine3 <- lm(Shannon ~ host_tribe + location, data = alphaframe)
lm.shannon.combine4 <- lm(Shannon ~ host_tribe + location + temperature , data = alphaframe)
lm.shannon.combine5 <- lm(Shannon ~ host_genus + location  , data = alphaframe)
lm.shannon.combine6 <- lm(Shannon ~ host_genus + location + temperature , data = alphaframe)

AIC(lm.shannon.combine1,lm.shannon.combine2,lm.shannon.combine3,lm.shannon.combine4,lm.shannon.combine5, lm.shannon.combine6 )
# subfamily + location better than 3 predictors


# for host genus tests
# Robustness Check against over fitting  (filter genus with ≥2 samples)
genus_to_keep <- alphaframe %>%
  count(host_genus) %>%
  filter(n >= 2) %>%
  pull(host_genus)

# Filter the dataframe
alphaframe.robust <- alphaframe %>%
  filter(host_genus %in% genus_to_keep)

lm.shannon.genus.robust <- lm(Shannon ~ host_genus  , data = alphaframe.robust)
summary(lm.shannon.genus.robust) # R-squared:  0.3848
Anova(lm.shannon.genus.robust, type=2)

lm.shannon.robust.combine <- lm(Shannon ~ host_genus + location  , data = alphaframe.robust)
summary(lm.shannon.robust.combine) # R-squared:  0.397 
Anova(lm.shannon.robust.combine, type=2)
# check colinearity
#table(alphaframe.robust$host_genus, alphaframe.robust$location)
#vif(lm.shannon.robust.combine) # variance inflation factor should be <3 
#alias(lm.shannon.robust.combine) # aliased coefficients in the model: location + host_genus

# combine with location
lm.robust.combine1 <- lm(Shannon ~ host_genus + location, data = alphaframe.robust)
lm.robust.combine2 <- lm(Shannon ~ host_genus + location + temperature, data = alphaframe.robust)
lm.robust.combine3 <- lm(Shannon ~ host_genus + location + kw, data = alphaframe.robust)
lm.robust.combine4 <- lm(Shannon ~ host_genus + location + storage , data = alphaframe.robust)

AIC(lm.robust.combine1,lm.robust.combine2,lm.robust.combine3,lm.robust.combine4 )
# only genus + location better than 3 predictor

summary(lm.robust.combine1)$adj.r.squared # 0.3969554
summary(lm.robust.combine2)$adj.r.squared # 0.3929983
summary(lm.robust.combine3)$adj.r.squared # 0.3968534
summary(lm.robust.combine4)$adj.r.squared # 0.3928543

Anova(lm.robust.combine1, type=2)
Anova(lm.robust.combine2, type=2) # temp not sign
Anova(lm.robust.combine3, type=2) # kw not sign
Anova(lm.robust.combine4, type=2) # storage not sign


# shannon best final model: host subfamily alone
summary(lm.shannon.subfam) 
summary(lm.shannon.subfam)$adj.r.squared # R-squared:   0.3088208 
Anova(lm.shannon.subfam, type=2)
plot(lm.shannon.subfam,2)
hist(lm.shannon.subfam$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon.subfam)) # Run Shapiro-Wilk test on residuals 

# shannon best final model: host genus alone
summary(lm.shannon.genus.robust) 
summary(lm.shannon.genus.robust)$adj.r.squared # R-squared:   0.391604 
Anova(lm.shannon.genus.robust, type=2)
plot(lm.shannon.genus.robust,2)
hist(lm.shannon.genus.robust$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon.genus.robust)) # Run Shapiro-Wilk test on residuals 


lm.subfamily.country <- lm(Shannon ~ host_subfamily + country , data = alphaframe)
summary(lm.subfamily.country) 
summary(lm.subfamily.country)$adj.r.squared # R-squared:   0.3076139 
Anova(lm.subfamily.country, type=2)
vif(lm.subfamily.country) # variance inflation factor should be <3 
alias(lm.subfamily.country) # all good
table(alphaframe$host_subfamily, alphaframe$country) # check colinearity
plot(lm.subfamily.country,2)
hist(lm.subfamily.country$residuals) # extract the residuals
shapiro.test(residuals(lm.subfamily.country)) # Run Shapiro-Wilk test on residuals 


lm.shannon.model <- lm(Shannon ~  host_subfamily  + location     , data=alphaframe)
summary(lm.shannon.model)
summary(lm.shannon.model)$adj.r.squared # 0.3230627
Anova(lm.shannon.model, type=2)
# location not significant when host_subfamily is in the model!!
vif(lm.shannon.model) # variance inflation factor should be <3
alias(lm.shannon.model) # all good
table(alphaframe$host_subfamily, alphaframe$location) # check colinearity
plot(lm.shannon.model,2)
hist(lm.shannon.model$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon.model)) # Run Shapiro-Wilk test on residuals 

# variance partitioning for lm model on Type II ANOVA
testme <- lm.shannon.model
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


lm.shannon.model2 <- lm(Shannon ~  host_tribe  + country     , data=alphaframe)
summary(lm.shannon.model2)
summary(lm.shannon.model2)$adj.r.squared # 0.3230627
Anova(lm.shannon.model2, type=2)
vif(lm.shannon.model2) # variance inflation factor should be <3
alias(lm.shannon.model2) # all good
table(alphaframe$host_tribe, alphaframe$country) # check colinearity
# variance partitioning for lm model on Type II ANOVA
testme <- lm.shannon.model2
anovaII_out <- Anova(testme, type = 2)
ss_res <- sum(residuals(testme)^2)
# Partial R² per term
partial_R2 <- anovaII_out$`Sum Sq` / (anovaII_out$`Sum Sq` + ss_res)
names(partial_R2) <- rownames(anovaII_out)
partial_R2

# Partial R2 or eta squared (library effectsize)
eta_squared(anovaII_out, partial = TRUE)


#### alpha Observed stats lm lmer 
lm.q0.family <- lm(Observed ~ host_family, data=alphaframe)
summary(lm.q0.family)
anova(lm.q0.family)

lm.q0.fam <- lm(Observed ~ host_subfamily , data = alphaframe)
summary(lm.q0.fam) # R-squared:  0.4632
anova(lm.q0.fam)

lm.q0.tribe <- lm(Observed ~ host_tribe , data = alphaframe)
summary(lm.q0.tribe) # R-squared:  0.6075
anova(lm.q0.tribe)

lm.q0.genus <- lm(Observed ~ host_genus  , data = alphaframe)
summary(lm.q0.genus) # R-squared:  0.6921
anova(lm.q0.genus)

lm.q0.country <- lm(Observed ~ country, data = alphaframe)
summary(lm.q0.country) # R-squared:  0.5425 
anova(lm.q0.country)

lm.q0.location <- lm(Observed ~ location, data = alphaframe)
summary(lm.q0.location) # Adjusted R-squared:  0.541
anova(lm.q0.location)

lm.q0.sublocation <- lm(Observed ~ sublocation, data = alphaframe)
summary(lm.q0.sublocation) # Adjusted R-squared:  0.2059 
anova(lm.q0.sublocation)

lm.q0.temp <- lm(Observed ~ temperature, data = alphaframe)
summary(lm.q0.temp) # Adjusted R-squared:  0.05995  
anova(lm.q0.temp)

lm.q0.year <- lm(Observed ~ year + host_genus, data = alphaframe)
summary(lm.q0.year) # Adjusted R-squared:  0.01247  
anova(lm.q0.year)

lm.q0.combine <- lm(Observed ~  host_subfamily + location + temperature, data = alphaframe)
summary(lm.q0.combine) # Adjusted R-squared:  0.01247  
Anova(lm.q0.combine,type = 2)

# Compare model fit. Lower AIC values or higher adjusted R-squared values indicate a better fit.
AIC(lm.q0.family, lm.q0.fam , lm.q0.tribe, lm.q0.genus, lm.q0.country, lm.q0.location, lm.q0.sublocation, lm.q0.temp, lm.q0.combine)
# Observed best AIC genus, subfamily, tribe
# But subfamily better F value

# Test individual subfamilies Shannon country

lm.shannon.coliadinae <- lm(Shannon ~ country, data = subset(alphaframe, host_subfamily == "Coliadinae"))
summary(lm.shannon.coliadinae)
anova(lm.shannon.coliadinae)
plot(lm.shannon.coliadinae,2)
plot(lm.shannon.coliadinae, which = 1)
hist(lm.shannon.coliadinae$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon.coliadinae)) # Run Shapiro-Wilk test on residuals

# Two sample tests
wilcox.test(Shannon ~  country, data = subset(alphaframe, host_subfamily == "Coliadinae"))
t.test(Shannon ~  country, data = subset(alphaframe, host_subfamily == "Coliadinae"))


lm.shannon.dismorphinae <- lm(Shannon ~ country, data = subset(alphaframe, host_subfamily == "Dismorphiinae"))
summary(lm.shannon.dismorphinae)
anova(lm.shannon.dismorphinae)
plot(lm.shannon.dismorphinae,2)
plot(lm.shannon.dismorphinae, which = 1)
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
plot(lm.shannon.heliconiinae,2)
plot(lm.shannon.heliconiinae, which = 1)
hist(lm.shannon.heliconiinae$residuals) # extract the residuals
shapiro.test(residuals(lm.shannon.heliconiinae)) # Run Shapiro-Wilk test on residuals
 



#### alpha split country Peru Germany---------------
alphaframe.Peru <- subset(alphaframe, country == "Peru")


# correlation temp and elevation
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


core.elevation.Peru <-
  ggplot(alphaframe.Peru, aes(x = elevation, y = genus_core)) +
  theme_line2() + 
  geom_point(shape = 17, size = 4, alpha = 0.8, aes(color=host_subfamily)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(label.y = 0.5, aes(label = after_stat(rr.label)),size=3) + stat_cor(label.y = 0.4, size=3, p.accuracy = 0.00001) +
  scale_color_manual(values=subfamily_colorfix) +
  labs(x="Elevation [m]",y="genus core [%]", fill="host subfamily") 

# Quick inspect alpha
alpha.elevation.Peru  <- plot_richness(Peru.species,x="elevation", measures=c("Shannon","Observed", "Simpson", "Fisher")) +
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
  
 
## For Shannon diversity use lm or lmer
lm.peru.subfam <- lm(Shannon ~ host_subfamily, data=alphaframe.Peru) # subfamily
summary(lm.peru.subfam)
lm.peru.genus <- lm(Shannon ~ host_genus , data=alphaframe.Peru) # genus
summary(lm.peru.genus)
lm.peru.elevation <- lm(Shannon ~ elevation, data=alphaframe.Peru) # elevation
summary(lm.peru.elevation)
lm.peru.temp <- lm(Shannon ~ temperature, data=alphaframe.Peru) # elevation
summary(lm.peru.temp)
anova(lm.peru.elevation)
lm.elevation_host <- lm(Shannon ~ host_subfamily+elevation, data=alphaframe.Peru) # integer
summary(lm.elevation_host)
Anova(lm.elevation_host,type = 2) # elevation not significant when accounting for host genus

lm.temp_host <- lm(Shannon ~ host_subfamily + temperature, data=alphaframe.Peru) # integer
summary(lm.temp_host)
Anova(lm.temp_host,type = 2)

lm.peru.all <- lm(Shannon ~ host_subfamily + host_genus + temperature + year, data=alphaframe.Peru) # integer
summary(lm.peru.all)

# Compare model fit. Lower AIC values or higher adjusted R-squared values indicate a better fit.
AIC(lm.peru.subfam, lm.peru.genus,lm.peru.elevation,lm.peru.temp, lm.elevation_host,lm.temp_host,lm.peru.all )
# best AIC subfamily + temp

summary(lm.peru.subfam)$adj.r.squared # 0.2302111
summary(lm.peru.genus)$adj.r.squared # 0.3044573 best  R2 genus (potential overfitting)
summary(lm.peru.elevation)$adj.r.squared # 0.158668
summary(lm.temp_host)$adj.r.squared # 0.2660948

# Shannon peru
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
  
core.genus.kw.Germany <- ggplot(alphaframe.Germany, aes(x = kw, y = genus_core)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 4, aes(color = host_subfamily), alpha = 0.8) +
  theme_line2() + theme(axis.text.y = element_text(face="italic") ) + # x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="kw",y="genus core", color = "host subfamily")

core.genus.sublocation.Germany <- ggplot(alphaframe.Germany, aes(x = sublocation, y = genus_core)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 4, aes(color = host_subfamily), alpha = 0.8) +
  theme_line2() + theme(axis.text.y = element_text(face="italic") ) + # x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  scale_colour_manual(values=subfamily_colorfix) +
  labs(x="sublocation",y="genus core", color = "host subfamily")

core.genus.sex.Germany <- ggplot(alphaframe.Germany, aes(x = host_subfamily, y = genus_core)) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 4, aes(color = sex), alpha = 0.8) +
  theme_line2() + theme(axis.text.y = element_text(face="italic") ) + # x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
  labs(x="host subfamily",y="genus core", color = "sex")



sample.rich.kw.Germany <-
  ggplot(alphaframe.Germany, aes(x = kw, y = Observed)) +
  theme_line2() + 
  geom_point(size = 4, alpha = 0.8, aes(color=host_genus)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(label.y = 0.0,formula = y ~ poly(x, 3), aes(label = after_stat(rr.label)),size=3) + stat_cor(label.y = 2,size=3, p.accuracy = 0.001) +
  ggtitle("sample.rich.kw.Germany") 
  facet_wrap(~ host_genus) #, scales = "free_x"

# Fits best to cubic polynomial fit (y ~ poly(x, 3))

## Fitting Observed ~ time (kw) to cubic polynomial model
lm.kw <- lm(Observed ~ kw, data=alphaframe.Germany) # integer
summary(lm.kw)
# Multiple R-squared:  0.1821,	Adjusted R-squared:  0.1716 
# F-statistic: 17.37 on 1 and 78 DF,  p-value: 7.9e-05

lm.kw.quad <- lm(Observed ~ poly(kw, 2), data = alphaframe.Germany)
summary(lm.kw.quad)
# Multiple R-squared:  0.1842,	Adjusted R-squared:  0.163 
# F-statistic: 8.695 on 2 and 77 DF,  p-value: 0.0003938
 
lm.kw.cubic <- lm(Observed ~ poly(kw, 3), data = alphaframe.Germany)
summary(lm.kw.cubic)
# Multiple R-squared:  0.4971,	Adjusted R-squared:  0.4772 
# F-statistic: 25.04 on 3 and 76 DF,  p-value: 2.278e-11
plot(lm.kw.cubic,2)
hist(lm.kw.cubic$residuals) # extract the residuals
shapiro.test(residuals(lm.kw.cubic))

# Compare model fit. Lower AIC values or higher adjusted R-squared values indicate a better fit.
AIC(lm.kw, lm.kw.quad, lm.kw.cubic)
summary(lm.kw)$adj.r.squared
summary(lm.kw.quad)$adj.r.squared
summary(lm.kw.cubic)$adj.r.squared
# cubic polynomial model fits best!
  

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



# Shannon calendar week
lm.shannon.kw <- lm(Shannon ~ kw, data=alphaframe.Germany.noAglais) # integer
summary(lm.shannon.kw)
anova(lm.shannon.kw) # kw significant

lm.shannon.host <- lm(Shannon ~ host_subfamily, data=alphaframe.Germany.noAglais) # integer
summary(lm.shannon.host)
anova(lm.shannon.host) # 

lm.shannon.temp <- lm(Shannon ~ temperature, data=alphaframe.Germany.noAglais) # integer
summary(lm.shannon.temp)
anova(lm.shannon.temp) # temp not sign

lm.shannon.polytemp <- lm(Shannon ~ poly(temperature, 2), data=alphaframe.Germany.noAglais) # integer
summary(lm.shannon.polytemp)
anova(lm.shannon.polytemp) # 

ggplot(alphaframe.Germany.noAglais, aes(x = temperature, y = Shannon)) +
  theme_line2() + 
  geom_point(size = 4, alpha = 0.8, aes(color=host_subfamily)) +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "black", linewidth = 0.5, se = TRUE) +
  #stat_regline_equation(label.y = 2.0,formula = y ~ poly(x, 2), aes(label = after_stat(rr.label)),size=3) +# stat_cor(label.y = 4,size=3, p.accuracy = 0.001) +
  scale_color_manual(values = subfamily_colorfix) +
  labs(x="temp",y="Shannon",
       color = "host subfamily") 



lm.shannon.all <- lm(Shannon ~ temperature + host_subfamily + kw, data=alphaframe.Germany.noAglais) # integer
summary(lm.shannon.all)
Anova(lm.shannon.all, type=2) #
vif(lm.shannon.all)

# variance partitioning for lm model on Type II ANOVA
testme <- lm.shannon.all
anovaII_out <- Anova(testme, type = 2)
ss_res <- sum(residuals(testme)^2)
# Partial R² per term
partial_R2 <- anovaII_out$`Sum Sq` / (anovaII_out$`Sum Sq` + ss_res)
names(partial_R2) <- rownames(anovaII_out)
partial_R2


lm.shannon.kw.host <- lm(Shannon ~ kw+host_subfamily, data=alphaframe.Germany.noAglais) # integer
summary(lm.shannon.kw.host)
Anova(lm.shannon.kw.host, type=2)
# Calendar week explains Shannon diversity independently of genus identity
vif(lm.shannon.kw.host) # variance inflation factor should be <3

AIC(lm.shannon.kw, lm.shannon.host, lm.shannon.temp, lm.shannon.all, lm.shannon.kw.host)
summary(lm.shannon.kw)$adj.r.squared
summary(lm.shannon.all)$adj.r.squared
summary(lm.shannon.kw.host)$adj.r.squared # 0.3230627


# variance partitioning for lm model on Type II ANOVA
testme <- lm.shannon.kw.host
anovaII_out <- Anova(testme, type = 2)
ss_res <- sum(residuals(testme)^2)
# Partial R² per term
partial_R2 <- anovaII_out$`Sum Sq` / (anovaII_out$`Sum Sq` + ss_res)
names(partial_R2) <- rownames(anovaII_out)
partial_R2

# Partial R2 or eta squared (library effectsize)
eta_squared(anovaII_out, partial = TRUE)



lm.shannon.kw.genus <- lm(Shannon ~ kw +host_genus , data=alphaframe.Germany.noAglais) # integer
summary(lm.shannon.kw.genus)
# kw:host_genus not significant, no additional effect of species turnover 
Anova(lm.shannon.kw.genus, type=2)
# No evidence that the effect of kw depends on host_genus (i.e., no genus turnover pattern)

AIC(lm.shannon.kw.host, lm.shannon.kw.genus)
# subfamily better

# variance partitioning for lm model on Type II ANOVA
testme <- lm.shannon.kw.genus
anovaII_out <- Anova(testme, type = 2)
ss_res <- sum(residuals(testme)^2)
# Partial R² per term
partial_R2 <- anovaII_out$`Sum Sq` / (anovaII_out$`Sum Sq` + ss_res)
names(partial_R2) <- rownames(anovaII_out)
partial_R2

pdf("plots_peru/02_sample_div_alpha_Shannon_02_country.pdf", width=8, height=6)
temp.elevation
core.elevation.Peru
alpha.elevation.Peru 
sample.shannon.elevation.Peru
sample.rich.elevation.Peru
sample.shannon.temp.Peru 
kw.temp
core.genus.kw.Germany
core.genus.sublocation.Germany 
core.genus.sex.Germany
sample.rich.kw.Germany 
shannon.kw.Germany.noAglais
shannon.kw.Germany.noAglais2
dev.off()



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
    guide = guide_legend(title = "Country", order = 1)
  ) +
  scale_fill_manual(
    values = subfamily_colorfix,
    guide = guide_legend(
      override.aes = list(shape = 21, size = 4, color = "black"),
      title = "Aglais (outgroup)",
      order = 2
    )
  ) +
  scale_color_manual(values = subfamily_colorfix) +
  labs(x="temperature [°C]",y="Shannon Index",
       color = "host subfamily") +
  coord_cartesian(ylim = c(0, 4), expand = TRUE,  default = FALSE,   clip = "on" ) +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

lm.shan.fit1 <- lm(Shannon ~  host_subfamily , data=alphaframe.noAglais)
lm.shan.fit2 <- lm(Shannon ~  temperature , data=alphaframe.noAglais) 
lm.shan.fit3 <- lm(Shannon ~  location , data=alphaframe.noAglais) 
lm.shan.fit4 <- lm(Shannon ~  host_subfamily + temperature , data=alphaframe.noAglais)
lm.shan.fit5 <- lm(Shannon ~  host_subfamily + location, data=alphaframe.noAglais)
lm.shan.fit6 <- lm(Shannon ~  host_subfamily + temperature + location, data=alphaframe.noAglais) 
# subfam + temp better than subfam + temp + location

AIC(lm.shan.fit1, lm.shan.fit2,lm.shan.fit3, lm.shan.fit4, lm.shan.fit5,lm.shan.fit6 )
# temp better predictor than location 
# final model subfamily * temp

## Shannon complete dataset from both countries 
lm.shan.temp.nA <- lm(Shannon ~  host_subfamily * temperature , data=alphaframe.noAglais) 
summary(lm.shan.temp.nA)
summary(lm.shan.temp.nA)$adj.r.squared # 0.4061391
Anova(lm.shan.temp.nA, type=2) # independent order
vif(lm.shan.temp.nA, type = "predictor") # use predictor for interaction models
alias(lm.shan.temp.nA)
plot(lm.shan.temp.nA,2)
hist(lm.shan.temp.nA$residuals) # extract the residuals
shapiro.test(residuals(lm.shan.temp.nA)) # Run Shapiro-Wilk test on residuals 
# F-statistic: 40.17 on 1 and 157 DF,  p-value: 2.354e-09

# Partial R2 variance partitioning for lm model on Type II ANOVA
testme <- lm.shan.temp.nA
anovaII_out <- Anova(testme, type = 2)
ss_res <- sum(residuals(testme)^2)
# Partial R² per term
partial_R2 <- anovaII_out$`Sum Sq` / (anovaII_out$`Sum Sq` + ss_res)
names(partial_R2) <- rownames(anovaII_out)
partial_R2

# Partial R2 or eta squared (library effectsize)
anovaII_out <- Anova(lm.shan.temp.nA, type=2) 
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
#vif(lm.temp.host.turnover, type = "predictor") # use predictor for interaction models
alias(lm.temp.host.turnover)
plot(lm.temp.host.turnover,2)
hist(lm.temp.host.turnover$residuals) # extract the residuals
shapiro.test(residuals(lm.temp.host.turnover)) # Run Shapiro-Wilk test on residuals 
kruskal.test(temperature ~ host_subfamily, data = alphaframe.noAglais)



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

# center/scale temperature
alphaframe$temp_c <- scale(alphaframe$temperature, center = TRUE, scale = FALSE)
lm.rich.tempc <- lm(Observed ~ host_subfamily * temp_c, data = alphaframe)
vif(lm.rich.tempc, type = "predictor")
summary(lm.rich.tempc)
Anova(lm.rich.tempc, type=2)

plot(alphaframe$Observed ~ alphaframe$temperature) 
lm.rich.temp <- lm(Observed ~     host_subfamily  * temperature   , data=alphaframe) # integer time
summary(lm.rich.temp)
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
anovaII_out <- Anova(lm.rich.temp, type=2) 
eta_squared(anovaII_out, partial = TRUE)


# Interaction effect on temp and subfamily
temp.q0_subfamily <- ggplot(alphaframe, aes(x = temperature, y = Observed, color=host_subfamily)) +
  theme_line2() + 
  geom_point(size=3, alpha = 0.8, aes(shape=country)) +
  geom_smooth(method = "lm", linewidth = 1, se = F) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = T) +
  labs(x="temperature [°C]",y="microbial richness (q=0)", title="", color="host subfamily" ) +
  #facet_wrap(~country) +
  scale_colour_manual(values=subfamily_colorfix) +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )




temp.q1 <- ggplot(alphaframe.noAglais, aes(x = temperature, y = hill_numbers )) +
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

core.temp <- ggplot(alphaframe, aes(x = genus_core , y = temperature  )) +
  theme_line2() + 
  geom_point(size = 4, alpha = 0.8, aes(color=host_subfamily)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(formula = y ~ poly(x, 2), label.y = 5, aes(label = after_stat(rr.label)),size=3) +
  stat_cor(method = "spearman", label.y = 0.6, size = 3, p.accuracy = 0.001) + # Spearman
  labs(x="genus core [%]",y="temperature [°C]", title="", color="host subfamily" ) +
  scale_colour_manual(values=subfamily_colorfix)


#### alpha shannon core 
shannon.core <- ggplot(alphaframe, aes(x = Shannon, y = genus_core)) +
  theme_line2() + 
  geom_point(size = 4, alpha = 0.8, aes(color=host_subfamily)) +
  geom_smooth(method = "auto", color = "black", linewidth = 0.5, se = TRUE) +
  stat_regline_equation(label.y = 0.4, aes(label = after_stat(rr.label)),size=3) + stat_cor(label.y = 0.6,size=3, p.accuracy = 0.001) +
  labs(x="Shannon Index",y="core genus", title="", color="host subfamily" ) +
  coord_cartesian(ylim = c(0, 1), expand = TRUE,  default = FALSE,   clip = "on" ) +
  scale_colour_manual(values=subfamily_colorfix)

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
  scale_x_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
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




pdf("plots_peru/02_sample_div_alpha_Shannon_03_core.pdf", width=8, height=6)
temp.shannon.noAglais
temp.shannon.noAglais2
temp.q0 
temp.q0_subfamily 
temp.q1
temp.q2
core.temp 
shannon.core
core.shannon.auto 
core.shannon 
core.shannon.temp 
dev.off()



sink("plots_peru/02_sample_div_alpha_Shannon_02.txt")
lm.results.summary
"Shannon host subfamily"
Anova(lm.shannon.subfam, type=2)
summary(lm.shannon.subfam)$adj.r.squared 
Anova(lm.subfamily.country, type=2)
summary(lm.subfamily.country)$adj.r.squared 
Anova(lm.shannon.model, type=2)
summary(lm.shannon.model)$adj.r.squared 
Anova(lm.shannon.model2, type=2)
summary(lm.shannon.model2)$adj.r.squared

"Host taxa"
anova(lm.shannon.pierinae)
anova(lm.shannon.coliadinae)
anova(lm.shannon.dismorphinae)
anova(lm.shannon.heliconiinae)

"Peru elevation"
Anova(lm.peru, type=2) 
summary(lm.peru)$adj.r.squared 

"Germany kw"
Anova(lm.shannon.kw.host)
summary(lm.shannon.kw.host)$adj.r.squared 
Anova(lm.shannon.kw.genus)
summary(lm.shannon.kw.genus)$adj.r.squared 

"Temperature"
Anova(lm.shan.temp.nA, type=2)
summary(lm.shan.temp.nA)$adj.r.squared 
Anova(lm.rich.temp, type=2) 
summary(lm.rich.temp)$adj.r.squared 
"Shannon Core "
anova(lm.quad)
summary(lm.quad)$adj.r.squared 
sink()


## 02 sample div beta NMDS -------------------
# datasets

# ASV-level
# data.ASV         #data.ASV.rel 
# sample.ASV       #sample.ASV.rel

# genus-level
# data.species     #data.species.rel 
# sample.species   #sample.species.rel 
# sample.filter    #sample.filter.rel 
# family-level
# sample.family    #sample.family.rel

sample.ASV.rel #LT removed (samples only) ASV level
sample.species.rel#LT removed (samples only) species level

#### beta div plots (NMDS)  -------------
sample.nmds <- ordinate(sample.species.rel, method="NMDS",distance = "bray", k=3, trymax=200)
sample.nmds.ASV <- ordinate(sample.ASV.rel, method="NMDS",distance = "bray", k=3, trymax=200)
# sample.nmds.family <- ordinate(sample.family.rel, method="NMDS",distance = "bray", k=3, trymax=200)

# Check goodness of fit with Shepards diagram
#stressplot(sample.nmds)
sample.nmds.stress <- round(sample.nmds$stress, 3)
sample.nmds.stress
# dimension : 2 Stress:     0.2818483 
# dimensions: 3 Stress:     0.2156491 # better 3 dimensions
# NMDS stress level 0.05	Excellent, 0.05–0.1 Good, 0.1–0.2 Fair, 0.2–0.3 Weak, > 0.3 Very poor

# NMDS
nmds.ASV <- plot_ordination(sample.ASV.rel,sample.nmds.ASV, color="host_subfamily", shape="country")+
  geom_point(size=4) + theme_grid() + stat_ellipse(aes(group = country))  +  ggtitle("sample.nmds.ASV"  )
nmds.species <- plot_ordination(sample.species.rel,sample.nmds, color="host_subfamily", shape="country")+
  geom_point(size=4) + theme_grid() + stat_ellipse(aes(group = factor(country))) +  ggtitle("sample.nmds.species" ) +
  annotate("text", x = Inf, y = -Inf, label = paste("Stress =", sample.nmds.stress), 
           hjust = 1.1, vjust = -0.5, size = 4)

# NMDS parameter
sample_data(sample.species.rel)$weight <- as.numeric(as.character(sample_data(sample.species.rel)$weight))
sample_data(sample.species.rel)$year <- as.factor(sample_data(sample.species.rel)$year)

sample.nmds.sublocation <- plot_ordination(sample.species.rel,sample.nmds, color="host_tribe", shape="study")+
  geom_point(size=4) + theme_grid() + facet_wrap(~sublocation)

#### beta div plots (NMDS) binder 
# combine NMDS sites with data frame to include custom alpha levels in figure
scores(sample.nmds)$sites
rownames(scores(sample.nmds)$sites) == row.names(sample_data(sample.species.rel)) # check if true before cbind
nmds.binder <- cbind(scores(sample.nmds)$sites,sample_data(sample.species.rel))

nmds.binder$samplesum <- sample_sums(sample.species)
nmds.binder$logsum <- log10(sample_sums(sample.species))
nmds.binder$rich <- estimate_richness(sample.species, measures=c("Observed", "Simpson", "Shannon"))
nmds.binder$actino <- sample_sums(subset_taxa(sample.species.rel, phylum=="Actinobacteria" ))
nmds.binder$entero <- sample_sums(subset_taxa(sample.species.rel, order=="Enterobacterales" ))
nmds.binder$familycore <- sample_sums(core.family.genus)
nmds.binder$genuscore <- sample_sums(core.genus.rel)
nmds.binder$LT7000 <- ifelse(sample_sums(sample.species) > 7000, 'HT', 'LT7000')
nmds.binder$LT5000 <- ifelse(sample_sums(sample.species) > 5000, 'HT', 'LT5000')
nmds.binder$LT3000 <- ifelse(sample_sums(sample.species) > 3000, 'HT', 'LT3000')

#str(nmds.binder)

NMDS.binder.wrap.logsum <-
  ggplot(nmds.binder , aes(x=NMDS1,y=NMDS2, shape=LT5000, color=rich$Observed, size=logsum ))+
  stat_ellipse(linewidth=1, alpha=0.8, aes(group = factor(country))) +
  geom_point( alpha=0.8)+
  theme_grid()+ #labs(title="all time points") +labs(color="sampling\ntime point") +
  theme(plot.title =element_text(size=10, face='bold')) +
  scale_colour_viridis_c(direction = -1) #+ facet_wrap(~country)

NMDS.binder.shannon.logsum <- ggplot(nmds.binder, aes(x=NMDS1,y=NMDS2,shape=LT5000, color=rich$Shannon, size=logsum))+
  stat_ellipse(linewidth=1, alpha=0.8, aes(group = factor(country))) +
  geom_point( alpha=0.8)+
  theme(plot.title =element_text(size=10, face='bold')) +
  scale_colour_viridis_c(direction = 1) +
  theme_grid() #+ labs(title="sample.species") +labs(color="Shannon")
  facet_wrap(~country)

NMDS.binder.wrap.shannon <-
  ggplot(nmds.binder , aes(x=NMDS1,y=NMDS2,color=rich$Shannon, shape=LT5000, size=actino))+
  #stat_ellipse(linewidth=1, alpha=0.5 ) +
  geom_point( alpha=0.8)+
  theme(plot.title =element_text(size=10, face='bold')) +
  scale_colour_viridis_c(direction = 1) + facet_wrap(~country) +
  theme_grid() #labs(title="all time points") +labs(color="sampling\ntime point") +

NMDS.binder.familycore <- ggplot(nmds.binder, aes(x=NMDS1,y=NMDS2,shape=country, color=genuscore, size=familycore))+
  stat_ellipse(linewidth=1, alpha=0.8, aes(group = factor(host_genus))) +
  geom_point( alpha=0.8)+
  theme(plot.title =element_text(size=10, face='bold')) +
  scale_colour_viridis_c(direction = -1) +
  theme_grid() + facet_wrap(~host_genus) + geom_label(aes(label = sampleID), size = 3)


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


pdf("plots_peru/02_sample_div_beta_NMDS_binder.pdf", width=8, height=8)
nmds.ASV
nmds.species 
NMDS.binder.wrap.logsum
NMDS.binder.shannon.logsum
NMDS.binder.wrap.shannon
NMDS.binder.familycore
NMDS.binder.tribe
NMDS.binder.tribe.wrap
NMDS.binder.genus.wrap
dev.off()



#### beta div plots (PCoA)  ------------------
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
sample.PCoA.location <- plot_ordination(sample.species.rel,sample.PCoA, color="host_tribe", shape = "study")+
  geom_point(size=4)+theme_grid() + facet_wrap(~location)
sample.PCoA.kw <- plot_ordination(sample.species.rel,sample.PCoA, color="kw", shape="country")+
  geom_point(size=4)+theme_grid() + stat_ellipse(aes(group = factor(country))) + facet_wrap(~country)

PCoA.host <- plot_ordination(sample.species.rel,sample.PCoA, shape = "country")+
  geom_point(size=4, aes(color = host_species))+theme_grid() + stat_ellipse(aes(group = country)) +
  scale_color_manual(values = species_colorfix) 


# Group unique genera together
# Export dataframe and count the number of samples per host_genus
sample_metadata_df <-  as.data.frame(sample_data(sample.species.rel))
table(sample_metadata_df$host_genus)
# Group genera with less than 1 specimen
low_rep_genera <- names(which(table(sample_metadata_df$host_genus) == 1))
sample_metadata_df$host_genus_grouped <- sample_metadata_df$host_genus
sample_metadata_df$host_genus_grouped[sample_metadata_df$host_genus %in% low_rep_genera] <- "unique genera"
# Update the sample data in the phyloseq object
sample_data(sample.species.rel)$host_genus_grouped <- sample_metadata_df$host_genus_grouped


PCoA.genus.wrap.label <- plot_ordination(sample.species.rel,sample.PCoA, color="host_genus", shape = "country")+
  geom_point(size=4)+theme_grid() + stat_ellipse(aes(group = host_genus))  + 
  facet_wrap(~host_genus_grouped) +
  scale_color_manual(values = genus_colorfix) + geom_label(aes(label = sampleID), size = 4)

PCoA.genus.wrap2 <- plot_ordination(sample.species.rel,sample.PCoA, color="host_subfamily", shape = "country")+
  geom_point(size=4)+theme_grid() + stat_ellipse(aes(group = host_genus))  + 
  facet_wrap(~host_genus_grouped) +
  scale_color_manual(values = subfamily_colorfix) #+ geom_label(aes(label = sampleID), size = 4)

PCoA.subfam.wrap <- plot_ordination(sample.species.rel,sample.PCoA, color="host_genus", shape = "country")+
  geom_point(size=4)+theme_grid() + stat_ellipse(aes(group = host_subfamily))  + 
  facet_wrap(~host_subfamily) +
  scale_color_manual(values = genus_colorfix) #  + geom_label(aes(label = sampleID), size = 3)

PCoA.tribe.wrap <- plot_ordination(sample.species.rel,sample.PCoA, color="host_genus", shape = "country")+
  geom_point(size=4)+theme_grid() + stat_ellipse(aes(group = host_tribe))  + 
  facet_wrap(~host_tribe) +
  scale_color_manual(values = genus_colorfix) #  + geom_label(aes(label = sampleID), size = 3)



pdf("plots_peru/02_sample_div_beta_PCoA.pdf", width=12, height=6)
PCoA.ASV
PCoA.species
sample.PCoA.genus
sample.PCoA.subfamily
sample.PCoA.location 
sample.PCoA.kw
PCoA.host
PCoA.genus.wrap.label
PCoA.genus.wrap2
PCoA.subfam.wrap
PCoA.tribe.wrap 
dev.off()

#### beta div plots (PCoA) binder

pcoa_scores <- sample.PCoA$vectors
pcoa_scores_df <- as.data.frame(pcoa_scores)
pcoa_scores_3 <- pcoa_scores_df[, 1:3]
pcoaframe.df <- as(sample_data(sample.species.rel),"data.frame")

pcoa.binder <- cbind(pcoa_scores_3,sample_data(sample.species.rel))
pcoa.binder$genus_core <- core.genus.rel.df$genus_core

axes_var <- sample.PCoA$values$Relative_eig * 100 # Export axis % variation 
xlab_pcoa1 <- paste0("PCoA1 (", round(axes_var[1], 1), "%)")
ylab_pcoa2 <- paste0("PCoA2 (", round(axes_var[2], 1), "%)")
ylab_pcoa3 <- paste0("PCoA3 (", round(axes_var[3], 1), "%)")

PCoA.host.overview  <- ggplot(pcoa.binder, aes(Axis.1, Axis.2)) +
  stat_ellipse(level = 0.9, linewidth = 0.8,  aes(group = host_subfamily, color = host_subfamily)) + #linetype = 5,
  geom_point(size = 4, aes(shape = country, color = host_subfamily), alpha = 0.8) +
  facet_wrap(~host_family) +
  scale_colour_manual(values = subfamily_colorfix) +
  labs(x = xlab_pcoa1, y = ylab_pcoa2, color = "host subfamily") +
  theme_grid() +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2) )

PCoA.axis13  <- ggplot(pcoa.binder, aes(Axis.1, Axis.3)) +
  stat_ellipse(level = 0.9, linewidth = 1, aes(group = host_subfamily, color = host_subfamily)) +
  geom_point(size = 5, aes(shape = country, color = host_subfamily), alpha = 0.8) +
  facet_wrap(~host_family) +
  scale_colour_manual(values = subfamily_colorfix) +
  labs(x = xlab_pcoa1, y = ylab_pcoa3) +
  theme_grid()

PCoA.axis23  <- ggplot(pcoa.binder, aes(Axis.2, Axis.3)) +
  stat_ellipse(level = 0.9, linewidth = 1, aes(group = host_subfamily, color = host_subfamily)) +
  geom_point(size = 5, aes(shape = country, color = host_subfamily), alpha = 0.8) +
  facet_wrap(~host_family) +
  scale_colour_manual(values = subfamily_colorfix) +
  labs(x = ylab_pcoa2, y = ylab_pcoa3) +
  theme_grid()

PCoA.host.overview.tribe  <- ggplot(pcoa.binder, aes(Axis.1, Axis.2)) +
  stat_ellipse(type = "t",level = 0.9, linewidth = 1, aes(group = host_tribe, color = host_tribe)) +
  geom_point(size = 4, aes(shape = country, color = host_tribe), alpha = 0.8) +
  facet_wrap(~host_family) +
  scale_colour_manual(values = tribe_colorfix) +
  labs(x = xlab_pcoa1, y = ylab_pcoa2) +
  theme_grid()

PCoA.subfamily.wrap  <- ggplot(pcoa.binder, aes(Axis.1, Axis.2)) +
  stat_ellipse(level = 0.9, linewidth = 0.7, aes(group = host_subfamily, color = host_tribe)) +
  geom_point(size = 5, aes(shape = country, color = host_tribe), alpha = 0.8) +
  facet_wrap(~host_subfamily) +
  scale_colour_manual(values = tribe_colorfix) +
  labs(x = xlab_pcoa1, y = ylab_pcoa2) +
  theme_grid()


PCoA.genus.wrap.tribe <- ggplot(pcoa.binder, aes(Axis.1, Axis.2), color = host_tribe)+
  geom_point(size=4, alpha = 0.8, aes(shape = country, color = host_tribe)) +
  stat_ellipse(aes(group = host_genus, color = host_tribe))  + 
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
  
pdf("plots_peru/02_sample_div_beta_PCoA_binder.pdf", width=12, height=6)
PCoA.host.overview
PCoA.axis13
PCoA.axis23
PCoA.host.overview.tribe
PCoA.subfamily.wrap
PCoA.genus.wrap.tribe
PCoA.host.country.wrap 
PCoA.host.country 
pcoa.binder.core 
dev.off()


#### betadisper betaframe------------
sample.vegdist <- vegan::vegdist(t(otu_table(sample.species.rel)), index="bray")
#sample.distance <- phyloseq::distance(otu_table(sample.species.rel), method = "bray") # identical outcome
all(rownames(nmds.binder) == names(betadisper(sample.vegdist, group = nmds.binder$host_subfamily)$distances)) # check if true


### New betaframe 
betaframe <- nmds.binder

# Put distances in data frame for plotting and linear model
betaframe$beta_subfamily <- betadisper(sample.vegdist, group=betaframe$host_subfamily)$distances
betaframe$beta_family <- betadisper(sample.vegdist, group=betaframe$host_family)$distances
betaframe$beta_genus <- betadisper(sample.vegdist, group=betaframe$host_genus)$distances
betaframe$beta_tribe <- betadisper(sample.vegdist, group=betaframe$host_tribe)$distances
betaframe$beta_country <- betadisper(sample.vegdist, group=betaframe$country)$distances
betaframe$beta_location <- betadisper(sample.vegdist, group=betaframe$location)$distances
betaframe$beta_sublocation <- betadisper(sample.vegdist, group=betaframe$sublocation)$distances
betaframe$host_subfamily_ordered <- factor(betaframe$host_subfamily, levels = subfamily_ordered)
betaframe$temperature <- betaframe$temp_week


# Simple base R plot
plot(betadisper(sample.vegdist, group = betaframe$host_subfamily))
plot(betadisper(sample.vegdist, group = betaframe$location))
plot(betadisper(sample.vegdist, group = betaframe$country))

betadisp.family <- ggplot(betaframe, aes(x=host_family, y=beta_family))+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  scale_colour_manual(values=subfamily_colorfix) +  theme_line2() +
  geom_smooth(method="lm", color="black", linewidth=0.5) +labs(y="beta dispersion", x='') #+

# Check normal distribution of data
hist(betaframe$beta_family) 
permutest(betadisper(sample.vegdist, betaframe$host_family), permutations = 999) 
# 5.8135    999  0.018 *

betadisp.subfamily.country <- ggplot(betaframe, aes(x=host_subfamily_ordered, y=beta_subfamily))+
  geom_violin( trim = T, adjust = 0.8, scale = "width" ) + # draw_quantiles = c(0.5),
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_tribe, shape=country),alpha = 0.8) +
  scale_colour_manual(values=tribe_colorfix) +  
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y="distance to centroid", x='', color ="host tribe") + #subtitle="beta dispersion", 
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

# Check normal distribution of data
hist(betaframe$beta_subfamily) 
# Tests if groups have different beta dispersion
permutest(betadisper(sample.vegdist, betaframe$host_subfamily), permutations = 9999) 
# Host subfamilies differ in beta dispersion, indicating variation in within-group community heterogeneity
# 8.2821   9999  1e-04 ***
 
# How to test betadisper group differences
anova(betadisper(sample.vegdist, group = betaframe$host_subfamily)) # permutest more robust, no normality assumption
TukeyHSD(betadisper(sample.vegdist, group = betaframe$host_subfamily))

# sort by subfamily
betaframe.sorted  <- betaframe  %>%
  group_by(host_subfamily) %>%  # Group by host_genus and country
  mutate(mean_beta = median(beta_subfamily, na.rm = TRUE)) %>%  # Calculate the mean of beta
  ungroup()

betadisp.subfamily.sorted <- ggplot(betaframe.sorted, aes(x=fct_reorder(host_subfamily, mean_beta), y=beta_subfamily))+
  #geom_boxplot() +
  #geom_violin(#draw_quantiles = c(0.5), trim = F,  scale = "width"
  #            ) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_tribe, shape=country),alpha = 0.8) +
  scale_colour_manual(values=tribe_colorfix) +  
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(subtitle="beta dispersion", y="distance to centroid", x='', color ="host tribe") +
  guides(
    color = guide_legend(order = 1), # manual legend order
    shape = guide_legend(order = 2)) 
  

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
  
betadisp.tribe <- ggplot(betaframe, aes(x=host_tribe, y=beta_tribe))+
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_tribe),alpha = 0.8) +
  scale_colour_manual(values=tribe_colorfix) +  theme_line2() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y="beta dispersion", x='')

# Check normal distribution of data
hist(betaframe$beta_tribe) 
permutest(betadisper(sample.vegdist, betaframe$host_tribe), permutations = 9999) 
# 15.509   9999  1e-04 ***

# sort by tribe
nmds.tribe.sorted  <- betaframe  %>%
  group_by(host_tribe) %>%  # Group by host_genus and country
  mutate(mean_beta = mean(beta_tribe, na.rm = TRUE)) %>%  # Calculate the mean of beta
  ungroup()

betadisp.tribe.sorted <- ggplot(nmds.tribe.sorted, aes(x=fct_reorder(host_tribe, mean_beta), y=beta_tribe))+
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  scale_colour_manual(values=subfamily_colorfix) +  theme_line2() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y="beta dispersion", x='')

betadisp.country <- ggplot(betaframe, aes(x=country, y=beta_country))+
  geom_boxplot(aes(group = country))+
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_subfamily),alpha = 0.8) +
  scale_colour_manual(values=subfamily_colorfix) +  theme_line2() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  #geom_smooth(method="lm", color="black", linewidth=0.5) +
  labs(y="beta distance", x='')

hist(betaframe$beta_country) 
permutest(betadisper(sample.vegdist, betaframe$country), permutations = 9999) 
#  1.2637   9999 0.2605 No variation in within-group community heterogeneity by country 
# Says nothing about group centroids (i.e., the actual composition of communities).


hist(betaframe$beta_location) 
permutest(betadisper(sample.vegdist, betaframe$location), permutations = 9999)
#  1.7072   9999 0.1674 No variation in within-group community heterogeneity by location 

location.select.df  <- betaframe  %>%
  group_by(location, country) %>%  # Group by location and country
  mutate(mean_beta = mean(beta_location, na.rm = TRUE)) %>%  # Calculate the mean of beta
  ungroup()

betadisp.location <- ggplot(location.select.df , aes(x=location, y=beta_location))+
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_tribe),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  facet_wrap(~ country, scales = "free_x") +
  scale_colour_manual(values=tribe_colorfix)


### Beta Dispersion genus 

# Simplify dataset / remove host_genus with only 1 count (as beta dispersion would be zero)
sample.species.rel.df <- as.data.frame(sample_data(sample.species.rel))

host_genus_counts <- sample.species.rel.df %>%
  group_by(host_genus) %>%
  summarise(n_samples = n())

# Filter out host genera that occur more than 2 times
genus_to_keep <- host_genus_counts %>%
  filter(n_samples >= 2) %>%
  pull(host_genus)

# Subset the phyloseq object to keep only the selected host genera
 genus.robust <- subset_samples(sample.species.rel, host_genus %in% genus_to_keep)


genus.vegdist <- vegan::vegdist(t(otu_table(genus.robust)), index="bray")
#genus.distance <- phyloseq::distance(otu_table(genus.robust ), method = "bray") # identical outcome
genus.robust.df  <- data.frame(sample_data(genus.robust))

# Tests if groups have different beta dispersion
permutest(betadisper(genus.vegdist, genus.robust.df$host_genus), permutations = 9999) 
# 6.3486   9999  1e-04 ***

# Put distances in data frame for plotting and linear model
genus.robust.df$beta <- betadisper(genus.vegdist, group=genus.robust.df$host_genus)$distances

genus.robust.mean  <- genus.robust.df  %>%
  group_by(host_genus, country) %>%  # Group by host_genus and country
  mutate(mean_beta = mean(beta, na.rm = TRUE)) %>%  # Calculate the mean of beta
  ungroup()

betadisp.genus   <- ggplot(genus.robust.mean , aes(x=fct_reorder(host_genus, mean_beta), y=beta))+
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_tribe, shape=country),alpha = 0.8) +
  scale_colour_manual(values=tribe_colorfix) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(subtitle="Within-Group Variation per Host Genus", y="distance to centroid", x='', color ="host tribe") +
  guides(
  color = guide_legend(order = 1), # manual legend order
  shape = guide_legend(order = 2) )


beta.lm.select <- lm(beta ~ host_genus, data=genus.robust.df)
summary(beta.lm.select)
# F-statistic: 6.872 on 16 and 157 DF,  p-value: 9.45e-12
Anova(beta.lm.select )

### betadisper sublocation 
#Beta Dispersion
# Simplify dataset / remove sublocation with only 1 count (as beta dispersion would be zero)
sublocation_counts <- sample.species.rel.df  %>% group_by(sublocation) %>%  summarise(n = n())

# Filter out host genera that occur more than 2 times
sublocation_to_keep <- sublocation_counts %>% filter(n > 1) %>% pull(sublocation)

# Subset the phyloseq object to keep only the selected host genera
sublocation.select.rel <- subset_samples(sample.species.rel, sublocation %in% sublocation_to_keep)

sublocation.vegdist <- vegan::vegdist(t(otu_table(sublocation.select.rel )), index="bray")
sublocation.select.df  <- as.data.frame(sample_data(sublocation.select.rel ))

# Tests if groups have different beta dispersion
permutest(betadisper(sublocation.vegdist, sublocation.select.df$sublocation), permutations = 999) 
# 2.7897   9999 0.0036 **


# Put distances in data frame for plotting and linear model
sublocation.select.df$beta <- betadisper(sublocation.vegdist, group=sublocation.select.df$sublocation)$distances
sublocation.select.df  <- sublocation.select.df  %>%
  group_by(sublocation, country) %>%  # Group by sublocation and country
  mutate(mean_beta = mean(beta, na.rm = TRUE)) %>%  # Calculate the mean of beta
  ungroup()

betadisp.sublocation  <- ggplot(sublocation.select.df , aes(x=fct_reorder(sublocation, mean_beta), y=beta))+
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0), size=4, aes(color = host_tribe),alpha = 0.8) +
  theme_line2()+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_colour_manual(values=tribe_colorfix) +
  facet_wrap(~ country, scales = "free_x") 


pdf("plots_peru/02_sample_div_beta_plot_betadisper.pdf", width=12, height=6)
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


#### beta bioenv envfit ---------------
# bioenv() is designed for numeric environmental variables only,  
# it calculates Euclidean distances and correlates them with your community dissimilarities.
numericvariables = which(sapply(sample_data(sample.species.rel), is.numeric))
numeric.df = data.frame(sample_data(sample.species.rel))[numericvariables]
bioenv(veganotu(sample.species.rel) ~  weight + kw + month + day + month_day + elevation + temp_week ,  numeric.df)
# elevation temp_week # with correlation  0.2380641  

# For plotting select only numerical variables for vector 
env_numerical<- c("elevation", "temp_week")
fit_to_nmds <- envfit(sample.nmds, numeric.df[, env_numerical], permutations = 999, na.rm = TRUE)
fit_to_nmds

# plot as vector if significant 
plot(sample.nmds, type = "n") 
points(sample.nmds, display = "sites", pch = 21, col = "black", bg = "lightblue")
plot(fit_to_nmds, p.max = 0.05, col = "red") # only significant arrows


#### Fit all environmental data with ordination (envfit for numeric and categorical)
all_metadata <- data.frame(sample_data(sample.species.rel))
#colnames(all_metadata) # select metadata of interest
env_metadata <- all_metadata[, c("host_subfamily", "host_tribe", "host_genus", "host_species", "sex", "location", "sublocation", "year", "kw",  "elevation", "country" , "temp_week" ,"storage")]
colSums(is.na(env_metadata)) # check which data has missing values (NAs)


# Fit environmental factors to ordination NMDS  (envfit)
fit_nmds <- envfit(sample.nmds, env_metadata, permutations = 999, na.rm = TRUE)
print(fit_nmds)
# does not allow to control for interaction or confounder effects like adonis2 but allows to show best R2 value
# ***VECTORS
# NMDS1    NMDS2     r2 Pr(>r)    
# elevation -0.23804 -0.97126 0.1658 0.0001 ***
# temp_week -0.74684 -0.66500 0.2372 0.0001 ***

# Goodness of fit:
#r2 Pr(>r)    
#host_subfamily 0.3035  1e-04 ***
#  host_tribe     0.4027  1e-04 ***
#  host_genus     0.6229  1e-04 ***
#  host_species   0.6545  1e-04 ***
#  sex            0.2116  1e-04 ***
#  location       0.2525  1e-04 ***
#  sublocation    0.3180  1e-04 ***
#  year           0.1845  1e-04 ***
#  country        0.2399  1e-04 ***


# Fit environmental factors to ordination PCoA with capscale
sample.capscale <- capscale(t(otu_table(sample.species.rel)) ~ 1, distance = "bray")
fit_PCoA <- envfit(sample.capscale, env_metadata, permutations = 999, na.rm = TRUE)
print(fit_PCoA )

#***VECTORS
#              MDS1     MDS2     r2 Pr(>r)    
#kw        -0.31061  0.95054 0.0265  0.120    
#elevation  0.02754 -0.99962 0.1368  0.001 ***
#temp_week  0.81164 -0.58415 0.1839  0.001 ***

#Goodness of fit:
#                   r2 Pr(>r)    
#host_subfamily 0.3626  0.001 ***
#host_tribe     0.4374  0.001 ***
#host_genus     0.6076  0.001 ***
#host_species   0.6421  0.001 ***
#sex            0.2155  0.001 ***
#location       0.2453  0.001 ***
#sublocation    0.3241  0.001 ***
#year           0.2632  0.001 ***
#country        0.2319  0.001 ***
#storage        0.2681  0.001 ***




#### PERMANOVA adonis2 --------------------------------
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

adonis.country  <- adonis2(phyloseq::distance(sample.species.rel, method = "bray")~ country     , data=betaframe,  permutations = 9999, by = "margin")
adonis.country
# country    1    3.592 0.05773 10.11  1e-04 *** # alone 6%

adonis.location  <- adonis2(phyloseq::distance(sample.species.rel, method = "bray")~ location     , data=betaframe,  permutations = 9999, by = "margin")
adonis.location
# location   3    4.809 0.07729 4.5511  1e-04 *** # alone 8% 

adonis.sublocation  <- adonis2(phyloseq::distance(sample.species.rel, method = "bray")~ sublocation     , data=betaframe,  permutations = 999, by = "margin")
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
# Adding environmental factors to species model increases R2 only slightly 


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


# check for colinearity
genus.robust.df$dummy <- seq_len(nrow(genus.robust.df)) # host_tribe + location  + elevation  + temp_week + kw
lm_dummy <- lm(dummy ~  host_tribe +  location + elevation +  temp_week + kw    , data = genus.robust.df )
#summary(lm_dummy)
vif(lm_dummy)   # Variance inflation factor VIF > 5–10 strong collinearity, consider removing variables
vif(lm_dummy, type = "predictor") # GVIF^(1/2Df) adjusted GVIF for categorical variables with multiple levels
alias(lm_dummy)
#                  GVIF Df GVIF^(1/(2*Df))
#host_tribe  74.231157 10        1.240307
#location   102.654525  3        2.163863
#elevation    8.472244  1        2.910712
#temp_week    9.849418  1        3.138378 max VIF 3.1 is still acceptable
#kw           2.028766  1        1.424347


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
# Extract marginal (partial) R2 values
marginal_r2 <- as.data.frame(adonis.robust.tribe)
r2_df_clean <- marginal_r2 |>
  subset(!(rownames(marginal_r2) %in% c("Residual", "Total")))

r2_df_clean <- data.frame(
  term = rownames(r2_df_clean),
  R2   = r2_df_clean$R2,
  p    = r2_df_clean$`Pr(>F)`)

r2_df_clean$term <- factor(r2_df_clean$term,
                           levels = c("host_tribe", "location", "kw", "elevation", "temp_week"),
                           labels = c("host tribe", "location", "calendar week", "elevation", "temperature"))

# Sort by explained variance
r2_df_clean$term <- factor(r2_df_clean$term,
  levels = r2_df_clean$term[order(r2_df_clean$R2, decreasing = TRUE)])

r2_df_clean$signif <- cut(r2_df_clean$p,
                          breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                          labels = c("***","**","*",""))

# Show effect size as marginal R2
marginalR2.plot <- ggplot(r2_df_clean, aes(x = term, y = R2)) +
  geom_col(fill = "grey70") +
  geom_text(aes(label = signif), 
            vjust = -0.5, # place above the bar
            size = 6) +    # adjust size
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.1))) + # add space above bars
  labs(
    x = NULL,
    y = "partial R² (marginal)"
  ) +
  theme_line2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  )


#### Variance partitioning beta diversity 
# decomposes shared vs unique fractions visually

community.matrix <- as(otu_table(genus.robust), "matrix")
if (taxa_are_rows(genus.robust)) {
  community.matrix <- t(community.matrix) }

# Hellinger transform
comm_hel <- decostand(community.matrix, method = "hellinger")
stopifnot(identical(rownames(comm_hel), rownames(genus.robust.df)))
genus.robust.df <- genus.robust.df[rownames(comm_hel), , drop = FALSE]

# check for colinearity
genus.robust.df$dummy <- seq_len(nrow(genus.robust.df))
lm_dummy <- lm(dummy ~  host_tribe +   temp_week + kw  + elevation + location  , data = genus.robust.df)
vif(lm_dummy)   # Variance inflation factor VIF > 5–10 strong collinearity, consider removing variables
vif(lm_dummy, type = "predictor") # GVIF^(1/2Df) adjusted GVIF for categorical variables with multiple levels
alias(lm_dummy)
# working combinations without colinearity:
# species + temp, kw,
# genus   + temp, kw, elev
# tribe   + temp, kw, elev, location
# subfam  + temp, kw, elev, location

# quick check
cor(genus.robust.df [, c("temp_week", "kw", "elevation")]) # quick check numbers
cor_matrix <- cor(genus.robust.df[, c("temp_week", "kw", "elevation")])
r2_matrix <- cor_matrix^2
r2_matrix


X_host     <- genus.robust.df[, "host_tribe", drop = FALSE]
X_location <- genus.robust.df[, "location", drop = FALSE]
#X_env  <- genus.robust.df[, c("temp_week"), drop = FALSE]
#X_env  <- genus.robust.df[, c("temp_week", "kw"), drop = FALSE]
#X_env  <- genus.robust.df[, c("temp_week", "kw", "elevation"), drop = FALSE]
X_env  <- genus.robust.df[, c("temp_week", "kw",  "location"), drop = FALSE]
#X_env  <- genus.robust.df[, c("temp_week", "kw",  "location" , "elevation"), drop = FALSE]

# List of predictor sets
predictors <- list(
  host = X_host,
  location = X_location,
  env = X_env
)

# Function to count unique values per column, works for factor, character, numeric
sapply(predictors, function(df) {
  sapply(df, function(col) length(unique(col)))
})

# varpart
vp1 <- varpart(
  comm_hel,
  X_host,
  X_env )
# warns if colinearity detected
vp1

#plot(vp1) # gives a Venn diagram-style plot

plot(vp1,
     digits = 1,
     bg = c("skyblue", "salmon"),
     Xnames = c("Host", "Environment"))

# subfamily 10% both5% env 4% 
# tribe 15%  both 5% env 4% tribe vs temp location = okay 
# genus 16% both 10 env 1%; genus 19% both 7% 1% elevation temp kw = okay
# species 19% both 10% 1%  

# only temp
# subfam 12% 4% 0%
# tribe 17% 4% 1%
# genus 22% 4% 0%
# species 24% 5% 1% temp kw

# varpart
vp2 <- varpart(
  comm_hel,
  X_host,
  X_location )

vp2

plot(vp2) # gives a Venn diagram-style plot

# check variance inflation, but less reliable than adjusted GVIF^(1/(2*Df))
vif.cca(rda(comm_hel ~ ., data = cbind(X_host, X_env)))

# Testing varpart fractions with condition approach
# Test host genus fraction controlling for environment
rda_host <- rda(comm_hel ~ host_tribe + Condition(kw + temp_week + location), 
                data = genus.robust.df)
anova(rda_host, permutations = 999)
RsquareAdj(rda_host)

# Test environment fraction controlling for host
rda_env <- rda(comm_hel ~ kw + temp_week + location +  Condition(host_tribe), 
               data = genus.robust.df)
anova(rda_env, permutations = 999)
RsquareAdj(rda_env)




sink("plots_peru/02_sample_div_beta_stats_varpart_ADONIS.txt")
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





## 03 taxa abundance core / endo / family  --------------

#### core abundance  
core.melt <- psmelt(core.genus.rel)

# order by subfamily_ordered
core.melt$host_subfamily_ordered <- factor(core.melt$host_subfamily, levels = rev(c(subfamily_ordered)))

# manually rebuild stripchart
core_stripchart.rebuild <- ggplot(core.melt, aes(x = host_subfamily_ordered, y = Abundance)) +
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
with(subset(core.melt, genus == "Acinetobacter"), 
     pairwise.wilcox.test(Abundance, host_subfamily, p.adjust.method = "BH"))

with(subset(core.melt, genus == "Apibacter"), 
     pairwise.wilcox.test(Abundance, host_subfamily, p.adjust.method = "BH", exact = FALSE))

with(subset(core.melt, genus == "Asaia"), 
     pairwise.wilcox.test(Abundance, host_subfamily, p.adjust.method = "BH"))

with(subset(core.melt, genus == "Bartonella"), 
     pairwise.wilcox.test(Abundance, host_subfamily, p.adjust.method = "BH", exact = FALSE))

with(subset(core.melt, genus == "Enterococcus"), 
     pairwise.wilcox.test(Abundance, host_subfamily, p.adjust.method = "BH", exact = FALSE))


# Run kruskal.test function per country 
core_stats_country <- run_kruskal(core.melt, "country")
core_stats_country

sig_genera <- core_stats_country %>% filter(p_adj < 0.05)
sig_genera

# Make label for figure
core_country_label <- core_stats_country %>%
  mutate(label = paste0("p = ", p_adj, " ", signif))

core_stripchart.country.stats <- ggplot(core.melt, aes(x = country, y = Abundance)) +
  geom_jitter(aes(color = country),width = 0.2, size = 3, alpha = 0.6) +  # points on top
  #geom_boxplot(outlier.shape = NA) +  # boxplot alpha = 0.8
  geom_violin() +
  stat_summary(fun = mean,  geom = "crossbar", width = 0.2,             
               color = "black",  linewidth = 0.5  ) +
  facet_wrap(~genus) + # , scales = "free"
  theme_grid() + theme(axis.text.x = element_text(size = 8)) +
  scale_colour_manual(values = country_colorfix) +
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +
  labs(x = "", y = "rel ab [%]", color = "country") +
  geom_text(
    data = core_country_label,
    aes(x = 1.5, # x-position for label 1.5 centered
        y = 0.95, # y-position above data points
        label = label
    ),
    inherit.aes = FALSE,
    size = 3  ) 

# manually rebuild stripchart
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


#### taxa_abundance ASV core -------------------

sample.ASV.core <- subset_taxa(sample.ASV.rel, genus %in% core.genus)
data.frame(sort(taxa_sums(sample.ASV.core),decreasing = TRUE))
#sample.ASV.core.top <- prune_taxa(taxa_sums(sample.ASV.core) > 0.04, sample.ASV.core) # cumulative abundance 4% (not clear)
#sample.ASV.core.top <- prune_taxa(taxa_sums(sample.ASV.core) / nsamples(sample.ASV.core) >= 0.01,  sample.ASV.core) # ≥ 1% relative abundance in entire dataset
#sample.ASV.core.top <- prune_taxa(taxa_sums(sample.ASV.core) / nsamples(sample.ASV.core) >= 0.00024,  sample.ASV.core) # same as cum 4% or 0.024% relative abundance in entire dataset
sample.ASV.core.top <- prune_taxa(apply(otu_table(sample.ASV.core), 1, max) >= 0.02, sample.ASV.core) # ≥ 2% relative abundance in at least one sample # better!

sample.ASV.core.top.melt <- psmelt(sample.ASV.core.top)

ASV.bubble.country <- ggplot(sample.ASV.core.top.melt,aes(fct_reorder(sampleID,collector),fct_reorder(OTU, genus, .desc=TRUE))) +
  geom_point(aes(size=Abundance, color= genus), alpha=0.8) +  
  theme_line2() +  
  theme(legend.position="right",legend.box="vertical", legend.text = element_text(size=9), legend.title= element_text(size=9), axis.text.y= element_text(size=9, face="italic"), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  facet_grid(~country, scales = "free", space = "free") +
  scale_color_manual(values=core_palette) +
  scale_size(name = "Rel. Abundance [%]") +
  guides(color = guide_legend(override.aes = list(size = 5))) + # bigger dots in legend 
  xlab("sample ID") + ylab("") 

ASV.bubble.core.country <- ggplot(sample.ASV.core.top.melt,aes(fct_reorder(sampleID,collector),fct_reorder(OTU, genus, .desc=TRUE))) +
  geom_point(aes(size=Abundance, color= genus), alpha=0.8) + #fill= genus), shape=21,color="black") + 
  theme_line2() +  
  facet_wrap(~genus+country, scales = "free", ncol = 4) +
  scale_color_manual(values=core_palette) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank() ) + 
  guides(color = guide_legend(override.aes = list(size = 5))) + # bigger dots in legend 
  labs(x="",y="", fill="genus core", size="rel. ab. [%]") 

ASV.bubble.core.subfamily <- ggplot(sample.ASV.core.top.melt,aes(fct_reorder(sampleID,collector),fct_reorder(OTU, genus, .desc=TRUE))) +
  geom_point(aes(size=Abundance, fill= genus), shape=21,color="black") + 
  theme_line2() +  
  facet_grid(~host_subfamily, scales = "free") +
  scale_size(name = "Rel. Abundance [%]") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank() ) + xlab("") + ylab("") +
  scale_fill_manual(values=core_palette) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # bigger dots in legend 


pdf("plots_peru/03_sample_taxa_abundance_core_ASV.pdf", width=8, height=6)
core_stripchart.rebuild
core_stripchart.country.stats
core_stripchart.rebuild.tribe
ASV.bubble.country
ASV.bubble.core.country
ASV.bubble.core.subfamily
dev.off()

sink("plots_peru/03_sample_taxa_abundance_core_ASV.txt")
"Kruskal.test core genus per subfamily"
data.frame(core_stats_subfam)
"Kruskal.test core genus per country"
data.frame(core_stats_country)
sink()


#### taxa_abundance Endosymbionts  ----------------

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

# Remove low abundances, set everything below 1% to zero
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
  theme_line2() + #theme(axis.text.x = element_text(angle = 60, hjust = 1, face="italic") ) + # x-axis angle
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + # y-axis 100%
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
  #facet_wrap(~country, scale = "free_y") +
  guides(fill = guide_legend(reverse = TRUE)) # reverse legend order

pdf("plots_peru/03_sample_taxa_abundance_endosymbiont.pdf", width=8, height=6)
endosymbiont.rel_stripchart
endosymbiont.rel_stripchart1
endo.bar2
endo.bar3
dev.off()


## 05 ggtree ASV alignment all taxa -----------------
# use reduced dataset sample.filter
##data.frame(sort(taxa_sums(sample.filter), decreasing = F)[1:100])

# percent of total reads remaining 
percent_retained <- (sum(otu_table(sample.filter)) / sum(otu_table(sample.species))) * 100
percent_retained # 98.0337% of total dataset

# Check if core taxa remain in sample.filter
sample.filter.rel <- transform_sample_counts(sample.filter, function(x) x/sum(x))
sample.filter.rel.g <- aggregate_taxa(sample.filter.rel, "genus")
data.frame(sort(taxa_sums(sample.filter.rel.g)))

filter.core.genus <- core_members(sample.filter.rel.g , detection = 0.01, prevalence = 20/100)
filter.core.genus.high <- core_members(sample.filter.rel.g, detection = 0.05, prevalence = 10/100) 
filter.core.genus.low <- core_members(sample.filter.rel.g, detection = 0.002, prevalence = 40/100) 

# sample.filter check taxa
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

# Collapse OTU table by genus and select most abundant ASV from ps object
sample.filter.genus <- tax_glom(sample.filter, taxrank = "genus", NArm=F) # glom by genus but keeps ASV name info
data.frame(sort(taxa_sums(sample.filter.genus)))
asv_names <- rownames(otu_table(sample.filter.genus))

# compare ASV name list both methods give same outcome
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

# Export ASV sequences
writeXStringSet(matched_all_seqs , file = "tree/ASV_export.fasta")
# Use to make tree SINA 1.2.12
# https://www.arb-silva.de/aligner/
# SSU, min identity with query 0.90, Number of neighbours 1
# tree building RAxML with GTR model and Gamma likelihood
# rename tree file in nwk and reimport into R


# Find phylum ASV name for Root tree
check_taxa <- subset_taxa(sample.filter.seq , phylum == "Bacteroidetes") # Fusobacteria Deinococcus-Thermus Verrucomicrobia
sort(data.frame(sum = taxa_sums(check_taxa)), decreasing = FALSE)
data.frame(
  genus = as.character(tax_table(check_taxa)[names(sort(taxa_sums(check_taxa), decreasing = T)), "genus"]),
  Taxa_Sum = sort(taxa_sums(check_taxa), decreasing = T))
OG <- "ASV99" # define outgroup ASV3819 ASV937

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

#SINA_tree2 <- read.tree("tree/arb-silva2.nwk")
#SINA_tree2.root <- root(SINA_tree2, outgroup = OG, resolve.root = TRUE) 
#sample.filter.SINA_tree2 <- merge_phyloseq(sample.filter.seq, SINA_tree2.root )
#SINA.tree2 <-  ggtree(sample.filter.SINA_tree2, layout="circular") + geom_tiplab(size = 3, align=TRUE, aes(label=genus, color=phylum)) + ggtitle("SINA RAxML tree")
#tree.root  <- drop.tip(SINA_tree2.root, setdiff(SINA_tree2.root$tip.label, asv_names))

# Choose best tree clean up tip labels
tree.root  <- drop.tip(SINA_tree.root, setdiff(SINA_tree.root$tip.label, asv_names))
sample.filter.tree <- merge_phyloseq(sample.filter.seq, tree.root )


# circular tree all genera
tree.all.circle <- ggtree(sample.filter.tree, layout="circular") + 
  geom_tiplab(size = 3, align=TRUE, aes(label=genus, color=phylum)) 

### 05 ggtree 02 core info 
tip_data <- data.frame(label = tree.root$tip.label)
tip_data$core <- ifelse(tip_data$label %in% core.genus.ASV.names, "Core", "Non-Core")
core_tips <- tip_data[tip_data$core == "Core", ]

# Ad phylum genus names for core tips
taxa_data <- as.data.frame(tax_table(sample.filter.tree))
core_tips$genus_name <- taxa_data$genus[match(core_tips$label, rownames(taxa_data))]
core_tips$phylum_name <- taxa_data$phylum[match(core_tips$label, rownames(taxa_data))]
core_tips

tree.core <- ggtree(tree.root, layout="circular") + 
  geom_tiplab(size = 3, align = TRUE, 
              aes(label = ifelse(label %in% core_tips$label, 
                                 core_tips$genus_name[match(label, core_tips$label)], ""))) # Only show core genera
  
# Select bee associated taxa 
select_bee_genus_list <- c("Gilliamella", "Snodgrassella", "Lactobacillus", "Bifidobacterium", "Frischella", "Bartonella", "Apibacter", "Apilactobacillus", "Commensalibacter")
bee_taxa <- taxa_data[rownames(taxa_data) %in% select_bee_genus_list | taxa_data[,"genus"] %in% select_bee_genus_list , ]

bee_tips <- data.frame(
  genus = bee_taxa[,"genus"],  # Extract genus names
  label = rownames(bee_taxa)  # Extract ASV/tip labels
)

tree.bee <- ggtree(sample.filter.tree, layout="circular") +
  geom_tree() +
  geom_tiplab(size = 3, align = TRUE, aes(label = ifelse(label %in% bee_tips$label, genus, ""))) +  # Only show core genera
  theme(legend.position = "right") 

tree.bee.asteriks <- ggtree(sample.filter.tree, layout="circular") +
  geom_tree() +
  geom_tiplab(size = 3, align = TRUE, aes(label = ifelse(label %in% bee_tips$label, "*", "")), color = "red") +  # Show * for matching tips
  theme(legend.position = "right")


pdf("plots_peru/05_ggtree_02_tree_core.pdf", width=12, height=6)
IQ.tree
SINA.tree
tree.all.circle
tree.core
tree.bee
tree.bee.asteriks
dev.off()


####  ggtree rounds tree  ---------------
taxonomy_df <- as.data.frame(tax_table(sample.filter.tree))
tree_tips <- phy_tree(sample.filter.tree)$tip.label
all(tree_tips %in% rownames(taxonomy_df)) # Check if all okay


taxonomy_df <- taxonomy_df[tree_tips, , drop = FALSE]
family_df <- data.frame(family = taxonomy_df$family)
rownames(family_df) <- rownames(taxonomy_df)
taxonomy_df$label <- rownames(taxonomy_df)
taxonomy_df$dummy <- seq_len(nrow(taxonomy_df))

# Adding total, relative and log abundance
otu_data <- otu_table(sample.filter.tree)
total_abundance <- rowSums(otu_data)
total_ab_df <- data.frame(total_ab = total_abundance)
taxonomy_df$total_ab <- total_ab_df[rownames(taxonomy_df), "total_ab"]

# ad rel ab (rough calculation without sample normalization)
relative_abundance <- total_abundance * 100/ (sum(total_abundance))
taxonomy_df$rel_ab <- relative_abundance[rownames(taxonomy_df)]

# ad rel ab normalized per sample (to account for different sampling depth)
sample.filter.tree.rel <- transform_sample_counts(sample.filter.tree, function(x) x/sum(x))
otu_data2 <- otu_table(sample.filter.tree.rel)
total_abundance2 <- rowSums(otu_data2) # cumulative rel abundance
relative_abundance2 <- total_abundance2 * 100/ (sum(total_abundance2))
taxonomy_df$rel_ab2 <- relative_abundance2[rownames(taxonomy_df)]

# ad log ab
log_total_abundance <- log10(total_abundance)
taxonomy_df$log_ab <- log_total_abundance[rownames(taxonomy_df)]

head(taxonomy_df)


sink("plots_peru/05_ggtree_03_layers.txt")
taxonomy_df
sink()

# phyla color by order for ggtree
phylum_order <- names(sort(tapply(taxonomy_df$total_ab, taxonomy_df$phylum, sum), decreasing = TRUE))
#phylumGradient <- colorRampPalette(brewer.pal(12, "Paired"))(length(phylum_order)) # gradient interpolation
#phylumPalette <- brewer.pal(12, "Paired")[seq_len(length(phylum_order))]  # use plain order of colors
#phylum_select = brewer.pal(n = 12, "Paired")[c(2, 4, 6, 8, 10)] # select only bright colors
phylumGradient <- colorRampPalette(brewer.pal(12, "Paired")[2:10])(length(phylum_order)) # gradient with selected bright colors
color_phyla <- setNames(phylumGradient, phylum_order)
# scale_fill_manual(values = color_phyla)



# prevalence sample.filter.tree
# Prevalence cutoff 10 reads minimum
prev_cutoff <- 20
prevalence_df <- data.frame(
  ASV = rownames(otu_table(sample.filter.tree)),
  prevalence = apply(otu_table(sample.filter.tree), 1, function(x) sum(x > prev_cutoff)) 
)
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
  prevalence = apply(otu_table(Satyrinae), 1, function(x) sum(x > prev_cutoff))
)
Satyrinae_prevalence_df$rel_prev <- Satyrinae_prevalence_df$prevalence / ncol(otu_table(Satyrinae))


Pierinae_prevalence_df <- data.frame(
  ASV = rownames(otu_table(Pierinae)),
  prevalence = apply(otu_table(Pierinae), 1, function(x) sum(x > prev_cutoff))
)
Pierinae_prevalence_df$rel_prev <- Pierinae_prevalence_df$prevalence / ncol(otu_table(Pierinae))

Heliconiinae_prevalence_df <- data.frame(
  ASV = rownames(otu_table(Heliconiinae )),
  prevalence = apply(otu_table(Heliconiinae ), 1, function(x) sum(x > prev_cutoff))
)
Heliconiinae_prevalence_df$rel_prev <- Heliconiinae_prevalence_df$prevalence / ncol(otu_table(Heliconiinae))

Coliadinae_prevalence_df <- data.frame(
  ASV = rownames(otu_table(Coliadinae)),
  prevalence = apply(otu_table(Coliadinae), 1, function(x) sum(x > prev_cutoff))
)
Coliadinae_prevalence_df$rel_prev <- Coliadinae_prevalence_df$prevalence / ncol(otu_table(Coliadinae))

Nymphalinae_prevalence_df <- data.frame(
  ASV = rownames(otu_table(Nymphalinae)),
  prevalence = apply(otu_table(Nymphalinae), 1, function(x) sum(x > prev_cutoff))
)
Nymphalinae_prevalence_df$rel_prev <- Nymphalinae_prevalence_df$prevalence / ncol(otu_table(Nymphalinae))


Dismorphiinae_prevalence_df <- data.frame(
  ASV = rownames(otu_table(Dismorphiinae)),
  prevalence = apply(otu_table(Dismorphiinae), 1, function(x) sum(x > prev_cutoff))
)
Dismorphiinae_prevalence_df$rel_prev <- Dismorphiinae_prevalence_df$prevalence / ncol(otu_table(Dismorphiinae))





combined_rel_prev.df <- data.frame(
  #ASV = Heliconiinae_prevalence_df$ASV,
  Satyrinae = Satyrinae_prevalence_df$rel_prev,
  #Dismorphiinae = Dismorphiinae_prevalence_df$rel_prev,
  Pierinae = Pierinae_prevalence_df$rel_prev,
  Heliconiinae = Heliconiinae_prevalence_df$rel_prev,
  Coliadinae = Coliadinae_prevalence_df$rel_prev,
  Nymphalinae = Nymphalinae_prevalence_df$rel_prev
)

rownames(combined_rel_prev.df) <- rownames(taxonomy_df)

prevcolors <- rev(brewer.pal(10, "RdBu")) #RdBu Spectral


# Core Tree fan including layers
fan.core <- ggtree(tree.root, layout="fan", size=0.2, open.angle=25)
fan.core2 <- ggtree(tree.root, layout="fan", size=0.2, open.angle=30) +
  geom_tiplab(aes(label = ""),  align = TRUE, offset = 0.01, linetype = "dotted")

fan.dots <- ggtree(tree.root, layout="fan", size=0.2, open.angle=25) +
  geom_point(aes(color = taxonomy_df$phylum),
             data = ~ .x[.x$isTip, ],  
             size = 2) +
             scale_color_manual(values = color_phyla)

fan.prev_all <- gheatmap(fan.core2, combined_rel_prev.df, 
                            width = 0.3,
                            offset = 0.03,
                            colnames_position = "bottom",
                            colnames_angle =90, font.size = 2.8,
                            hjust = 1,
                            colnames_offset_y = 0) +
  #scale_fill_distiller(palette="RdYlBu")  #RdYlBu, #RdBu
  scale_fill_viridis_c(option = "inferno", name = "Prevalence") 
  
  
 

fan.core.bar <- fan.prev_all    + new_scale_fill() +
  geom_fruit(
    data=taxonomy_df, 
    geom=geom_bar, 
    mapping=aes(y=label, x=rel_ab, fill = phylum),
    pwidth=0.5,
    offset = 0.44, # 0.2
    orientation="y", 
    stat="identity",
    axis.params = list(
      axis = "x",
      #title = "rel",
      #title.size = 3,
      text.size = 2.8,
      hjust = 0.5, # 0.5
      vjust = 1,
      line.color = "grey"
              ),
    grid.params = list(size = 0.2, color = "grey")) + scale_fill_manual(values = color_phyla)
 

fan.core.label <- fan.core.bar  + geom_tiplab(size = 3, align = T, offset = 0.19, #nudge_y = 0.1, # offset 0.13
            aes(label = ifelse(label %in% core_tips$label, 
                               core_tips$genus_name[match(label, core_tips$label)], "")), linetype = NA)   # Only show core genera
  theme(legend.position = "right") 

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
  
  


fan.core.annotation <- fan.core.bee + annotate(
  "text", 
  x = 0,                 # Centered on the circular layout
  y = 0,                 # Adjust distance from the center
  label = "rel. ab [%]", # Axis label text
  size = 2.8,              # Text size
  angle = 0,            # Adjust angle for readability
  hjust = -4,
  vjust = 2.5
) +
  theme(
    legend.key.size = unit(0.6, "lines"),  # smaller legend (usually 1 line)
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.5),     # place legend closer to plot
    legend.justification = c(0, 0.5)           # Right-center of the legend box
  ) +
  theme(
    plot.margin = margin(0, 40, 0, 0)  # increase right margin for legend
  ) 

# figure border highlight 
# +  theme(panel.border = element_rect(color = "blue", fill = NA, linewidth = 1)) +
# +  theme(plot.background = element_rect(color = "red", fill = NA, linewidth = 1)) 
  
  

             

#### ggtree square tree  -----
#square.core <- ggtree(tree.root, size=0.2)
head(taxonomy_df)
all(tree.root$tip.label %in% rownames(taxonomy_df)) #check if labels are okay
all(tree.root$nodes %in% rownames(taxonomy_df)) #check if labels are okay

# use ps instead, but does not work later with bar axis
#square.core <- ggtree(sample.filter.tree) +
#  geom_tiplab(size = 2, hjust = 1,  align=T,
#              offset = 0.02,
#              aes(label=genus, color=phylum)) +
#  scale_x_continuous(expand = c(0.4, 0)) 


# Or merge tree with tax data
# Ad a row label before merge
colnames(taxa_data) <- paste0("tree_", colnames(taxa_data)) # rename taxa to avoid name conflict with taxonomy_df
taxa_data2 <- taxa_data %>% rownames_to_column("label")
# Ad taxdata to tree
square.tree <- ggtree(tree.root, size=0.2) %<+% taxa_data2# + xlim(NA, 0.3) # ad taxa_data2 dataframe for taxa color 
square.tree2 <- ggtree(tree.root, size=0.2)  

# Ad tax info as tip label
square.core <- square.tree + geom_tiplab(aes(label=tree_genus, color=tree_phylum
                                   ),   linetype = NA, # NA  "dotted" 
                      align=T,size=2.5, hjust = 1, offset = 0.08, #family='mono'#, linetype = "dotted" , linesize = .7, offset = 0.035 0.005
                      show.legend = FALSE, lineend = "round",  lineheight = 0.9
                      ) + scale_color_manual(values = color_phyla)




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


square.core.bar <- square.prev_all     + new_scale_fill() +
  geom_fruit(
    data=taxonomy_df, 
    geom=geom_bar, 
    mapping=aes(y=label, x=rel_ab, fill=phylum), #, fill = phylum
    pwidth=0.5, # 0.5 50% space for geom bar
    offset = 0.45, # 0.25 0.18 0.31
    orientation="y", 
    stat="identity",
    axis.params = list(
      axis = "x",
      title = NULL,
      text.size = 3,
      hjust = 0.5, 
      vjust = 1.5,
      line.color = "black"
      ),
    grid.params = list(size = 0.2, color = "black")
  ) +  # coord_cartesian(clip = "off") + 
  scale_fill_manual(values = color_phyla) +
  labs(x = "Rel. Abundance [%]") +
  theme(axis.title.x = element_text(
      size = 10,         # Adjust title font size
      vjust = -1,        # Lower the title below the axis
      hjust = 0.85          # Align the title to the far right (under bars)
        ) ) # Adjust tick size 


square.label <- square.core.bar  + geom_tiplab(size = 2.5, align = T, offset = 0.19, #nudge_y = 0.3 , # 0.13 0.085
                                           aes(label = ifelse(label %in% core_tips$label, 
                                                              core_tips$genus_name[match(label, core_tips$label)], "")), linetype = NA)  # Only show core genera
square.bee.core <- square.label + 
  geom_tiplab(size = 5, align = TRUE, offset = 0.08, #nudge_y = -0.4,
              aes(label = ifelse(label %in% bee_tips$label, "\u2022", "")),
              linetype = NA,
              family = "mono")   # Only show * for matching tips


square.point <- square.bee.core + geom_point(aes(color = taxonomy_df$phylum, x = x + 0.005),
                                              data = ~ .x[.x$isTip, ],  
                                              size = 2, show.legend = FALSE) +
                                              scale_color_manual(values = color_phyla)

square.annotation <-  square.point +   geom_treescale(fontsize=3, linesize=0.1, x=0, y=-2) +
  theme(legend.position.inside=c(1.3, 0.5),
        #legend.background=element_rect(fill=NA),
        #legend.title=element_blank(),
        #legend.text=element_text(size=12),
        #legend.spacing.y = unit(0.2, "cm"),
  ) 

# convert inch to cm -> 16 / 2.54 
pdf("plots_peru/05_ggtree_03_layers.pdf", width= 20 / 2.54 , height= 16 / 2.54)
fan.core.label
fan.core.annotation
square.annotation 
dev.off()



## 06 Country Color Chord Alluvial  --------------

# remove host_genus with only 1 count
core.genus.rel.df <- as.data.frame(sample_data(core.genus.rel))
host_counts <- core.genus.rel.df  %>% group_by(host_subfamily) %>%  summarise(n = n())

# Filter out host genera that occur more than 2 times
hosts_to_keep <- host_counts %>% filter(n > 2) %>% pull(host_subfamily)

# Subset the phyloseq object to keep only the selected host genera
core.filtered <- subset_samples(core.genus.rel, host_subfamily %in% hosts_to_keep)

# Network country
countrymat <- otu_table(merge_samples(core.filtered,group="country"))
plotweb(data.frame(t(otu_table(countrymat))))
#chordDiagram(countrymat) # works inverse taxa country
matcountry = t(otu_table(countrymat))
#chordDiagram(matcountry) # does not work

# Make long format for chord diagram and alluvial plot 
df_long <- as.data.frame(as.table(as.matrix(matcountry))) 
colnames(df_long) <- c("Taxa", "Sample", "Abundance")
#chordDiagram(df_long) # country taxa works

# customn color
countryColors <- brewer.pal(n = min(12, nrow(countrymat)), name = "Dark2") # Colors for country Dark2
names(countryColors) <- c(rownames(countrymat))
sectorColorsCountry <- c(core_palette, countryColors)

# Network connection dprime
dfun(t(matcountry)) # degree distribution, how connected nodes are (number of partners)
# dprime (d′) Normalized Specialization Index. It ranges from 0 (no specialization) to 1 (perfect specialist).

H2fun(t(matcountry), H2_integer = F) 
# 𝐻′₂ Standardized Specialization Index (generalization to the entire interaction web)


### 06 Country Alluvial plot 

# Sort dataframe like in plotweb
# Compute the ratio of abundance per taxon across countries 
mat_frame <- as.data.frame(matcountry) # mat

taxa_abundance <- mat_frame %>%
  mutate(Ratio = Peru / (Germany + Peru))

taxa_abundance <- taxa_abundance %>%
  arrange(Ratio)

sorted_taxa <- rownames(taxa_abundance)

df_sorted <- df_long %>% # sort dataframe
mutate(
  Taxa = factor(Taxa, levels = sorted_taxa),  # Sort Taxa based on sorted_taxa
  Sample = factor(Sample, levels = sort(levels(Sample)))  # Alphabetical order for countries
  ) %>%
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
    color = "black"
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
  scale_fill_manual(values=core_palette) +
  scale_x_discrete(limits = c("Taxa", "Sample")) +
  theme_heat() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),panel.border = element_blank(), legend.position = "none") 


# normalize for rel abundance
df_sorted_norm <- df_sorted %>%
  group_by(Sample, Taxa) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  complete(Sample, Taxa, fill = list(RelAbundance = 0)) %>%  # Ensure all taxa appear in all groups
  filter(RelAbundance > 0)  # Remove any zero values to clean the visualization

AlluTaxaCountry_dual <- ggplot(df_sorted_norm,
                           aes(x = Sample, stratum = Taxa, alluvium = Taxa,
                               y = RelAbundance, fill = Taxa, label = Taxa)) +
  geom_flow(alpha = 0.5) +  # Smooth transitions between host groups
  geom_stratum(alpha = 0.8, color = "black") +  # Host group categories
  geom_text(stat = "stratum", size = 3) +  # Add taxa names
  theme_line2() +
  scale_fill_manual(values=core_palette) +
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

core.melt2 <- core.melt
core.melt2$genus <- factor(core.melt$genus, levels = sorted_taxa)

# apply taxa order to label
core_country_label$genus <- factor(core_country_label$genus,
                                   levels = sorted_taxa)

core_stripchart.only.country.sorted <- ggplot(core.melt2, aes(x = country, y = Abundance)) +
  geom_jitter(aes(color = country),width = 0.2, size = 3, alpha = 0.6) +  # points on top
  geom_violin() +
  stat_summary(fun = mean,  geom = "crossbar", width = 0.2,             
               color = "black",  linewidth = 0.5  ) +
  facet_wrap(~ genus) +
  theme_grid() +
  theme(axis.text.x = element_text(size = 8)) +
  scale_colour_manual(values = country_colorfix) +
  scale_y_continuous(labels = scales::label_percent(scale = 100)) +
  labs(x = "", y = "rel ab [%]", color = "country") +
  geom_text(
    data = core_country_label,
    aes(x = 1.5, # x-position for label 1.5 centered
        y = 0.95, # y-position above data points
        label = label
    ),
    inherit.aes = FALSE,
    size = 3  ) 


pdf("plots_peru/06_CountryAlluvial_color.pdf", width=8, height=8)
Allu_sorted
AlluTaxaCountry
AlluTaxaCountryreverse
AlluTaxaCountry2
AlluTaxaCountry_dual
AlluCountry
core_stripchart.only.country.sorted
dev.off()



## 07 Final Figures ggarrange -----

empty.plot <- ggplot() + theme_void()

betadisp.subfamily.theme <- betadisp.subfamily.ordered  + theme(
  legend.key.size = unit(0.9, "lines") )

PCoA.genus.wrap.tribe_noleg <- PCoA.genus.wrap.tribe + theme(legend.position = "none")

core.genus.boxplot.subfam.no.text <- core.genus.boxplot.subfam + theme(axis.text.y = element_blank()) 


fig1.arrange <- ggarrange(empty.plot , alpha.shannon2,  PCoA.host.country ,    labels = c('a', 'b', 'c'),
                          common.legend = T, legend = "right", ncol = 3,   nrow = 1, align = "h", widths = c(1, 1, 1) )

fig1.arrangeb <- ggarrange(empty.plot , alpha.shannon2,  alpha.shannon.tribe2 ,    labels = c('a', 'b', 'c'),
                          common.legend = T, legend = "right", ncol = 3,   nrow = 1, align = "h", widths = c(1, 1, 1) )


fig2.betadisp.tribe.wrap <- ggarrange( betadisp.subfamily.theme, PCoA.host.overview.tribe ,   labels = c('a', 'b'),
                                       common.legend = T, legend = "right", ncol = 2,   nrow = 1, align = "h", widths = c(1.2 ,2) )

# Fig 3 gen.core.abundance patchwork
gen.core.abundance <- (order.core + theme(legend.position = "right"))+ 
  (sample.order.top.bar.flip + theme(legend.position = "right") ) +
  (order.fam.country.mirrored.genus + theme(legend.position = "none") ) +
  plot_layout(ncol = 3, guides = "collect", widths = c(0.5,0.8, 2)) &  # , guides = "collect"
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14, hjust = 0),
        plot.tag.position = c(0.01, 0.98),
        plot.margin = margin(3, 3, 3, 3),
        legend.key.size = unit(0.9, "lines") )      # key box size

# Fig 5 net core
net.core.comp <- 
  (AlluTaxaCountryreverse) +
  (core.genus.abundance.subfam + theme(legend.position = "none") ) +
  (core.genus.boxplot.subfam.no.text + theme(legend.position = "none") ) +
  (core.shannon+ theme(legend.position = "none") ) +
  plot_layout(ncol = 4, widths = c(0.9, 0.7, 0.7, 1.4),nrow=1, height = c(1,1)) &  
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14, hjust = 0),
        plot.tag.position = c(0.01, 0.98),
        plot.margin = margin(4, 4, 4, 4))  # styling tags



# supplement S8 stripchart arrange
stripchart.arrange <- ggarrange( core.genus.boxplot , core_stripchart.rebuild, labels = c('a', 'b'),
           common.legend = T, legend = "none", ncol = 2,   nrow = 1, align = "h", widths = c(1, 2) )


# Supplement family core
fam.core.abundance2 <- ggarrange(family.core,  core.family.abundance, core.family.boxplot, labels = c('a', 'b', 'c'),
                                 common.legend = F, legend = "right", ncol = 3,   nrow = 1, align = "h" )


# supplement S3 PCoA genus tribe
betadisp.tribe.wrap2 <- ggarrange(betadisp.genus, PCoA.genus.wrap.tribe_noleg, marginalR2.plot,   labels = c('a', 'b', 'c'),
                                  common.legend = F, ncol = 3,   nrow = 1, align = "h", widths = c(1,1.2,0.4) ) # , legend = "right"


# supplement S7 Endosymbiont
endosymbiont.arrange <- ggarrange(endo.bar3, endo.bar2,    labels = c('a', 'b'),
          common.legend = T, legend = "right", ncol = 2,   nrow = 1, align = "h", widths = c(1,2) )



### 07 Final Figure patchwork

# supplement S1 alpha.country patchwork
plot_alpha.country <- (alpha.shannon.tribe2 + theme(legend.position = "none"))+ 
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
plot_alpha.elevation <- (
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
patch_venn_core_noncore_plot <- (left_col | right_col) +
  plot_layout(widths = c(0.7,1.2)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 14, hjust = 0),
    #plot.tag.position = c(0.02, 0.98),
    #plot.margin = margin(3, 1, 3, 1),
    legend.key.size = unit(0.8, "lines"),      # key box size legend
  )



### ggsave Figures -------------

# figure dimensions
# 1/3 page  20 x 7 cm
# 2/3 page  20 x 20 cm
# full page 20 x 25 cm
# 1 column  10 x 12 cm

ggsave("plots_peru/Fig1.png", plot = fig1.arrange, width = 30, height = 9, units = "cm", dpi = 300)
ggsave("plots_peru/Fig1.pdf", plot = fig1.arrange, width = 30, height = 9, units = "cm", dpi = 300)

ggsave("plots_peru/Fig2.png", plot = fig2.betadisp.tribe.wrap , width = 30, height = 9, units = "cm", dpi = 300)
ggsave("plots_peru/Fig2.pdf", plot = fig2.betadisp.tribe.wrap , width = 30, height = 9, units = "cm", dpi = 300)

ggsave("plots_peru/Fig3.png", plot = gen.core.abundance  , width = 30, height = 10, units = "cm", dpi = 300)
ggsave("plots_peru/Fig3.pdf", plot = gen.core.abundance  , width = 30, height = 10, units = "cm", dpi = 300)

ggsave("plots_peru/Fig4.png", plot = fan.core.annotation, width = 20, height = 16, units = "cm", dpi = 300)
ggsave("plots_peru/Fig4.pdf", plot = fan.core.annotation, width = 20, height = 16, units = "cm", dpi = 300)

ggsave("plots_peru/Fig5.png", plot = net.core.comp, width = 30, height = 10, units = "cm", dpi = 300)
ggsave("plots_peru/Fig5.pdf", plot = net.core.comp, width = 30, height = 10, units = "cm", dpi = 300)


### ggsave Supplemental Figures -----

ggsave("plots_peru/FigS1_alpha.country.png", plot = plot_alpha.country , width = 24, height = 14, units = "cm", dpi = 300)
ggsave("plots_peru/FigS1_alpha.country.pdf", plot = plot_alpha.country , width = 24, height = 14, units = "cm", dpi = 300)

ggsave("plots_peru/FigS2_alphaelevation.png", plot = plot_alpha.elevation , width = 24, height = 20, units = "cm", dpi = 300)
ggsave("plots_peru/FigS2_alphaelevation.pdf", plot = plot_alpha.elevation , width = 24, height = 20, units = "cm", dpi = 300)

ggsave("plots_peru/FigS3_beta_genus.png", plot = betadisp.tribe.wrap2 , width = 30, height = 13, units = "cm", dpi = 300)
ggsave("plots_peru/FigS3_beta_genus.pdf", plot = betadisp.tribe.wrap2 , width = 30, height = 13, units = "cm", dpi = 300)

ggsave("plots_peru/FigS4_corevenn_alternative.png", plot = patch_venn_core_noncore_plot , width = 34, height = 22, units = "cm", dpi = 300)
ggsave("plots_peru/FigS4_corevenn_alternative.pdf", plot = patch_venn_core_noncore_plot , width = 34, height = 22, units = "cm", dpi = 300)

ggsave("plots_peru/FigS5_core_square.png", plot = square.annotation , width = 25, height = 25, units = "cm", dpi = 300)
ggsave("plots_peru/FigS5_core_square.pdf", plot = square.annotation , width = 25, height = 25, units = "cm", dpi = 300)

ggsave("plots_peru/FigS6_country_core.png", plot = core_stripchart.only.country.sorted, width = 24, height = 15, units = "cm", dpi = 300)
ggsave("plots_peru/FigS6_country_core.pdf", plot = core_stripchart.only.country.sorted, width = 24, height = 15, units = "cm", dpi = 300)

ggsave("plots_peru/FigS7_bubble.country.png", plot = ASV.bubble.core.country , width = 25, height = 30, units = "cm", dpi = 300)
ggsave("plots_peru/FigS7_bubble.country.pdf", plot = ASV.bubble.core.country , width = 25, height = 30, units = "cm", dpi = 300)

ggsave("plots_peru/FigS8_endo.png", plot = endosymbiont.arrange , width = 25, height = 12, units = "cm", dpi = 300)
ggsave("plots_peru/FigS8_endo.pdf", plot = endosymbiont.arrange , width = 25, height = 12, units = "cm", dpi = 300)

ggsave("plots_peru/FigS9_corestrip.png", plot = stripchart.arrange , width = 24, height = 12, units = "cm", dpi = 300)
ggsave("plots_peru/FigS9_corestrip.pdf", plot = stripchart.arrange , width = 24, height = 12, units = "cm", dpi = 300)

ggsave("plots_peru/FigSx_core-fam_genus.png", plot = fam.core.abundance2, width = 35, height = 15, units = "cm", dpi = 300)
ggsave("plots_peru/FigSx_core-fam_genus.pdf", plot = fam.core.abundance2, width = 35, height = 15, units = "cm", dpi = 300)


# optional: Run extended workflow for extra figures
if (file.exists("R_16S_AW_extension.R")) {
  message("Running extended workflow...")
  source("R_16S_AW_extension.R", echo = TRUE, local = TRUE)
}

message("Pipeline finished at ", Sys.time())
