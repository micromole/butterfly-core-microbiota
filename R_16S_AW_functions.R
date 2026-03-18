### Load custom themes and functions ----
theme_line2 <- function(base_size = 11) {
  theme_linedraw()%+replace%
    theme(
      panel.grid.major =  element_blank(),
      panel.grid.minor =  element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(colour = "black", size=10, margin = margin(b = 3)),
      #panel.border = element_rect(colour = "black", fill = NA),
      axis.text= element_text(color="black",size=10),
      axis.ticks = element_line(color = "black")
    ) }

theme_heat <- function(base_size = 11) {
  theme_linedraw()%+replace%
    theme(
      panel.grid.major =  element_blank(),
      panel.grid.minor =  element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      strip.text = element_text(colour = "black", size=10, margin = margin(b = 3)),
      #panel.border = element_rect(colour = "black", fill = NA),
      axis.text= element_text(color="black",size=10),
      axis.ticks = element_line(color = "black"),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank()
    ) }


theme_grid <- function(base_size = 11) {
  theme_linedraw()%+replace%
    theme(
      panel.grid.major =  element_line(linetype=1,color="grey90"),
      panel.grid.minor =  element_line(linetype=1,color="grey90"),
      strip.background = element_blank(),
      strip.text = element_text(colour = "black", size=10, margin = margin(b = 3)),
      #panel.border = element_rect(colour = "black", fill = NA),
      axis.text= element_text(color="black",size=10),
      axis.ticks = element_line(color = "black")
    ) }


##function "replace_tax_prefixes"
replace_tax_prefixes <- function(phyloseq){
  tmp_tax_table <- apply(tax_table(phyloseq), c(1, 2),function(y) gsub("^\\w:","",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub("_spc_.*","_spc",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub("_"," ",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub(";$","",y))
  tax_table(phyloseq)<- tmp_tax_table
  return(phyloseq)
}

##function "propagate_incomplete_taxonomy"
propagate_incomplete_taxonomy <- function(phyloseq){
  taxranks <- colnames(tax_table(phyloseq))
  for (i in 2:length(taxranks)){
    tax_table(phyloseq)[tax_table(phyloseq)[,taxranks[i]]=="",taxranks[i]]<-paste(tax_table(phyloseq)[tax_table(phyloseq)[,taxranks[i]]=="",taxranks[i-1]],"_spc",sep="")
  }
  return(phyloseq)
}

##function "calc_prevalence"
calc_prevalence<- function(ps, rank = "Phylum") {
  
  # Calculate prevalence per taxon
  prev <- apply(
    X = otu_table(ps),
    MARGIN = ifelse(taxa_are_rows(ps), 1, 2),
    FUN = function(x) sum(x > 0)
  )
  
  # Combine with taxonomy and abundance
  prevdf <- data.frame(
    Prevalence = prev,
    TotalAbundance = taxa_sums(ps),
    tax_table(ps)
  )
  
  # Summarize by taxonomic rank
  rank_summary <- plyr::ddply(
    prevdf,
    rank,
    function(df) cbind(
      MeanPrevalence = mean(df$Prevalence),
      SumPrevalence = sum(df$Prevalence)
    )
  )
  
  return(list(
    taxa_table = prevdf,
    rank_summary = rank_summary
  ))
}

## Label samples with low throughput with LT (function)
label_low_throughput <- function(phyloseq, threshold){
  sample_names(phyloseq)[sample_sums(phyloseq)<threshold]<-paste(sample_names(phyloseq)[sample_sums(phyloseq)<threshold],"|LT",threshold,sep="")
  return(phyloseq)
}


### Test for variables that correlate best with community matrix
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

# Run kruskal.test output as table: run_kruskal(core.melt, "host_subfamily")
run_kruskal <- function(data, group_var) {
  data %>%
    group_by(genus) %>%
    summarise(
      test = list(kruskal.test(reformulate(group_var, "Abundance")))
    ) %>%
    mutate(
      statistic = sapply(test, \(x) x$statistic),
      df        = sapply(test, \(x) x$parameter),
      p_value   = sapply(test, \(x) x$p.value),
      p_adj     = p.adjust(p_value, method = "BH")
    ) %>%
    select(-test) %>%
    mutate(
      statistic = round(statistic, 3),
      p_value = sprintf("%.3f", p_value),
      p_adj   = sprintf("%.3f", p_adj),
      signif = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
}