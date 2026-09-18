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


## function "replace_tax_prefixes" by Alexander Keller
replace_tax_prefixes <- function(phyloseq){
  tmp_tax_table <- apply(tax_table(phyloseq), c(1, 2),function(y) gsub("^\\w:","",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub("_spc_.*","_spc",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub("_"," ",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub(";$","",y))
  tax_table(phyloseq)<- tmp_tax_table
  return(phyloseq)
}

## function "propagate_incomplete_taxonomy" by Alexander Keller
propagate_incomplete_taxonomy <- function(phyloseq){
  taxranks <- colnames(tax_table(phyloseq))
  for (i in 2:length(taxranks)){
    tax_table(phyloseq)[tax_table(phyloseq)[,taxranks[i]]=="",taxranks[i]]<-paste(tax_table(phyloseq)[tax_table(phyloseq)[,taxranks[i]]=="",taxranks[i-1]],"_spc",sep="")
  }
  return(phyloseq)
}

## function "calc_prevalence" for prevalence abundance plots
calc_prevalence <- function(ps, tax_rank = "phylum") {
  
  # 1. Prevalence per taxon
  prev <- apply(
    X = otu_table(ps),
    MARGIN = if (taxa_are_rows(ps)) 1 else 2,
    FUN = function(x) sum(x > 0)
  )
  
  # 2. Build dataframe (ASV level)
  prevdf <- data.frame(
    Prevalence = prev,
    TotalAbundance = taxa_sums(ps),
    as.data.frame(tax_table(ps))
  )
  
  # 3. Subset to taxa present
  prevdf_filtered <- subset(
    prevdf,
    prevdf[[tax_rank]] %in% get_taxa_unique(ps, tax_rank)
  )
  
  return(prevdf_filtered)
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

# function "plate_frame" to show read abundance per extraction plate as heatmap
plate_frame <- function(ps) {
  title <- deparse(substitute(ps))  # get object name
  df <- cbind(
    as.data.frame(sample_data(ps)),
    reads = sample_sums(ps)
  )
  df %>%
    separate(well, into = c("Row", "Column"), sep = "(?<=\\D)(?=\\d)", convert = TRUE) %>%
    arrange(Row, Column) %>%
    mutate(Row = factor(Row, levels = c("H","G","F","E","D","C","B","A"))) %>%
    ggplot(aes(x = Column, y = Row, fill = reads, label = reads)) +
    geom_tile() +
    geom_text(size = 2) +
    facet_wrap(~plate) +
    scale_fill_viridis_c(direction = -1) +
    labs(title = title, x = "Column", y = "Row") +
    theme_minimal()
}

# get partial R2 (variance partitioning for lm model on Type II ANOVA)
get_partial_R2 <- function(model) {
  anovaII <- Anova(model, type = 2)
  ss_res <- sum(residuals(model)^2)
  partial_R2 <- anovaII$`Sum Sq` / 
    (anovaII$`Sum Sq` + ss_res)
  
  data.frame(
    term = rownames(anovaII),
    Df = anovaII$Df,
    F = anovaII$`F value`,
    p = anovaII$`Pr(>F)`,
    partial_R2 = partial_R2,
    row.names = NULL   )
}



# Run kruskal.test output as table: run_kruskal(core.melt, "host_subfamily")
run_kruskal <- function(data, group_var) {
  data %>%
    group_by(genus) %>%
    summarise(
      test = list(kruskal.test(reformulate(group_var, "Abundance"))),
      highest_in = names(which.max(tapply(
        rank(Abundance),
        .data[[group_var]],
        mean,
        na.rm = TRUE
      )))
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
      #p_value = sprintf("%.3f", p_value),
      #p_adj   = sprintf("%.3f", p_adj),
      signif = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
}
