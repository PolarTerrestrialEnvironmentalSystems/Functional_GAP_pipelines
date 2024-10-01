library("plyr")
library("dplyr")
library("tidyverse")
library("ggplot2")
library("tidyr")
library("data.table")

setwd("set-to-the-path")

# load prepared ncbi tax reference file
ncbi_df <- read.csv("ncbi_df.csv")

# If you do not have, ncbi file preparation is the following: 
# you need a taxdump from ncbi website
# taxdir="taxdump"
# 
# # taxon name (eg.scientific name, blast name..., Not nodes!)
# ncbi_names=getnames(taxdir)
# 
# # read rankedlineage.dmp. 
# # The ranks are species, genus, family, order, class, phylum, kingdom, and superkingdom.
# rankedlineage=as.data.frame(read_tsv(paste0(taxdir, "/rankedlineage.dmp"),
#                                      col_names = c("taxid", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"),
#                                      col_types = ("i-c-c-c-c-c-c-c-c-c-")))
# # read full name lineage
# fullnamelineage=as.data.frame(read_tsv(paste0(taxdir, "/fullnamelineage.dmp"),
#                                        col_names = c("taxid", "tax_name", "fullnamelineage"),
#                                        col_types = ("i-c-c-")))
# 
# ncbi_df=merge(fullnamelineage[c("taxid", "fullnamelineage")], rankedlineage, by = "taxid")

# read pydamage all result
pydamage <- read.csv("all_pydamage_result.csv", header = F)
pydamage <- pydamage[-1, ]

# fixed all the columns
fixed_colname <- c("sample_name", "contig_id","predicted_accuracy","null_model_p0","null_model_p0_stdev","damage_model_p",
                   "damage_model_p_stdev","damage_model_pmin","damage_model_pmin_stdev","damage_model_pmax","damage_model_pmax_stdev",
                   "pvalue","qvalue","RMSE","nb_reads_aligned","coverage",
                   "reflen","CtoT-0","CtoT-1","CtoT-2","CtoT-3",
                   "CtoT-4","CtoT-5","CtoT-6","CtoT-7","CtoT-8",
                   "CtoT-9","CtoT-10","CtoT-11","CtoT-12","CtoT-13",
                   "CtoT-14","CtoT-15","CtoT-16","CtoT-17","CtoT-18",
                   "CtoT-19")

colnames(pydamage) <- fixed_colname

# these rows should be discarded
pydamage <- pydamage[-which(pydamage$sample_name == "reference"),]

# check the file
dim(pydamage)

# load kraken file from contigs
kraken <- read.csv("all_contig.kraken", sep=" ", header = F)
unique(kraken$sample_name)

# fix the colnames to merge
colnames(kraken) <- c("sample_name","contig_id","taxid")

# get rid of extension from the name
kraken$sample_name <- gsub("_conf0_contig\\.kraken", "", kraken$sample_name)

# check the ncbi reference file
colnames(ncbi_df)

# merge kraken and ncbi file
classified_ncbi=merge(kraken, ncbi_df, by = "taxid", all.x=T)

# check both files
unique(kraken$sample_name)
unique(pydamage$sample_name)


#Find NA If there is###

genus_taxa <- classified_ncbi %>% select(taxid,sample_name,contig_id,genus, kingdom)
NA_values <- which(is.na(genus_taxa$genus))
genus_taxa2 <- genus_taxa[-NA_values,]

# Merge sample_name and contig_id into a new column
classified_ncbi2 <- classified_ncbi %>%
  mutate(merged_name = paste(sample_name, contig_id, sep = "_merged_"))

# pydamage <- pydamage %>%
#   mutate(merged_name = paste(sample_name, contig_id, sep = "_merged_"))

# which(pydamage$sample_name == "JK081L-2_S16_JK081L-2_S2")

pydamage_all=merge(classified_ncbi2, pydamage, by = "merged_name", all.x=T)


pydamage_all=pydamage_all[-which(is.na(pydamage_all$sample_name.y)),]
pydamage_all$nb_reads_aligned = (as.numeric(pydamage_all$nb_reads_aligned))


# save the final file as R data
save(pydamage_all, file = "pydamage_all_result_taxon.RData")


### General Trend ###
#setwd("set-to-the-path")
#load("pydamage_all_result_taxon_wide_age_numeric.RData")

ka_filter       <- 0
reads_filter    <- 0
accuracy_filter <- 0.5
pval_filter     <- 0
qval_filter     <- 0
taxa <- "Eukaryota"
len_filter <- 600
# -------------------------------------------------------------------------------------------------
pydamage_all_numeric$reflen = as.numeric(pydamage_all_numeric$reflen)
pydamage_eukaryota <- pydamage_all_numeric %>% filter(superkingdom == "Eukaryota" & reflen >= len_filter)

pydamage_eukaryota_selected <- pydamage_eukaryota %>%
  #filter(tax_name %in% taxalist_eukaryota_selected_list$tax_name) %>%
  # filter(tax_name == taxa) %>%
  filter(ka > ka_filter & nb_reads_aligned > reads_filter) %>%
  filter(predicted_accuracy >= accuracy_filter & pvalue >= pval_filter & qvalue >= qval_filter)

eukaryota_selected_ntaxa  <- length(unique(pydamage_eukaryota_selected$tax_name))
eukaryota_selected_nreads <- sum(pydamage_eukaryota_selected$nb_reads_aligned)

# eukaryota_selected_taxa_df <- data.frame(Taxon = sort(unique(pydamage_eukaryota_selected$tax_name)))
# write.csv(eukaryota_selected_taxa_df, file = paste0(datapath,"Pydamage/pydamage_eukaryota_selected_taxalist_ka_",ka_filter,"_reads_",reads_filter,"_accuracy_",accuracy_filter,"_pval_",pval_filter,"_qval_",qval_filter,".csv"), row.names = FALSE)

# -------------------------------------------------------------------------------------------------

pydamage_eukaryota_selected_long <- pydamage_eukaryota_selected %>%
  select(starts_with("CtoT")) %>%
  pivot_longer(starts_with("CtoT"),
               names_to = "position",
               values_to = "freq") %>%
  arrange(position) %>%
  separate(col = "position",
           into = c("pos_par1", "pos_par2"),
           sep = "-",
           remove = FALSE,
           extra = "drop") %>%
  mutate(pos_par2 = as.numeric(pos_par2))


pydamage_eukaryota_position_label <- unique(pydamage_eukaryota_selected_long[,c(1,3)])
pydamage_eukaryota_position_label <- pydamage_eukaryota_position_label %>% arrange(pos_par2) %>% rename(pos_label = position)

pydamage_eukaryota_selected_ct_df <- pydamage_eukaryota_selected %>% select(starts_with("CtoT"))

pydamage_eukaryota_selected_ct_trend <- data.frame(pydamage_eukaryota_position_label, mean = colMeans(pydamage_eukaryota_selected_ct_df), median = apply(pydamage_eukaryota_selected_ct_df, 2, median))

pydamage_eukaryota_mean <- ggplot(data = pydamage_eukaryota_selected_long, aes(x = pos_par2, y = freq)) +
  geom_point(color = "grey80",size = 0.6) +
  geom_line(data = pydamage_eukaryota_selected_ct_trend, aes(x = pos_par2, y = mean), color = "red") +
  scale_x_continuous(breaks = pydamage_eukaryota_position_label$pos_par2, labels = pydamage_eukaryota_position_label$pos_par2) +
  scale_y_continuous(breaks = seq(0,range(pydamage_eukaryota_selected_long$freq)[2],0.1), labels = seq(0,range(pydamage_eukaryota_selected_long$freq)[2],0.1)) +
  labs(title = taxa,
       subtitle = paste0(taxa,"|", "ka > ",ka_filter," | reads > ", reads_filter," | predicted accuracy >= ",accuracy_filter," | p-value >= ",pval_filter," | q-value >= ",qval_filter," | n_taxa = ",eukaryota_selected_ntaxa," | n_reads = ",eukaryota_selected_nreads,"\ntrendline = mean")) +
  labs(y="C>T Frequency", x="position read") +
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA), panel.grid.minor = element_blank()) +
  theme(axis.title.x=element_text(vjust = -2),
        axis.title.y=element_text(vjust = 4)) +
  theme(plot.title=element_text(size=16),
        plot.subtitle=element_text(size=12)) +
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.spacing.x = unit(1, "mm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.ticks.length = unit(1, "mm"))

#pydamage_eukaryota_mean


ggsave(pydamage_eukaryota_mean, file = paste0("pydamage_Viridiplantae_ka_",taxa, "_",ka_filter,"_nb_reads_",reads_filter,"_accuracy_",accuracy_filter,"_pval_",pval_filter,"_qval_",qval_filter,"_len_",len_filter,"_mean.png"), width = 30, height = 16, units = "cm", dpi = 300)


### Regression Analysis ###
setwd("set-to-the-path")
load("pydamage_all_result_taxon_wide_age_numeric.RData")

  # settings:
  
  accuracy_filter <- 0.5
  reflen_filter   <- 1000
  
  # -------------------------------------------------------------------------------------------------

  taxa <- "Eukaryota"
  
  # -------------------------------------------------------------------------------------------------
  pydamage_all_numeric$reflen = as.numeric(pydamage_all_numeric$reflen)
  
  pydamage_eukaryota <- pydamage_all_numeric %>% filter(superkingdom == "Eukaryota")
  
  pydamage_eukaryota <- pydamage_eukaryota %>% mutate(contig_num = 1)
  
  pydamage_eukaryota_selected <- pydamage_eukaryota %>%
    #filter(tax_name %in% taxalist_eukaryota_selected_list$tax_name) %>%
    # filter(tax_name == taxa) %>%
    filter(predicted_accuracy >= accuracy_filter & reflen >= reflen_filter)
  
  eukaryota_selected_ntaxa  <- length(unique(pydamage_eukaryota_selected$tax_name))
  eukaryota_selected_nreads <- sum(pydamage_eukaryota_selected$nb_reads_aligned)
  
  # save(pydamage_eukaryota_selected, file = paste0(datapath,"Pydamage/pydamage_eukaryota_selected_acc06.RData"))
  # save(pydamage_eukaryota_selected, file = paste0(datapath2,"Pydamage/pydamage_eukaryota_selected","_ka_",ka_filter,"_reads_",reads_filter,"_acc_",accuracy_filter,"_pos1_",pos1_filter,"_pval_",pval_filter,"_qval_",qval_filter,".RData"))
  
  # -------------------------------------------------------------------------------------------------
  
  pydamage_eukaryota_selected_taxsum <- pydamage_eukaryota_selected %>%
    select(tax_name, nb_reads_aligned, age, ka) %>%
    group_by(ka, tax_name) %>%
    mutate(sumcount = sum(nb_reads_aligned)) %>%
    ungroup() %>%
    group_by(tax_name) %>%
    mutate(taxa_occ = n_distinct(ka)) %>%
    ungroup() %>%
    arrange(-sumcount) %>%
    select(-nb_reads_aligned) %>%
    distinct()
  
  # -------------------------------------------------------------------------------------------------
  
  pydamage_eukaryota_selected_taxon_df <- pydamage_eukaryota_selected %>%
    select(tax_name, `CtoT-0`, age, ka, nb_reads_aligned) %>% rename(pos1 = `CtoT-0`) %>%
    mutate(age = as.numeric(age))
  
  cor.test(pydamage_eukaryota_selected_taxon_df$age, pydamage_eukaryota_selected_taxon_df$pos1)
  
  pydamage_eukaryota_selected_taxon_df_15ka <- pydamage_eukaryota_selected_taxon_df %>% filter(age <= 15000)
  
  corcoeff <- cor.test(pydamage_eukaryota_selected_taxon_df_15ka$age, pydamage_eukaryota_selected_taxon_df_15ka$pos1)
  cor.test(pydamage_eukaryota_selected_taxon_df_15ka$age, pydamage_eukaryota_selected_taxon_df_15ka$pos1)
  
  quantile(pydamage_eukaryota_selected_taxon_df_15ka$pos1, probs = seq(0,1,0.05))
  
  linreg <- lm(pos1 ~ age, data = pydamage_eukaryota_selected_taxon_df_15ka)
  
  pydamage_eukaryota_selected_taxon_df_15ka_trendline <- data.frame(age = sort(unique(pydamage_eukaryota_selected_taxon_df_15ka$age)))
  
  pydamage_eukaryota_selected_taxon_df_15ka_trendline <- data.frame(age = sort(unique(pydamage_eukaryota_selected_taxon_df_15ka$age)),
                                                                   predict(linreg, newdata = pydamage_eukaryota_selected_taxon_df_15ka_trendline, interval = 'confidence'))
  
  # -------------------------------------------------------------------------------------------------
  
  pydamage_eukaryota_ct1_age <- ggplot(data = pydamage_eukaryota_selected_taxon_df_15ka, aes(x = age, y = pos1)) +
    geom_point(color = "grey80", shape = 3, size = 0.6) +
    
    geom_ribbon(data = pydamage_eukaryota_selected_taxon_df_15ka_trendline, aes(y = fit, ymin = lwr, ymax = upr), fill = "red", alpha = 0.1) +
    geom_line(data = pydamage_eukaryota_selected_taxon_df_15ka_trendline, aes(x = age, y = fit), color = "red", linewidth = 0.5) +
    
    scale_x_continuous(limits = c(0,15000), breaks = seq(0,25000,5000), labels = seq(0,25000,5000)) +
    scale_y_continuous(limits = c(0,0.5)) +
    labs(title = taxa,
         subtitle = paste0(taxa,"  |  pipeline = pyDamage","  |  n_reads = ",sum(pydamage_eukaryota_selected_taxon_df_15ka$nb_reads_aligned),"  |  r = ",round(corcoeff$estimate, digits=2), ifelse(corcoeff$p.value < 0.001, "  p-value < 0.001", paste0("  p-value = ",round(corcoeff$p.value, digits=3))))) +
    labs(y="C to T substitution frequency", x="age") +
    theme_bw() +
    theme(panel.grid.major=element_line(colour=NA), panel.grid.minor = element_blank()) +
    theme(axis.title.x=element_text(vjust = -2),
          axis.title.y=element_text(vjust = 4)) +
    theme(plot.title=element_text(size=18),
          plot.subtitle=element_text(size=16)) +
    theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 16),
          panel.spacing.x = unit(1, "mm"),
          panel.border = element_rect(fill = NA, colour = "black"),
          axis.ticks.length = unit(1, "mm"))
  
  
  pydamage_eukaryota_ct1_age
  
  
  ggsave(pydamage_eukaryota_ct1_age, file = paste0("pydamage_eukaryota_",taxa,"_accuracy_",accuracy_filter,"_reflength_",reflen_filter,"_ct1_age.png"), width = 30, height = 16, units = "cm", dpi = 300)
  
  # -------------------------------------------------------------------------------------------------
