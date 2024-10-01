# Author: Ugur Cabuk, ugur.cabuk@awi.de
library("dplyr")
library("tidyverse")
library("tidypaleo")
library("ggplot2")
library("readr")
library("tidyr")
library("vegan")
library("data.table")

setwd("")

analyze_kegg_diversity <- function(eggnog_file, gene_abundance_file, output_file) {

#eggnog input file
eggnog_input <- fread(eggnog_file)

# gene abundance file
quant_merge <- fread(gene_abundance_file, header = TRUE)
colnames(quant_merge)[1] <- "gene_id"
  
colnames(quant_merge) <- c("gene_id",18715,17948,6400,9190,11536,13746,22982,21753,20630,19593,3023,16946,76,724,1858,4381,5168,15599,14669,7562,6890,8513,20059,22518,264,1212,2504,17222,19192,21068,9816,12992,10428,11791,14020,14315,15123,16084,10091,13930,16541,13588,"gene_length")
quant_merge$gene_length <- NULL


#change the col name of egnnog
colnames(eggnog_input)[1] <- "gene_id"

quant_eggnog_merge <- merge(eggnog_input,quant_merge, by="gene_id",all=T)
quant_eggnog_merge_2 <-quant_eggnog_merge[-which(is.na(quant_eggnog_merge$eggNOG_OGs)),]

quant_eggnog_merge_long <- quant_eggnog_merge_2 %>%
  pivot_longer(
    cols = `18715`:`13588`,
    names_to = "Ages",
    values_to = "TPM_value"
)

ko_abund <- quant_eggnog_merge_long %>%
   group_by(Ages) %>%
   summarise(SUM_CPM=sum(TPM_value))

write.csv(ko_abund, file = paste0(output_file, "_eggnog_annotation.csv"))

# discard non taxon ids
quant_eggnog_merge_2 <-quant_eggnog_merge[-which(is.na(quant_eggnog_merge$eggNOG_OGs)),]

# keep eukaryotes or prokaryotes, here should be changed !
quant_eggnog_merge_3 = quant_eggnog_merge_2[grepl("Eukaryota", quant_eggnog_merge_2$eggNOG_OGs),]


quant_eggnog_merge_long <- quant_eggnog_merge_3 %>%
  pivot_longer(
    cols = `18715`:`13588`,
    names_to = "Ages",
    values_to = "TPM_value"
)


ko_abund <- quant_eggnog_merge_long %>%
   group_by(Ages) %>%
   summarise(SUM_CPM=sum(TPM_value))

write.csv(ko_abund, file = paste0(output_file, "_eukaryota_annotation.csv"))

#if one gene has two KO IDs, discard them.
quant_eggnog_merge_4 <- quant_eggnog_merge_3[!grepl(",", quant_eggnog_merge_3$KEGG_ko),]

quant_eggnog_merge_4 <- quant_eggnog_merge_4 %>%
  mutate(KEGG_ko= strsplit(as.character(KEGG_ko),',', perl=TRUE)) %>%
  unnest(KEGG_ko)
quant_eggnog_merge_4$KEGG_ko<-  gsub('ko:','', quant_eggnog_merge_4$KEGG_ko)

#if one gene has "-" in KO col, discard them.
quant_eggnog_merge_4 <- quant_eggnog_merge_4[-which(quant_eggnog_merge_4$KEGG_ko == "-"),]

#Write the statistics to the files
gene_number_ko <- length((quant_eggnog_merge_4$gene_id))
uniq_ko_ids <- length(unique(quant_eggnog_merge_4$KEGG_ko))


write.table(data.frame(Gene_number_KO_assignment = gene_number_ko, Uniq_KO_assignment = uniq_ko_ids),
              file = paste0(output_file, "_eukaryotes_KO_gene_counts.txt"), sep = "\t", row.names = FALSE)

which(is.na(quant_eggnog_merge_4$KEGG_ko))

#wide to long format
quant_eggnog_merge_long <- quant_eggnog_merge_3 %>%
  pivot_longer(
    cols = `18715`:`13588`,
    names_to = "Ages",
    values_to = "TPM_value"
)

ko_abund <- quant_eggnog_merge_long %>%
   group_by(Ages,KEGG_ko) %>%
   summarise(SUM_CPM=sum(TPM_value))

write.csv(ko_abund, file = paste0(output_file, "_eukaryotes_KOs_abundance.csv"))

### PROKARYOTES ###

# keep eukaryotes or prokaryotes
quant_eggnog_merge_3 = quant_eggnog_merge_2[!grepl("Eukaryota", quant_eggnog_merge_2$eggNOG_OGs),]

quant_eggnog_merge_long <- quant_eggnog_merge_3 %>%
  pivot_longer(
    cols = `18715`:`13588`,
    names_to = "Ages",
    values_to = "TPM_value"
)

ko_abund <- quant_eggnog_merge_long %>%
   group_by(Ages) %>%
   summarise(SUM_CPM=sum(TPM_value))

write.csv(ko_abund, file = paste0(output_file, "_prokaryota_annotation.csv"))

#if one gene has two KO IDs, discard them.
quant_eggnog_merge_4 <- quant_eggnog_merge_3[!grepl(",", quant_eggnog_merge_3$KEGG_ko),]

quant_eggnog_merge_4 <- quant_eggnog_merge_4 %>%
  mutate(KEGG_ko= strsplit(as.character(KEGG_ko),',', perl=TRUE)) %>%
  unnest(KEGG_ko)
quant_eggnog_merge_4$KEGG_ko<-  gsub('ko:','', quant_eggnog_merge_4$KEGG_ko)

#if one gene has "-" in KO col, discard them.
quant_eggnog_merge_4 <- quant_eggnog_merge_4[-which(quant_eggnog_merge_4$KEGG_ko == "-"),]

gene_number_ko <- length((quant_eggnog_merge_4$gene_id))
uniq_ko_ids <- length(unique(quant_eggnog_merge_4$KEGG_ko))


write.table(data.frame(Gene_number_KO_assignment = gene_number_ko, Uniq_KO_assignment = uniq_ko_ids),
              file = paste0(output_file, "_prokaryotes_KO_gene_counts.txt"), sep = "\t", row.names = FALSE)

which(is.na(quant_eggnog_merge_4$KEGG_ko))

#wide to long format
quant_eggnog_merge_long <- quant_eggnog_merge_3 %>%
  pivot_longer(
    cols = `18715`:`13588`,
    names_to = "Ages",
    values_to = "TPM_value"
)

# ALL gene
all_quant_eggnog_merge_long <- quant_eggnog_merge_3 %>%
quant_eggnog_merge_long <- quant_eggnog_merge_3 %>%
  pivot_longer(
    cols = `18715`:`13588`,
    names_to = "Ages",
    values_to = "TPM_value"
)

#ROW_KO_COUNTS
write.csv(quant_eggnog_merge_long, file = paste0(output_file, ".csv"))

ko_abund <- quant_eggnog_merge_long %>%
   group_by(Ages,KEGG_ko) %>%
   summarise(SUM_CPM=sum(TPM_value))

write.csv(ko_abund, file = paste0(output_file, "_prokaryotes_KOs_abundance.csv"))
}


# Normalized
# PROKGAP
prokgap <- analyze_kegg_diversity("non_redundant_PROKGAP_protein_eggNOG.emapper.annotations",
                       "prokGAP_all_lake_lama_gene_quant.CPM.6m.fixed.tsv",
                       "prokgap_6m_CPM_eggnog_abundance")
# EUKGAP
eukgap <- analyze_kegg_diversity("non_redundant_EUKGAP_protein_eggNOG.emapper.annotations",
                       "eukGAP_all_lake_lama_gene_quant.CPM.6m.fixed.tsv",
                       "eukgap_6m_CPM_eggnog_abundance")

# PRECLASS_PROKGAP
preclass_euk <- analyze_kegg_diversity("non_redundant_preclass_prokGAP_protein_eggNOG.emapper.annotations",
                       "preclass_prokGAP_all_lake_lama_gene_quant.CPM.6m.fixed.tsv",
                       "preclass_prokgap_6m_CPM_eggnog_abundance")

# PRECLASS_EUKGAP
preclass_prok <- analyze_kegg_diversity("non_redundant_preclass_eukGAP_protein_eggNOG.emapper.annotations",
                       "preclass_eukGAP_all_lake_lama_gene_quant.CPM.6m.fixed.tsv",
                       "preclass_eukgap_6m_CPM_eggnog_abundance")


# RAW COUNT

# PROKGAP
prokgap <- analyze_kegg_diversity("non_redundant_PROKGAP_protein_eggNOG.emapper.annotations",
                       "prokGAP_all_lake_lama_gene_quant.raw_count.fixed.tsv",
                       "prokgap_raw_count_eggnog_abundance")
 EUKGAP
eukgap <- analyze_kegg_diversity("non_redundant_EUKGAP_protein_eggNOG.emapper.annotations",
                       "eukGAP_all_lake_lama_gene_quant.raw_count.fixed.tsv",
                       "eukgap_raw_count_eggnog_abundance")

 PRECLASS_PROKGAP
preclass_euk <- analyze_kegg_diversity("non_redundant_preclass_prokGAP_protein_eggNOG.emapper.annotations",
                       "preclass_prokGAP_all_lake_lama_gene_quant.raw_count.fixed.tsv",
                       "preclass_prokgap_raw_count_eggnog_abundance")

 PRECLASS_EUKGAP
preclass_prok <- analyze_kegg_diversity("non_redundant_preclass_eukGAP_protein_eggNOG.emapper.annotations",
                       "preclass_eukGAP_all_lake_lama_gene_quant.raw_count.fixed.tsv",
                       "preclass_eukgap_raw_count_eggnog_abundance")
