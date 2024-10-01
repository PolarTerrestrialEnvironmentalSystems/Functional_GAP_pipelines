library("plyr")
library("dplyr")
library("tidyverse")
library("ggplot2")
library("data.table")

# setwd(set-to-your-environment)

taxonomy_abund <- function(eggnog_file, gene_abundance_file, taxonomy_file, output_file) {
eggnog_input <- fread(eggnog_file)

quant_merge <- fread(gene_abundance_file, header = T)
colnames(quant_merge)[1] <- "gene_id"

# manually filled with ages
  
#colnames(quant_merge) <- c("gene_id","18716","17948","6379","9157","11508","13743","22974","21740","20620","19590","3025","16946","53","712","1817","4322","5149","15596","14667","7569","6869","8491","20053","22505","260","1195","2482","17212","19192","21057","9805","12977","10402","11769","14320","15125","16079","18408","10086","13959","16535","13553","gene_length")
colnames(quant_merge) <- c("gene_id",18715,17948,6400,9190,11536,13746,22982,21753,20630,19593,3023,16946,76,724,1858,4381,5168,15599,14669,7562,6890,8513,20059,22518,264,1212,2504,17222,19192,21068,9816,12992,10428,11791,14020,14315,15123,16084,10091,13930,16541,13588,"gene_length")
quant_merge$gene_length <- NULL

#change the col name of egnnog
colnames(eggnog_input)[1] <- "gene_id"
#merge the eggnogg annotation and read counts
quant_eggnog_merge <- merge(eggnog_input,quant_merge, by="gene_id",all=T)
quant_eggnog_merge_2= quant_eggnog_merge

#keep only eukaryotes
quant_eggnog_merge_2 <-quant_eggnog_merge[-which(is.na(quant_eggnog_merge$eggNOG_OGs)),]
#quant_eggnog_merge_2 = quant_eggnog_merge_2[grepl("Eukaryota", quant_eggnog_merge_2$eggNOG_OGs),]
gene_id_eggnog_based <- quant_eggnog_merge_2$gene_id


### TAXONOMY from NCBI ####
taxon_NCBI <- fread(taxonomy_file, header = F)
# colnames(taxon_NCBI_2)
taxon_NCBI_2 <- separate(data = taxon_NCBI, col = V5, into = c("species", "genus","family","order","class","phylum","kingdom","superkingdom"), sep = ";")
taxon_NCBI_3 <- taxon_NCBI_2 %>% select(c("V1", "V4", "species", "genus","family","order","class","phylum","kingdom","superkingdom"))
taxon_NCBI_3$species <- as.character(taxon_NCBI_3$species)
taxon_NCBI_3[taxon_NCBI_3 == ""] <- "Unclassified"
taxon_NCBI_3$species <- sub("^uc_.*|unknown", "Unclassified", taxon_NCBI_3$species)
taxon_NCBI_3$genus <- sub("^uc_.*|unknown", "Unclassified", taxon_NCBI_3$genus)
taxon_NCBI_3$family <- sub("^uc_.*|unknown", "Unclassified", taxon_NCBI_3$family)
taxon_NCBI_3$order <- sub("^uc_.*|unknown", "Unclassified", taxon_NCBI_3$order)
taxon_NCBI_3$class <- sub("^uc_.*|unknown", "Unclassified", taxon_NCBI_3$class)
taxon_NCBI_3$phylum <- sub("^uc_.*|unknown", "Unclassified", taxon_NCBI_3$phylum)
taxon_NCBI_3$kingdom <- sub("^uc_.*|unknown", "Unclassified", taxon_NCBI_3$kingdom)
taxon_NCBI_3$superkingdom <- sub("^uc_.*|unknown", "Eukaryota", taxon_NCBI_3$superkingdom)
taxon_NCBI_4 <- replace(taxon_NCBI_3, is.na(taxon_NCBI_3), "Unclassified")
taxon_NCBI_4$superkingdom <- sub("^Unclassif.*", "Eukaryota", taxon_NCBI_4$superkingdom)
colnames(taxon_NCBI_4)[1] <- "gene_id"

quant_eggnog_merge_taxon <- merge(quant_eggnog_merge_2,taxon_NCBI_4, by="gene_id",all.x=T)
quant_eggnog_merge_taxon <- quant_eggnog_merge_taxon[grepl("Streptophyta", quant_eggnog_merge_taxon$phylum),]
quant_eggnog_merge_3 = quant_eggnog_merge_taxon


#if one gene has two KO IDs, discard them.
quant_eggnog_merge_4 <- quant_eggnog_merge_3[!grepl(",", quant_eggnog_merge_3$KEGG_ko),]
quant_eggnog_merge_4 <- quant_eggnog_merge_4 %>%
  mutate(KEGG_ko= strsplit(as.character(KEGG_ko),',', perl=TRUE)) %>%
  unnest(KEGG_ko)
quant_eggnog_merge_4$KEGG_ko<-  gsub('ko:','', quant_eggnog_merge_4$KEGG_ko)

#if one gene has "-" in KO col, discard them.
quant_eggnog_merge_4 <- quant_eggnog_merge_4[-which(quant_eggnog_merge_4$KEGG_ko == "-"),]

which(is.na(quant_eggnog_merge_4$KEGG_ko))

quant_eggnog_merge_long <- quant_eggnog_merge_4 %>%
  pivot_longer(
    cols = `18716`:`13553`,
    names_to = "Ages",
    values_to = "TPM_value"
  )

# ALL gene
all_quant_eggnog_merge_long <- quant_eggnog_merge_3 %>%
  pivot_longer(
    cols = `18716`:`13553`,
    names_to = "Ages",
    values_to = "TPM_value"
  )

# ALL gene taxonomy family level
all_gene_abund <- all_quant_eggnog_merge_long %>%
  group_by(Ages,family) %>%
  summarise(SUM_TPM = sum(TPM_value))

gene_abund_rel <- all_gene_abund %>% group_by(Ages) %>% mutate(rel_abund_family = (SUM_TPM / sum(SUM_TPM))*100)
write.csv(gene_abund_rel, file = paste0(output_file, ".csv"))
}


taxonomy_abund("out.eggnog/non_redundant_PROKGAP_protein_eggNOG.emapper.annotations",
                       "prokGAP_all_lake_lama_gene_quant.CPM.fixed.tsv","eggnog_eukaryote_OG_proteins_prokGAP_taxonomy.tsv"
                       "plant_family_prokgap)


taxonomy_abund("out.eggnog/non_redundant_EUKGAP_protein_eggNOG.emapper.annotations",
                       "eukGAP_all_lake_lama_gene_quant.CPM.fixed.tsv","eggnog_eukaryote_OG_proteins_eukGAP_taxonomy.tsv"
                       "plant_family_eukgap)


taxonomy_abund("out.eggnog/non_redundant_preclass_prokGAP_protein_eggNOG.emapper.annotations",
                       "preclass_prokGAP_all_lake_lama_gene_quant.CPM.fixed.tsv","eggnog_eukaryote_OG_proteins_preclass_prokGAP_taxonomy.tsv"
                       "plant_family_preclass_prokgap)

taxonomy_abund("out.eggnog/non_redundant_preclass_eukGAP_protein_eggNOG.emapper.annotations",
                       "preclass_eukGAP_all_lake_lama_gene_quant.CPM.fixed.tsv","eggnog_eukaryote_OG_proteins_preclass_eukGAP_taxonomy.tsv"
                       "plant_family_preclass_eukgap)
