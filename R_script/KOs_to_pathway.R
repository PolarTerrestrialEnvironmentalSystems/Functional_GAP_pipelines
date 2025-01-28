#=============================================================================#
# R-script 
# author: Ugur Cabuk
# Description: R functional kegg and taxonomy annotation of eukaryotic proteins
#=============================================================================#
rm(list = ls())

library("dplyr")
library("tidyverse")
library("tidypaleo")
library("ggplot2")
library("readr")
library("tidyr")
library("vegan")
library("data.table")

setwd("set-pathway-to-KEGG-required files")
# Reference file from KEGG Website
kegg2map <- read.delim("kegg_ko_to_map.txt",header = F)
path <- read.delim("kegg_pathway_mod_draft.tsv",header = T)
ko_name <- read.delim("ko_numbers_reference_list.tsv",row.names = NULL,header = T)

ko_name2 <- ko_name[!duplicated(ko_name$KO_IDS), ]

kegg2map2 <- kegg2map[!grepl("path:ko", kegg2map$V1),]

colnames(kegg2map2) <- c("pathway_id","KO_IDs")
colnames(path)

final_kegg_db <- merge(kegg2map2,path, by="pathway_id", all=T)

colnames(ko_name) <- c("V2","V1")
colnames(final_kegg_db)[2] <- "V2"

final_kegg_db$V2<-  gsub('ko:','', final_kegg_db$V2)

colnames(ko_name) <- c("KO_IDS","V2","V3")
colnames(final_kegg_db)[2] <- "KO_IDS"

kegg_ref_db <- merge(ko_name2,final_kegg_db, by="KO_IDS", all=T)
colnames(kegg_ref_db)[1:3] <- c("KEGG_ko","Gene_name","Description")
kegg_ref_db_2 <- kegg_ref_db[-which(is.na(kegg_ref_db$Gene_name)),]


### ANALYZE THE DATA for prokGAP ###

setwd("set-pathway-to-input-files")

# load the function
analyze_kegg_diversity <- function(eggnog_file, gene_abundance_file,metadata, output_file) {
  #eggnog input file
  
  eggnog_input <- fread(eggnog_file)
  
  # gene abundance file
  quant_merge <- fread(gene_abundance_file, header = TRUE)
  colnames(quant_merge)[1] <- "gene_id"
  
  #colnames(quant_merge) <- c("gene_id","18716","17948","6379","9157","11508","13743","22974","21740","20620","19590","3025","16946","53","712","1817","4322","5149","15596","14667","7569","6869","8491","20053","22505","260","1195","2482","17212","19192","21057","9805","12977","10402","11769","14320","15125","16079","18408","10086","13959","16535","13553","gene_length")
  #colnames(quant_merge) <- c("gene_id",18715,17948,6400,9190,11536,13746,22982,21753,20630,19593,3023,16946,76,724,1858,4381,5168,15599,14669,7562,6890,8513,20059,22518,264,1212,2504,17222,19192,21068,9816,12992,10428,11791,14020,14315,15123,16084,10091,13930,16541,13588,"gene_length")
  
  metadata <- fread(metadata)
  rename_map <- setNames(metadata$Age, metadata$Name)
  colnames(quant_merge) <- ifelse(colnames(quant_merge) %in% names(rename_map), rename_map[colnames(quant_merge)], colnames(quant_merge))
  
  age_columns <- intersect(colnames(quant_merge), rename_map)
  age_columns <- age_columns[order(as.numeric(age_columns))]
  
  #change the col name of egnnog
  colnames(eggnog_input)[1] <- "gene_id"
  quant_merge$gene_length <- NULL
  
  
  #change the col name of egnnog
  colnames(eggnog_input)[1] <- "gene_id"
  
  quant_eggnog_merge <- merge(eggnog_input,quant_merge, by="gene_id",all=T)
  quant_eggnog_merge_2 <-quant_eggnog_merge[-which(is.na(quant_eggnog_merge$eggNOG_OGs)),]
  
  # keep eukaryotes or prokaryotes
  quant_eggnog_merge_3 = quant_eggnog_merge_2[grepl("Streptophyta", quant_eggnog_merge_2$eggNOG_OGs),]
  
  #if one gene has two KO IDs, discard them.
  quant_eggnog_merge_4 <- quant_eggnog_merge_3[!grepl(",", quant_eggnog_merge_3$KEGG_ko),]
  
  quant_eggnog_merge_4 <- quant_eggnog_merge_4 %>%
    mutate(KEGG_ko= strsplit(as.character(KEGG_ko),',', perl=TRUE)) %>%
    unnest(KEGG_ko)
  quant_eggnog_merge_4$KEGG_ko<-  gsub('ko:','', quant_eggnog_merge_4$KEGG_ko)
  
  #if one gene has "-" in KO col, discard them.
  quant_eggnog_merge_4 <- quant_eggnog_merge_4[-which(quant_eggnog_merge_4$KEGG_ko == "-"),]
  
  
  which(is.na(quant_eggnog_merge_4$KEGG_ko))
  
  #wide to long format
  #quant_eggnog_merge_long <- quant_eggnog_merge_4 %>%
  #  pivot_longer(
  #    cols = `18716`:`13553`,
  #    names_to = "Ages",
  #    values_to = "TPM_value"
  #  )
  
  quant_eggnog_merge_long <- quant_eggnog_merge_4 %>%
    pivot_longer(
      cols = all_of(age_columns),
      names_to = "Ages",
      values_to = "count"
    )
  
  
  #merge kegg reference file with your prepared eggnog annotations
  kegg_ref_db_3 <- kegg_ref_db_2[,1:2]
  
  #  unfunctional line/// kegg_ref_db_4 <- kegg_ref_db_3[which(!duplicated(kegg_ref_db_3$KEGG_ko)), ]
  
  final_df <- merge(quant_eggnog_merge_long,kegg_ref_db_2, by="KEGG_ko", all=T)
  
  # IF any NAs in AGEs, will be discarded
  final_df <- final_df[-which(is.na(final_df$Ages)),]
  
  
  final_df$level_1[is.na(final_df$pathway_id)] <- "No assign level 1"
  final_df$level_2[is.na(final_df$pathway_id)] <- "No assign level 2"
  final_df$level_3[is.na(final_df$pathway_id)] <- "No assign level 3"
  
  #asıl script bu
  final_df_2 <- final_df[-which(is.na(final_df$level_2)),]
  
  #IF any NAs in KEGG level_3, will be discarded
  #final_df_2 <- final_df[-which(is.na(final_df$level_3)),]
  
  # make the ages numberic
  final_df_2$Ages = as.numeric(final_df_2$Ages)
  
  
  # These are the pathways belongs to prokaryotes
  
  prokaryote_paths <- c("Microbial metabolism in diverse environments","Degradation of aromatic compounds","Carbon fixation pathways in prokaryotes",
                        "Methane metabolism","Primary bile acid biosynthesis","Secondary bile acid biosynthesis","Steroid hormone biosynthesis",
                        "Taurine and hypotaurine metabolism","O-Antigen repeat unit biosynthesis","O-Antigen nucleotide sugar biosynthesis",
                        "Peptidoglycan biosynthesis","Teichoic acid biosynthesis","Lipoarabinomannan (LAM) biosynthesis","Exopolysaccharide biosynthesis","Pinene, camphor and geraniol degradation",
                        "Type I polyketide structures","Biosynthesis of ansamycins","Biosynthesis of enediyne antibiotics","Nonribosomal peptide structures","Polyketide sugar unit biosynthesis",
                        "Tetracycline biosynthesis","Biosynthesis of siderophore group nonribosomal peptides","Biosynthesis of vancomycin group antibiotics",
                        "Caffeine metabolism","Monobactam biosynthesis","Clavulanic acid biosynthesis","Novobiocin biosynthesis","Phenazine biosynthesis",
                        "Neomycin, kanamycin and gentamicin biosynthesis","Streptomycin biosynthesis","Acarbose and validamycin biosynthesis","Benzoate degradation",
                        "Benzoate degradation","Fluorobenzoate degradation","Chloroalkane and chloroalkene degradation","Chlorocyclohexane and chlorobenzene degradation","Toluene degradation",
                        "Xylene degradation","Nitrotoluene degradation","Ethylbenzene degradation","Styrene degradation","Atrazine degradation","Caprolactam degradation","Bisphenol degradation",
                        "Dioxin degradation","Naphthalene degradation","Polycyclic aromatic hydrocarbon degradation","Furfural degradation","Steroid degradation","Metabolism of xenobiotics by cytochrome P450",
                        "Drug metabolism - cytochrome P450","Drug metabolism - other enzymes","ABC transporters","Bacterial secretion system","Two-component system",
                        "Biofilm formation - Vibrio cholerae","Biofilm formation - Pseudomonas aeruginosa","Biofilm formation - Escherichia coli","Bacterial chemotaxis","Flagellar assembly","Regulation of actin cytoskeleton", "Quorum sensing",
                        "beta-Lactam resistance","Chemical carcinogenesis - DNA adducts","Chemical carcinogenesis - receptor activation","Cationic antimicrobial peptide","Vancomycin resistance","Lipopolysaccharide biosynthesis","Carbapenem biosynthesis",
                        "Antifolate resistance","Endocrine resistance","Platinum drug resistance","EGFR tyrosine kinase inhibitor resistance","Pathways in cancer","MicroRNAs in cancer","Transcriptional misregulation in cancer",
                        "Proteoglycans in cancer","Viral carcinogenesis","Central carbon metabolism in cancer","Choline metabolism in cancer","PD-L1 expression and PD-1 checkpoint pathway in cancer","Cell cycle - Caulobacter",
                        "Signaling pathways regulating pluripotency of stem cells","Prodigiosin biosynthesis","Penicillin and cephalosporin biosynthesis")
  
  # prokaryote level 3 paths will be discarded
  final_df_eukaryote_paths <- filter(final_df_2, !grepl(paste(prokaryote_paths, collapse="|"), level_3))
  
  final_df_eukaryote_paths$Ages = as.numeric(final_df_eukaryote_paths$Ages)
  final_df_eukaryote_paths <- final_df_eukaryote_paths[order(final_df_eukaryote_paths$Ages),]
  final_df_eukaryote_paths = final_df_eukaryote_paths[!grepl("No assign level 3", final_df_eukaryote_paths$level_3),]
  
  ko_abund <- final_df_eukaryote_paths %>%
    group_by(Ages,KEGG_ko,level_3) %>%
    summarise(SUM_COUNT = sum(count))
  
  
  #ko_abund_rel <- ko_abund %>%
  #        group_by(Ages) %>%
  #        mutate(rel_abund_fam = (SUM_TPM / sum(SUM_TPM))*100)
  
  write.csv(ko_abund, file = paste0(output_file, ".csv"))
  
}

# INPUT FILES
# Call the function for Normalized_Gene_Count
# PROKGAP
prokgap <- analyze_kegg_diversity("non_redundant_PROKGAP_protein.emapper.annotations",
                                  "numreads_sample_fixed.tsv","metadata.tsv",
                                  "prokgap_KO_pathway_CPM_eggnog_abundance")
# EUKGAP
eukgap <- analyze_kegg_diversity("non_redundant_EUKGAP_protein_eggNOG.emapper.annotations",
                                 "eukGAP_all_lake_lama_gene_quant.CPM.6m.fixed.tsv","metadata.tsv"
                                 "eukgap_KO_pathway_CPM_eggnog_abundance")

# PRECLASS_PROKGAP
preclass_prok <- analyze_kegg_diversity("non_redundant_preclass_prokGAP_protein_eggNOG.emapper.annotations",
                                        "preclass_prokGAP_all_lake_lama_gene_quant.CPM.6m.fixed.tsv","metadata.tsv"
                                        "preclass_prokgap_KO_pathway_CPM_eggnog_abundance")

# PRECLASS_EUKGAP
preclass_euk <- analyze_kegg_diversity("non_redundant_preclass_eukGAP_protein_eggNOG.emapper.annotations",
                                       "preclass_eukGAP_all_lake_lama_gene_quant.CPM.6m.fixed.tsv","metadata.tsv"
                                       "preclass_eukgap_KO_pathway_CPM_eggnog_abundance")


# INPUT FILES; If you want to work on raw count instead of normalized count in the catalog
# # RAW COUNT
# # PROKGAP
# prokgap <- analyze_kegg_diversity("non_redundant_PROKGAP_protein_eggNOG.emapper.annotations",
#                                   "prokGAP_all_lake_lama_gene_quant.raw_count.fixed.tsv","metadata.tsv"
#                                   ""prokgap_KO_pathway_raw_eggnog_abundance"")
# # EUKGAP
# eukgap <- analyze_kegg_diversity("non_redundant_EUKGAP_protein_eggNOG.emapper.annotations",
#                                  "eukGAP_all_lake_lama_gene_quant.raw_count.fixed.tsv","metadata.tsv"
#                                  "eukgap_KO_pathway_raw_eggnog_abundance")
# 
# # PRECLASS_PROKGAP
# preclass_euk <- analyze_kegg_diversity("non_redundant_preclass_prokGAP_protein_eggNOG.emapper.annotations",
#                                        "preclass_prokGAP_all_lake_lama_gene_quant.raw_count.fixed.tsv","metadata.tsv"
#                                        "preclass_prokgap_KO_pathway_raw_eggnog_abundance")
# 
# # PRECLASS_EUKGAP
# preclass_prok <- analyze_kegg_diversity("non_redundant_preclass_eukGAP_protein_eggNOG.emapper.annotations",
#                                         "preclass_eukGAP_all_lake_lama_gene_quant.raw_count.fixed.tsv","metadata.tsv"
#                                         "preclass_eukgap_KO_pathway_raw_eggnog_abundance")
