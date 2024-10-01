setwd("")
library(data.table)
library(ggvenn)
library(gridExtra)
# Eukaryota
euk_prokgap <- fread("prokGAP_final_df_row_counts_lama.csv")
euk_eukgap <- fread("eukGAP_final_df_row_KO_counts_lama.csv")
euk_preclass_prokgap <- fread("preclass_prokGAP_final_df_row_KO_counts_lama.csv")
euk_preclass_eukgap <- fread("preclass_eukGAP_final_df_row_KO_counts_lama.csv")

uniq_prokgap <- data.frame(unique(euk_prokgap$KEGG_ko))
uniq_prokgap$approach <- "prokgap"
colnames(uniq_prokgap)[1] <- "KEGG_Ko"

uniq_eukgap <- data.frame(unique(euk_eukgap$KEGG_ko))
uniq_eukgap$approach <- "eukgap"
colnames(uniq_eukgap)[1] <- "KEGG_Ko"

uniq_preclass_prokgap <- data.frame(unique(euk_preclass_prokgap$KEGG_ko))
uniq_preclass_prokgap$approach <- "preclass_prokgrap"
colnames(uniq_preclass_prokgap)[1] <- "KEGG_Ko"

uniq_preclass_eukgap <- data.frame(unique(euk_preclass_eukgap$KEGG_ko))
uniq_preclass_eukgap$approach <- "preclass_eukgap"
colnames(uniq_preclass_eukgap)[1] <- "KEGG_Ko"

# Combine all unique KEGG_ko data
euk_all_data <- rbind(uniq_preclass_eukgap, uniq_preclass_prokgap, uniq_prokgap, uniq_eukgap)
# count it
length(unique(euk_all_data$KEGG_Ko))

EUK = list(EukGAP = uniq_eukgap$KEGG_Ko, ProkGAP = uniq_prokgap$KEGG_Ko, 'Pre-class EukGAP'=uniq_preclass_eukgap$KEGG_Ko, 'Pre-class ProkGAP'= uniq_preclass_prokgap$KEGG_Ko)

# Prokaryotes 

# Eukaryota
prok_prokgap <- fread("prokgap_row_KO_count_prokaryotes.csv")
prok_eukgap <- fread("eukgap_row_KO_count_prokaryotes.csv")
prok_preclass_prokgap <- fread("preclass_prokgap_row_KO_count_prokaryotes.csv")
prok_preclass_eukgap <- fread("preclass_eukgap_row_KO_count_prokaryotes.csv")

uniq_prokgap <- data.frame(unique(prok_prokgap$KEGG_ko))
uniq_prokgap$approach <- "prokgap"
colnames(uniq_prokgap)[1] <- "KEGG_Ko"

uniq_eukgap <- data.frame(unique(prok_eukgap$KEGG_ko))
uniq_eukgap$approach <- "eukgap"
colnames(uniq_eukgap)[1] <- "KEGG_Ko"

uniq_preclass_prokgap <- data.frame(unique(prok_preclass_prokgap$KEGG_ko))
uniq_preclass_prokgap$approach <- "preclass_prokgrap"
colnames(uniq_preclass_prokgap)[1] <- "KEGG_Ko"

uniq_preclass_eukgap <- data.frame(unique(prok_preclass_eukgap$KEGG_ko))
uniq_preclass_eukgap$approach <- "preclass_eukgap"
colnames(uniq_preclass_eukgap)[1] <- "KEGG_Ko"

# Combine all unique KEGG_ko data
prok_all_data <- rbind(uniq_preclass_eukgap, uniq_preclass_prokgap, uniq_prokgap, uniq_eukgap)
# count it
length(unique(prok_all_data$KEGG_Ko))

PROK = list(EukGAP = uniq_eukgap$KEGG_Ko, ProkGAP = uniq_prokgap$KEGG_Ko, 'Pre-class EukGAP'=uniq_preclass_eukgap$KEGG_Ko, 'Pre-class ProkGAP'= uniq_preclass_prokgap$KEGG_Ko)


# color for legend

custom_palette <- c("gold","dodgerblue3","darkorchid3","forestgreen")

euk_venn <- ggvenn(
  EUK, 
  fill_color = custom_palette,
  stroke_size = 0, set_name_size = 3
) + scale_x_continuous(expand = expansion(mult = .2)) + labs(title = "Eukaryotic KOs") + theme(plot.title = element_text(hjust = 0.5))

prok_venn <- ggvenn(
  PROK, 
  fill_color = custom_palette,
  stroke_size = 0, set_name_size = 3
) + scale_x_continuous(expand = expansion(mult = .2)) + labs(title = "Prokaryotic KOs") + theme(plot.title = element_text(hjust = 0.5))

# Extract approach names
approach_names <- c(rep("EukGAP", length(EUK$EukGAP)),
                    rep("ProkGAP", length(EUK$ProkGAP)),
                    rep("Pre-class EukGAP", length(EUK$`Pre-class EukGAP`)),
                    rep("Pre-class ProkGAP", length(EUK$`Pre-class ProkGAP`)))
# Combine all lists into one column
combined_list <- c(EUK$EukGAP, EUK$ProkGAP, EUK$`Pre-class EukGAP`, EUK$`Pre-class ProkGAP`)

# Create a dataframe with two columns
euk_combined_df <- data.frame(Gene = combined_list, Approach = approach_names)

# View the dataframe
print(euk_combined_df)
uniq_all_pipeline_euk <- unique(combined_df$Gene)
length(uniq_all_pipeline_euk)

approach_names <- c(rep("EukGAP", length(PROK$EukGAP)),
                    rep("ProkGAP", length(PROK$ProkGAP)),
                    rep("Pre-class EukGAP", length(PROK$`Pre-class EukGAP`)),
                    rep("Pre-class ProkGAP", length(PROK$`Pre-class ProkGAP`)))
# Combine all lists into one column
combined_list <- c(PROK$EukGAP, PROK$ProkGAP, PROK$`Pre-class EukGAP`, PROK$`Pre-class ProkGAP`)

# Create a dataframe with two columns
prok_combined_df <- data.frame(Gene = combined_list, Approach = approach_names)

# View the dataframe
print(prok_combined_df)

uniq_all_pipeline_euk <- unique(euk_combined_df$Gene)
uniq_all_pipeline_prok <- unique(prok_combined_df$Gene)

length(uniq_all_pipeline_prok)

mixed_euk_prok = list("Eukaryotic KOs" = uniq_all_pipeline_euk, "Prokaryotic KOs" = uniq_all_pipeline_prok)

mixed_venn <- ggvenn(
  mixed_euk_prok, 
  fill_color = custom_palette,
  stroke_size = 0, set_name_size = 4
) + scale_x_continuous(expand = expansion(mult = .2)) + labs(title = "")

# Merge all plots
# Adjust the figures
heights <- c(2, 2, 2) 
scale_factor_euk_prok <- 2.5
layout_matrix <- rbind(
  c(1, 2),     
  c(3, 3)       
)

merged_plot <- arrangeGrob(
  euk_venn, prok_venn, mixed_venn,
  ncol = 2,
  layout_matrix = layout_matrix
  # heights = heights * c(scale_factor_euk_prok, scale_factor_euk_prok, 1)
)

# Display the merged plot
merged_plot_2 <- grid.arrange(merged_plot)
