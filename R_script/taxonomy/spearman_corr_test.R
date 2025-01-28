library(corrplot)
library(psych)
library(ggplot2)
# getwd()
setwd("")
# Read the CSV files
prokgap <- read.csv("plant_family_raw_count.csv")
prokgap$approach <- "ProkGAP"
prokgap$X <- NULL

eukgap <- read.csv("plant_family_eukgap_raw_count.csv")
eukgap$approach <- "EukGAP"
eukgap$X <- NULL

preclass_prokgap <- read.csv("plant_family_preclass_prokgap_raw_count.csv")
preclass_prokgap$approach <- "Pre-class ProkGAP"
preclass_prokgap$X <- NULL

preclass_eukgap <- read.csv("plant_family_preclass_eukgap_raw_count.csv")
preclass_eukgap$approach <- "Pre-class EukGAP"
preclass_eukgap$X <- NULL

# needs to be generated from kraken taxonomy.
kraken <- read.table("kraken_family_abundance.tsv")
kraken$approach <- "Kraken"

# Rename columns
colnames(kraken) <- c("family", "Ages", "SUM_COUNT", "rel_abund_family", "approach")

# Combine the data frames
combined_data <- rbind(prokgap,eukgap,preclass_prokgap,preclass_eukgap,kraken)
combined_data <- rbind(combined_data)

# Selected taxa
selected_taxa <- c("Rosaceae", "Ranunculaceae", "Poaceae", "Pinaceae", "Salicaceae","Betulaceae", "Ericaceae", "Brassicaceae", "Asteraceae")

tax_abund <- combined_data %>%
  group_by(Ages,family,approach) %>%
  summarise(SUM_COUNT = sum(SUM_COUNT))

tax_abund_rel <- tax_abund %>%
  group_by(Ages,approach) %>%
  mutate(rel_abund_fam = (SUM_COUNT / sum(SUM_COUNT))*100)

tax_abund_rel$rel_abund_fam <- as.numeric(tax_abund_rel$rel_abund_fam)
tax_abund_rel$family <- factor(tax_abund_rel$family, levels = unique(tax_abund_rel$family))

tax_abund_rel_fil <- tax_abund_rel[tax_abund_rel$family %in% selected_taxa, ]
tax_abund_rel = tax_abund_rel_fil

# Reorder levels of the approach factor
tax_abund_rel$approach <- factor(tax_abund_rel$approach, levels = c("Kraken", "ProkGAP", "EukGAP", "Pre-class ProkGAP","Pre-class EukGAP"))  # Add other approaches as needed
# Adjusting facet_wrap and y-axis label
tax_abund_rel_subset <- tax_abund_rel %>%
  filter(approach %in% c("Kraken", "ProkGAP"))

# Filter out rows where rel_abund_family is 0
tax_abund_rel_subset <- tax_abund_rel_subset %>%
  filter(rel_abund_fam != 0)

# only prokgap vs kraken
subset_prokgap_kraken <- subset(tax_abund_rel, approach %in% c("Kraken", "ProkGAP"))

new_df <- subset_prokgap_kraken %>% 
  filter(approach %in% c("ProkGAP", "Kraken")) %>% 
  select(Ages, family, approach,rel_abund_fam)

tax_abund_wide <- reshape2::dcast(new_df, Ages ~ family + approach, value.var = "rel_abund_fam")

# Replace NA values with 0
tax_abund_wide[is.na(tax_abund_wide)] <- 0

colnames(tax_abund_wide) <- gsub("_", " ", colnames(tax_abund_wide))


correlation_matrix <- cor(tax_abund_wide[, -1], method ="spearman")  # Excluding the first column (Ages)

# p-values for correlations
cor_pmat <- cor.mtest(tax_abund_wide[, -1])$p  # Excluding the first column (Ages)
View(cor_pmat)

# Order variables by family names
variable_order <- order(colnames(tax_abund_wide)[-1])
correlation_matrix_ordered <- correlation_matrix[variable_order, variable_order]
cor_pmat_ordered <- cor_pmat[variable_order, variable_order]

png("correlation_matrix_new8.png", width = 3200, height = 3200, units = "px", res = 300)
corrplot(correlation_matrix_ordered, method = "number", type = "upper", order = "original", addrect = 8, tl.col = "black", insig='pch', tl.srt = 45, diag = FALSE, mar=c(0,0,2,0)+0.1)
dev.off()
