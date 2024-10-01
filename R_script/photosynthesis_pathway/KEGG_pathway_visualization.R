# Bar plot and heatmap

setwd("set-to-the-path")

prokgap <- fread("prokgap_KO_pathway_CPM_eggnog_abundance.csv")
eukgap <- fread("eukgap_KO_pathway_CPM_eggnog_abundance.csv")
preclass_eukgap  <- fread("preclass_eukgap_KO_pathway_CPM_eggnog_abundance.csv")
preclass_prokgap <- fread("preclass_prokgap_KO_pathway_CPM_eggnog_abundance.csv")

# subset
prokgap$approach <- "ProkGAP"
eukgap$approach <- "EukGAP"
preclass_prokgap$approach <- "Pre-class ProkGAP"
preclass_eukgap$approach <- "Pre-class EukGAP"

all_photosynthesis <- rbind(prokgap,eukgap,
                            preclass_prokgap,preclass_eukgap)

ko_abund <- all_photosynthesis %>%
  group_by(Ages,level_3,approach, KEGG_ko) %>%
  summarise(SUM_TPM = sum(SUM_TPM))

# custom palette and levels for approach
palette <- c("ProkGAP" = "dodgerblue3", 
             "EukGAP" = "gold", 
             "Pre-class ProkGAP" = "darkorchid3", 
             "Pre-class EukGAP" = "forestgreen")

approach_levels <- c("ProkGAP", "EukGAP", "Pre-class ProkGAP", "Pre-class EukGAP")

photo_subset=subset(ko_abund, level_3 == "Photosynthesis")

photo_subset_filtered <- photo_subset %>% group_by(Ages) %>%
  mutate(value = case_when(
    approach == "ProkGAP" ~ SUM_TPM,
    approach == "EukGAP" ~ SUM_TPM,
    approach == "Pre-class ProkGAP" ~ SUM_TPM,
    approach == "Pre-class EukGAP" ~ SUM_TPM,
    TRUE ~ SUM_TPM  # Keep original value for other cases
  ))

x=subset(photo_subset_filtered, level_3 == "Photosynthesis")

photo_subset <- photo_subset_filtered %>%
  group_by(Ages,level_3,approach) %>%
  summarise(SUM_TPM = sum(SUM_TPM))

####
new_age_depth <- as.data.frame(unique(photo_subset$Ages))
colnames(new_age_depth) <- "old_age"
photo_subset2 = merge(photo_subset, Ages, by= "Ages")
####

# Filter data for photosynthesis level 3
photosynthesis <- ggplot(photo_subset2, aes(x = factor(New_Ages), y = log(SUM_TPM), fill = approach)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(x = "Calibrated years BP", y = "log(Normalized Gene Count)", fill = "Approach", title = "") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "top",
    # legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 15, hjust = 0.5),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm")  # Adjust legend key size
  ) +
  scale_fill_manual(values = palette, limits = approach_levels) +  # Set palette and limits
  guides(fill = guide_legend(title = "prediction pipelines", position = "bottom"))  # Customize legend title and position

# Print the plot
print(photosynthesis)
ggsave("streptophyta_kegg3_photosynthesis.png", photosynthesis, width = 8, height = 5, dpi = 300)

### KOs in the pathway ###
photo_subset <- photo_subset_filtered %>%
  group_by(Ages,level_3,approach, KEGG_ko) %>%
  summarise(SUM_TPM = sum(SUM_TPM))
unique(photo_subset$Ages)

new_age_depth <- as.data.frame(unique(photo_subset$Ages))
colnames(new_age_depth) <- "old_age"
photo_subset2 = merge(photo_subset, Ages, by= "Ages")

new_order <- c("ProkGAP", "EukGAP", "Pre-class ProkGAP","Pre-class EukGAP")

photo_subset2$approach <- factor(photo_subset2$approach, levels = new_order)
subset_test <- subset(photo_subset2,approach == "Pre-class ProkGAP")
length(unique(photo_subset2$KEGG_ko))
heatmap <- ggplot(photo_subset2,aes(x = factor(New_Ages), y = KEGG_ko , fill = log(SUM_TPM))) +
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "steelblue", na.value = "white", name = "log(Normalized gene count)") + # Include color for 0 values
  labs(x="Calibrated years BP",y = "Photosynthesis KOs", title = "") +
  facet_wrap(~approach, ncol = 2) +
  theme_minimal(base_size = 5) +
  custom_theme

print(heatmap)
ggsave("Photosynthesis_KOs_heatmap_streptophyta.png", plot = heatmap, width = 9, height = 9, units = "in", dpi = 300)
