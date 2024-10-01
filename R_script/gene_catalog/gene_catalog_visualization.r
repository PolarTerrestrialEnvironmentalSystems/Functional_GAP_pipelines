setwd("set-to-the-path")

library(patchwork)
library(egg)
library(ggh4x)
lama_abundance_stats <- fread("gene_catalog_report.csv", header = TRUE)

lama_abundance_stats_df <- as.data.frame(lama_abundance_stats)

# Normalization Visualization

long_CPM <- melt(setDT(lama_abundance_stats_df), id.vars = c("Sample","Ages","Total_reads","Description",
                                                         "ProkGAP","EukGAP",
                                                         "pre-class EukGAP",
                                                         "pre-class ProkGAP"), variable.name = "Values")

# Check whether all samples are in

length(unique(long_CPM$Ages))

# need to be reorder


neworder <- c("Non_redundant_gene_catalog_abundance_no_annotation",
              "eggNOG annotation","Eukaryota_OGs_Annotation","Prokaryota_OGs_Annotation","Eukaryota_KOs_Annotation","Prokaryota_KOs_Annotation")

long2 <- arrange(transform(long_CPM,
                           Description=factor(Description,levels=neworder)),Description)
						   
						   
# rename the column names

long2$Description <- as.character(long2$Description)
long2$Description[which(long2$Description == "Non_redundant_gene_catalog_abundance_no_annotation")] <- "Non redundant gene catalog"
long2$Description[which(long2$Description == "Prokaryota_OGs_Annotation")] <- "Prokaryotic gene subset"
long2$Description[which(long2$Description == "Eukaryota_OGs_Annotation")] <- "Eukaryotic gene subset "
long2$Description[which(long2$Description == "Prokaryota_KOs_Annotation")] <- "Prokaryotic KOs"
long2$Description[which(long2$Description == "Eukaryota_KOs_Annotation")] <- "Eukaryotic KOs"

long2$Description <- factor(long2$Description, levels = unique(long2$Description))


# Here two things happen; column names fixed and added another column to separate Ages.

long2_epoch <- long2 %>%
  mutate(epoch = factor(
    ifelse(Ages >= 53 & Ages <= 13959, "53 to 13959", 
           ifelse(Ages >= 14220 & Ages <= 24000, "14320 to 22974", NA)),
    levels = c("53 to 13959", "14320 to 22974")
  )) %>%mutate(Values = gsub("_CPM", "", Values)) %>% 
  # Remove "_CPM" from the Values colum
  mutate(Values = gsub("preclass_ProkGAP", "pre-class ProkGAP", Values)) %>%
  mutate(Values = gsub("preclass_EukGAP", "pre-class EukGAP", Values))
  
# Order


long2_epoch$Values <- factor(long2_epoch$Values, 
                       levels = c("ProkGAP", "EukGAP", "pre-class ProkGAP", "pre-class EukGAP"))
					   


# Multipy with number of genes in the catalog

long2_epoch_filtered <- long2_epoch %>%
  mutate(value = case_when(
    Values == "ProkGAP" ~ value * 6568483 / 1000000,
    Values == "EukGAP" ~ value * 4631120/ 1000000,
    Values == "pre-class ProkGAP" ~ value * 5930831/ 1000000,
    Values == "pre-class EukGAP" ~ value * 380105/ 1000000,
    # Add more conditions for other values if needed
    TRUE ~ value  # for values not matching any condition, keep them unchanged
  ))
  
  
# pre-defined custom them

custom_theme <- theme_minimal() +
  theme(
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 10),
    legend.position = "bottom",
    legend.key.size = unit(1.2, "lines"),  # Adjust legend key size
    legend.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt"),  # Adjust legend margin
    legend.box = "none",axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),
    # legend.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 11), # Adjust strip text size, face, and color
    strip.placement = "outside"  # Move strip label to to
  )
  

long2_epoch_filtered$Values <- as.character(long2_epoch_filtered$Values)

long2_epoch_filtered$Values[which(long2_epoch_filtered$Values == "pre-class ProkGAP")] <- "Pre-class ProkGAP"
long2_epoch_filtered$Values[which(long2_epoch_filtered$Values == "pre-class EukGAP")] <- "Pre-class EukGAP"

long2_epoch_filtered$Values <- factor(long2_epoch_filtered$Values)
# again ordered that you want
new_order <- c("ProkGAP", "EukGAP", "Pre-class ProkGAP","Pre-class EukGAP")
long2_epoch_filtered$Values <- factor(long2_epoch_filtered$Values, levels = new_order)

# Color palette
palette <- c("ProkGAP" = "dodgerblue3", 
             "EukGAP" = "gold", 
             "Pre-class ProkGAP" = "darkorchid3", 
             "Pre-class EukGAP" = "forestgreen")

plot_data <- ggplot(long2_epoch_filtered, aes(x = factor(Ages), y = value, group = Values, fill = Values)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  facet_grid2(Description~., scales = "free", space = "free_x", switch = "x",labeller = label_wrap_gen(width = 10),
              strip = strip_themed(
                size = "variable", # Shrinks 2nd layer
                text_x = elem_list_text(face = c(1, 2)), # 2nd layer in bold
                by_layer_x = TRUE
              ))+
  labs(x = "Calibrated years BP", y = "Normalized gene/protein count", fill = "prediction pipelines") +
  # theme() +  # Rotate x-axis labels # Rotate x-axis labels
  scale_fill_manual(values = palette) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  expand_limits(y = 0) + custom_theme  # Apply the custom theme

plot_data

ggsave("gene_abundance_catalog_comparison.png", plot_data, width = 9, height = 9, units = "in", dpi = 300)
ggsave("gene_abundance_catalog_comparison.pdf", plot_data, width = 9, height = 9, units = "in", dpi = 300)


# Raw Count Visualization

long_count <- melt(setDT(lama_abundance_stats_df), id.vars = c("Sample","Ages","Total_reads","Description",
                                                         "ProkGAP_CPM","EukGAP_CPM",
                                                         "preclass_EukGAP_CPM",
                                                         "preclass_ProkGAP_CPM"), variable.name = "Values")


# Convert Description to character
long_count$Description <- as.character(long_count$Description)
long_count$Description[which(long_count$Description == "Non_redundant_gene_catalog_abundance_no_annotation")] <- "Non redundant gene catalog"
long_count$Description[which(long_count$Description == "Prokaryota_OGs_Annotation")] <- "Prokaryotic gene subset"
long_count$Description[which(long_count$Description == "Eukaryota_OGs_Annotation")] <- "Eukaryotic gene subset "
long_count$Description[which(long_count$Description == "Prokaryota_KOs_Annotation")] <- "Prokaryotic KOs"
long_count$Description[which(long_count$Description == "Eukaryota_KOs_Annotation")] <- "Eukaryotic KOs"

# Convert Description back to factor with desired levels
long_count$Description <- factor(long_count$Description, levels = unique(long_count$Description))

long_count$Values <- as.character(long_count$Values)

long_count$Values[which(long_count$Values == "pre-class ProkGAP")] <- "Pre-class ProkGAP"
long_count$Values[which(long_count$Values == "pre-class EukGAP")] <- "Pre-class EukGAP"

long_count$Values <- factor(long_count$Values)

new_order <- c("ProkGAP", "EukGAP", "Pre-class ProkGAP","Pre-class EukGAP")

long_count$Values <- factor(long_count$Values, levels = new_order)

palette <- c("ProkGAP" = "dodgerblue3", 
             "EukGAP" = "gold", 
             "Pre-class ProkGAP" = "darkorchid3", 
             "Pre-class EukGAP" = "forestgreen")

plot_count <- ggplot(long_count, aes(x = factor(Ages), y = value, group = Values, fill = Values)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  facet_grid2(Description~., scales = "free", space = "free_x", switch = "x",labeller = label_wrap_gen(width = 10),
              strip = strip_themed(
                size = "variable", # Shrinks 2nd layer
                text_x = elem_list_text(face = c(1, 2)), # 2nd layer in bold
                by_layer_x = TRUE
              ))+
  labs(x = "Calibrated years BP", y = "Normalized gene/protein count", fill = "prediction pipelines") +
  
  # theme() +  # Rotate x-axis labels # Rotate x-axis labels
  
  scale_fill_manual(values = palette) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  expand_limits(y = 0) + custom_theme  # Apply the custom theme
plot_count
ggsave("gene_abundance_raw_pseudocount_comparison.png", plot_count, width = 9, height = 9, units = "in", dpi = 300)
ggsave("gene_abundance_raw_pseudocount_comparison.png", plot_count, width = 9, height = 9, units = "in", dpi = 300)


### BOX PLOT ###

subset_data <- subset(long2_epoch_filtered, Description %in% c("eggNOG annotation", "Eukaryotic gene subset ", "Prokaryotic gene subset", "Eukaryotic KOs", "Prokaryotic KOs"))

# Define a custom theme with modifications
custom_theme <- theme_minimal() +
  theme(
    # strip.text = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()# Remove legend box
  )


long2_epoch <- subset_data %>%
  mutate(epoch = factor(
    ifelse(Ages >= 76 & Ages <= 11700, "Holocene", 
           ifelse(Ages > 11700   & Ages <= 24000, "Glacial", NA)),
    levels = c("Holocene", "Glacial")
  ))
subset_data = long2_epoch
# Assuming 'subset_data' is your data frame
palette <- c("Holocene" = "red2", 
             "Glacial" = "lightblue")
library(ggpubr)

# Filter data for eggNOG annotation
eggnog_data <- subset_data[subset_data$Description == "eggNOG annotation", ]
p_eggnog_data <- ggplot(long2_epoch, aes(x = factor(epoch), y = value,fill=epoch,group=epoch)) + 
  geom_boxplot() +
  labs(x = "", y = "Normalized Gene Abundance") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + custom_theme + theme(legend.position = "none") +
  ggtitle("") +
  facet_nested_wrap(vars(Values,Description), scales = "free_y")+
  # facet_wrap2(Values~Description, scales = "free") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  stat_compare_means(
    method = "t.test",
    # label = "p.signif", # Increase label size for better visibility
    # label = "p.signif", size=5,  # Use 'p.signif' for significant markers
    comparisons = list(c("Holocene", "Glacial")),  # Specify groups for comparison
    hide.ns = F,  # Hide non-significant comparisons
    vjust = 1.2,
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", ""))# Customize significance symbols
  )


# Filter data for Eukaryotic gene subset
euk_data <- subset_data[subset_data$Description == "Eukaryotic gene subset ", ]
max(euk_data$value)
p_euk_data <- ggplot(euk_data, aes(x = factor(epoch), y = value,fill=epoch,group=epoch)) + 
  geom_boxplot() +
  labs(x = "", y = "") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))  + custom_theme + theme(legend.position = "none") +
  ggtitle("Eukaryotic Gene Subset")  +
  facet_grid(~Values, scales = "free", space = "free")+
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  stat_compare_means(
    method = "t.test",
    # label = "p.signif", size=5,  # Use 'p.signif' for significant markers
    comparisons = list(c("Holocene", "Glacial")),  # Specify groups for comparison
    hide.ns = F,  # Hide non-significant comparisons
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", ""))# Customize significance symbols
  )


# Filter data for Prokaryotic gene subset
prok_data <- subset_data[subset_data$Description == "Prokaryotic gene subset", ]
p_prok_data <- ggplot(prok_data, aes(x = factor(epoch), y = value,fill=epoch,group=epoch)) + 
  geom_boxplot() +
  labs(x = "", y = "") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))  + custom_theme + theme(legend.position = "none") +
  ggtitle("Prokaryotic Gene Subset")  +
  facet_grid(~Values, scales = "free", space = "free")+
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  stat_compare_means(
    method = "t.test",
    # label = "p.signif", size=5,  # Use 'p.signif' for significant markers
    comparisons = list(c("Holocene", "Glacial")),  # Specify groups for comparison
    hide.ns = F,  # Hide non-significant comparisons
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", ""))# Customize significance symbols
  )


# Filter data for Eukaryotic KOs
euk_kos_data <- subset_data[subset_data$Description == "Eukaryotic KOs", ]
p_euk_kos_data <- ggplot(euk_kos_data, aes(x = factor(epoch), y = value,fill=epoch,group=epoch)) + 
  geom_boxplot() +
  labs(x = "", y = "") +
  scale_fill_manual(values = palette) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))  + custom_theme + 
  ggtitle("Eukaryotic KOs") +
  facet_grid(~Values, scales = "free", space = "free")+
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  stat_compare_means(
    method = "t.test",
    # label = "p.signif", size=5,  # Use 'p.signif' for significant markers
    comparisons = list(c("Holocene", "Glacial")),  # Specify groups for comparison
    hide.ns = F,  # Hide non-significant comparisons
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", ""))# Customize significance symbols
  )


# Filter data for Prokaryotic KOs
prok_kos_data <- subset_data[subset_data$Description == "Prokaryotic KOs", ]
p_prok_kos_data <- ggplot(subset(prok_kos_data), aes(x = factor(epoch), y = value,fill=epoch,group=epoch)) + 
  geom_boxplot() +
  labs(x = "", y = "") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))  + custom_theme + 
  ggtitle("Prokaryotic KOs") +
  # facet_wrap(~ Values, scales = "free_y", nrow = 1) +  # Organize by 'Values' and set number of rows
  facet_grid(~Values, scales = "free", space = "free")+
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  stat_compare_means(
    method = "t.test",
    # label = "p.signif", size=5,  # Use 'p.signif' for significant markers
    comparisons = list(c("Holocene", "Glacial")),  # Specify groups for comparison
    hide.ns = F,  # Hide non-significant comparisons
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", ""))# Customize significance symbols
  )



composite_plot2 <- p_eggnog_data + p_euk_data + p_prok_data + p_euk_kos_data + p_prok_kos_data
all_merged <- composite_plot2 +plot_layout(guides = 'collect',nrow = 1) & theme(legend.position = 'bottom') 

print(all_merged)

ggsave("gene_abundance_plot_v4.png", p_eggnog_data, width = 10, height = 6, units = "in", dpi = 300)




### all dataset plot ###

p_all_eggnog_data<- ggplot(long2_epoch, aes(x = "", y = value ,fill=Values)) + 
  geom_boxplot(fatten=3) +
  facet_wrap(~Description, ncol=6,scale= "free_y",labeller = label_wrap_gen(width = 10)
  )+
  theme_minimal(base_size = 15) +
  labs(x = "", y = "Normalized gene/protein count", fill = "prediction pipelines") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  # Rotate x-axis labels and adjust size
        axis.text.y = element_text(size = 14),  # Adjust size of y-axis labels
        axis.title = element_text(size = 14),  # Adjust size of axis titles
        strip.text = element_text(size = 14),  # Adjust size of facet labels
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        legend.position = "bottom") + # Remove legend
  scale_fill_manual(values = pipeline_palette) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  # scale_y_continuous(labels = scales::comma_format(scale = 1e-3)) +  # Display Y-axis labels in K (thousands)
  expand_limits(y = 0)  # En
# p_contig_data

ggsave("assembly_box_plot_publication_plot_v8.png", p_all_eggnog_data, width = 12, height = 6, units = "in", dpi = 300)
