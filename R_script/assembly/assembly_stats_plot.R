library(ggpubr)
library(patchwork)
library(ggh4x)
library(plyr)
library(dplyr)
library(tidyverse)
### quast report ###

setwd("")

lama_quast_report <- fread("lake_lama_quast_report.csv", header = TRUE)

# Define the relevant metrics
metrics <- c("Contig Number", "N50", "Largest Contig", "Raw Data read counts", "Filtered and error corrected total read count")

# Extract sample names
samples <- lama_quast_report$Sample_name

# Reshape the data
plot_data <- melt(lama_quast_report, id.vars = "AGE", measure.vars = metrics)

# Rename columns to remove spaces and make them consistent
colnames(plot_data) <- c("Age", "Metric", "Value")

# Order the facets
plot_data$Metric <- factor(plot_data$Metric, levels = c(
  "Raw Data read counts",
  "Total read count after deduplication",
  "Filtered and error corrected total read count",
  "Contig Number",
  "Largest Contig",
  "N50"
))

# Define units for each metric
units <- c("count", "bp", "bp", "bp", "bp")

# Combine metrics and units

# Combine metrics and units
metric_labels <- paste(metrics, "(", units, ")", sep = " ")

# Custom labeller function to add metric labels for each facet
my_labeller <- function(variable, value) {
  return(paste(value, metric_labels[as.numeric(value)], sep = " "))
}

# Create a function to generate a gradient of grey colors
generate_grey_gradient <- function(num_colors) {
  return(colorRampPalette(c("black", "grey90"))(num_colors))
}

# Define your data
# Assuming plot_data contains your data

# Define the number of metrics
num_metrics <- 5

# Generate a gradient of grey colors
metric_colors <- generate_grey_gradient(num_metrics)


# Function to add units to facet labels based on Metric and specify order
add_units_to_facet_labels <- function(metric) {
  if (metric == "Filtered and error corrected total read count") {
    return("Preprocessed read count")
  } else if (metric == "Contig Number") {
    return("Contig number count")
  } else if (metric == "Largest Contig") {
    return("Largest contig (bp)")
  } else if (metric == "N50") {
    return("N50 (bp)")
  } else {
    return("Raw read count")  # Default label without units
  }
}

# Apply add_units_to_facet_labels function to modify facet labels and specify order
plot_data$Facet_Label <- sapply(plot_data$Metric, add_units_to_facet_labels)

# Specify the desired order of facet labels
facet_order <- c("Raw read count", "Preprocessed read count", "Contig number count", 
                 "Largest contig (bp)", "N50 (bp)")

# Convert Facet_Label to factor with specified order
plot_data$Facet_Label <- factor(plot_data$Facet_Label, levels = facet_order)

# Plotting
plot <- ggplot(plot_data, aes(x = factor(Age), y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", color = NA) + # NA removes the border around bars
  facet_grid2(Facet_Label~., scales = "free_y", space = "free_x", switch = "y",labeller = label_wrap_gen(),
              strip = strip_themed(
                size = "variable", # Shrinks 2nd layer
                text_x = elem_list_text(face = c(1, 2)), # 2nd layer in bold
                by_layer_x = TRUE
              ))+
  theme_minimal(base_size = 12) +
  labs(x = "Calibrated years BP", y = NULL) +  # Remove y-axis label
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  # Rotate x-axis labels and adjust size
        axis.text.y = element_text(size = 10),  # Adjust size of y-axis labels
        axis.title = element_text(size = 12),  # Adjust size of axis titles
        strip.text = element_text(size = 12),  # Adjust size of facet labels
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        legend.position = "none") + # Remove legend
  scale_fill_manual(values = metric_colors) +  # Assign gradient of grey colors
  # scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_y_continuous(labels = scales::comma_format(scale = 1e-3)) +  # Display Y-axis labels in K (thousands)
  expand_limits(y = 0)  # Ensure that y-axis starts at 0
print(plot)

ggsave("publication_assembly_plot_final.pdf", plot, width = 9, height = 11, units = "in", dpi = 300)
ggsave("publication_assembly_plot_final.png", plot, width = 9, height = 11, units = "in", dpi = 300)

### BLANKS SAMPLES ####

blank_report <- fread("blank_raw_data_stats.csv", header = TRUE)

# Subset the data for LB and EB samples
lb_data <- subset(blank_report, Sample == "LB")
eb_data <- subset(blank_report, Sample == "EB")

# Define a custom theme with modifications
custom_theme <- theme_minimal() +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()# Remove legend box
  )

p_lb_boxplot <- ggplot(blank_report, aes(x = "", y = `Raw data total read count`, fill = Sample)) +
  # geom_boxplot() +
  geom_boxplot(width = 0.4, position=position_dodge(width=1)) + 
  labs(x = "", y = "Raw read count", fill = "blanks") +
  theme_minimal() + scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) + custom_theme +
  facet_grid(~Sample, scales = "free", space = "free") +
  theme( # Rotate x-axis labels and adjust size
    axis.text.y = element_text(size = 10),  # Adjust size of y-axis labels
    axis.title = element_text(size = 10),  # Adjust size of axis titles
    # strip.text = element_text(size = 12),  # Adjust size of facet labels
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),legend.position = "none") + # Remove legend)+
  # Remove minor gridlines
  # Remove legend
  # scale_fill_manual(values = metric_colors) +  # Assign gradient of grey colors
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_y_continuous(labels = scales::comma_format(scale = 1e-3)) + 
  expand_limits(y = 0)  # En

p_eb_boxplot <- ggplot(blank_report, aes(x = "", y = `Total read counts after filter`, fill = Sample)) +
  geom_boxplot(width = 0.4, position=position_dodge(width=1)) + 
  labs(x = "", y = "Preprocessed read count", fill = "blanks") +
  theme_minimal() + scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) + custom_theme +
  facet_grid(~Sample, scales = "free", space = "free") +
  theme( # Rotate x-axis labels and adjust size
    axis.text.y = element_text(size = 10),  # Adjust size of y-axis labels
    axis.title = element_text(size = 10),  # Adjust size of axis titles
    # Adjust size of facet labels
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),legend.position = "bottom") + # Remove minor gridlines
  # Remove legend
  # scale_fill_manual(values = Sample) +  # Assign gradient of grey colors
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_y_continuous(labels = scales::comma_format(scale = 1e-3)) + 
  expand_limits(y = 0)  # En

# Arrange plots side by side
boxplot_combined <- p_lb_boxplot + p_eb_boxplot + plot_layout(guide = "collect") +
  plot_layout(nrow = 1)  & theme(legend.position = 'bottom')

# Save the combined plot to a file
ggsave("lb_eb_raw_boxplot.pdf", boxplot_combined, width = 8, height = 8, dpi = 300)

#### BOX PLOT Categorical ####


plot_data_epoch <- plot_data %>%
  mutate(epoch = factor(
    ifelse(Age >= 76 & Age <= 12992, "Holocene", 
           ifelse(Age > 13588 & Age <= 24000, "Glacial", NA))))

# Reorder the 'epoch' factor levels
plot_data_epoch$epoch <- factor(plot_data_epoch$epoch, 
                                levels = c("Holocene", "Glacial"))


# Define a custom theme with modifications
custom_theme <- theme_minimal() +
  theme(
    # strip.text = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()# Remove legend box
  )


# Assuming 'subset_data' is your data frame
palette <- c("Holocene" = "red2", 
             "Glacial" = "lightblue")


contig_data <- plot_data_epoch[plot_data_epoch$Metric == "Contig Number", ]
p_contig_data<- ggplot(contig_data, aes(x = factor(epoch), y = Value,fill=epoch,group=epoch)) + 
  geom_boxplot() +
  labs(x = "", y = "count") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + custom_theme + 
  ggtitle("Contig Number")+
  stat_compare_means(
    method = "t.test",
    # label = "p.signif", size=5,  # Use 'p.signif' for significant markers
    comparisons = list(c("Holocene", "Glacial")),  # Specify groups for comparison
    hide.ns = F,  # Hide non-significant comparisons
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", ""))  # Customize significance symbols
  )

N50_data <- plot_data_epoch[plot_data_epoch$Metric == "N50", ]
p_N50<- ggplot(N50_data, aes(x = factor(epoch), y = Value,fill=epoch,group=epoch)) + 
  geom_boxplot() +
  labs(x = "", y = "bp") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + custom_theme + 
  ggtitle("N50 (bp)")+
  stat_compare_means(
    method = "t.test",
    # label = "p.signif", size=5,  # Use 'p.signif' for significant markers
    comparisons = list(c("Holocene", "Glacial")),  # Specify groups for comparison
    hide.ns = F,  # Hide non-significant comparisons
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", ""))  # Customize significance symbols
  )


Largest_Contig_data <- plot_data_epoch[plot_data_epoch$Metric == "Largest Contig", ]
p_Largest_Contig_data<- ggplot(Largest_Contig_data, aes(x = factor(epoch), y = Value,fill=epoch,group=epoch)) + 
  geom_boxplot() +
  labs(x = "", y = "bp") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + custom_theme +
  ggtitle("Largest Contig (bp)")+
  stat_compare_means(
    method = "t.test",
    # label = "p.signif", size=5,  # Use 'p.signif' for significant markers
    comparisons = list(c("Holocene", "Glacial")), # Specify groups for comparison
    hide.ns = F,  # Hide non-significant comparisons
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", ""))  # Customize significance symbols
  )


filtered_data <- plot_data_epoch[plot_data_epoch$Metric == "Filtered and error corrected total read count", ]
p_filtered_data <- ggplot(filtered_data, aes(x = factor(epoch), y = Value,fill=epoch,group=epoch)) + 
  geom_boxplot() +
  labs(x = "", y = "count") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + custom_theme + theme(legend.position = "none") +
  ggtitle("Preprocessed read count")+
  stat_compare_means(
    method = "t.test",
    # label = "p.signif", size=5,  # Use 'p.signif' for significant markers
    comparisons = list(c("Holocene", "Glacial")),  # Specify groups for comparison
    hide.ns = F,  # Hide non-significant comparisons
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", ""))  # Customize significance symbols
  )



raw_read_data <- plot_data_epoch[plot_data_epoch$Metric == "Raw Data read counts", ]
p_raw_read_data <- ggplot(raw_read_data, aes(x = factor(epoch), y = Value,fill=epoch,group=epoch)) + 
  geom_boxplot() +
  labs(x = "", y = "count") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + custom_theme + theme(legend.position = "none") +
  ggtitle("Raw read count")+
  # facet_grid(~Values, scales = "free", space = "free")
  stat_compare_means(
    method = "t.test",
    # label = "p.signif", size=5,  # Use 'p.signif' for significant markers
    comparisons = list(c("Holocene", "Glacial")),  # Specify groups for comparison
    hide.ns = F,  # Hide non-significant comparisons
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", ""))  # Customize significance symbols
  )



composite_plot2 <- p_raw_read_data + p_filtered_data + p_contig_data + p_Largest_Contig_data + p_N50  # Plot for other categories

# Save or display the composite plots
all_merged <- composite_plot2 +plot_layout(guides = 'collect',nrow = 1) & theme(legend.position = 'bottom') 

print(all_merged)

ggsave("assembly_box_plot_publication_plot_v5.png", all_merged, width = 22, height = 8, units = "in", dpi = 300)



### BOX PLOT WHOLE DATASET ###

# here also we add how many reads map to the contigs in whole dataset
# lama_quast_report

# Define the relevant metrics
metrics <- c("Contig Number", "N50", "Largest Contig", "Raw Data read counts", "Filtered and error corrected total read count","Number of read count mapping to contigs")

# Extract sample names
samples <- lama_quast_report$Sample_name

# Reshape the data
box_plot_data <- melt(lama_quast_report, id.vars = "AGE", measure.vars = metrics)

# Rename columns to remove spaces and make them consistent
colnames(box_plot_data) <- c("Age", "Metric", "Value")

# Order the facets
box_plot_data$Metric <- factor(box_plot_data$Metric, levels = c(
  "Raw Data read counts",
  "Total read count after deduplication",
  "Filtered and error corrected total read count",
  "Contig Number",
  "Number of read count mapping to contigs",
  "Largest Contig",
  "N50"
))

# Define units for each metric
units <- c("count", "bp", "bp","count", "bp", "bp")

# Combine metrics and units

# Combine metrics and units
metric_labels <- paste(metrics, "(", units, ")", sep = " ")

num_metrics <- 6
# Generate a gradient of grey colors
metric_colors <- generate_grey_gradient(num_metrics)

p_assembly_data <- ggplot(box_plot_data, aes(x = "", y = Value ,fill=Metric)) + 
  geom_boxplot(fatten=3) +
  facet_wrap(~Metric, ncol=6,scale= "free_y",labeller = label_wrap_gen(width = 15))+
  theme_minimal(base_size = 12) +
  labs(x = "Total Dataset", y = NULL) +  # Remove y-axis label
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  # Rotate x-axis labels and adjust size
        axis.text.y = element_text(size = 14),  # Adjust size of y-axis labels
        axis.title = element_text(size = 12),  # Adjust size of axis titles
        strip.text = element_text(size = 12),  # Adjust size of facet labels
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        legend.position = "none") + # Remove legend
  scale_fill_manual(values = metric_colors) +  # Assign gradient of grey colors
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_y_continuous(labels = scales::comma_format(scale = 1e-3)) +  # Display Y-axis labels in K (thousands)
  expand_limits(y = 0)  # En


ggsave("assembly_box_plot_no_split_publication_plot.png", p_assembly_data , width = 12, height = 4, units = "in", dpi = 300)
ggsave("assembly_box_plot_no_split_publication_plot.pdf", p_assembly_data , width = 12, height = 4, units = "in", dpi = 300)
