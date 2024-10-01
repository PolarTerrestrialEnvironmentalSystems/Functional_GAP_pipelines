library(GGally)
library(tidyr)
load("eukaryotic_KOs_richness.RData")
load("assembly_stats.Rdata")

largest_contig <- plot_data[which(plot_data$Metric == "N50"),]
colnames(largest_contig)[1] <- "Ages"
contig_number <- plot_data[which(plot_data$Metric == "Contig Number"),]
colnames(contig_number)[1] <- "Ages"
raw_data <- plot_data[which(plot_data$Metric == "Raw Data read counts"),]
colnames(raw_data)[1] <- "Ages"

# combined_result
combined_result$Ages <- as.integer(as.numeric(gsub("\\.$", "", combined_result$Ages)))
all_pipeline_selected <- combined_result %>% select("qD", "Ages","approach")

all_pipeline_wide <- pivot_wider(
  all_pipeline_selected,
  names_from = approach,
  values_from = qD
)

# ProkGAP case
eukaryotic_KOs <- all_pipeline_wide %>% select("ProkGAP", "Ages")
eukaryotic_KOs <- eukaryotic_KOs[order(eukaryotic_KOs$Ages), ]

merged_data <- inner_join(eukaryotic_KOs, contig_number, by = "Ages")
colnames(merged_data)[1] <- "Eukaryotic KO Richness"
merged_data$Metric <- NULL
merged_data$Facet_Label <- NULL
colnames(merged_data)[3] <- "Contig number (N)"
merged_data = inner_join(merged_data, largest_contig, by = "Ages")
merged_data$Facet_Label <- NULL
colnames(merged_data)[5] <- "Median contig length (bp)"
merged_data$Metric <- NULL
merged_data = inner_join(merged_data, raw_data, by = "Ages")
merged_data$Facet_Label <- NULL
colnames(merged_data)[6] <- "Short read count (N)"
merged_data$Metric <- NULL
merged_data$variable <- NULL
# ggpairs(merged_data)

merged_data <- merged_data[, c(2, 1,3,4,5)] 
merged_data$`Eukaryotic KO Richness` <- as.integer(merged_data$`Eukaryotic KO Richness`)

printVar <- function(x, y) {
  vals <- cor.test(x, y, method = "spearman")[c("estimate", "p.value")]
  names(vals) <- c("rho", "p")
  paste(names(vals), signif(unlist(vals), 2), collapse = "\n")
}

my_fn <- function(data, mapping, ...) {
  xData <- eval_data_col(data, mapping$x)
  yData <- eval_data_col(data, mapping$y)
  colorData <- eval_data_col(data, mapping$colour)
  if (is.null(colorData)) {
    colorData <- rep(1, length(xData))
  }
  byGroup <- by(data.frame(xData, yData), colorData, function(i) printVar(i[, 1], i[, 2]))
  byGroup <- data.frame(col = names(byGroup), label = as.character(byGroup))
  byGroup$x <- 0.5
  byGroup$y <- seq(0.8 - 0.3, 0.2, length.out = nrow(byGroup))  
  # correlation
  mainCor <- printVar(xData, yData)
  p <- ggplot(data = data, mapping = mapping) +
    annotate(x = 0.5, y = 0.8, label = mainCor, geom = "text", size = 3) +
    theme(legend.position = "bottom") +
    theme_void() + ylim(c(0, 1))
  p
}
p <- ggpairs(merged_data, columns = 1:ncol(merged_data),
             columnLabels = gsub('.', ' ', colnames(merged_data), fixed = TRUE), 
             labeller = label_wrap_gen(10), upper = list(continuous = wrap("cor", method = "spearman"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        # axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
# check the corr
print(p)

## for all pipelines

merged_data <- inner_join(all_pipeline_wide, contig_number, by = "Ages")
colnames(merged_data)[7] <- "Contig number (N)"
merged_data$Metric <- NULL
merged_data = inner_join(merged_data, largest_contig, by = "Ages")
colnames(merged_data)[8] <- "Median contig length (bp)"
merged_data$Metric <- NULL
merged_data = inner_join(merged_data, raw_data, by = "Ages")
colnames(merged_data)[9] <- "Short read count (N)"
merged_data$Metric <- NULL

# Convert columns to integers
merged_data <- merged_data %>%
  mutate(
    Ages = as.integer(Ages),
    ProkGAP = as.integer(ProkGAP),
    EukGAP = as.integer(EukGAP),
    `Pre-class ProkGAP` = as.integer(`Pre-class ProkGAP`),
    `Pre-class EukGAP` = as.integer(`Pre-class EukGAP`),
  )
# check
str(merged_data)

p <- ggpairs(merged_data, columns = 1:ncol(merged_data),
             columnLabels = gsub('.', ' ', colnames(merged_data), fixed = TRUE), 
             labeller = label_wrap_gen(10), upper = list(continuous = wrap("cor", method = "spearman"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        # axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

ggsave("KOs_richness_all_pipeline_correlation_plot_v1.png", p, width = 9, height = 9, units = "in",dpi = 300)
