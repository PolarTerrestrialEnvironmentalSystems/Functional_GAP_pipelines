setwd("")
load("abundance_analysis_results.RData")
results_list

Ages <- c(76,264,724,1212,1858,2504,3023,4381,5168,6400,6890,7562,8513,9190,9816,10091,10428,11536,11791,12992,13588,13746,13930,14020,14315,14669,15123,15599,16084,16541,16946,17222,17948,18715,19192,19593,20059,20630,21068,21753,22518,22982)

results_list[[1]]$DataInfo
prokgap_result <- results_list[[1]]$iNextEst$size_based
prokgap_result = subset(prokgap_result, m=="46564")
prokgap_result$approach <- "ProkGAP"
prokgap_result$Ages <- Ages

eukgap_result <- results_list[[2]]$iNextEst$size_based
eukgap_result = subset(eukgap_result, m=="63945")
eukgap_result$approach <- "EukGAP"
eukgap_result$Ages <- Ages

preclass_prokgap_result <- results_list[[3]]$iNextEst$size_based
preclass_prokgap_result = subset(preclass_prokgap_result, m=="41075")
preclass_prokgap_result$approach <- "Pre-class ProkGAP"
preclass_prokgap_result$Ages <- Ages

preclass_eukgap_result <- results_list[[4]]$iNextEst$size_based
preclass_eukgap_result = subset(preclass_eukgap_result, m=="5628")
preclass_eukgap_result$approach <- "Pre-class EukGAP"
preclass_eukgap_result$Ages <- Ages

combined_result <- rbind(prokgap_result,eukgap_result,preclass_prokgap_result,preclass_eukgap_result)

palette <- c("ProkGAP" = "dodgerblue3", 
             "EukGAP" = "gold", 
             "Pre-class ProkGAP" = "darkorchid3", 
             "Pre-class EukGAP" = "forestgreen")

p <- ggplot(subset(combined_result), aes(x = factor(Ages), y = qD, fill = approach)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  labs(x = "Age / cal yr BP", y = "Value", fill = "Variable") +
  theme_minimal()+ scale_fill_manual(values = palette)
print(p)

prokgap_assemblage <- results_list[[1]]$DataInfo$Assemblage
prokgap_s_obs <- results_list[[1]]$DataInfo$S.obs
prokgap_df <- data.frame(Ages = Ages, Assemblage = prokgap_assemblage, S.obs = prokgap_s_obs)
prokgap_df$approach <- "ProkGAP"

eukgap_assemblage <- results_list[[2]]$DataInfo$Assemblage
eukgap_s_obs <- results_list[[2]]$DataInfo$S.obs
eukgap_df <- data.frame(Ages = Ages, Assemblage = eukgap_assemblage, S.obs = eukgap_s_obs)
eukgap_df$approach <- "EukGAP"

preclass_prokgap_assemblage <- results_list[[3]]$DataInfo$Assemblage
preclass_prokgap_s_obs <- results_list[[3]]$DataInfo$S.obs
preclass_prokgap_df <- data.frame(Ages = Ages, Assemblage = preclass_prokgap_assemblage, S.obs = preclass_prokgap_s_obs)
preclass_prokgap_df$approach <- "Pre-class ProkGAP"

preclass_eukgap_assemblage <- results_list[[4]]$DataInfo$Assemblage
preclass_eukgap_s_obs <- results_list[[4]]$DataInfo$S.obs
preclass_eukgap_df <- data.frame(Ages = Ages, Assemblage = preclass_eukgap_assemblage, S.obs = preclass_eukgap_s_obs)
preclass_eukgap_df$approach <- "Pre-class EukGAP"

combined_result_no_resample <- rbind(prokgap_df,eukgap_df,preclass_prokgap_df,preclass_eukgap_df)

palette <- c("ProkGAP" = "dodgerblue3", 
             "EukGAP" = "gold", 
             "Pre-class ProkGAP" = "darkorchid3", 
             "Pre-class EukGAP" = "forestgreen")

p <- ggplot(subset(combined_result_no_resample), aes(x = factor(Ages), y = S.obs,fill = approach)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  labs(x = "Age / cal yr BP", y = "Value", fill = "Variable") +
  theme_minimal()+ scale_fill_manual(values = palette)
print(p)

combined_result_no_resample$method <- "Original"
combined_result$method <- "Resampled"

selected_columns <- c("Assemblage", "Ages", "qD", "approach", "method")

new_df <- combined_result[, selected_columns, drop = FALSE]
colnames(new_df)[3] <- "S.obs"

merged_all_result <- rbind(combined_result_no_resample,new_df)

new_order <- c("ProkGAP", "EukGAP", "Pre-class ProkGAP","Pre-class EukGAP")

merged_all_result$approach <- factor(merged_all_result$approach, levels = new_order)

merged_all_result$approach <- as.character(merged_all_result$approach)

merged_all_result$approach[which(merged_all_result$approach == "prokGAP")] <- "ProkGAP"

merged_all_result$approach <- factor(merged_all_result$approach)

plot <- ggplot(merged_all_result, aes(x = factor(Ages), y = S.obs, fill = approach)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = NA, size = 3.8) +
  labs(x = "Calibrated years BP", y = "(KO) richness", title = "") +
  facet_wrap(~ method, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.position = "top",
    # legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm")
  ) +
  scale_fill_manual(values = palette) +
  guides(fill = guide_legend(title = "prediction pipelines",position = "bottom")) 
plot <- plot +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14), # Rotate x-axis labels for better readability
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.spacing = unit(1, "lines") # Increase space between panels
  )


ggsave("euk_prok_observed_merged_plot.png", plot, width = 8, height = 6, dpi = 300)
ggsave("euk_prok_observed_merged_plot.pdf", plot, width = 8, height = 6, dpi = 300)
tiff("test.tiff", units="in", width=5, height=5, res=300)
dev.off()
print(plot)


### box plot ###

# Define a custom theme with modifications for publication
custom_theme <- theme_minimal() +
  theme(
    strip.text = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

merged_all_result_epoch$approach <- as.character(merged_all_result_epoch$approach)

merged_all_result_epoch$approach[which(merged_all_result_epoch$approach == "prokGAP")] <- "ProkGAP"

merged_all_result_epoch$approach <- factor(merged_all_result_epoch$approach)

new_order <- c("ProkGAP", "EukGAP", "Pre-class ProkGAP","Pre-class EukGAP")

merged_all_result_epoch$approach <- factor(merged_all_result_epoch$approach, levels = new_order)

### BOX PLOT ### 

p_all_kos_data <- ggplot(subset(merged_all_result_epoch), aes(x = "", y = S.obs,fill=approach,group=approach)) + 
  geom_boxplot() +
  labs(x = "", y = "KO Richness",fill = "prediction pipelines") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))  + custom_theme + 
  ggtitle("") +
  facet_grid(~method, scales = "free", space = "free") + theme(legend.position = "bottom")
p_all_kos_data
ggsave("euk_prok_observed_boxplot_nosplit.png", p_all_kos_data, width = 18, height = 4, dpi = 300)
