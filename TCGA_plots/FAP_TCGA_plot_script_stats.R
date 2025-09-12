# ===============================
# FAP Expression Analysis Script
# ===============================

# --- Load libraries ---
suppressPackageStartupMessages({
  library(TCGAplot)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggsci)
  library(rstudioapi)
})

# --- Choose working directory interactively ---
wd <- rstudioapi::selectDirectory(caption = "Choose output folder")
setwd(wd)
message("Working directory set to: ", wd)

# --- Ask user what to run ---
choice <- menu(c("Full analysis (plot + stats)", "Only collect stats"),
               title = "Choose analysis mode:")

# ==========================================================
# 1. Generate plot_data and ratios (always needed for stats)
# ==========================================================
p <- pan_boxplot("FAP", palette = "npg", legend = "right", method = "wilcox.test")
plot_data <- p$data

ratios <- plot_data %>%
  group_by(Cancer, Group) %>%
  summarise(mean_expr = mean(FAP, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = mean_expr) %>%
  mutate(ratio = Tumor / Normal)

# ==========================================================
# 2. If full analysis: filtering, ordering, and plotting
# ==========================================================
if (choice == 1) {
  # Filter for valid cancers (both Tumor & Normal available)
  group_counts <- plot_data %>%
    group_by(Cancer, Group) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = n, values_fill = 0)
  
  valid_cancers <- group_counts %>%
    filter(Tumor > 0 & Normal > 0) %>%
    pull(Cancer)
  
  plot_data <- plot_data %>% filter(Cancer %in% valid_cancers)
  ordered_cancers <- ratios %>%
    arrange(desc(ratio)) %>%
    pull(Cancer)
  ordered_cancers <- ordered_cancers[ordered_cancers %in% valid_cancers]
  plot_data$Cancer <- factor(plot_data$Cancer, levels = ordered_cancers)
  
  # Compute annotation labels
  pval_annot <- plot_data %>%
    group_by(Cancer) %>%
    summarise(
      p.value = tryCatch(
        wilcox.test(FAP[Group == "Tumor"], FAP[Group == "Normal"])$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    )
  
  annot_df <- left_join(
    pval_annot,
    ratios %>% select(Cancer, ratio),
    by = "Cancer"
  )
  
  annot_df$label <- ifelse(
    is.na(annot_df$ratio) | is.na(annot_df$p.value),
    "n/a",
    paste0(
      "ratio = ", sprintf("%.2f", annot_df$ratio), "\n",
      "p = ", sprintf("%.2g", annot_df$p.value)
    )
  )
  
  y_positions <- plot_data %>%
    group_by(Cancer) %>%
    summarise(y = max(FAP, na.rm = TRUE) * 1.18)
  annot_df <- left_join(annot_df, y_positions, by = "Cancer")
  
  # Create the plot
  plot_obj <- ggplot(plot_data, aes(x = Cancer, y = FAP, fill = Group)) +
    geom_boxplot() +
    scale_fill_npg() +
    geom_text(
      data = annot_df,
      aes(x = Cancer, y = y, label = label),
      inherit.aes = FALSE,
      size = 3.2,
      vjust = 0
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    labs(
      title = "FAP Expression in Tumor vs Normal (Ordered by Tumor/Normal Ratio)",
      x = "Cancer Type",
      y = "Expression"
    )
  
  # Save plot
  ggsave("FAP_plot.pdf", plot_obj, width = 10, height = 6)
  message("Figure saved: FAP_plot.pdf")
}

# ==========================================================
# 3. Stats collection (always run)
# ==========================================================
sample_counts <- plot_data %>%
  group_by(Cancer, Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = n, values_fill = 0)

expr_summary <- plot_data %>%
  group_by(Cancer, Group) %>%
  summarise(
    n = n(),
    mean_expr = mean(FAP, na.rm = TRUE),
    median_expr = median(FAP, na.rm = TRUE),
    sd_expr = sd(FAP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(n, mean_expr, median_expr, sd_expr),
    values_fill = NA
  )

test_results <- plot_data %>%
  group_by(Cancer) %>%
  summarise(
    test = "Wilcoxon rank-sum test (unpaired, two-sided)",
    p.value = tryCatch(
      wilcox.test(FAP[Group == "Tumor"], FAP[Group == "Normal"])$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

stats_table <- ratios %>%
  left_join(expr_summary, by = "Cancer") %>%
  left_join(test_results, by = "Cancer")

# Save stats table
write.csv(stats_table, "FAP_TCGA_GTEx_stats.csv", row.names = FALSE)
message("Statistics table saved: FAP_TCGA_GTEx_stats.csv")

# --- End of Script ---
