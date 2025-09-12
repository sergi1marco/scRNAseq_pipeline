# Load required libraries
library(TCGAplot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(scales)
library(ggrepel)
library(viridis)

# Generate the pan_boxplot for FAP and extract data
# If you encounter an \"unexpected string constant\" error here,
# please ensure the double quotes around \"FAP\" are standard straight quotes (")
# and not \"smart quotes\" or other non-standard characters.
# Also, ensure the file is saved with UTF-8 encoding.
p <- pan_boxplot("FAP", palette="npg", legend="right", method="wilcox.test")
plot_data <- p$data

# Calculate tumor/normal mean expression ratio for each cancer type
ratios <- plot_data %>%
  group_by(Cancer, Group) %>%
  summarise(mean_expr = mean(FAP, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = mean_expr) %>%
  mutate(ratio = Tumor / Normal)

# Order cancer types by decreasing tumor/normal ratio
ordered_cancers <- ratios %>%
  arrange(desc(ratio)) %>%
  pull(Cancer)

# Filter for cancer types with both Tumor and Normal samples
group_counts <- plot_data %>%
  group_by(Cancer, Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = n, values_fill = 0)

valid_cancers <- group_counts %>%
  filter(Tumor > 0 & Normal > 0) %>%
  pull(Cancer)

plot_data <- plot_data %>% filter(Cancer %in% valid_cancers)
ordered_cancers <- ordered_cancers[ordered_cancers %in% valid_cancers]

# Reorder Cancer factor in plot_data
plot_data$Cancer <- factor(plot_data$Cancer, levels = ordered_cancers)

# Calculate Wilcoxon p-values for each cancer type
# Calculate Wilcoxon p-values for each cancer type
pval_annot <- plot_data %>%
  group_by(Cancer) %>%
  summarise(
    p.value = tryCatch(
      wilcox.test(FAP[Group == "Tumor"], FAP[Group == "Normal"], exact = FALSE)$p.value, # <<< Added exact = FALSE
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

# Merge ratio info for annotation
annot_df <- left_join(
  pval_annot,
  ratios %>% select(Cancer, ratio),
  by = "Cancer"
)

# Format ratio and p-value labels for display
annot_df$label <- ifelse(
  is.na(annot_df$ratio) | is.na(annot_df$p.value),
  "n/a",
  paste0(
    "ratio = ", sprintf("%.2f", annot_df$ratio), "\n",
    "p = ", sprintf("%.2g", annot_df$p.value)
  )
)

# Calculate y-position for annotation labels
y_positions <- plot_data %>%
  group_by(Cancer) %>%
  summarise(y = max(FAP, na.rm = TRUE) * 1.18)
annot_df <- left_join(annot_df, y_positions, by = "Cancer")

# Box plot code (remains unchanged)
ggplot(plot_data, aes(x = Cancer, y = FAP, fill = Group)) +
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

# SCATTER PLOT CODE STARTS HERE
filtered_ratios <- ratios %>% filter(Cancer %in% valid_cancers)

mean_tumor_expr <- plot_data %>%
  filter(Group == "Tumor") %>%
  group_by(Cancer) %>%
  summarise(
    mean_tumor_FAP = mean(FAP, na.rm = TRUE),
    sd_tumor_FAP = sd(FAP, na.rm = TRUE),
    n_tumor_FAP = n(),
    .groups = "drop"
  ) %>%
  mutate(
    se_tumor_FAP = sd_tumor_FAP / sqrt(n_tumor_FAP)
  )

scatter_data <- left_join(filtered_ratios, mean_tumor_expr, by = "Cancer") %>%
  drop_na(ratio, mean_tumor_FAP) %>%
  filter(mean_tumor_FAP > 0)

# Automatically set axis ranges based on data
# Calculate min/max for error bars to ensure full visibility
min_x_for_err <- min(scatter_data$ratio - scatter_data$se_tumor_FAP, na.rm = TRUE)
max_x_for_err <- max(scatter_data$ratio + scatter_data$se_tumor_FAP, na.rm = TRUE)
min_y_for_err <- min(scatter_data$mean_tumor_FAP - scatter_data$se_tumor_FAP, na.rm = TRUE)
max_y_for_err <- max(scatter_data$mean_tumor_FAP + scatter_data$se_tumor_FAP, na.rm = TRUE)

# Apply a multiplicative buffer for log scales
# Ensure min values are strictly positive before buffering
buffer_x_min <- ifelse(min_x_for_err > 0, min_x_for_err * 0.8, min(scatter_data$ratio[scatter_data$ratio > 0], na.rm = TRUE) * 0.8)
buffer_x_max <- max_x_for_err * 1.2
buffer_y_min <- ifelse(min_y_for_err > 0, min_y_for_err * 0.5, min(scatter_data$mean_tumor_FAP[scatter_data$mean_tumor_FAP > 0], na.rm = TRUE) * 0.5)
buffer_y_max <- max_y_for_err * 2.0


# Ask user for X-intercept for vertical line
message("\nTo define the central vertical line for quadrants,")
message("please enter an X-value (e.g., 1 for no ratio change, leave blank for no line):")
x_intercept_input <- readline(prompt="Enter X-intercept value: ")
x_intercept <- as.numeric(x_intercept_input)

# Ask user for Y-intercept for horizontal line
message("\nTo define the central horizontal line for quadrants,")
message("please enter a Y-value (must be > 0 for log scale, e.g., 1 for baseline expression, leave blank for no line):")
y_intercept_input <- readline(prompt="Enter Y-intercept value: ")
y_intercept <- as.numeric(y_intercept_input)

# Generate the scatter plot and assign it to an object

scatter_plot_obj <- ggplot(scatter_data, aes(x = ratio, y = mean_tumor_FAP)) +
  # Make error bars black and draw them *behind* the dots
  geom_errorbar(
    aes(
      ymin = mean_tumor_FAP - se_tumor_FAP,
      ymax = mean_tumor_FAP + se_tumor_FAP,
      group = Cancer
    ),
    color = "black",
    width = 0.05
  ) +
  # Add black stroke to dots and change color mapping to fill
  geom_point(aes(fill = log10(ratio)), color = "black", shape = 21, size = 5.5, stroke = 1.2) +
  # Use geom_text_repel for intelligent label placement
  ggrepel::geom_text_repel(aes(label = Cancer),
                           size = 4.5,
                           box.padding = 0.8,
                           point.padding = 0.8,
                           segment.color = 'grey50',
                           max.overlaps = Inf
  ) +
  # Use scale_fill_distiller as fill is now mapped to log10(ratio)
  scale_fill_distiller(palette = "BrBG",
                       type = "div",
                       name = "Tumor/Normal Ratio",
                       breaks = c(log10(0.5), log10(1), log10(5), log10(10)), # <<< Added specific breaks (log10 of desired values)
                       labels = function(x) { 10^x }
  ) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA), # Added panel border
    axis.line = element_line(colour = "black"),
    legend.position = "right"
  ) +
  scale_x_log10(
    limits = c(buffer_x_min, buffer_x_max),
    breaks = c(0.5, 1, 5, 10), # <<< MODIFIED: Set specific breaks for X-axis
    labels = scales::number_format()
  ) +
  scale_y_log10(
    limits = c(buffer_y_min, buffer_y_max),
    breaks = c(0.5, 1, 5), # <<< MODIFIED: Set specific breaks for Y-axis
    labels = scales::number_format()
  ) +
  annotation_logticks(sides = "bl") +
  coord_fixed(ratio = 1)

# Add vertical line for quadrant division if x_intercept is provided and valid
if (!is.na(x_intercept) && x_intercept > 0) {
  scatter_plot_obj <- scatter_plot_obj +
    geom_vline(xintercept = x_intercept, linetype = "11", color = "black", linewidth = 0.8)
}

# Add horizontal line for quadrant division if y_intercept is provided and valid
if (!is.na(y_intercept) && y_intercept > 0) {
  scatter_plot_obj <- scatter_plot_obj +
    geom_hline(yintercept = y_intercept, linetype = "11", color = "black", linewidth = 0.8)
}

# Print the plot to display it
print(scatter_plot_obj)

# Ask user for directory path to save the PDF
message("Please enter the full directory path to save the PDF file (e.g., C:/Users/YourName/Documents/Plots/ or /Users/YourName/Documents/Plots/):\")")
save_directory <- readline(prompt="Directory Path: ")

# Define the filename
output_filename <- "FAP_Tumor_Scatter_Plot.pdf"

# Combine the directory and filename
pdf_path <- file.path(save_directory, output_filename)

# Save the plot
if (dir.exists(save_directory)) {
  ggsave(
    filename = pdf_path,
    plot = scatter_plot_obj,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300
  )
  message(paste("Graph saved successfully to:", pdf_path))
} else {
  message("Error: The specified directory does not exist. Please check the path and try again.")
}
