# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(grid) # for unit() in margin()
})

# Parse user selection input allowing commas and hyphen ranges, robustly
parse_selection_input <- function(input_str, max_idx) {
  input_str <- gsub(" ", "", input_str)  # Remove spaces
  parts <- unlist(strsplit(input_str, ","))
  selected <- integer(0)
  for (p in parts) {
    if (grepl("-", p)) {
      bounds <- as.integer(unlist(strsplit(p, "-")))
      if (length(bounds) == 2 && all(!is.na(bounds))) {
        if (bounds[1] <= bounds[2]) {
          selected <- c(selected, seq(bounds[1], bounds[2]))
        } else {
          message(sprintf("Ignoring invalid range: %s", p))
        }
      } else {
        message(sprintf("Skipping invalid range: %s", p))
      }
    } else {
      val <- as.integer(p)
      if (!is.na(val)) {
        selected <- c(selected, val)
      } else {
        message(sprintf("Skipping invalid value: %s", p))
      }
    }
  }
  selected <- unique(selected)
  selected <- selected[selected >= 1 & selected <= max_idx]
  return(sort(selected))
}

# Main function to plot percent positive bar plot by cell type(s)
plot_gene_positive_bar <- function(seurat_obj, gene_name) {
  if (!gene_name %in% rownames(seurat_obj)) {
    stop(sprintf("Gene '%s' not found in the Seurat object.", gene_name))
  }
  
  celltype_col <- "scATOMIC_prediction"
  DefaultAssay(seurat_obj) <- "RNA"
  expr_vec <- Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "data")[gene_name, ]
  is_pos <- as.numeric(expr_vec > 0)
  meta <- seurat_obj@meta.data
  df <- data.frame(celltype = as.character(meta[[celltype_col]]), positive = is_pos)
  
  pct_df <- df %>%
    group_by(celltype) %>%
    summarize(
      n_total = n(),
      n_pos = sum(positive, na.rm = TRUE),
      pct_pos = 100 * n_pos / n_total,
      .groups = "drop"
    ) %>%
    filter(n_pos > 0) %>%
    arrange(desc(pct_pos))
  
  if (nrow(pct_df) == 0) {
    message(sprintf("No cells express '%s'. Nothing to plot.", gene_name))
    return(invisible(NULL))
  }
  
  cat("Available populations:\n")
  for (i in seq_len(nrow(pct_df))) {
    cat(sprintf("%d: %s (%.1f%% positive)\n", i, pct_df$celltype[i], pct_df$pct_pos[i]))
  }
  
  user_input <- readline(prompt = "Enter numbers of populations to include (e.g. 1,3-5): ")
  selected_indices <- parse_selection_input(user_input, nrow(pct_df))
  
  if (length(selected_indices) == 0) {
    message("No valid populations selected. Plotting all.")
    selected_indices <- seq_len(nrow(pct_df))
  }
  
  selected_df <- pct_df[selected_indices, , drop = FALSE]
  
  colors <- rep("#D3D3D3", nrow(selected_df))
  colors[1] <- "#F1948A"
  selected_df$color <- colors
  
  # --- Existing code for percentage plot ---
  p_bar <- ggplot(selected_df, aes(x = reorder(celltype, -pct_pos), y = pct_pos, fill = color)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = sprintf("%.1f%%", pct_pos)), vjust = -0.5, size = 3) +
    scale_fill_identity() +
    labs(title = NULL, x = NULL, y = "FAP+ cells (percentage of total)") +
    ylim(0, max(100, ceiling(max(selected_df$pct_pos) / 10) * 10)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.3),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
      axis.title.x = element_blank()
    ) +
    coord_cartesian(clip = "off")
  
  pdf_width <- max(7, 0.7 * nrow(selected_df))
  pdf_file <- sprintf("PercentPositive_%s_bySelectedCellTypes.pdf", gsub("[^A-Za-z0-9._-]", "_", gene_name))
  grDevices::pdf(pdf_file, width = pdf_width, height = 4.5)
  print(p_bar)
  grDevices::dev.off()
  
  message(sprintf("Saved filtered bar plot to '%s'.", pdf_file))
  
  # --- NEW CODE for raw counts plot ---
  p_counts <- ggplot(selected_df, aes(x = reorder(celltype, -n_pos), y = n_pos, fill = color)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = n_pos), vjust = -0.5, size = 3) + # Use n_pos for the label
    scale_fill_identity() +
    labs(title = NULL, x = NULL, y = "Number of FAP+ cells") +
    ylim(0, max(selected_df$n_pos) * 1.1) + # Set y-axis limit based on max count
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.3),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
      axis.title.x = element_blank()
    ) +
    coord_cartesian(clip = "off")
  
  pdf_counts_file <- sprintf("RawCounts_%s_bySelectedCellTypes.pdf", gsub("[^A-Za-z0-9._-]", "_", gene_name))
  grDevices::pdf(pdf_counts_file, width = pdf_width, height = 4.5)
  print(p_counts)
  grDevices::dev.off()
  
  message(sprintf("Saved raw counts bar plot to '%s'.", pdf_counts_file))
}

# Example interactive call
gene_name <- readline(prompt = "Enter gene symbol to plot: ")
plot_gene_positive_bar(seurat_object, gene_name)