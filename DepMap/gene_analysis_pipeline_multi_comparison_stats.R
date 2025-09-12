# === Gene Expression Analysis in Cancer vs Normal Cell Lines (Multiple Genes) ===

# Load required packages
required_packages <- c("tidyverse", "ggplot2", "ggpubr", "biomaRt", "ExperimentHub", "AnnotationHub", "ComplexHeatmap", "RColorBrewer", "viridis", "grid", "rstatix") # ADDED 'rstatix'
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(biomaRt)
  library(ExperimentHub)
  library(AnnotationHub)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(viridis)
  library(grid)
  library(rstatix) # ADDED 'rstatix'
})

# Set working directory
working_dir <- "C:/Users/Sergi Marco/OneDrive - University of Glasgow/SMARCO/AVACTA project/Databases/DepMap/correlations"
if (!dir.exists(working_dir)) dir.create(working_dir, recursive = TRUE)
setwd(working_dir) # Set working directory initially

# Create output directory
output_dir <- file.path(working_dir, "Gene_Analysis_Outputs_MultiGene") # Changed output directory
if (!dir.exists(output_dir)) dir.create(output_dir)

# Prompt for genes (multiple genes)
cat("\n=== Gene Expression Analysis ===\n")
gene_symbols_input <- readline(prompt = "Enter gene symbols to analyze, separated by commas (e.g., FAP, CD44, TP53): ")
gene_symbols <- str_split(gene_symbols_input, pattern = ",\\s*")[[1]] %>% toupper() %>% unique()

if (length(gene_symbols) == 0) {
  stop("No gene symbols entered. Please provide at least one gene symbol.")
}

# Load DepMap data
cat("\nLoading DepMap data...\n")
eh <- ExperimentHub(cache = working_dir)
depmap_data <- query(eh, "depmap")
depmap_df <- as.data.frame(mcols(depmap_data))

# Load latest TPM expression dataset
tpm_idx <- grep("TPM_", depmap_df$title)
if (length(tpm_idx) == 0) stop("No TPM expression data found in DepMap")
latest_tpm_idx <- which.max(as.Date(depmap_df$rdatadateadded[tpm_idx]))
latest_tpm <- eh[[names(depmap_data)[tpm_idx[latest_tpm_idx]]]]
cat("Loaded expression data:", depmap_df$title[tpm_idx[latest_tpm_idx]], "\n")
tpm_df <- as.data.frame(latest_tpm)

# Load latest metadata
meta_idx <- grep("metadata_", depmap_df$title)
if (length(meta_idx) == 0) stop("No metadata found in DepMap")
latest_meta_idx <- which.max(as.Date(depmap_df$rdatadateadded[meta_idx]))
metadata <- eh[[names(depmap_data)[meta_idx[latest_meta_idx]]]]
metadata_df <- as.data.frame(metadata)

# Initialize a list to store expression data for all genes
all_genes_expr_data <- list()
validated_gene_symbols <- c()

# Loop through each gene to get info and expression
cat("\nGetting gene information and expression for provided genes...\n")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

for (g_symbol in gene_symbols) {
  cat(paste0("Processing gene: ", g_symbol, "\n"))
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "external_synonym"),
    filters = "hgnc_symbol",
    values = g_symbol,
    mart = ensembl
  )
  
  if (nrow(gene_info) == 0) {
    # Retry using synonym
    gene_info <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol", "external_synonym"),
      filters = "external_synonym",
      values = g_symbol,
      mart = ensembl
    )
    if (nrow(gene_info) == 0) {
      cat(paste("Warning: Gene", g_symbol, "not found in Ensembl by symbol or synonym. Skipping.\n"))
      next
    }
    cat("Mapped input symbol to known Ensembl synonym:", gene_info$hgnc_symbol[1], "\n")
    current_gene_symbol <- gene_info$hgnc_symbol[1]
  } else {
    current_gene_symbol <- gene_info$hgnc_symbol[1]
  }
  
  ensembl_id <- gene_info$ensembl_gene_id[1]
  
  # Extract expression
  gene_expr_data <- NULL
  if ("gene_name" %in% colnames(tpm_df)) {
    gene_expr_data <- tpm_df %>%
      filter(gene_name == current_gene_symbol) %>%
      dplyr::select(depmap_id, expression = rna_expression)
  } else if (any(grepl(ensembl_id, colnames(tpm_df)))) {
    gene_col <- grep(ensembl_id, colnames(tpm_df), value = TRUE)[1]
    gene_expr_data <- tpm_df %>%
      dplyr::select(depmap_id, expression = all_of(gene_col))
  }
  
  if (is.null(gene_expr_data) || nrow(gene_expr_data) == 0 || all(is.na(gene_expr_data$expression))) {
    cat(paste("Warning: No expression data found for gene", current_gene_symbol, ". Skipping.\n"))
    next
  }
  
  # Rename the 'expression' column to the gene symbol for merging
  colnames(gene_expr_data)[colnames(gene_expr_data) == "expression"] <- current_gene_symbol
  all_genes_expr_data[[current_gene_symbol]] <- gene_expr_data
  validated_gene_symbols <- c(validated_gene_symbols, current_gene_symbol)
}

if (length(all_genes_expr_data) == 0) {
  stop("No valid gene expression data found for any of the provided genes. Exiting.")
}

# Merge all gene expression data into a single dataframe
# Start with the first gene's data, then join others
merged_expr_df <- all_genes_expr_data[[1]]
if (length(all_genes_expr_data) > 1) {
  for (i in 2:length(all_genes_expr_data)) {
    merged_expr_df <- left_join(merged_expr_df, all_genes_expr_data[[i]], by = "depmap_id")
  }
}

# Merge with metadata and classify
analysis_data_multi_gene <- merged_expr_df %>%
  left_join(metadata_df, by = "depmap_id") %>%
  dplyr::select(depmap_id, cell_line_name, primary_disease, lineage, all_of(validated_gene_symbols)) %>% # Select all validated gene columns
  mutate(
    cell_type = ifelse(primary_disease == "Non-Cancerous", "Normal", "Cancer"),
    subgroup = lineage
  ) %>%
  mutate(
    primary_disease = ifelse(primary_disease == "Non-Cancerous", "zNon-Cancerous", primary_disease)
  ) %>%
  drop_na(all_of(validated_gene_symbols)) # Drop rows with NA in any gene expression column

# Individual gene analysis (optional, but good for understanding)
# This part is largely the same, but now it loops for each gene.
# You might want to remove it or modify it to save individual plots only if needed,
# as the main goal is the combined heatmap.
cat("\nPerforming individual gene analysis (histograms, summaries)...\n")
for (g_symbol in validated_gene_symbols) {
  current_gene_data <- analysis_data_multi_gene %>%
    dplyr::select(depmap_id, cell_line_name, primary_disease, lineage, expression = all_of(g_symbol), cell_type, subgroup)
  
  cancer_data_single <- current_gene_data %>% filter(cell_type == "Cancer")
  normal_data_single <- current_gene_data %>% filter(cell_type == "Normal")
  
  if (nrow(cancer_data_single) > 0) {
    cat(paste0("\n=== Analysis for ", g_symbol, " in Cancer cell lines ===\n"))
    cat("Number of cell lines:", nrow(cancer_data_single), "\n")
    print(summary(cancer_data_single$expression))
    p_cancer <- ggplot(cancer_data_single, aes(x = expression)) +
      geom_histogram(fill = "firebrick", bins = nclass.Sturges(cancer_data_single$expression)) +
      labs(title = paste(g_symbol, "Expression in Cancer Cell Lines"),
           x = "Expression (TPM)", y = "Number of Cell Lines") +
      theme_minimal()
    ggsave(file.path(output_dir, paste0(g_symbol, "_Cancer_expression.pdf")), p_cancer, width = 8, height = 6)
  } else {
    cat(paste0("\nNo Cancer cell lines found for gene ", g_symbol, "\n"))
  }
  
  if (nrow(normal_data_single) > 0) {
    cat(paste0("\n=== Analysis for ", g_symbol, " in Normal cell lines ===\n"))
    cat("Number of cell lines:", nrow(normal_data_single), "\n")
    print(summary(normal_data_single$expression))
    p_normal <- ggplot(normal_data_single, aes(x = expression)) +
      geom_histogram(fill = "steelblue", bins = nclass.Sturges(normal_data_single$expression)) +
      labs(title = paste(g_symbol, "Expression in Normal Cell Lines"),
           x = "Expression (TPM)", y = "Number of Cell Lines") +
      theme_minimal()
    ggsave(file.path(output_dir, paste0(g_symbol, "_Normal_expression.pdf")), p_normal, width = 8, height = 6)
  } else {
    cat(paste0("\nNo Normal cell lines found for gene ", g_symbol, "\n"))
  }
  
  # Comparative plot + stats for each gene (optional)
  if (nrow(cancer_data_single) > 0 && nrow(normal_data_single) > 0) {
    combined_single_gene <- bind_rows(
      mutate(cancer_data_single, group = "Cancer"),
      mutate(normal_data_single, group = "Normal")
    )
    
    cat(paste0("\n=== Comparative Analysis for ", g_symbol, " ===\n"))
    comp_plot_single_gene <- ggplot(combined_single_gene, aes(x = group, y = expression, fill = group)) +
      geom_boxplot() +
      scale_fill_manual(values = c("Cancer" = "firebrick", "Normal" = "steelblue")) +
      labs(title = paste(g_symbol, "Expression: Cancer vs Normal"),
           x = "Cell Type", y = "Expression (TPM)") +
      theme_minimal() +
      stat_compare_means(method = "t.test", label = "p.format")
    
    ggsave(file.path(output_dir, paste0(g_symbol, "_cancer_vs_normal.pdf")), comp_plot_single_gene, width = 8, height = 6)
    
    t_test_single <- t.test(expression ~ group, data = combined_single_gene)
    wilcox_single <- wilcox.test(expression ~ group, data = combined_single_gene)
    
    cat("\nT-test results:\n")
    print(t_test_single)
    cat("\nWilcoxon test results:\n")
    print(wilcox_single)
  }
}

# --- START: STATISTICAL COMPARISON ACROSS PRIMARY DISEASE GROUPS ---
cat("\nPerforming statistical comparisons for each gene across PRIMARY DISEASE groups...\n")

# DEBUG: Check if validated_gene_symbols is populated
cat("\nGenes selected for statistical testing:", paste(validated_gene_symbols, collapse = ", "), "\n")
if (length(validated_gene_symbols) == 0) {
  cat("WARNING: No validated gene symbols found for statistical analysis. Check previous steps.\n")
}

# DEBUG: Check initial state of analysis_data_multi_gene
cat("\nDimensions of analysis_data_multi_gene:", dim(analysis_data_multi_gene), "\n")
cat("Head of analysis_data_multi_gene:\n")
print(head(analysis_data_multi_gene))
cat("Table of primary_disease in analysis_data_multi_gene (before filtering for stats):\n") # Changed label
print(table(analysis_data_multi_gene$primary_disease))


# --- NEW: Filter out primary_disease groups with too few samples for statistical analysis ---
min_samples_per_group_for_stats <- 3 # You can change this value if needed (e.g., to 5)
cat(paste0("\nFiltering primary disease groups to ensure at least ", min_samples_per_group_for_stats, " samples for statistical analysis...\n"))

# Calculate counts per primary_disease group
disease_group_counts <- analysis_data_multi_gene %>%
  count(primary_disease) %>%
  filter(n >= min_samples_per_group_for_stats)

# Get list of primary_disease groups that meet the criteria
valid_disease_groups <- disease_group_counts$primary_disease

# Filter analysis_data_multi_gene to include only these valid groups for statistics
analysis_data_multi_gene_filtered_for_stats <- analysis_data_multi_gene %>%
  filter(primary_disease %in% valid_disease_groups) %>%
  # Optionally, drop unused factor levels after filtering for cleaner output
  mutate(primary_disease = droplevels(as.factor(primary_disease)))


if (nrow(analysis_data_multi_gene_filtered_for_stats) == 0) {
  stop(paste0("After filtering, no cell lines remain with primary disease groups having at least ", min_samples_per_group_for_stats, " samples. Adjust 'min_samples_per_group_for_stats' or check data."))
}
cat(paste0("Number of cell lines after filtering for min samples per group (for stats): ", nrow(analysis_data_multi_gene_filtered_for_stats), "\n"))
cat("Primary Disease groups included in statistical analysis:\n")
print(table(analysis_data_multi_gene_filtered_for_stats$primary_disease))

# Initialize a list to store results for each gene
gene_comparison_results <- list()

# Loop through each gene to perform statistical tests against primary_disease
for (gene in validated_gene_symbols) {
  cat(paste0("\n--- Debugging Gene: ", gene, " ---\n"))
  
  # Create a temporary data frame for the current gene's expression and primary_disease
  # NOW using 'analysis_data_multi_gene_filtered_for_stats'
  current_gene_data <- analysis_data_multi_gene_filtered_for_stats %>%
    dplyr::select(Expression = !!sym(gene), DiseaseGroup = primary_disease) %>%
    drop_na(Expression, DiseaseGroup) # Ensure no NAs in current gene's expression or primary_disease
  
  # DEBUG: Check current_gene_data after NA removal
  cat(paste0("  Rows in current_gene_data for ", gene, " after NA removal: ", nrow(current_gene_data), "\n"))
  if (nrow(current_gene_data) == 0) {
    cat(paste0("  Skipping gene '", gene, "' as current_gene_data is empty after NA removal.\n"))
    gene_comparison_results[[gene]] <- "No data after NA removal."
    next
  }
  
  # Convert DiseaseGroup to a factor (important for statistical tests)
  current_gene_data$DiseaseGroup <- as.factor(current_gene_data$DiseaseGroup)
  
  # DEBUG: Check DiseaseGroup distribution for the current gene
  cat(paste0("  DiseaseGroup distribution for ", gene, ":\n")) # Changed
  print(table(current_gene_data$DiseaseGroup)) # Changed
  cat(paste0("  Unique DiseaseGroups: ", length(unique(current_gene_data$DiseaseGroup)), "\n")) # Changed
  cat(paste0("  Min samples per DiseaseGroup: ", min(table(current_gene_data$DiseaseGroup)), "\n")) # Changed
  
  
  # Check if there are enough data points and DiseaseGroups for the test
  # This condition now uses 'min_samples_per_group_for_stats' and should ideally not trigger often
  if (length(unique(current_gene_data$DiseaseGroup)) < 2 || min(table(current_gene_data$DiseaseGroup)) < min_samples_per_group_for_stats) { # CHANGED to use new variable
    cat(paste0("  Skipping gene '", gene, "' due to insufficient data or DiseaseGroups (less than 2 unique DiseaseGroups or less than ", min_samples_per_group_for_stats, " observations in a DiseaseGroup after NA removal, despite initial filtering).\n")) # Changed
    gene_comparison_results[[gene]] <- "Insufficient data for comparison."
    next
  }
  
  # Option 1: Kruskal-Wallis Test (Non-parametric, robust to non-normal data)
  # CHANGED: formula to use 'DiseaseGroup'
  kw_test_result <- kruskal.test(Expression ~ DiseaseGroup, data = current_gene_data)
  cat(paste0("  Kruskal-Wallis p-value for ", gene, ": ", round(kw_test_result$p.value, 6), "\n"))
  
  
  # If Kruskal-Wallis is significant, perform Dunn's post-hoc test
  if (kw_test_result$p.value < 0.05) { # Using 0.05 as a significance threshold for overall test
    # Dunn's test (with Benjamini-Hochberg FDR correction for multiple comparisons)
    # CHANGED: formula to use 'DiseaseGroup'
    dunn_posthoc <- dunn_test(Expression ~ DiseaseGroup, data = current_gene_data, p.adjust.method = "BH")
    gene_comparison_results[[gene]] <- list(
      Overall_KruskalWallis = kw_test_result,
      Dunn_PostHoc = dunn_posthoc
    )
    cat(paste0("  Gene '", gene, "': Significant Kruskal-Wallis (p=", round(kw_test_result$p.value, 4), "). Dunn's post-hoc performed.\n"))
    cat("  Dunn's post-hoc results:\n")
    print(dunn_posthoc) # DEBUG: Print post-hoc results directly
    
  } else {
    gene_comparison_results[[gene]] <- list(
      Overall_KruskalWallis = kw_test_result,
      Message = paste0("Kruskal-Wallis not significant (p=", round(kw_test_result$p.value, 4), "). No post-hoc needed.")
    )
    cat(paste0("  Gene '", gene, "': Kruskal-Wallis not significant (p=", round(kw_test_result$p.value, 4), ").\n"))
  }
}

cat("\nStatistical comparisons complete. Summarizing results...\n")

# DEBUG: Check the contents of gene_comparison_results
cat("\nContents of gene_comparison_results (first few entries):\n")
print(head(gene_comparison_results))


# --- NEW: Calculate and save Mean Expression per Primary Disease Group ---
cat("\nCalculating mean expression per gene per primary disease group...\n")
mean_expression_df_temp <- analysis_data_multi_gene_filtered_for_stats %>% # Use the filtered data consistent with stats
  pivot_longer(
    cols = all_of(validated_gene_symbols),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  group_by(Gene, primary_disease) %>%
  summarise(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = 'drop') %>%
  arrange(Gene, primary_disease)

write.csv(mean_expression_df_temp, file.path(output_dir, "Mean_Expression_by_PrimaryDisease_Group.csv"), row.names = FALSE)
cat(paste0("Mean expression per primary disease group saved to: ", file.path(output_dir, "Mean_Expression_by_PrimaryDisease_Group.csv"), "\n"))


# Summarize and save results
kruskal_summary_df <- data.frame(
  Gene = character(),
  Kruskal_P_Value = numeric(),
  stringsAsFactors = FALSE
)

dunn_results_df_temp <- data.frame() # Use a temporary name to avoid conflict with read later

for (gene in names(gene_comparison_results)) {
  res <- gene_comparison_results[[gene]]
  
  # Ensure it's a list with KruskalWallis results, not a skip message
  if (is.list(res) && "Overall_KruskalWallis" %in% names(res)) {
    kw_p_value <- res$Overall_KruskalWallis$p.value
    kruskal_summary_df <- rbind(kruskal_summary_df, data.frame(Gene = gene, Kruskal_P_Value = kw_p_value))
    
    if ("Dunn_PostHoc" %in% names(res) && !is.null(res$Dunn_PostHoc) && inherits(res$Dunn_PostHoc, "data.frame")) { # Ensure it's a dataframe
      dunn_df_temp <- res$Dunn_PostHoc %>%
        as.data.frame() %>%
        mutate(Gene = gene, .before = 1) # Add gene column at the beginning
      dunn_results_df_temp <- rbind(dunn_results_df_temp, dunn_df_temp)
    } else if ("Dunn_PostHoc" %in% names(res) && !inherits(res$Dunn_PostHoc, "data.frame")) {
      cat(paste0("DEBUG: Dunn_PostHoc for gene '", gene, "' exists but is not a data.frame. Type: ", class(res$Dunn_PostHoc)[1], "\n"))
    }
  } else {
    cat(paste0("DEBUG: Gene '", gene, "' result skipped from summary (not a full result list).\n"))
  }
}

# DEBUG: Check final dimensions of summary dataframes
cat("\nFinal dimensions of kruskal_summary_df:", dim(kruskal_summary_df), "\n")
cat("Final dimensions of dunn_results_df_temp:", dim(dunn_results_df_temp), "\n")

# Order by p-value for overall Kruskal-Wallis test
kruskal_summary_df <- kruskal_summary_df %>% arrange(Kruskal_P_Value)

# Save the results - CHANGED FILENAMES
write.csv(kruskal_summary_df, file.path(output_dir, "Kruskal_Wallis_Overall_Gene_Comparison_Results_PrimaryDisease.csv"), row.names = FALSE)
cat(paste0("\nOverall Kruskal-Wallis results saved to: ", file.path(output_dir, "Kruskal_Wallis_Overall_Gene_Comparison_Results_PrimaryDisease.csv"), "\n"))

if (nrow(dunn_results_df_temp) > 0) {
  write.csv(dunn_results_df_temp, file.path(output_dir, "Dunn_PostHoc_Pairwise_Comparison_Results_PrimaryDisease.csv"), row.names = FALSE)
  cat(paste0("Dunn's post-hoc results saved to: ", file.path(output_dir, "Dunn_PostHoc_Pairwise_Comparison_Results_PrimaryDisease.csv"), "\n"))
} else {
  cat("\nNo significant Kruskal-Wallis results found, so no Dunn's post-hoc file generated.\n")
}

# --- Optional: Visualize a significant gene with a boxplot (REMOVED as requested) ---

# --- END: STATISTICAL COMPARISON ACROSS PRIMARY DISEASE GROUPS ---

# --- further statistical analysis (Post-processing of generated CSVs) ---
cat("\n--- Starting Post-analysis: Identifying Significantly Higher Expression ---\n")

# 1. No need to load libraries again if they are already loaded globally at the top
# 2. No need to set working directory again. It's already correctly set to 'working_dir'
#    and files are accessed via 'output_dir' path.
# setwd(output_dir) # REMOVED: This line is redundant and can be confusing.

# 3. Load the generated CSV files
# Check if Mean_Expression_by_PrimaryDisease_Group.csv exists before reading
mean_expr_file_path <- file.path(output_dir, "Mean_Expression_by_PrimaryDisease_Group.csv")
if (file.exists(mean_expr_file_path)) {
  mean_expression_df <- read.csv(mean_expr_file_path)
  cat(paste0("Loaded: ", mean_expr_file_path, "\n"))
} else {
  stop(paste0("Error: Mean_Expression_by_PrimaryDisease_Group.csv not found at ", mean_expr_file_path, ".\n",
              "Please ensure the main script ran successfully and generated this file."))
}


# Check if Dunn_PostHoc_Pairwise_Comparison_Results_PrimaryDisease.csv exists before reading
dunn_results_file_path <- file.path(output_dir, "Dunn_PostHoc_Pairwise_Comparison_Results_PrimaryDisease.csv")
if (file.exists(dunn_results_file_path)) {
  dunn_results_df <- read.csv(dunn_results_file_path)
  cat(paste0("Loaded: ", dunn_results_file_path, "\n"))
} else {
  cat(paste0("Warning: Dunn_PostHoc_Pairwise_Comparison_Results_PrimaryDisease.csv not found at ", dunn_results_file_path, ".\n",
             "This file is only generated if at least one gene had a significant Kruskal-Wallis test.\n",
             "Skipping 'Significantly Higher Expression' analysis as pairwise p-values are not available.\n"))
  # Create an empty dataframe to avoid errors later if this file is missing
  significantly_higher_expression <- data.frame() # Initialize as empty
}


# 4. Define your significance threshold (e.g., 0.05 for p.adj)
significance_threshold <- 0.05

# 5. Filter for statistically significant pairwise comparisons and combine with mean expression data
# Only proceed if dunn_results_df has data (i.e., if the file was loaded)
if (exists("dunn_results_df") && nrow(dunn_results_df) > 0) {
  significantly_higher_expression <- dunn_results_df %>%
    # Filter for comparisons that are statistically significant
    filter(p.adj < significance_threshold) %>%
    
    # Join with the mean expression data for group1
    left_join(mean_expression_df, by = c("Gene" = "Gene", "group1" = "primary_disease")) %>%
    rename(Mean_Expression_group1 = Mean_Expression) %>%
    
    # Join with the mean expression data for group2
    left_join(mean_expression_df, by = c("Gene" = "Gene", "group2" = "primary_disease")) %>%
    rename(Mean_Expression_group2 = Mean_Expression) %>%
    
    # Determine which group has the significantly higher expression
    mutate(
      Higher_Expression_Group = case_when(
        Mean_Expression_group1 > Mean_Expression_group2 ~ group1,
        Mean_Expression_group2 > Mean_Expression_group1 ~ group2,
        TRUE ~ NA_character_ # Should ideally not happen for significant differences
      ),
      Lower_Expression_Group = case_when(
        Mean_Expression_group1 > Mean_Expression_group2 ~ group2,
        Mean_Expression_group2 > Mean_Expression_group1 ~ group1,
        TRUE ~ NA_character_
      )
    ) %>%
    
    # Re-evaluate Mean_Higher and Mean_Lower more accurately using rowwise for dynamic selection
    rowwise() %>% # Apply the following row by row
    mutate(
      Mean_Higher = ifelse(Higher_Expression_Group == group1, Mean_Expression_group1, Mean_Expression_group2),
      Mean_Lower = ifelse(Lower_Expression_Group == group1, Mean_Expression_group1, Mean_Expression_group2)
    ) %>%
    ungroup() %>% # Remove rowwise grouping
    
    # Select and arrange relevant columns for the final summary
    dplyr::select(
      Gene,
      Higher_Expression_Group,
      Lower_Expression_Group,
      Mean_Higher,
      Mean_Lower,
      p.adj,
      p.adj.signif,
      statistic # Include statistic for context if desired
    ) %>%
    
    # Filter out rows where Higher_Expression_Group couldn't be determined (e.g., if means are exactly equal or NA)
    drop_na(Higher_Expression_Group) %>%
    
    # Sort by gene and then by significance for easier review
    arrange(Gene, p.adj)
  
  # 6. Print the summary and save to a new CSV file
  cat("\n--- Summary of Significantly Higher Gene Expression per Primary Disease Group ---\n")
  if (nrow(significantly_higher_expression) > 0) {
    print(significantly_higher_expression)
    output_filename_higher <- file.path(output_dir, "Significantly_Higher_Expression_Summary_PrimaryDisease.csv")
    write.csv(significantly_higher_expression, output_filename_higher, row.names = FALSE)
    cat(paste0("\nSummary saved to: ", output_filename_higher, "\n"))
  } else {
    cat("No significant 'higher' expression differences found at the chosen threshold (p.adj < ", significance_threshold, ").\n")
  }
  
} else {
  cat("\n'Significantly Higher Expression' analysis skipped because Dunn's post-hoc results were not available or empty.\n")
  # Initialize an empty dataframe to signify no results
  significantly_higher_expression <- data.frame(
    Gene = character(),
    Higher_Expression_Group = character(),
    Lower_Expression_Group = character(),
    Mean_Higher = numeric(),
    Mean_Lower = numeric(),
    p.adj = numeric(),
    p.adj.signif = character(),
    statistic = numeric(),
    stringsAsFactors = FALSE
  )
}


# --- End of further analysis --- #

# Prepare data for combined heatmap
# The `analysis_data_multi_gene` now contains all cell lines and all gene expressions.
# We need to make it into a matrix where rows are cell lines and columns are genes.

# FIRST, ensure cell_line_name is not NA. If it is, use depmap_id as a fallback.
# This prevents make.unique from creating NAs if original cell_line_name was NA.
analysis_data_multi_gene_cleaned <- analysis_data_multi_gene %>%
  mutate(
    cell_line_name_for_rownames = ifelse(is.na(cell_line_name) | cell_line_name == "", depmap_id, cell_line_name)
  )

# THEN, apply make.unique to the cleaned row names
analysis_data_for_heatmap <- analysis_data_multi_gene_cleaned %>%
  mutate(cell_line_name_unique = make.unique(as.character(cell_line_name_for_rownames))) %>%
  dplyr::select(cell_line_name_unique, primary_disease, lineage, cell_type, subgroup, all_of(validated_gene_symbols))

# Now create the expression matrix
expr_mat_combined <- as.matrix(analysis_data_for_heatmap %>% dplyr::select(all_of(validated_gene_symbols)))
rownames(expr_mat_combined) <- analysis_data_for_heatmap$cell_line_name_unique

# Annotation data for rows (cell lines)
annotation_row_combined <- analysis_data_for_heatmap %>%
  dplyr::select(cell_type, subgroup, primary_disease)
rownames(annotation_row_combined) <- analysis_data_for_heatmap$cell_line_name_unique

# --- Rest of your code remains the same from here ---

# Set cell_type, primary_disease, and subgroup as factors with desired order
annotation_row_combined$cell_type <- factor(annotation_row_combined$cell_type, levels = c("Cancer", "Normal"))
annotation_row_combined$primary_disease <- factor(annotation_row_combined$primary_disease, levels = sort(unique(annotation_row_combined$primary_disease)))
annotation_row_combined$subgroup <- factor(annotation_row_combined$subgroup, levels = sort(unique(annotation_row_combined$subgroup)))

# Order the matrix and annotation data by cell_type and primary_disease
ordered_idx_combined <- order(annotation_row_combined$cell_type, annotation_row_combined$primary_disease)
annotation_row_combined <- annotation_row_combined[ordered_idx_combined, ]
expr_mat_combined <- expr_mat_combined[ordered_idx_combined, , drop = FALSE]

# Convert to character for row_split
annotation_row_combined$cell_type <- as.character(annotation_row_combined$cell_type)
annotation_row_combined$subgroup <- as.character(annotation_row_combined$subgroup)
annotation_row_combined$primary_disease <- as.character(annotation_row_combined$primary_disease) # This line was added in previous discussion

# Set up row_split for the heatmap
row_split_combined <- annotation_row_combined$primary_disease
# Define annotation colors
ann_colors_combined <- list(
  cell_type = c(Cancer = "firebrick", Normal = "steelblue"),
  subgroup = setNames(
    colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(unique(annotation_row_combined$subgroup))),
    unique(annotation_row_combined$subgroup)
  ),
  primary_disease = setNames(
    viridis_pal(option = "D")(length(unique(annotation_row_combined$primary_disease))),
    unique(annotation_row_combined$primary_disease)
  )
)

# Row annotation for the heatmap
row_anno_combined <- rowAnnotation(
  CellType = anno_simple(annotation_row_combined$cell_type, col = ann_colors_combined$cell_type),
  # Disease = anno_simple(annotation_row_combined$primary_disease, col = ann_colors_combined$primary_disease), # Commented out
  Subgroup = anno_simple(annotation_row_combined$subgroup, col = ann_colors_combined$subgroup),
  annotation_width = unit(c(0.2, 0.5), "cm"), # Adjusted to 2 columns
  gap = unit(2, "mm") # This gap is between CellType and Subgroup
)

# Heatmap for multiple genes
cat("\nGenerating combined heatmap...\n")

# Draw to RStudio plot window first
cat("\nShowing combined heatmap in RStudio Plots window...\n")
ht_combined <- Heatmap(
  expr_mat_combined,
  name = "Expression (TPM)",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  left_annotation = row_anno_combined,
  row_split = row_split_combined,
  cluster_row_slices = FALSE,
  column_title = NULL, # Title removed as per your request
  col = colorRampPalette(RColorBrewer::brewer.pal(11, "YlOrRd"))(100),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),
  width = ncol(expr_mat_combined) * unit(10, "mm")
  # To control the gap between subgroup and heatmap, add 'gap = unit(X, "mm")' here
  # E.g., gap = unit(5, "mm")
)
draw(ht_combined, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

# Now, open the PDF device and draw to it
pdf(file.path(output_dir, paste0("Combined_Gene_Expression_Heatmap_", paste(validated_gene_symbols, collapse = "_"), ".pdf")),
    width = 4 + 0.5 * length(validated_gene_symbols),
    height = 0.03 * nrow(expr_mat_combined) + 5)

draw(ht_combined, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

# Save session info
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info_multigene.txt"))
cat("\n=== Analysis complete ===\n")
cat("All outputs saved to:", output_dir, "\n")

# --- START: CORRELATION ANALYSIS FOR FAP AND KRT19 ---
cat("\nPerforming pairwise correlation analysis for FAP and KRT19...\n")

# Define the two genes for comparison
gene_A <- "FAP"
gene_B <- "KRT19"

# Check if both genes exist in the data
if (!all(c(gene_A, gene_B) %in% validated_gene_symbols)) {
  stop(paste0("Error: Both '", gene_A, "' and '", gene_B, "' must be in the input gene list for correlation analysis."))
}

# Use the data that has already been cleaned and merged
correlation_data <- analysis_data_multi_gene %>%
  dplyr::select(primary_disease, all_of(c(gene_A, gene_B))) %>%
  drop_na() # Remove any rows with missing expression data for these two genes

# Initialize a list to store results
correlation_results_by_disease <- list()
all_correlation_data <- data.frame()

# Loop through each unique cancer type (primary_disease)
unique_diseases <- unique(correlation_data$primary_disease)
for (disease in unique_diseases) {
  # Skip "zNon-Cancerous" group
  if (disease == "zNon-Cancerous") {
    next
  }
  
  disease_data <- correlation_data %>%
    filter(primary_disease == disease)
  
  # Check if there are enough samples for correlation (at least 3)
  if (nrow(disease_data) < 3) {
    cat(paste0("Skipping '", disease, "' due to insufficient samples (", nrow(disease_data), ").\n"))
    next
  }
  
  # Perform Spearman correlation test
  cor_test_result <- cor.test(
    x = disease_data[[gene_A]],
    y = disease_data[[gene_B]],
    method = "spearman"
  )
  
  # Store the result
  correlation_results_by_disease[[disease]] <- list(
    correlation = cor_test_result$estimate,
    p_value = cor_test_result$p.value,
    count = nrow(disease_data)
  )
  
  # Prepare data for plotting later (optional, but good for visualization)
  plot_df <- disease_data %>%
    mutate(
      primary_disease = disease,
      label = paste0(
        "R = ", round(cor_test_result$estimate, 2),
        ", p = ", format.pval(cor_test_result$p.value, digits = 2)
      )
    )
  all_correlation_data <- bind_rows(all_correlation_data, plot_df)
}

# Convert results to a data frame for easy viewing and saving
correlation_summary_df <- tibble(
  Primary_Disease = names(correlation_results_by_disease),
  Correlation_R = map_dbl(correlation_results_by_disease, "correlation"),
  P_Value = map_dbl(correlation_results_by_disease, "p_value"),
  Sample_Count = map_dbl(correlation_results_by_disease, "count")
) %>%
  arrange(P_Value) %>% # Sort by significance
  mutate(Significance = case_when(
    P_Value < 0.001 ~ "***",
    P_Value < 0.01 ~ "**",
    P_Value < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

# Save the summary table to a CSV file
output_filename_cor <- file.path(output_dir, paste0(gene_A, "_vs_", gene_B, "_Correlation_by_Disease.csv"))
write.csv(correlation_summary_df, output_filename_cor, row.names = FALSE)
cat(paste0("\nCorrelation summary saved to: ", output_filename_cor, "\n"))

# Print the summary
cat("\n--- Summary of Correlation between", gene_A, "and", gene_B, "by Cancer Type ---\n")
print(correlation_summary_df)

# Visualize the correlations using a scatter plot for each disease type
if (nrow(all_correlation_data) > 0) {
  cor_plot <- ggplot(all_correlation_data, aes(x = .data[[gene_A]], y = .data[[gene_B]])) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "darkred") +
    geom_text(aes(x = -Inf, y = Inf, label = label),
              hjust = -0.1, vjust = 1.2,
              size = 3, fontface = "italic") +
    facet_wrap(~ primary_disease, scales = "free") +
    labs(
      title = paste0("Expression Correlation of ", gene_A, " vs ", gene_B, " by Cancer Type"),
      x = paste0(gene_A, " Expression (TPM)"),
      y = paste0(gene_B, " Expression (TPM)")
    ) +
    theme_minimal()
  
  output_filename_plot <- file.path(output_dir, paste0(gene_A, "_vs_", gene_B, "_Correlation_Plot.pdf"))
  ggsave(output_filename_plot, cor_plot, width = 12, height = 9)
  cat(paste0("\nCorrelation plot saved to: ", output_filename_plot, "\n"))
}

# --- END: CORRELATION ANALYSIS FOR FAP AND KRT19 ---
