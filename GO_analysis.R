# -----------------------------
# Global options
# -----------------------------
options(future.globals.maxSize = 8000 * 1024^2)
options(scipen = 100)

# -----------------------------
# 0. Load packages
# -----------------------------
required_cran <- c("tidyverse", "rstudioapi", "httr", "jsonlite")
for(pkg in required_cran){
  if(!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "GO.db")
for(pkg in bioc_pkgs){
  if(!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# -----------------------------
# 1. Pick CSV interactively
# -----------------------------
data_file <- rstudioapi::selectFile(
  caption = "Select cytokine CSV file", 
  filter = "CSV Files (*.csv)"
)
data <- read.csv(data_file, stringsAsFactors = FALSE)

# Rename first column as Cytokine
colnames(data)[1] <- "Cytokine"

# Detect treatment columns automatically (all except Cytokine)
treatment_cols <- colnames(data)[-1]

# -----------------------------
# 2. Enumerate columns and ask user to choose
# -----------------------------
cat("Available treatments:\n")
for(i in seq_along(treatment_cols)){
  cat(i, ":", treatment_cols[i], "\n")
}

selected_numbers <- readline(
  prompt = "Enter the number(s) of the treatment(s) to analyze, separated by commas: "
)
selected_numbers <- as.numeric(unlist(strsplit(selected_numbers, ",")))
selected_treatments <- treatment_cols[selected_numbers]
cat("Selected treatments:", paste(selected_treatments, collapse = ", "), "\n")

# -----------------------------
# 3. Prepare fold-change data
# -----------------------------
fc_data <- dplyr::select(data, Cytokine, dplyr::all_of(selected_treatments))
rownames(fc_data) <- fc_data$Cytokine
fc_data <- dplyr::select(fc_data, -Cytokine)

# -----------------------------
# 4. Direction-specific GO enrichment
# -----------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

for(trt in selected_treatments){
  
  # Check if saved results exist
  res_file <- paste0(trt, "_GO_results.rds")
  
  if(file.exists(res_file)){
    cat("Saved GO results found for", trt, "- skipping enrichment calculation.\n")
    enrichment_res <- readRDS(res_file)
    ego_up <- enrichment_res$ego_up
    ego_down <- enrichment_res$ego_down
    up_genes <- enrichment_res$up_genes
    down_genes <- enrichment_res$down_genes
  } else {
    cat("=== Treatment:", trt, "===\n")
    
    up_genes <- rownames(fc_data)[fc_data[[trt]] > 1]
    down_genes <- rownames(fc_data)[fc_data[[trt]] < 1]
    
    pdf(file = paste0(trt, "_GO_Enrichment.pdf"), width = 10, height = 8)
    
    ego_up <- ego_down <- NULL
    
    # -----------------------------
    # GO enrichment: UP
    # -----------------------------
    if(length(up_genes) > 0){
      ego_up <- tryCatch({
        enrichGO(
          gene = up_genes,
          OrgDb = org.Hs.eg.db,
          keyType = "SYMBOL",
          ont = "BP",
          pAdjustMethod = "BH",
          qvalueCutoff = 0.05,
          readable = TRUE
        )
      }, error = function(e) NULL)
      
      if(!is.null(ego_up) && nrow(ego_up@result) > 0 && any(ego_up@result$p.adjust < 0.05)){
        p_go_up <- dotplot(ego_up, showCategory = 20, font.size = 8) + 
          ggtitle(paste(trt, "- Upregulated Cytokines GO")) +
          theme(plot.title = element_text(size = 12))
        print(p_go_up)
        cat("GO UP: Found", nrow(ego_up@result), "enriched terms\n")
      } else {
        plot.new()
        text(0.5, 0.5, paste(trt, "- No Significant GO Terms\nfor Upregulated Genes"), cex=1.2)
        cat("GO UP: No significant terms found\n")
      }
    } else {
      plot.new()
      text(0.5, 0.5, paste(trt, "- No Upregulated Genes for GO"), cex=1.2)
      cat("GO UP: No genes to analyze\n")
    }
    
    # -----------------------------
    # GO enrichment: DOWN
    # -----------------------------
    if(length(down_genes) > 0){
      ego_down <- tryCatch({
        enrichGO(
          gene = down_genes,
          OrgDb = org.Hs.eg.db,
          keyType = "SYMBOL",
          ont = "BP",
          pAdjustMethod = "BH",
          qvalueCutoff = 0.05,
          readable = TRUE
        )
      }, error = function(e) NULL)
      
      if(!is.null(ego_down) && nrow(ego_down@result) > 0 && any(ego_down@result$p.adjust < 0.05)){
        p_go_down <- dotplot(ego_down, showCategory = 20, font.size = 8) + 
          ggtitle(paste(trt, "- Downregulated Cytokines GO")) +
          theme(plot.title = element_text(size = 12))
        print(p_go_down)
        cat("GO DOWN: Found", nrow(ego_down@result), "enriched terms\n")
      } else {
        plot.new()
        text(0.5, 0.5, paste(trt, "- No Significant GO Terms\nfor Downregulated Genes"), cex=1.2)
        cat("GO DOWN: No significant terms found\n")
      }
    } else {
      plot.new()
      text(0.5, 0.5, paste(trt, "- No Downregulated Genes for GO"), cex=1.2)
      cat("GO DOWN: No genes to analyze\n")
    }
    
    dev.off()
    cat("PDF saved:", paste0(trt, "_GO_Enrichment.pdf\n"))
    
    # -----------------------------
    # SAVE enrichment results for later use
    # -----------------------------
    saveRDS(
      list(ego_up = ego_up, ego_down = ego_down,
           up_genes = up_genes, down_genes = down_genes),
      file = paste0(trt, "_GO_results.rds")
    )
    cat("GO enrichment results saved:", paste0(trt, "_GO_results.rds\n"))
  }
}

# -----------------------------
# 5. Interactive compressed GO UP dotplot (module 4 aesthetics, compressed)
# -----------------------------
trt <- readline(prompt = "Enter the treatment name to plot compressed GO dotplot: ")
res_file <- paste0(trt, "_GO_results.rds")

if(file.exists(res_file)){
  enrichment_res <- readRDS(res_file)
  ego_up <- enrichment_res$ego_up
  up_genes <- enrichment_res$up_genes
  
  if(!is.null(ego_up) && length(up_genes) > 0){
    gene_ratio_cutoff <- as.numeric(readline(prompt = "Enter GeneRatio cutoff for compressed GO UP dotplot (e.g., 0.31): "))
    max_terms <- as.numeric(readline(prompt = "Enter the maximum number of GO terms to display (e.g., 20): "))
    
    ego_up@result$GeneRatioNum <- sapply(strsplit(as.character(ego_up@result$GeneRatio), "/"),
                                         function(x) as.numeric(x[1])/as.numeric(x[2]))
    ego_up@result <- ego_up@result[ego_up@result$GeneRatioNum > gene_ratio_cutoff, ]
    
    if(nrow(ego_up@result) > 0){
      ego_up@result$Description <- factor(
        ego_up@result$Description, 
        levels = ego_up@result$Description[order(ego_up@result$GeneRatioNum, decreasing = TRUE)]
      )
      
      show_n <- min(max_terms, nrow(ego_up@result))
      
      pdf(file = paste0(trt, "_GO_UP_dotplot_cutoff", gene_ratio_cutoff, "_compressed.pdf"),
          width = 6, height = max(2.5, show_n * 0.15))
      
      p_go_up <- dotplot(ego_up, showCategory = show_n, font.size = 10) +
        ggtitle(paste(trt, "- Upregulated Cytokines GO (GeneRatio >", gene_ratio_cutoff, ")")) +
        theme_minimal() +
        coord_flip() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(margin = margin(r = 2), size = 8),
          plot.title = element_text(size = 12),
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
        )
      
      print(p_go_up)
      dev.off()
      
      cat(trt, "- Compressed GO UP dotplot saved with GeneRatio >", gene_ratio_cutoff,
          "(top", show_n, "terms)\n")
    } else {
      cat(trt, "- No GO terms with GeneRatio >", gene_ratio_cutoff, "\n")
    }
  } else {
    cat(trt, "- No upregulated genes or GO results available\n")
  }
} else {
  cat("File", res_file, "not found. Please run module 4 first to generate GO results.\n")
}
