library(tidyverse)
library(rstatix)
library(ggalluvial)
library(gt)
library(spatstat.geom)
library(vroom)
library(progressr)
library(readr)

# Set the working directory to the location of the script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --- Check for saved RDS file or Load and Process Data ---
rds_file <- "allCells.rds"
if (file.exists(rds_file)) {
  cat("RDS file found. Loading data from", rds_file, "...\n")
  allCells <- readRDS(rds_file)
} else {
  cat("RDS file not found. Processing raw data...\n")
  
  # --- File Paths and Data Loading ---
  batch1_file <- "230920_batch1_multiplex_data.tsv"
  batch2_file <- "230921_batch2_multiplex_data.tsv"
  cohort_file <- "subtype_and_sample_type.csv"
  
  # Load and process Batch 1 data with a progress bar
  Batch1Data <- vroom(batch1_file, progress = TRUE) %>%
    select(ImageName = Name,
           FOV = Image,
           Phenotype = `ObjectInfoNuclei - LabelName`,
           Region = `ObjectInfoNuclei - ROIName`,
           x = CellX,
           y = CellY) %>%
    filter(!str_detect(FOV,"027.+PLEX_5|055.+PLEX_1|046.+53140")) %>%
    filter(Phenotype != "Class 2", Region != "Clear") %>%
    mutate(ImageName = str_replace(ImageName, "LIVER_RB045-UB000651", "LIVER_RB045-2-UB000651"),
           ImageName = str_replace(ImageName, "(RB[0-9]{2,3})(-[A-Z])(.+)", "\\1-1\\2"),
           PatientID = str_extract(ImageName, "RB[0-9]{2,3}-[12]"),
           PatientID = str_replace(PatientID, "(RB)([0-9]{2})-", "\\10\\2-"),
           Phenotype = str_remove_all(Phenotype, "\\(Opal [0-9]{3}\\)| "),
           FOV = str_replace_all(FOV, "\\\\", "/") %>% basename())
  
  # Load and process Batch 2 data with a progress bar
  Batch2Data <- vroom(batch2_file, progress = TRUE) %>%
    select(ImageName = Name,
           FOV = Image,
           Phenotype = `ObjectInfo - LabelName`,
           Region = `ObjectInfo - ROIName`,
           x = CellX,
           y = CellY) %>%
    filter(!str_detect(FOV,"1042.+23038|660931.+20912")) %>%
    filter(Phenotype != "Class 2", Region != "Clear") %>%
    mutate(PatientID = str_extract(ImageName, "LIVER_([0-9]+)-", group = 1),
           Phenotype = str_remove_all(Phenotype, "\\(Opal [0-9]{3}\\)| "))
  
  # Combine data from both batches
  allCells <- rbind(Batch1Data, Batch2Data) %>%
    mutate(PatientID = case_when(PatientID == "61434" ~ "614634",
                                 PatientID =="59640" ~ "596640",
                                 PatientID =="660931" ~ "66093",
                                 .default = PatientID))
  
  # Add cohort data
  cohortData <- read_csv(cohort_file) %>%
    mutate(PatientID = ifelse(str_detect(PatientID, "RB") & !str_detect(PatientID, "-"),
                              paste0(PatientID, "-1"),
                              PatientID))
  
  allCells <- allCells %>% left_join(cohortData)
  
  # Save the combined data as an RDS file for faster loading next time
  saveRDS(allCells, rds_file)
  cat("Data processed and saved to", rds_file, "for future use.\n\n")
}

# --- Interactive Phenotype Selection ---

# Get unique phenotypes for user to choose from
all_phenotypes <- unique(allCells$Phenotype)
# Exclude "Class 2" and "Clear" if they somehow made it this far
all_phenotypes <- all_phenotypes[all_phenotypes != "Class 2" & all_phenotypes != "Clear"]

cat("--- Available Phenotypes (Reference) ---\n")
for (i in 1:length(all_phenotypes)) {
  cat(sprintf("%d: %s\n", i, all_phenotypes[i]))
}

# Get user input for reference phenotype and validate it
valid_choice <- FALSE
while (!valid_choice) {
  user_input <- readline(prompt="Enter the number for the phenotype you want to use as REFERENCE: ")
  choice_index <- as.integer(user_input)
  
  if (!is.na(choice_index) && choice_index >= 1 && choice_index <= length(all_phenotypes)) {
    reference_phenotype <- all_phenotypes[choice_index]
    cat("You have selected:", reference_phenotype, "\n\n")
    valid_choice <- TRUE
  } else {
    cat("Invalid selection. Please enter a number from the list.\n")
  }
}

# Get a new list of phenotypes to choose from for plotting
phenotypes_for_plotting <- all_phenotypes[all_phenotypes != reference_phenotype]
cat("--- Available Phenotypes (for plotting against reference) ---\n")
for (i in 1:length(phenotypes_for_plotting)) {
  cat(sprintf("%d: %s\n", i, phenotypes_for_plotting[i]))
}

# Get user input for phenotypes to plot and validate it
valid_plot_choices <- FALSE
while (!valid_plot_choices) {
  user_input_plot <- readline(prompt="Enter the numbers (e.g., 1, 3, 5) for the phenotypes you want to plot: ")
  choice_indices_plot <- as.integer(strsplit(gsub(" ", "", user_input_plot), ",")[[1]])
  
  # Check if all inputs are valid numbers and within range
  if (all(!is.na(choice_indices_plot)) && all(choice_indices_plot >= 1) && all(choice_indices_plot <= length(phenotypes_for_plotting))) {
    phenotypes_to_analyze <- phenotypes_for_plotting[choice_indices_plot]
    cat("You have selected to plot:", paste(phenotypes_to_analyze, collapse = ", "), "\n\n")
    valid_plot_choices <- TRUE
  } else {
    cat("Invalid selection. Please enter one or more numbers from the list, separated by commas.\n")
  }
}

# --- Closest Neighbor Analysis Functions ---

# Define the function to get distance data for a given FOV and two phenotypes
getDistanceData <- function(data, Phenotype1, Phenotype2, k = 1) {
  idColString = ifelse(k == 1, "which", paste0("which.", k))
  
  window <- owin(xrange = c(min(data$x), max(data$x)),
                 yrange = c(min(data$y), max(data$y)))
  
  cellType1 <- data %>% filter(Phenotype == Phenotype1) %>% drop_na()
  cellType2 <- data %>% filter(Phenotype == Phenotype2) %>% drop_na()
  
  if (nrow(cellType1) == 0 | nrow(cellType2) == 0) {
    return(NULL)
  }
  
  distData <- nncross(ppp(cellType1$x, cellType1$y, window),
                      ppp(cellType2$x, cellType2$y, window),
                      k = k) %>%
    rename_with(~str_remove(., paste0("\\.", k))) %>%
    mutate(target = Phenotype2)
  
  distData <- distData %>%
    mutate(from.x = cellType1$x,
           from.y = cellType1$y,
           to.x = cellType2$x[which],
           to.y = cellType2$y[which]) %>%
    cbind(cellType1 %>% select(Region, PatientID, subtype, sampleType = `sample type`))
  
  distData
}

# --- Perform the analysis and generate the plot ---

# Filter and nest the data for analysis using the selected reference phenotype
nestedCells <- allCells %>%
  filter(`sample type` == "resection") %>%
  # Note: Retained the original fix for FAP since it was specifically requested
  mutate(Phenotype = ifelse(str_detect(Phenotype, "FAP"), "FAP", Phenotype)) %>%
  nest(.by = FOV)

# Map over the nested data and phenotypes to calculate distances, filtering out NULL results
handlers(handler_progress())
distanceDataframe <- with_progress({
  p <- progressor(steps = nrow(nestedCells))
  
  nestedCells %>%
    pmap_dfr(function(FOV, data) {
      p()
      map_dfr(phenotypes_to_analyze, function(phenotype) {
        dist_data <- data %>% getDistanceData(reference_phenotype, phenotype)
        if (!is.null(dist_data)) {
          dist_data %>% mutate(FOV = FOV)
        }
      })
    })
})

# --- Visualizing the groups ---

# Generate the density plot
distanceDataframe %>%
  ggplot(aes(x = dist, color = target)) +
  geom_density() +
  scale_x_log10() +
  labs(x = "Distance (log10)", y = "Density", color = "Target Phenotype") +
  theme_minimal() +
  ggtitle(paste("Density Plot of Distances from", reference_phenotype))

# --- Statistical Analysis and Saving Results ---

# Aggregate distances per PatientID, FOV, and target phenotype
aggDistances <- distanceDataframe %>%
  group_by(PatientID, FOV, target) %>%
  summarise(median_dist = median(dist, na.rm = TRUE), .groups = "drop")

cat("\n--- Kruskal-Wallis Test on aggregated medians ---\n")
kruskal_test_result <- kruskal_test(aggDistances, median_dist ~ target)
print(kruskal_test_result)
write_csv(kruskal_test_result, "kruskal_test_results_aggregated.csv")
cat("\nKruskal-Wallis results saved to kruskal_test_results_aggregated.csv\n\n")

# --- Statistical Analysis and Saving Results ---

# Aggregate distances per PatientID, FOV, and target phenotype
aggDistances <- distanceDataframe %>%
  group_by(PatientID, FOV, target) %>%
  summarise(median_dist = median(dist, na.rm = TRUE), .groups = "drop")

cat("\n--- Kruskal-Wallis Test on aggregated medians ---\n")
kruskal_test_result <- kruskal_test(aggDistances, median_dist ~ target)
print(kruskal_test_result)
write_csv(kruskal_test_result, "kruskal_test_results_aggregated.csv")
cat("\nKruskal-Wallis results saved to kruskal_test_results_aggregated.csv\n\n")

# If the Kruskal-Wallis test is significant, perform Dunn's test
if (kruskal_test_result$p < 0.05) {
  
  cat("--- Dunn's Post-Hoc Test on aggregated medians (all pairwise comparisons) ---\n")
  
  # Dunn’s test (pairwise Wilcoxon) across all pairs
  dunn_test_result <- pairwise_wilcox_test(
    aggDistances, 
    median_dist ~ target, 
    p.adjust.method = "bonferroni"
  )
  
  # Effect sizes for all pairs
  effect_size_result <- wilcox_effsize(
    aggDistances, 
    median_dist ~ target
  ) %>% 
    select(group1, group2, effsize = effsize)
  
  dunn_test_result <- dunn_test_result %>%
    left_join(effect_size_result, by = c("group1", "group2"))
  
  # Save results
  write_csv(dunn_test_result, "dunn_test_results_aggregated.csv")
  print(dunn_test_result)
  cat("\nDunn's test results (all pairwise comparisons) saved to dunn_test_results_aggregated.csv\n\n")
  
} else {
  cat("The Kruskal-Wallis test was not significant (p >= 0.05), so no post-hoc test is needed.\n\n")
}
