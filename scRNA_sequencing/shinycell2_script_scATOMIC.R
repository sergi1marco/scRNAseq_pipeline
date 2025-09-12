# ===========================================================
# Interactive ShinyCell2 app generator
# Opens in default browser + prepares plot folder
# ===========================================================

# 1. Install/load required packages
if (!requireNamespace("ShinyCell2", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("the-ouyang-lab/shinycell2")
}
if (!requireNamespace("rstudioapi", quietly = TRUE)) {
  install.packages("rstudioapi")
}

library(Seurat)
library(ShinyCell2)
library(shiny)
library(rstudioapi)

# 2. Prompt user to select Seurat object file (.rds)
seurat_path <- rstudioapi::selectFile(
  caption = "Select your Seurat object (.rds)",
  label = "Choose file",
  filter = "rds",
  existing = TRUE
)

if (is.null(seurat_path) || seurat_path == "") {
  stop("❌ No file selected. Exiting.")
}

message("✅ Selected file: ", seurat_path)

# Load object
seu <- readRDS(seurat_path)

# 3. Create ShinyCell2 config
scConf <- createConfig(seu)

# 4. Set scATOMIC_prediction as default metadata if available
if ("scATOMIC_prediction" %in% scConf$meta) {
  scConf$meta <- c("scATOMIC_prediction",
                   setdiff(scConf$meta, "scATOMIC_prediction"))
}

# 5. Define output directory for Shiny app + plots
base_dir <- dirname(seurat_path)
app_dir <- file.path(base_dir, "shinyApp")
plot_dir <- file.path(base_dir, "shinycell2_plots")

if (!dir.exists(app_dir)) dir.create(app_dir)
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# 6. Build ShinyCell2 app
makeShinyFiles(seu, scConf, shiny.prefix = "myApp", shiny.dir = app_dir)
makeShinyCodes(shiny.title = "Single-Cell App with scATOMIC",
               shiny.prefix = "myApp",
               shiny.dir = app_dir)

message("📂 A folder for saving plots has been created at: ", plot_dir)
message("👉 Use the download buttons in the Shiny app to save plots there.")

# 7. Launch the app in system default browser
options(shiny.launch.browser = TRUE)
shiny::runApp(app_dir)
