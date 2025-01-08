# List of required packages
packages <- c("shiny", "shinythemes", "DT", "dplyr", "readr", 
              "ggplot2", "ggplotify", "pheatmap", "ggbeeswarm", "ggrepel", 
              "pathfindR", "tidyr", "shinyjs", "plotly", 
              "htmlwidgets", "grid", "zip", "shinyFiles")

# Function to check and install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(lapply(packages, install_if_missing))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of Bioconductor packages
bioc_packages <- c("DESeq2", "SummarizedExperiment")

# Function to check and install missing Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Install BioVis from GitHub if not already installed
if (!requireNamespace("BioVis", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("LMilazzo/BioVis")
}

cat("All required packages are installed.\n")
