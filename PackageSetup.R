# List of required packages
packages <- c("shiny", "shinythemes", "DT", "dplyr", "readr", 
              "ggplot2", "ggplotify", "patchwork", "matrixStats", "circlize", 
              "colourpicker", "ggbeeswarm", "ggrepel", 
              "pathfindR", "tidyr", "shinyjs", "plotly", 
              "htmlwidgets", "png", "leaflet", "devtools", "shiny")

# Function to check and install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = FALSE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(lapply(packages, install_if_missing))

if (!requireNamespace("BiocManager", quietly = FALSE)){
  install.packages("BiocManager")
}

bioc_packages <- c("DESeq2", "ComplexHeatmap", "SummarizedExperiment")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = FALSE)) {
    BiocManager::install(pkg)
  }
}

if (!requireNamespace("BioVis", quietly = FALSE)){
  devtools::install_github("LMilazzo/BioVis")
}

cat("All required packages are installed.\n")
