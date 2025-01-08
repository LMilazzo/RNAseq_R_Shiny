
runLibs <- function() {
  # For shiny running
  library(shiny)
  library(shinythemes)
  library(shinyjs)
  library(shinyFiles)
  library(DT)
  library(htmlwidgets)
  
  # For data
  library(dplyr)
  library(tidyr)
  library(readr)
  library(DESeq2)
  library(SummarizedExperiment)
  library(BioVis)
  library(pathfindR)
  
  # For plotting
  library(ggplot2)
  library(ggplotify)
  library(pheatmap)
  library(ggbeeswarm)
  library(ggrepel)
  library(plotly)
  library(grid)
  
  # Additional
  library(zip)
}
