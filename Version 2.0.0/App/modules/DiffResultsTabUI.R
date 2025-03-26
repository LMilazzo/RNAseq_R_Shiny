diffResultsTab <- function(id){
  
  ns <- NS(id)
  
  fluidPage(
    
    uiOutput(ns("ContrastSelectionBox")),
   
    tabsetPanel(
    
      # # Results / Summary TabPanel
      tabPanel(
        title = "Results / Summary",
        ResultsSummaryUI(ns("Results-Summary"))
      ),

      #Principle Components Analysis
      tabPanel(
        title = "PCA",
        ResultsPcaUI(ns("Results-PCA"))
      ),
      #
      # #Correlation Heatmap
      tabPanel(
       title = "Correlation",
       ResultsCorrelationUI(ns("Results-Corr"))
      ),
      #
      # # Gene Counts
      tabPanel(
        title = "View Gene Counts",
        ResultsGeneCountsUI(ns("Results-Count"))
      ),
      
      # Volcano Plot
      tabPanel(
        title = "Volcano Plot",
        ResultsVolcanoUI(ns("Results-Volcano"))
      )
      
    )
   
  ) 
  
}