diffIntermediateTab <- function(id){
  
  ns <- NS(id)
  
  fluidPage(
    
    fluidRow(
      actionButton(
        ns('RunDeseq'), 
        'Run',
        style = "background-color: #4CAF50; color: black; 
          border-color: white; background-image: none;
          margin-top: 0px ; margin-bottom: 0px; margin-right: 5%;"
      ),
      style = "text-align: center; margin-top: 10px; margin-bottom: 10px"
    ),
    
    tabsetPanel(

      #Raw Data Preview----
      tabPanel(
        "Raw Data Preview",
        
        fluidRow(
          column(
            width = 9,
            div(style = "max-width: 100%; overflow-x: auto;", DTOutput(ns("RawDataTable")))
          ),
          
          column(
            width = 3,
            div(
              h4("Search"),
              textInput(ns("GenePreviewSummary"), "", width = "100%"),
              DTOutput(ns("SearchedGene")),
              style = "border: 2px solid black; padding: 5px; padding-bottom: 15px;
              background: #131517; border-radius: 15px; text-align: center;"
            )
          ),
          style = "margin-top: 10px;"
        )
      ),
      
      #Meta Data Preview----
      tabPanel(
        "Meta Data",
        
        fluidRow(
          column(
            width = 12,
            div(style = "max-width: 100%; overflow-x: auto;", DTOutput(ns("MetaDataDumpPreview")))
          )
        )
        
      )
    #----
    )
  )
}