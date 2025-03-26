DiffResultsServer <- function(id, Data, Uploads){
  moduleServer(id, function(input, output, session){
    
    ns <- session$ns
    
    #Data ~~~~~~~~~
    # DiffResults <- list(
    #    "contrasts" = list, can be null
    #    "ddsc" = Reactive, can be null
    #    "Results" = Reactive
    #    "vstObject" = Reactive
    #    "vstCounts" = Reactive
    #    "NormCounts" = Reactive
    #    "status" = Reactive 0, 1
    # )
    #~~~~~~~~~~~~~~~
    #Uploads ~~~~~~~~~
    # DiffUploads <- list(
    #   "RawCounts" = Reactive,
    #   "GeneNames" = Reactive,
    #   "MetaData" = Reactive,
    #   "PrevResults" = Reactive list,
    #   "status" = Reactive 0, 1, 2
    # )
    # #~~~~~~~~~~~~~~~
    
    
    #     CONTRAST SELECTION
    #     (only shown it contrasts available)
    #_______________________________
    output$ContrastSelectionBox <- renderUI({
  
      if(is.null(Data$contrasts())){return()}
      
      selectInput(
        ns("NewContrast"), "Contrast",
        choices = names(Data$contrasts()),
      )
      
    })
    observeEvent(input$NewContrast, {
      
      if(is.null(Data$contrasts())){return()}
      
      new <- Data$contrasts()[[input$NewContrast]]
      
      Data$Results( 
        as.data.frame(
          results(Data$ddsc(), 
                  contrast = new
                  )))
      
      showNotification("Contrast Changed, all data will be updated.", type = "message", duration = 1.5)
      
    })
    
    
    
    
    #Module Servers
    ResultsSummaryServer("Results-Summary", Data, Uploads)
    ResultsPcaServer("Results-PCA", Data, Uploads)
    ResultsCorrelationServer("Results-Corr", Data, Uploads)
    ResultsGeneCountsServer("Results-Count", Data, Uploads)
    ResultsVolcanoServer("Results-Volcano", Data, Uploads)
    
  })
}