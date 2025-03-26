
server <- function(input, output, session){
  
  #Reload Application Logic
  observeEvent(input$ReloadApp, {runjs("location.reload();")})

  
  #       APPLICATION STAGE 1 RAT DEG
  #________________________________________________
  
  
  #Differential Gene Expression Uploads Data
  # DiffUploads <- list(
  #   "RawCounts" = Reactive,
  #   "GeneNames" = Reactive,
  #   "MetaData" = Reactive,
  #   "PrevResults" = Reactive list,
  #   "status" = Reactive 0, 1, 2
  # )
  #~~~~~~~~~~~~~~~
  DiffUploads <- DiffUploadsServer("Diff-Uploads")
  
  #Differential Gene Expression Preview Data ( returns the calculated results)
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
  DiffResults <- DiffIntermediateServer("Diff-Intermediate", 
    DiffUploads #Data from the uploaded files to be run through DESeq2
  )
  
  #The UI Elements that make up a bulk of the analysis and app
  DiffResultsServer("Diff-Results", 
    DiffResults, #Results of DESeq2 or the uploaded files
    DiffUploads #Data from the uploaded files used for things like gene names and if the raw data is needed
  )
  
  #Validate Diff Display For Current Status
  observe({
    invalidateLater(500)
    
    #Case when there is any data at all
    if(DiffUploads$status() != 0 ){
      shinyjs::hide("Diff-Uploads")
    }
    
    #Case when data has been uploaded but deseq not run
    if(DiffUploads$status() == 1 && DiffResults$status() == 0){
      shinyjs::show("Diff-Intermediate")
    }
    
    #Case when there is results data 
    if(DiffUploads$status() > 0 && DiffResults$status() > 0){
      shinyjs::hide("Diff-Intermediate")
      shinyjs::show("Diff-Results")
    }
    
    
  })
  
  
  #       APPLICATION STAGE 2 RAT PATHS
  #________________________________________________
  
  PathUploads <- PathUploadsServer("Path-Uploads")
  
}