DiffUploadsServer <- function(id){
  moduleServer(id, function(input, output, session){
    
    ns <- session$ns
    
    status <- reactiveVal(0)
    
    #       NEW DATA
    #____________________________
    
    RCV.Counts.Raw <- reactiveVal(NULL)
    RCV.GeneNames <- reactiveVal(NULL)
    RCV.MetaData <- reactiveVal(NULL)
    
    #User Uploads new experiment data modal
    observeEvent(input$UploadRawDiffData, {
      showModal(
        modalDialog(
          div( 
            h5("Counts Matrix"),
            p("Accepted File Types:  .csv   .tsv "),
            p("Format: ( gene_id, gene_name, .samples... ) "),
            p("The gene id column should contain a unique identifier for every row"),
            
            fileInput(ns("Diff_Raw_Counts"), "Counts Matrix"),
            
            h5("Sample Conditions Table"),
            p("Accepted File Types:  .csv   .tsv "),
            p("Format: (sample, conditions...) "),
            p("Sample column should include sample names that are found in your experiment data"),
            
            fileInput(ns("Diff_Meta_Data"), "Sample Conditions")
          ),
          footer = div(
            div(
              checkboxInput(ns("CustomFiltering"), "Custom Filtering", value = FALSE),
              style = "font-size: 16px; padding: 2px; margin: 5px;"
            ),
            actionButton(ns("LoadRawDiffData"), "Finish"),
            style = "display: inline-flex; justify-content: center; align-items: center;
            text-align: center; position: relative;"
          ), easyClose = TRUE
        )
      )
    })
    
    #User Confirms new Experiment data
    observeEvent(input$LoadRawDiffData, {
      removeModal()
      
      #         COUNTS DATA
      #______________________________
      
      #Asserting File Type
      RCV.Counts.Raw(AssertFile(input$Diff_Raw_Counts$datapath))
      if(!is.data.frame(RCV.Counts.Raw()) || !all(c("gene_name", "gene_id") %in% colnames(RCV.Counts.Raw()))){
        showErrorModal("(Raw Counts) Invalid File Type Error or Missing Data") 
        return()
      }
      
      #Collect Row Sums
      sums <- rowSums(RCV.Counts.Raw() %>% select(-gene_name, -gene_id))

      #Custom Filtering Modal
      if(input$CustomFiltering){
        showModal(
          modalDialog(
            h4("Would you like to edit the filtering process?"),
            
            numericInput(ns("FilterTarget"), "Minimum Row Sum: ", value = 10, min = 0),
            
            div(
              renderText(paste0("Rows Removed: ", sum(sums <= input$FilterTarget))),
              style = "color: red;"
            ),
            
            footer = actionButton(ns("FilterRawData"), "Continue")
          )
        )
      } 
      #Custom Filtering Skipped
      else{
        removeModal()
        x <- readCountsDF(RCV.Counts.Raw(), 0)
        if(!is.list(x)){
          return()
        }
        RCV.Counts.Raw(x$counts)
        RCV.GeneNames(x$gene_names)
        status(1)
      }  
      
      #          META DATA
      #_______________________________
      RCV.MetaData(AssertFile(input$Diff_Meta_Data$datapath))
      if(!is.data.frame(RCV.MetaData())){
        showErrorModal("(Sample Conditions) Invalid File Type Error") 
        return()
      }
      
      
      RCV.MetaData(readMetaDataDF(RCV.MetaData()))
      
    })
    
    #Custom Filtering Logic
    observeEvent(input$FilterRawData, {
      removeModal()
      x <- readCountsDF(RCV.Counts.Raw(), input$FilterTarget)
      if(!is.list(x)){
        return()
      }
      RCV.Counts.Raw(x$counts)
      RCV.GeneNames(x$gene_names)
      status(1)
    })
    
    #     OLD DATA
    #____________________________
    
    RCV.Results <- reactiveVal(NULL)
    
    #User Uploads old data
    observeEvent(input$UploadPrevDiffData ,{
      showModal(
        modalDialog(
          div( 
            h5("DEG Experiment"),
            p("Accepted File Types:  .csv   .tsv "),
            p("Format: ( gene_id, gene_name, log2FoldChange, padj, .samples... ) "),
            p("- The gene id column should contain a unique identifier for every row"),
            p("- Sample counts may be included and are recommended in reviewing an experiment for all features"),
            p("- Samples should each have their own column with the sample name proceeding a '.' (standard dot/period)"),
            
            fileInput(ns("Diff_Prev_Results"), "DEG Results"),
            
            h5("Sample Conditions Table"),
            p("Accepted File Types:  .csv   .tsv "),
            p("Format: (sample, conditions...) "),
            p("Sample column should include sample names that are found in your experiment data"),
            
            fileInput(ns("Diff_Meta_Data"), "Sample Conditions")
            
          ),footer = actionButton(ns("LoadPrevDiffData"), "Finish"), easyClose = TRUE
        )
      )
    })
    
    #User Confirms prev Experiment data
    observeEvent(input$LoadPrevDiffData, {
      removeModal()
      
      #     RESULTS DATA
      #______________________________
      RCV.Results(AssertFile(input$Diff_Prev_Results$datapath))
      if(!is.data.frame(RCV.Results())){
        showErrorModal("(Results) Invalid File Type Error") 
        return()
      }
      
      #          META DATA
      #_______________________________
      RCV.MetaData(AssertFile(input$Diff_Meta_Data$datapath))
      if(!is.data.frame(RCV.MetaData())){
        showErrorModal("(Sample Conditions) Invalid File Type Error") 
        return()
      }
      RCV.MetaData(readMetaDataDF(RCV.MetaData()))
    
      x <- readDegDF(RCV.Results(), RCV.MetaData())
      if(!is.list(x)){
        return()
      }
      RCV.GeneNames(x$GeneNames)
      RCV.Results(x)
      status(2)
    })
    
    
    #     RETURN
    #____________________________
    
    return(
      
      list(
        "RawCounts" = RCV.Counts.Raw,
        "GeneNames" = RCV.GeneNames,
        "MetaData" = RCV.MetaData,
        "PrevResults" = RCV.Results,
        "status" = status
      )
      
    )
    
  })
}