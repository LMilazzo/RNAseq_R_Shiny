DiffIntermediateServer <- function(id, Data){
  moduleServer(id, function(input, output, session){
    
    ns <- session$ns
  
    #     Experiment Setup
    #____________________________
    
    formula.text <- reactiveVal("")
    numOfInteractions <- reactiveVal(0)
    
    #Reload the formula
    observe({
      invalidateLater(1000)
      
      mainterms <- paste(input$MainTerms, collapse = " + ")
      
      if((numOfInteractions() > 0)){
        interactionterms <- sapply(1:numOfInteractions(), function(i) {
          # Check the value of each checkbox dynamically based on its ID
          
          paste0( input[[paste0("Interaction", i, "Left")]], ":",  
                  input[[paste0("Interaction", i, "Right")]] )
          
        })
        interactionterms <- paste(interactionterms, collapse = " + ")
        if(length(input$MainTerms) > 0){
          interactionterms <- paste0(" + ", interactionterms)
        }
      }else{interactionterms <- ""}
      
      
      formula.text(paste0(mainterms, interactionterms))
      
    })
    
    #Text output of the formula
    output$FormulaText <- renderText({
      
      paste0("~ ", formula.text())
      
    })
    
    #Main term selection
    output$MainTermSelection <- renderUI({
      
      options <- colnames(Data$MetaData())
      
      checkboxGroupInput(
        ns("MainTerms"), "", choices = options, selected = options[1]
      )
      
    })
    
    #Remove Interaction
    observeEvent(input$RemoveInteraction,{
      if(numOfInteractions() > 0){
        numOfInteractions(numOfInteractions()-1)
      }
    })
    #Add Interaction
    observeEvent(input$AddInteraction,{numOfInteractions(numOfInteractions()+1)})
    
    #Interaction Term selection
    output$InteractionTermSelection <- renderUI({
      
      req(numOfInteractions() > 0)
      
      options <- colnames(Data$MetaData())
      
      interactionTags <- lapply(1:numOfInteractions(), function(i){
        
        fluidRow(
          column( 
            width = 6,
            
            div(
              selectInput(
                ns(paste0("Interaction", i, "Left")),
                "Variable 1: ",
                choices = options,
                selected = input[[paste0("Interaction", i, "Left")]]
              ),
              style = "display: flex; align-items: center; justify-content: center; padding: 20px;"
            )
            
          ),
          column(
            width = 6,
            div(
              selectInput(
                ns(paste0("Interaction", i, "Right")),
                "Variable 2:",
                choices = options,
                selected = input[[paste0("Interaction", i, "Right")]]
              ),
              style = " display: flex; align-items: center; justify-content: center; padding: 20px"
            )
          )
        )
        
      })
      
      tagList(interactionTags)
      
    })
    
    
    #     RAW COUNTS
    #____________________________
    
    # Data table for raw counts
    output$RawDataTable <- renderDT({
      
      if(is.null(Data$RawCounts())){return()}
      
      Data$RawCount()
      
    })
    
    # Data table for searched gene summary
    output$SearchedGene <- renderDT({
      
      index <- tolower(Data$GeneNames()) == tolower(input$GenePreviewSummary)
      
      gene <- Data$RawCounts()[index,]
      
      if(is.null(gene) || nrow(gene) != 1){ return() }
      print(gene)
      
      datatable(
        data.frame(
          Gene = input$GenePreviewSummary,
          Sum = sum(gene),
          Highest = max(gene),
          Mean = sum(gene) / ncol(gene),
          Lowest = min(gene)
        ) %>% t(),
        colnames = NULL,
        options = list(dom = "t")
      )
      
    })
    
    #     Meta Data
    #____________________________
    
    # Data table with a bulk dump of the meta data
    output$MetaDataDumpPreview <- renderDT({
      
      datatable(
        Data$MetaData(),
        options = list(dom = "ft")
      )
    
    })
    
    
    #     Run DESeq2
    #____________________________
    
    #RESULTS DATA
    RCV.CONTRASTS <- reactiveVal(NULL)
    RCV.DDSC <- reactiveVal(NULL)
    RCV.RESULTS.DDSC <- reactiveVal(NULL)
    RCV.VST_OBJ <- reactiveVal(NULL)
    RCV.VST_COUNTS <- reactiveVal(NULL)
    RCV.NORM_COUNTS <- reactiveVal(NULL)
    status <- reactiveVal(0)
    
    #Create formula modal
    observeEvent(input$RunDeseq,{
      
      showModal({
        modalDialog(
          #Stuff----
          div(
            h3("Design Formula"),
            textOutput(ns("FormulaText")),
            style = "text-align: center; font-size: 16px; border: 1px solid white;
            padding-bottom: 5px; border-radius: 10px;"
          ),
          
          fluidRow(
            column(
              width = 4,
              div(
                h4("Main Terms"),
                uiOutput(ns("MainTermSelection")),
                style = " text-align: center; border:
              margin: 5px; padding 5px; border: 1px solid #4CAF50; margin: 5px;
                border-radius: 10px"
              )
            ),
            column(
              width = 8,
              div(
                h4("Interaction Terms"),
                actionButton(ns("RemoveInteraction"), "Remove", 
                style="background-color: #FF6F61; font-size: 12px; 
                color: black; border-color: black; background-image: none; margin-right: 7%"),
                
                actionButton(ns("AddInteraction"), " Add ",
                style="background-color: #4CAF50; font-size: 12px;
                color: black; border-color: black; background-image: none; margin-right: 5%"),
                
                uiOutput(ns("InteractionTermSelection")),
                
                style = "text-align: center; margin: 5px; padding 5px;"
              )
            )
          ),
          #----
          size = "l",
          footer = tagList(
            div(
              
              actionButton(ns("ConfirmDeseq"), "Confirm", 
              style="background-color: #4CAF50; color: black; 
              border-color: black; background-image: none; "),
              
              tags$div(
                class = "tip-mark", ' ? ',
                tags$div(class = "tool-tip", 
                HTML("<p>The design formula in DESeq2 specifies the experimental conditions and factors to be tested for differential expression.
                  Include variables from your metadata tthat should be includes in the analysis. For more complex designs you may add interaction effect
                  between variables.</p><p>In the case that your run fails upon submission it is suggested that you check the DESeq2 library documentation for the displayed error.</p>")),
                style = "margin-left: 25px;"
              ),
              
              style = "justify-content: bottom; align-items: bottom; display: inline-flex;"
            )
          ),easyClose = TRUE
        )
      })
      
    })
    
    #Run logic
    observeEvent(input$ConfirmDeseq, {
      removeModal()
      
      #md <- Data$MetaData()
      #data <- Data$RawCounts()
      formula <- paste0("~", formula.text()) %>% as.formula()
      
      
      # results <- list(
      #   
      #   "contrast_list" = contrast_list,
      #   "ddsc" = DESeqDataSet,
      #   "results_ddsc" = as.data.frame(results(DESeqDataSet, contrast = default[[1]])),
      #   "vst_obj" = vst(DESeqDataSet, blind = TRUE, nsub = 50),
      #   "vst_counts" = as.data.frame(assay(vst(DESeqDataSet, blind = TRUE, nsub = 50))),
      #   "normalized_counts" = as.data.frame(counts(DESeqDataSet, normalized = TRUE))
      #   
      # )
      Diff_Results <- runDeseq(Data$MetaData(), Data$RawCounts(), formula)
      
      RCV.CONTRASTS(Diff_Results$contrast_list)
      RCV.DDSC(Diff_Results$ddsc)
      RCV.RESULTS.DDSC(Diff_Results$results_ddsc) 
      RCV.VST_OBJ(Diff_Results$vst_obj)
      RCV.VST_COUNTS(Diff_Results$vst_counts) 
      RCV.NORM_COUNTS(Diff_Results$normalized_counts) 
      status(1)
      
    })
  
    #     RETURN
    #____________________
    
    observe({
      invalidateLater(500)
      if(Data$status() == 2){
        RCV.CONTRASTS(NULL)
        RCV.DDSC(NULL)
        RCV.RESULTS.DDSC(Data$PrevResults()$Results)
        RCV.VST_OBJ(Data$PrevResults()$vst_obj)
        RCV.VST_COUNTS(Data$PrevResults()$vst_counts)
        RCV.NORM_COUNTS(Data$PrevResults()$NormalizedCounts)
        status(2)
      }
    })
    
    return(
      list(
        "contrasts" = RCV.CONTRASTS,
        "ddsc" = RCV.DDSC,
        "Results" = RCV.RESULTS.DDSC,
        "vstObject" = RCV.VST_OBJ,
        "vstCounts" = RCV.VST_COUNTS,
        "NormCounts" = RCV.NORM_COUNTS,
        "status" = status
      )
    )
   
    
  })
}