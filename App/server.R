#_________________________________Shiny Server__________________________________________
server <- function(input, output) {
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_____________Reactive Values______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  raw_counts <- reactiveVal(NULL)
  gene_names <- reactiveVal(NULL)
  filtered_counts <- reactiveVal(NULL)
  nFiltered <- reactiveVal(0)
  metaData <- reactiveVal(NULL)
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_______________Observables________________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-----Observe UI Event------# ----> merged_gene_counts_uploaded_file
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  observeEvent(input$merged_gene_counts_uploaded_file,
    {
      #Upon submission of UI object [ merged_gene_counts_uploaded_file ]      
      
      #Goal:
      #             |#----->> gene_names
      #   input ----|
      #             |#----->> raw_counts

      #Expect as return c(raw_counts, gene_names)           
      func_return <- Set_Clean_Counts(input$merged_gene_counts_uploaded_file$datapath)
      
      #Error PG1.0
      if(length(func_return) == 1){
        output$errorMessagesPG1.0 <- renderUI({
          tags$p(style = "color: red;","Error: Raw Counts upload must be .tsv file")
        })
        return()
      }else{output$errorMessagesPG1.0 <- NULL}
      
      raw_counts(data.frame(func_return[1]))
      
      gene_names(data.frame(func_return[2]))
      
    }             
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-----Observe UI Event------# ----> meta_data_conditions_uploaded_file
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  observeEvent(input$meta_data_conditions_uploaded_file,
    {
      
      req(input$meta_data_conditions_uploaded_file)
      
      metaDataFilePath  <- input$meta_data_conditions_uploaded_file$datapath
      
      #Error PG1.1
      if( !grepl('.csv',metaDataFilePath , fixed=TRUE) ){
        output$errorMessagesPG1.1 <- renderUI({
          tags$p(style = "color: red;","Error: Meta Data Table upload expected a .csv file")
        })
        return()
      }else{output$errorMessagesPG1.1 <- NULL}
      
      md <- read.csv(metaDataFilePath)
      
      #Error PG1.2
      if( ncol(md) != 2 ){
        output$errorMessagesPG1.2 <- renderUI({
          tags$p(style = "color: red;","Error: Incorrect .csv format")
        })
        return()
      }else(output$errorMessagesPG1.2 <- NULL)
      
      cond <- factor(md[,2])
      
      coldata <- data.frame(cond)
      
      rownames(coldata) <- md[,1]
      
      metaData(coldata)
    }
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #---Observe Value Event-----# ----> !is.null(raw_counts)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  observeEvent(raw_counts(),
    {
      
      if(is.null(raw_counts())){return()}
      
      func_return <- filterCounts(raw_counts())
      
      filtered_counts(data.frame(func_return[1]))
      
      nFiltered(func_return[2])
      
    }
  )
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######____________Output $ Objects______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #------------errorMessages(pg#)----------#
  #~~~~~~~~~~~~~~~Text Value~~~~~~~~~~~~~~~#
  #For displaying Error messegaes
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #--------raw_counts_PreviewTable---------#
  #~~~~~~~~~~~~~~Data Table~~~~~~~~~~~~~~~~#
  #A Data Table preview of what is inside the uploaded counts table 
  #after formating correctly
    #Expected Format:
    #axis.names| Gene Name | Sample1 | Sample 2 
    #Gene_ID_1 |    x      |    x    |    x
    #Gene_ID_2 |    x      |    x    |    x
  output$raw_counts_PreviewTable <- renderDT(
    {
      
      req(input$merged_gene_counts_uploaded_file)
      
      df <- raw_counts()
      
      if(!is.null(df)){
        
        df$gene <- row.names(df)
        
        #Some editing to display the gene names as well
        df <- left_join(df, gene_names())
        
        df <-df %>% tibble::column_to_rownames('gene')
        
        data.frame(df[, c((ncol(df)), 1:(ncol(df)-1))])
        
      }
      
    }, rownames = TRUE, options = list(pageLength=5)
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #--------sample_conditions_Preview-------#
  #~~~~~~~~~~~~~~~~~Table~~~~~~~~~~~~~~~~~~#
  #A Data Table preview of what is inside the uploaded counts table 
  #after formating correctly
  #Expected Format:
  #axis.names| Gene Name | Sample1 | Sample 2 
  #Gene_ID_1 |    x      |    x    |    x
  #Gene_ID_2 |    x      |    x    |    x
  output$sample_conditions_PreviewTable <- renderTable(
    {
      if(!is.null(metaData())){
        
        t(metaData())
        
      }
    }
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-------filtered_counts_PreviewTable-----#
  #~~~~~~~~~~~~~~Data Table~~~~~~~~~~~~~~~~#
  #A Data Table preview of what is inside the uploaded counts table 
  #after filtering correctly for low count genes
    #Expected Format:
    #axis.names| Gene Name | Sample1 | Sample 2 
    #Gene_ID_1 |    x      |    x    |    x
    #Gene_ID_2 |    x      |    x    |    x
  output$filtered_counts_PreviewTable <- renderDT(
    {
      req(!is.null(filtered_counts))
      
      #                       subject   setting
      df <- makePreviewTable(filtered_counts(),"Full")
      
      if(is.null(df)){
        #Nothing is displayed if the preview table function comes back empty
      }else{
        df$gene <- row.names(df)
        
        #Some editing to display the gene names as well
        df <- left_join(df, gene_names())
        
        df <-df %>% tibble::column_to_rownames('gene')
        
        data.frame(df[, c((ncol(df)), 1:(ncol(df)-1))])
      }
      
    }
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #----------------nFiltered---------------#
  #~~~~~~~~~~~~~~~~Text Value~~~~~~~~~~~~~~#
  #The number of row filtered from filterCounts
  output$nFiltered <- renderText(
    {
      message <- paste("There were ", nFiltered(), "rows/genes omited due to row sums < 10")
      message
    }
  )
  
}#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X

