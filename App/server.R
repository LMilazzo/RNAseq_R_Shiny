#_________________________________Shiny Server__________________________________________
server <- function(input, output) {
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_____________Reactive Values______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  raw_counts <- reactiveVal(NULL)
  gene_names <- reactiveVal(NULL)
  filtered_counts <- reactiveVal(NULL)
  metaData <- reactiveVal(NULL)
  ddsc <- reactiveVal(NULL)
  
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
      
      filtered_counts(data.frame(func_return))
    }
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #---Observe Value Event-----# ----> raw_counts & metaData -> ddsc
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  observeEvent(c(filtered_counts(), metaData()),
    {
      print("trigger")
      if(is.null(filtered_counts())){return()}
      if(is.null(metaData())){return()}
      print("will run")
      
      datamatrix <- as.matrix(filtered_counts())
      
      coldata <- metaData()
      
      cond <- factor(metaData()[,1])
      
      print(head(datamatrix))
      print(coldata)
      print(cond)
      
      d <- DESeqDataSetFromMatrix(countData = datamatrix, 
                                     colData = coldata, 
                                     design = ~ cond)
      
      showModal(modalDialog("Running DESeq...", footer=NULL))
      
      ddsc(DESeq(d))
      
      showModal(modalDialog("Done!!", easyClose=TRUE, footer=NULL))
      
      Sys.sleep(1)
      
      removeModal()
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
  output$sample_conditions_PreviewTable <- renderTable(
    {
      if(!is.null(metaData())){
        
        t(metaData())
        
      }
    }
  )
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #----------------nFiltered---------------#
  #~~~~~~~~~~~~~~~~Text Value~~~~~~~~~~~~~~#
  #The number of row filtered from filterCounts
  output$dataInfo <- renderUI(
    {
      
      if(is.null(raw_counts())){
        return(tags$p(style = "color: red;","Upload Data"))
      }
      
      r <- raw_counts()
      f <- filtered_counts()
      diff <- nrow(r) - nrow(f)
      #diff <- tags$p(style = "color: red; margin: 0; padding: 0;",diff)
      
      Top5RowSum <- head(r[order(-rowSums(r)),], 5)
      Top5RowSum$gene <- row.names(Top5RowSum)
      Top5RowSum <- left_join(Top5RowSum, gene_names())
      Top5RowSum <- Top5RowSum %>% tibble::column_to_rownames('gene_name')
      Top5RowSum <- datatable(Top5RowSum[,1:ncol(Top5RowSum)-1], options=list(pageLength=5, dom="t"))
      
      div(
        
        h1("Data Set Info"),
        
        h3("# Rows Filtered"),
        
        p(span(diff, style="color: red;"), " rows were removed from the data set with row sums < 10"),
        
        h3("Top 5 Genes by Row Sum"),
        
        Top5RowSum
        
        )
    }
  )
  
}#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X

