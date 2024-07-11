#_________________________________Shiny Server__________________________________________
server <- function(input, output) {
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_____________Reactive Values______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  run_question <- reactiveVal(NULL)
  raw_counts <- reactiveVal(NULL)
  gene_names <- reactiveVal(NULL)
  filtered_counts <- reactiveVal(NULL)
  metaData <- reactiveVal(NULL)
  ddsc <- reactiveVal(NULL)
  up <- reactiveVal(NULL)
  down <- reactiveVal(NULL)
  noR <- reactiveVal(NULL)
  vst_counts <- reactiveVal(NULL)
  
  
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
        showModal(modalDialog(
          tags$p(style = "color: red;","Error: Raw Counts upload must be a valid file type (.tsv, .csv)"), easyClose = TRUE, footer=NULL))
        return()
      }
      
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
      if( !grepl('.csv',metaDataFilePath , fixed=TRUE) && !grepl('.tsv',metaDataFilePath , fixed=TRUE)){
        showModal(modalDialog(
          tags$p(style = "color: red;","Error: Meta Data Table upload expected a .csv or .tsv file"), easyClose = TRUE, footer=NULL))
        return()
      }
      
      if(grepl('.tsv',metaDataFilePath , fixed=TRUE)){
        md <- read.csv(metaDataFilePath, sep="\t")
      }else{
        md <- read.csv(metaDataFilePath)
      }
      
      #Error PG1.2
      if( ncol(md) != 2 ){
        showModal(modalDialog(
          tags$p(style = "color: red;","Error: Incorrect file format"), easyClose = TRUE, footer=NULL))
        return()
      }
      
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

      r <- raw_counts()
      f <- filtered_counts()
      diff <- nrow(r) - nrow(f)

      message <- div(p(span(diff, style="color: red;"), " rows were removed from the data set with row sums < 10"))
    
      showModal(modalDialog(message, easyClose = TRUE, footer=NULL))
    }
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #---Observe input Event-----# ----> raw_counts & metaData 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#  Performs DESeq 
  observeEvent(input$run_DESeq2,
    {
      if(is.null(filtered_counts()) || nrow(filtered_counts()) <= 0 ){
        message <- p(span("Merged counts table required", style="color: red;"))
        showModal(modalDialog(message, easyClose=TRUE, footer=NULL))
        return()
      }
      if(is.null(metaData())){
        message <- p(span("Meta Data table required", style="color: red;" ))
        showModal(modalDialog(message, easyClose=TRUE, footer=NULL))
        return()
      }
      
      Sys.sleep(.5)
      showModal(modalDialog(
          div(p("Would you like to generate a DEG dataset with the given raw counts & meta data"), 
          p("This action will lock you out of uploading an already generated DEG table and will run DESeq2 with default parameters for your data"))
          ,
          footer=tagList(
            actionButton('cancel', 'Cancel'),
            actionButton('run', 'Run')
          )
        )
      )
    }
  )
  observeEvent(input$cancel,{
    run_question(NULL)
    removeModal()
  })
  observeEvent(input$run,{
    run_question(TRUE)
    removeModal()
  })
  observe({
    decision <- run_question()
    
    if(!is.null(decision)){
      if(decision){
        removeModal()
        
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
        
        removeModal()
       
        Cvst <- vst(ddsc(), blind=TRUE)
        vst_counts(assay(Cvst))
        
      }else{
        return()
      }
    }else{
      return()
    }
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #---Observe input Event-----# ----> DEG_analysis_data
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#  Takes in the DEG format data
  observeEvent(input$DEG_analysis_data,
   {
 
     filePath <- input$DEG_analysis_data$datapath
     
     if( !grepl('.tsv', filePath, fixed=TRUE) && !grepl('.csv', filePath, fixed=TRUE)){
      showModal(modalDialog(
        tags$p(style = "color: red;","Error: Meta Data Table upload expected a .csv or .tsv file"), easyClose = TRUE, footer=NULL)
        )
      return()
     }
      
     if(grepl('.tsv', filePath, fixed=TRUE)){
       data <- read.csv(filePath, sep="\t")
     }else{
       data <- read.csv(filePath)
     }
     
     i_swear_you_better_be_the_right_format <- c("gene","gene_name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")

     colNames <- colnames(data)
     
     if(ncol(data) < 8){
       showModal(modalDialog(
         tags$p(style = "color: red;","This file does not have the required data"), easyClose = TRUE, footer=NULL)
       )
       return()
     }
     if(FALSE %in% (colNames[1:8] == i_swear_you_better_be_the_right_format)){
       showModal(modalDialog(
         tags$p(style = "color: red;","This file may not have the correct format or data"), easyClose = TRUE, footer=NULL)
       )
       return()
     }
     
     if(ncol(data) < 9){
       showModal(modalDialog(
         tags$p(style = "color: red;","The vst counts are not present in this file"), easyClose = TRUE, footer=NULL)
       )
       return()
     }
     
     sample_columns <- colNames[9:ncol(data)]
     
     if(length(sample_columns) < 2){
       showModal(modalDialog(
         tags$p(style = "color: red;","There is only one sample column present in your file"), easyClose = TRUE, footer=NULL)
       )
       return()
     }
     
     gene_names(data[,2:1])
     
     deg <- data[1:8]
     
     deg <- deg %>% select(-gene_name)
     
     deg <- deg %>% tibble::column_to_rownames('gene')
     
     ddsc(deg)
     
     stuff <- data.frame(data[,9:ncol(data)])
     stuff <- stuff %>% mutate(gene = gene_names()$gene)
     stuff <- stuff %>% tibble::column_to_rownames('gene')
     
     vst_counts(stuff)
    }
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #---Observe Value Event-----# ----> ddsc()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#  splits set by given cutoff values
  observeEvent(list(ddsc(),input$cutOffs),
   {
     
     if(is.null(ddsc())){
       return()
     }
     if(is.data.frame(ddsc())){
       func_return <- splitByExpr(data.frame(ddsc()), gene_names(), input$cutOffs )

       up(data.frame(func_return[1]))
       down(data.frame(func_return[2]))
       noR(data.frame(func_return[3]))
     }else{
       
       func_return <- splitByExpr( data.frame(results(ddsc())), gene_names(), input$cutOffs )
       
       up(data.frame(func_return[1]))
       down(data.frame(func_return[2]))
       noR(data.frame(func_return[3]))
     }
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
        
        df <- df %>% tibble::column_to_rownames('gene')
        
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
  #-----DESeq_Expression_Analysis_Tables---#
  #~~~~~~~~~~~~~~~~~HTML UI~~~~~~~~~~~~~~~~#
  output$DESeq_Expression_Analysis_Tables <- renderUI(
    {
      
      return_ui <- NULL
      
      if(is.null(up())){
        return_ui <- h1(span("No Data Yet",style="color: red;"))
        return_ui
      }else
      if(is.null(down())){
        return_ui <- h1(span("No Data Yet",style="color: red;"))
        return_ui
      }else
      if(is.null(noR())){
        return_ui <- h1(span("No Data Yet",style="color: red;"))
        return_ui
      }else{

        div(
          
          div(h1("Up Regulated Genes: ",  nrow(up() %>% filter(padj < as.numeric(input$pvalue) ) ) )), 
          
          datatable(up() %>% filter(padj < as.numeric(input$pvalue) ), rownames=FALSE, option=list(pageLength=7)),
          
          div(h1("Down Regulated Genes: ", nrow(down() %>% filter(padj < as.numeric(input$pvalue) ) ) )),
          
          datatable(down()%>% filter(padj < as.numeric(input$pvalue) ), rownames=FALSE, option=list(pageLength=7)),
          
          div(h1("Genes with low diffrential expression: ", nrow(noR() %>% filter(padj < as.numeric(input$pvalue) ) ) )),
          
          datatable(noR()  %>% filter(padj < as.numeric(input$pvalue) ), rownames=FALSE, option=list(pageLength=7))
        )
      
      }
      
      
      
    }
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-----pg2_display_upload_DEG_data--------#
  #~~~~~~~~~~~~~Conditional Ui~~~~~~~~~~~~~#
  # Decides whether to display the upload form for the DEG data based on whether 
  # the user has uploaded counts and metadata to perform DESeq2
  #@return either an upload form or nothing
  
  output$pg2_display_upload_DEG_data <- renderUI(
    {
      if(is.null(ddsc())){
        widget <- fileInput("DEG_analysis_data", "Differential Expression data .tsv/.csv")
        widget
      }else{
        return()
      }
    }
  )
}#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X

