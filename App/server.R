#_________________________________Shiny Server__________________________________________
server <- function(input, output) {
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_____________Reactive Values______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  raw_counts <- reactiveVal(NULL)
  gene_names <- reactiveVal(NULL)
  filtered_counts <- reactiveVal(NULL)
  nFiltered <- reactiveVal(0)
  
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
      
      if(length(func_return) == 1){
        output$errorMessagesPG1.0 <- renderText({
          func_return
        })
      }else{
        
        output$errorMessages <- NULL
        
        raw_counts(data.frame(func_return[1]))
        
        gene_names(data.frame(func_return[2]))
        
      }
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
    req(input$raw_counts_matrix_Preview)
    
    #                       subject         setting
    df <- makePreviewTable(raw_counts(),input$raw_counts_matrix_Preview)
    
    if(is.null(df)){
      #Nothing is displayed if the preview table function comes back empty
    }else{
      df$gene <- row.names(df)
      
      #Some editing to display the gene names as well
      df <- left_join(df, gene_names())
      
      df <-df %>% tibble::column_to_rownames('gene')
      
      data.frame(df[, c((ncol(df)), 1:(ncol(df)-1))])
      
      }
    
    }, rownames = TRUE, caption= "Raw Counts"
  )

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #------geneID_geneName_PreviewTable------#
  #~~~~~~~~~~~~~~Data Table~~~~~~~~~~~~~~~~#
  #A Data Table preview of all gene ids and their gene names
    #Expected Format:
    #axis.names|    gene    | gene_name
    #   1      | Gene_ID_1  |   Actin    
    #   2      | Gene_ID_2  |   COOR4_7    
  output$geneID_geneName_PreviewTable <- renderDT(
    {
      req(input$merged_gene_counts_uploaded_file)
      req(input$geneID_geneName_Preview)
      
      #                 subject
      makePreviewTable(gene_names(),input$geneID_geneName_Preview)
      
    }, rownames = TRUE, caption= "Gene IDs & Gene name"
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

