#_________________________________Shiny Server__________________________________________
server <- function(input, output) {
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_____________Reactive Values______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  run_question <- reactiveVal(NULL) #Determines when to run DESeq2
  kill_deg <- reactiveVal(NULL) #if DESeq2 has been run the DEG upload no longer available
  
  raw_counts <- reactiveVal(NULL)
  gene_names <- reactiveVal(NULL)
  filtered_counts <- reactiveVal(NULL)
  metaData <- reactiveVal(NULL)
  
  ddsc <- reactiveVal(NULL) #The Object DESeq2 
  results_ddsc <- reactiveVal(NULL) #Just a dataframe of results
  
  #Lists of regulated genes -0.5, 0.5 FoldChange cut off
  up <- reactiveVal(NULL)
  down <- reactiveVal(NULL)
  noR <- reactiveVal(NULL) #No regulation
  
  vst_counts_VST_OBJECT <- reactiveVal(NULL) #variance stabalized counts OBJECT
  vst_counts <- reactiveVal(NULL)#variance stabalized counts

  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_______________Observables________________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-----Observe UI Event------# ----> merged_gene_counts_uploaded_file
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  observeEvent(input$merged_gene_counts_uploaded_file,{
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
      
    })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #---Observe Value Event-----# ----> !is.null(raw_counts)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  observeEvent(raw_counts(),{
      
      if(is.null(raw_counts())){return()}
  
      func_return <- filterCounts(raw_counts())

      filtered_counts(data.frame(func_return))

      r <- raw_counts()
      f <- filtered_counts()
      diff <- nrow(r) - nrow(f)

      message <- div(p(span(diff, style="color: red;"), " rows were removed from the data set with row sums < 10"))
    
      showModal(modalDialog(message, easyClose = TRUE, footer=NULL))
    })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-----Observe UI Event------# ----> meta_data_conditions_uploaded_file
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  observeEvent(input$meta_data_conditions_uploaded_file,{
      
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
    })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #---Observe input Event-----# ----> run_DESeq2
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  #Performs DESeq with standard 1 condition design
  #sets results_ddsc() to the result matrix of the DESeq2 object
  #If the app is used with a raw counts matrix and DESeq2 is run ddsc() will be
  #a DESeqDataSet NOT a results matrix like when using DEG results 
  observeEvent(input$run_DESeq2,{
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
          div(p("Would you like to generate a DEG dataset with the given raw counts & meta data."), 
          p("This action will lock you out of uploading an already generated DEG table and will run DESeq2 with default parameters for your data."),
          p("If you have previously uploaded a DEG file it will be scrapped from the application, all tables will display the new DESeq2 data, visuals and graphs will also be updated accordingly."))
          ,
          footer=tagList(
            actionButton('cancel', 'Cancel'),
            actionButton('run', 'Run')
          )
        )
      )
    }) #Observe to see if everything is good
  observeEvent(input$cancel,{
    run_question(NULL)
    removeModal()
  }) #if cancel sumbit run DESeq2
  observeEvent(input$run,{
    run_question(TRUE)
    kill_deg(TRUE)
    removeModal()
  })#if confirm sumbit run DESeq2
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
        results_ddsc(data.frame(results(ddsc())))
        
        showModal(modalDialog("Done!!", easyClose=TRUE, footer=NULL))
        
        removeModal()
       
        Cvst <- vst(ddsc(), blind=TRUE)
        vst_counts(data.frame(assay(Cvst)))
        vst_counts_VST_OBJECT(Cvst)
        
      }else{
        return()
      }
    }else{
      return()
    }
    
  }) #Run
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #---Observe input Event-----# ----> DEG_analysis_data
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  #Takes in the DEG format data and puts it a format that will be processed by the 
  #rest of the app
  #vst_counts() is set here
  #ddsc() in the case of an uploaded DEG file is just a datafram == to results_ddsc()
  #results_ddsc() is the same as ddsc() if DEG is used rather than DESeq2
  observeEvent(input$DEG_analysis_data,{
 
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
    
     colNames <- colnames(data)
     
     if(ncol(data) < 8){
       showModal(modalDialog(
         tags$p(style = "color: red;","This file does not have the required data"), easyClose = TRUE, footer=NULL)
       )
       return()
     }
     
     format_check <- c("gene","gene_name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
     
     if(FALSE %in% (colNames[1:8] == format_check)){
       showModal(modalDialog(
         tags$p(style = "color: red;","This file may not have the correct format or data"), easyClose = TRUE, footer=NULL)
       )
       return()
     }
     
     if(ncol(data) < 9){
       showModal(modalDialog(
         tags$p(style = "color: red;","There doesn't seem to be any variance stabalized counts in your file, without them you may lose access to some features"), easyClose = TRUE, footer=NULL)
       )
       
       gene_names(data[,2:1])
       
       deg <- data[1:8]
       
       deg <- deg %>% select(-gene_name)
       
       deg <- deg %>% tibble::column_to_rownames('gene')
       
       results_ddsc(deg)
       
       return()
       
     }else{
       
       gene_names(data[,2:1])
       
       deg <- data[1:8]
       
       deg <- deg %>% select(-gene_name)
       
       deg <- deg %>% tibble::column_to_rownames('gene')
       
       results_ddsc(deg)
       
       sample_columns <- colNames[9:ncol(data)]
       stuff <- data.frame(data[,9:ncol(data)])
       stuff <- stuff %>% mutate(gene = gene_names()$gene)
       stuff <- stuff %>% tibble::column_to_rownames('gene')
       
       vst_counts(stuff)
       
       if(is.null(metaData())){
         showModal(modalDialog(
           tags$p(style = "color: red;","For a full analysis it is recomended to return to page one and submit a proper meta data file"), easyClose = TRUE, footer=NULL))
       }
       return()
    }
  })
  #If metaData and the DEG analysis info has been uploaded than ddsc() object DESeq2 dataset, 
  #and a DESeqTransform object will take the place of vst_counts_VST_OBJECT
  #A DESeqTransform is similar in that you can see the assays, in this case the counts they uploaded
  #(they have not been transformed like in the vst object as uploading a DEG you assumed responsibility for
  #normilization, due to this the size factors are not available for display), 
  #you may use plotPCA to plot principle components which is the primary use of the vst_counts
  # ALL other visuals volcano and counts plots do not require this specific object.
  observeEvent(list(metaData(), input$DEG_analysis_data), {
    
    if(is.null(input$DEG_analysis_data)){return()}
    if(is.null(metaData())){return()}
    if(is.null(vst_counts())){return()}

    vstcounts <- vst_counts()    
    
    datamatrix <- as.matrix(round(vstcounts))
    
    coldata <- metaData()
    
    cond <- factor(metaData()[,1])
    
    print(head(datamatrix))
    print(coldata)
    print(cond)
    
    dds <- DESeqDataSetFromMatrix(countData = datamatrix, colData = coldata, design = ~ cond)
    
    ddsc(dds)
    
    se <- SummarizedExperiment(assays = list(counts = vstcounts), colData = colData(dds))
    
    vsd <- DESeqTransform(se)
    
    vst_counts_VST_OBJECT(vsd)
    
  })
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #---Observe Value Event-----# ----> ddsc()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  #splits set by given cutoff values
  #Default:
  #down regulated  < -0.5
  #up regulated > 0.5
  #Anything else is classified as no regulation
  #up(), down(), noR() are all set and reset here upon inputs
  observeEvent(list(results_ddsc(),input$cutOffs),{
     
     if(is.null(results_ddsc())){
       return()
     }
     func_return <- splitByExpr(results_ddsc(), gene_names(), input$cutOffs )

     up(data.frame(func_return[1]))
     down(data.frame(func_return[2]))
     noR(data.frame(func_return[3]))
     
   })
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######____________Output $ Objects______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#______________________________Page 1____________________________#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #--------raw_counts_PreviewTable---------#
  #~~~~~~~~~~~~~~Data Table~~~~~~~~~~~~~~~~#
  #A Data Table preview of what is inside the uploaded counts table 
  #after formating correctly
    #Expected Format:
    #axis.names| Gene Name | Sample1 | Sample 2 
    #Gene_ID_1 |    x      |    x    |    x
    #Gene_ID_2 |    x      |    x    |    x
  output$raw_counts_PreviewTable <- renderDT({
      
      req(input$merged_gene_counts_uploaded_file)
      
      df <- raw_counts()
      
      if(!is.null(df)){
        
        df$gene <- row.names(df)
        
        #Some editing to display the gene names as well
        df <- left_join(df, gene_names())
        
        df <- df %>% tibble::column_to_rownames('gene')
        
        data.frame(df[, c((ncol(df)), 1:(ncol(df)-1))])
        
      }
      
    }, rownames = TRUE, options = list(pageLength=5))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #--------sample_conditions_Preview-------#
  #~~~~~~~~~~~~~~~~~Table~~~~~~~~~~~~~~~~~~#
  #Just prints the samples and their conditions 
  #May be helpful to ensure your file was read correctly
  output$sample_conditions_PreviewTable <- renderTable({
      if(!is.null(metaData())){
        
        t(metaData())
        
      }
    })
  
#______________________________Page 2____________________________#
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-----pg2_display_upload_DEG_data--------#
  #~~~~~~~~~~~~~Conditional UI~~~~~~~~~~~~~#
  # Decides whether to display the upload form for the DEG data based on whether 
  # the user has uploaded counts and metadata to perform DESeq2
  #@return either an upload form or nothing
  output$pg2_display_upload_DEG_data <- renderUI({
      
    if(is.null(kill_deg())){
        widget <- fileInput("DEG_analysis_data", "Differential Expression data .tsv/.csv")
        widget
      }else{
        return()
      }
    })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-----DESeq_Expression_Analysis_Tables---#
  #~~~~~~~~~~~~~~~~~HTML UI~~~~~~~~~~~~~~~~#
  #A ui div element featuring a series of interactable tables that show gene data
  output$DESeq_Expression_Analysis_Tables <- renderUI({
      
      return_ui <- NULL
      
      if(is.null(up())){
        return_ui <- p(span("No Data Yet",style="color: red;"))
        return_ui
      }else
      if(is.null(down())){
        return_ui <- p(span("No Data Yet",style="color: red;"))
        return_ui
      }else
      if(is.null(noR())){
        return_ui <- p(span("No Data Yet",style="color: red;"))
        return_ui
      }else{
      
        p <- as.numeric(input$pvalue)
        
        t1 <- up() %>% select(input$display_col, 'padj')
        t2 <- down() %>% select(input$display_col, 'padj')
        t3 <- noR() %>% select(input$display_col, 'padj')
        
        
        t1 <- t1 %>% filter(padj < as.numeric(input$pvalue))
        t2 <- t2 %>% filter(padj < as.numeric(input$pvalue))
        t3 <- t3 %>% filter(padj < as.numeric(input$pvalue))
        
        div(
          
          div(h2("Up Regulated Genes: ",  nrow(t1) )), 
          
          datatable(t1, rownames=FALSE, option=list(pageLength=7)),
          
          div(h2("Down Regulated Genes: ", nrow(t2) )),
          
          datatable(t2, rownames=FALSE, option=list(pageLength=7)),
          
          div(h2("Genes Between Selected Cutoffs: ", nrow(t3) )),
          
          datatable(t3, rownames=FALSE, option=list(pageLength=7))
        )
      
      }
      
      
      
    })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #--------DEG_Distribution_Histogram------#
  #~~~~~~~~~~~~~~~~~Histogram~~~~~~~~~~~~~~#
  #A simple histogram for a visual display of the distribution or regulation
  output$DEG_Distribution_Histogram <- renderPlot({
      if(is.null(results_ddsc())){return()}
      data <- results_ddsc() %>% filter(padj < as.numeric(input$pvalue))
      data <- data %>% select('log2FoldChange')
      
      ggplot(data, aes(x=log2FoldChange)) + 
        geom_histogram(bins=50, fill='darkslategrey') + 
        xlab("Log2 Fold Change") + 
        ylab("") + 
        theme_minimal() + 
        theme(panel.grid.major.x = element_blank(), 
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.line.x = element_blank(),
              axis.line.y = element_blank(),
              plot.margin = margin(75, 15, 0, -5, "pt"),
              axis.text.x = element_text(color = "black", size = 15),
              axis.title.x = element_text(size = 15),) 
      
    })
  #UI for histogram display and title 
  #(Prevents the blank plot box from showing when no data)
  output$Distribution_Histogram_ui <- renderUI({
      if(is.null(results_ddsc())){return()}
      
      div(
        
        HTML("<h4>Distribution of Fold Change</h4>")
        
        ,
        
        plotOutput('DEG_Distribution_Histogram')
      
        )
    })
  
  #______________________________Page 3____________________________#
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-------vst_counts_PreviewTable----------#
  #~~~~~~~~~~~~~~~~dataTable~~~~~~~~~~~~~~~#
  # A DataTable with the variance stabilized transformation of the counts
  output$vst_counts_PreviewTable <- renderUI({
    
    return_ui <- NULL
    
    if(is.null(vst_counts())){
      return_ui <- p(span("No Data Yet",style="color: red;"))
      return_ui
    }else{
    
      vst <- vst_counts() %>% round(digits = 4)
      
      vst <- vst %>% mutate(gene = rownames(vst_counts()))
      
      vst <- left_join(vst, gene_names()) %>% select(-gene)
      
      pvals <- results_ddsc() %>% mutate(gene = rownames(results_ddsc())) %>%
        left_join(gene_names()) %>% select(gene_name, padj)
      
      vst <- left_join(vst, pvals, relationship="many-to-many") 
      
      vst <- vst %>% filter(padj < as.numeric(input$pvaluePg3))
      
      vst <- vst[order(vst$padj),c((ncol(vst) - 1), ncol(vst), 1:(ncol(vst) - 2))]
      
      return_ui <- datatable(vst, rownames=FALSE, options = list(pageLength=15))
      
      return_ui
      
    }
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #----------Extra_vst_count_info----------# #Size factors only right now
  #~~~~~~~~~~~~~~~~UI output~~~~~~~~~~~~~~~#
  output$Extra_vst_count_info <- renderUI({
    
    if(is.null(vst_counts())){
      return_ui <- p(span("No Data Yet",style="color: red;"))
      return_ui
    }else{
     
      sF <- div()
      if(!is.null(kill_deg())){
        size <- as.matrix(sizeFactors(ddsc()) %>% round(digits = 3))
        print(size)
        sF <- div(h4("Size Factors"), renderTable({t(size)}))
      }else{print("doody butt")}
      
      div(sF)
    }
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #----------------PCA_PLOTS---------------# 
  #~~~~~~~~~~~~~~~~UI output~~~~~~~~~~~~~~~#
  output$PCA_PLOTS <- renderUI({
      #TODO
  })
  
  
  
  
  
}#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X

