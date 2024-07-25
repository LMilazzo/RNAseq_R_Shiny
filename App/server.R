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
  normalized_counts <- reactiveVal(NULL) #normalized counts from the DESeq
  
  #Lists of regulated genes -0.5, 0.5 FoldChange cut off
  up <- reactiveVal(NULL)
  down <- reactiveVal(NULL)
  noR <- reactiveVal(NULL) #No regulation
  
  vst_Obj <- reactiveVal(NULL) #variance stabalized counts OBJECT
  vst_counts <- reactiveVal(NULL)#variance stabalized counts
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_______________Observables________________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-----Observe UI Event------# ----> merged_gene_counts_uploaded_file
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Process:
  # Reads the uploaded file with Set_Clean_Counts() function returning a list with 
  # the raw counts matrix and a dataframe of gene names and ids
  # Takes out gene names and ids and sets gene_names() rownames = gene_id, col1 = gene names
  # Sets raw_counts() colnames = samples, rownames = gene ids
  # Sets filtered_counts() by calling filterCounts() on the raw_counts() to filter any rows with
  # sum < 10
  # 
  # Error:
  # Will return an error if neither a .csv or .tsv were uploaded
  # Will return an error if either gene names or gene ids columns were not present
  # Will return an error if there is a duplicate value in the gene_id column
  # 
  # Notes: 
  # Gene ids CAN be arbitrary, but MUST be unique
  # Will display the number of rows filtered from the data
  observeEvent(input$merged_gene_counts_uploaded_file,{
      
      #Expect as return c(raw_counts, gene_names)           
      func_return <- Set_Clean_Counts(input$merged_gene_counts_uploaded_file$datapath)
      
      if(length(func_return) == 1){
        showModal(modalDialog(
          div(
            span(
              div(
                p("Possible Errors: "),
                p(" - Raw Counts upload was not valid file type (.tsv, .csv)"),
                p(" - Raw Counts upload did not include required 'gene_name' and 'gene_id' columns"),
                p(" - There were duplicate identifiers in the gene_id column (ids be abitrary but must be unique)")
              ),
              style="color: red;")
          ), 
        easyClose = TRUE, footer=NULL))
        return()
      }
      
      raw_counts(data.frame(func_return[1]))##################################################Reactive Assignment
      gene_names(data.frame(func_return[2]))##################################################Reactive Assignment
      
      filtered_counts( filterCounts(raw_counts()))##################################################Reactive Assignment
      
      diff <- nrow( raw_counts() ) - nrow( filtered_counts() )
      
      showModal(modalDialog(div(p(span(diff, style="color: red;"), 
          " rows were removed from the data set with row sums < 10")), 
        easyClose = TRUE, footer=NULL))
    })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-----Observe UI Event------# ----> meta_data_conditions_uploaded_file
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Process:
  # Will set the path to the datapath of the uploaded file
  # Will read the data from the file
  # Sets the metaData() dataframe to be a single column with factored conditions
  # that werre present in column 2 of the uploaded file
  # rownames for metaData() will be the sample names present in col 1 of the uploaded file
  # 
  # Error:
  # Will return an error if the upload is not a .csv or .tsv file
  # Will return an error if the uploaded file does not have at least 2 columns to work with
  # 
  # Notes:
  # More columns are fine however only the first and second are considered
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
      if( ncol(md) < 2 ){
        showModal(modalDialog(
          tags$p(style = "color: red;","Error: Incorrect file format"), easyClose = TRUE, footer=NULL))
        return()
      }
      
      coldata <- data.frame(md[,2:ncol(md)])
      
      for (i in ncol(coldata)) {
        coldata[,i] <- as.factor(coldata[,i])
      }
      
      rownames(coldata) <- md[,1]
      
      colN <- colnames(md)[-1]
      
      colnames(coldata) <- colN
      
      print(coldata)
      
      metaData(coldata)##################################################Reactive Assignment
    })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #---Observe input Event-----# ----> run_DESeq2
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  # Process:
  # Determines whether upon input of the run differential expression button the 
  # required data has been uploaded
  # 
  # Error:
  # will return an error if the filtered_counts() object is null or as 0 rows
  # will return an error if the metaData() object is null
  # 
  # Notes:
  # Displays a modal with information and warnings about running DESeq
  # Displays a cancel or run button in response to the information above
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
    })
  # Cancels the modal and does not run DESeq if cancel button is selected
  observeEvent(input$cancel,{
    run_question(NULL)##################################################Reactive Assignment
    removeModal()
  }) 
  # Sets the run_question() object to true if selected 
  # This is used so that if run is selected multiple times with the same data DESeq is not
  # activated again and the user has to wait Also so when input$run is observed in the future
  # DESeq2 does not activate again by mistake
  # Sets kill_DEG() to TRUE to remove the option on page 2 for uploaded DEG results and starting 
  # from there
  observeEvent(input$run,{
    run_question(TRUE)##################################################Reactive Assignment
    kill_deg(TRUE)##################################################Reactive Assignment
    removeModal()
  })
  # Runs DESeq 2 upon input of input$run with mostly default parameters and single 
  # condition design.
  # Sets the ddsc() object to be a DESeqDataSet class object with the data from the experiment
  # Sets the results_ddsc() to be a dataframe consisting of the output of running results() on 
  # the ddsc() object. This should just be a dataframe with the DEG results
  # Sets vst_counts() to be the assay() data.frame of the variance stabilized transformation of the 
  # counts matrix
  # Sets vst_Obj() to be the variance stabilized counts object obtained from the vst() funciton
  observe({
    decision <- run_question()
    
    if(!is.null(decision)){
      if(decision){
        removeModal()
        
        coldata <- metaData()
        
        condition <- coldata$condition
        
        print(head(filtered_counts()))
        print(coldata)
        print(condition)
        
        ddsc <- DESeqDataSetFromMatrix(countData = round(filtered_counts()), 
                                    colData = coldata, 
                                    design = ~ condition)
        
        showModal(modalDialog("Running DESeq...", footer=NULL))
        
        d <- DESeq(ddsc)
        
        ddsc(d)##################################################Reactive Assignment
        
        results_ddsc(data.frame(results(ddsc())))##################################################Reactive Assignment
        
        showModal(modalDialog("Done!!", easyClose=TRUE, footer=NULL))
        
        removeModal()
       
        Cvst <- vst(ddsc(), blind=TRUE)
        normalized_counts(data.frame(counts(ddsc(), normalized=TRUE)))##################################################Reactive Assignment
        vst_counts(data.frame(assay(Cvst)))##################################################Reactive Assignment
        vst_Obj(Cvst)##################################################Reactive Assignment
        
      }else{
        return()
      }
    }else{
      return()
    }
    
  }) 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #---Observe input Event-----# ----> DEG_analysis_data
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  # Process:
  # On page 2 if DEG results are uploaded they will be read
  # The DEG results will be put in the results_ddsc() this should be the first 8 cols
  # of any uploaded DEG results.
  # The gene_names() object is set to the first and second columns of the DEG results dataframe
  # results_ddsc() is  set to be the DEG results -gene_name columns and with gene_id as rownames
  # Sets vst_counts() to any columns after the 8th which are designated to user variance stabalized
  # transformation counts. The counts put in vst_counts() with this method are expected to be normalized
  # vst() is not called on them.
  # 
  # Error:
  # Will return an error if.tsv or .csv are not in the file path
  # Will return an error if there are less than 8 columns indicated the required 
  # data isnt present
  # Will return an error if the header format does not match the expected format
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
     
     possible_columns <- c("gene_id","gene_name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
     
     common_columns <- intersect(possible_columns, colNames)

     #Required stuff
     #Set geneID
     if(! "gene_id" %in% common_columns){
       showModal(modalDialog(
         tags$p(style = "color: red;","Error: a gene id column is required with header gene_id"), easyClose = TRUE, footer=NULL)
       )
       return()
     }else{df <- data.frame(data %>% select(gene_id)) }
     #Set gene names
     if(! "gene_name" %in% common_columns){
       showModal(modalDialog(
         tags$p(style = "color: red;","Error: a gene name column is required with header gene_name"), easyClose = TRUE, footer=NULL)
       )
       return()
     }else{df <- df %>% mutate(gene_name = data$gene_name)}
     #set fold change
     if(! "log2FoldChange" %in% common_columns){
       showModal(modalDialog(
         tags$p(style = "color: red;","Error: a fold change column is required with header log2FoldChange"), easyClose = TRUE, footer=NULL)
       )
       return()
     }else{df <- df %>% mutate(log2FoldChange = data$log2FoldChange)}
     #Set padj
     if(! "padj" %in% common_columns){
       showModal(modalDialog(
         tags$p(style = "color: red;","Error: a padj column is required with header padj"), easyClose = TRUE, footer=NULL)
       )
       return()
     }else{df <- df %>% mutate(padj = data$padj)}
     
     #Optional stuff
     if("baseMean" %in% common_columns){df <- df %>% mutate(baseMean  = data$baseMean)}
     if("lfcSE" %in% common_columns){df <- df %>% mutate(lfcSE = data$lfcSE)}
     if("stat" %in% common_columns){df <- df %>% mutate(stat = data$stat)}
     if("pvalue" %in% common_columns){df <- df %>% mutate(pvalue = data$pvalue)}
    
     
     gene_names(df[,2:1])##################################################Reactive Assignment
    
     df <- df %>% tibble::column_to_rownames('gene_id')
     df <- df %>% select(-gene_name)
     
     results_ddsc(df)##################################################Reactive Assignment
     
     sample_columns <- data[,grep("^\\.", colnames(data), value = TRUE)]
     
     if(ncol(sample_columns) == 0){
       showModal(modalDialog(
         tags$p(style = "color: red;","Some features not available because of missing sample data"), easyClose = TRUE, footer=NULL))
       return()
     }
     
     if(! nrow(sample_columns) == nrow(gene_names())){
       showModal(modalDialog(
         tags$p(style = "color: red;","Error: number of rows in the counts data != rows in gene names list"), easyClose = TRUE, footer=NULL)
       )
       return()
     }
    
     sample_counts <- sample_columns %>% 
       mutate(gene_id = gene_names()$gene_id)%>%
       tibble::column_to_rownames('gene_id')
     
     normalized_counts(sample_counts)##################################################Reactive Assignment
     
     if(is.null(metaData())){
       showModal(modalDialog(
         tags$p(style = "color: red;","For a full analysis it is recomended to return to page one and submit a proper meta data file"), easyClose = TRUE, footer=NULL))
     }
  })
  #If metaData and the DEG analysis info has been uploaded than ddsc() object DESeq2 dataset, 
  #and a DESeqTransform object will take the place of vst_counts_VST_OBJECT
  #A DESeqTransform is similar in that you can see the assays, in this case the counts they uploaded
  #(they have not been transformed like in the vst object as uploading a DEG you assumed responsibility for
  #normilization, due to this the size factors are not available for display), 
  #you may use plotPCA to plot principle components which is the primary use of the vst_counts
  # ALL other visuals volcano and counts plots do not require this specific object.
  observeEvent(list(metaData(), input$DEG_analysis_data),{
    
    if(is.null(input$DEG_analysis_data)){return()}
    if(is.null(metaData())){return()}
    if(is.null(normalized_counts())){return()}

   
    if(nrow(metaData()) != ncol(normalized_counts())){
      showModal(modalDialog(
        tags$p(style = "color: red;","Error: Number of samples in meta data is not the same as the number of samples in the counts data"), easyClose = TRUE, footer=NULL))
      return()
    }
    
    sample_names_in_metaData <- rownames(metaData())
  
    counts <- round(normalized_counts())    
    
    colnames(counts) <- sample_names_in_metaData
    
    v <- varianceStabilizingTransformation(as.matrix(counts), blind=TRUE)
    
    vst_counts(data.frame(v))##################################################Reactive Assignment

    se <- SummarizedExperiment(assays = list(counts = v), colData = metaData())

    vsd <- DESeqTransform(se)
    
    vst_Obj(vsd)##################################################Reactive Assignment
   
  })

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #---Observe Value Event-----# ----> results_ddsc()
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
        
        df$gene_id <- row.names(df)
        
        #Some editing to display the gene names as well
        df <- left_join(df, gene_names())
        
        df <- df %>% tibble::column_to_rownames('gene_id')
        
        data.frame(df[, c((ncol(df)), 1:(ncol(df)-1))])
        
      }
      
    }, rownames = TRUE, options = list(pageLength=5))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #--------sample_conditions_Preview-------#
  #~~~~~~~~~~~~~~~~~Table~~~~~~~~~~~~~~~~~~#
  #Just prints the samples and their conditions 
  #May be helpful to ensure your file was read correctly
  output$sample_conditions_PreviewTable <- renderTable({
      if(!is.null(metaData())){t(metaData())}
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
        widget <- fileInput("DEG_analysis_data", 
          div(
            p("Differential Expression data .tsv/.csv"),
            p("Required data is: an abitrary but unique gene_id column, a gene_name column, a log2FoldChange column, and an adjusted P-value columns"),
            p("If you upload sample counts you MUST denote sample colums with sample names starting with ' . ' eg: ' .T96.s1 '"),
            p("Sample names in meta data file should match the order of the sample columns in the DEG file, names do not have to match and will be overrided by the ones in the meta data file"),
            p("Sample counts must be a variance stabalized transformation of your raw counts")
          )
        )
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
        
        present_columns <- colnames(up())
        display_columns  <- intersect(present_columns, input$display_col)
     
        t1 <- up() %>% select('gene_name', display_columns, 'padj')
        t2 <- down() %>% select('gene_name', display_columns, 'padj')
        t3 <- noR() %>% select('gene_name', display_columns, 'padj')
        
        
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
              axis.text.y = element_text(color = "black", size = 15),
              axis.title.x = element_text(size = 15)) 
      
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
  output$normalized_counts_PreviewTable <- renderUI({
    
    return_ui <- NULL
    
    if(is.null(normalized_counts())){
      return_ui <- p(span("No Data Yet",style="color: red;"))
      return_ui
    }else{
    
      normalized <- normalized_counts() %>% round(digits = 4)
      
      normalized <- normalized %>% mutate(gene_id = rownames(normalized_counts()))
      
      
      normalized <- left_join(normalized, gene_names()) %>% select(-gene_id)
      
      pvals <- results_ddsc() %>% mutate(gene_id = rownames(results_ddsc())) %>%
        left_join(gene_names()) %>% select(gene_name, padj)
      
      normalized <- left_join(normalized, pvals, relationship="many-to-many") 
      
      normalized <- normalized %>% filter(padj < as.numeric(input$pvaluePg3))
      
      normalized <- normalized[order(normalized$padj),c((ncol(normalized) - 1), ncol(normalized), 1:(ncol(normalized) - 2))]
      
      return_ui <- datatable(normalized, rownames=FALSE, options = list(pageLength=15))
      
      return_ui
      
    }
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #----------Extra_vst_count_info----------# #Size factors only right now
  #~~~~~~~~~~~~~~~~UI output~~~~~~~~~~~~~~~#
  output$Extra_normalized_count_info <- renderUI({
    
    if(is.null(normalized_counts())){
      return_ui <- p(span("No Data Yet",style="color: red;"))
      return_ui
    }else{
     
      sF <- div()
      if(!is.null(kill_deg())){
        size <- as.matrix(sizeFactors(ddsc()) %>% round(digits = 3))
        print(size)
        sF <- div(h4("Size Factors"), renderTable({t(size)}))
      }else{print("No Size Factors for pre-normalized data")}
      
      div(sF)
    }
    
  })
  
#______________________________Page 4____________________________#
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-------principle_component_plots--------#
  #~~~~~~~~~~~~~~~~plot~~~~~~~~~~~~~~~~~~~~#
  # A ggplot using the pcaPlot data
  output$pca_plot1 <- renderPlot({
    
    if(is.null(vst_Obj())){return()}
    
    if(is.null(input$pca_cond)){return()}
    
    pca_data <- plotPCA(vst_Obj(), intgroup=input$pca_cond, ntop=input$pca_nrow, returnData=TRUE)
    
    percentVar <- round(100 * attr(pca_data, "percentVar"))
    
    aesthetics <- aes(x = pca_data$PC1, y = pca_data$PC2 )
    
    x <- length(input$pca_cond)
   
    if(x >= 1){
      aesthetics$colour = as.name(input$pca_cond[1])
    }
    if(x >= 2){
      aesthetics$shape = as.name(input$pca_cond[2])
    }
    if(x >= 3){
      aesthetics$size = as.name(input$pca_cond[3])
    }
    
    pca_plot <- ggplot(pca_data,aesthetics)+
      
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      
      theme_minimal() +

      labs(title = input$title_pca_plot,
           subtitle = input$subtitle_pca_plot,
           caption = input$caption_pca_plot,
           color = if (x >= 1) input$pca_cond[1] else NULL,
           shape = if (x >= 2) input$pca_cond[2] else NULL,
           size = if (x >= 3) input$pca_cond[3] else NULL
           ) +

      theme(plot.margin = margin(10, 10, 10, 10, "pt"),
            axis.title.x = element_text(color='black', size = 20, margin = margin(15, 15, 15, 15, "pt")),
            axis.title.y = element_text(color='black', size=20, margin = margin(15, 15, 15, 15, "pt")),
            axis.text.y = element_text(color='black', size=15),
            axis.text.x = element_text(color='black', size=15),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.line.x = element_line(color='grey'),
            axis.line.y = element_line(color='grey'),

            legend.text = element_text(color='black', size=20),
            legend.title = element_text(color='black', size=20),
            legend.margin = margin(15, 15, 15, 15, "pt"),

            plot.title = element_text(color='black', size=30, margin=margin(20,20,5,10,"pt")),
            plot.subtitle = element_text(color='black', size=20, margin=margin(5,5,15,10,"pt")),
            plot.caption = element_text(color='black', size=15, margin=margin(10, 10, 10, 10, "pt")),

            aspect.ratio = 1) +
        
      coord_fixed()

    #Change size scale according to present number of variables
    if(x >= 3){
      pca_plot <- pca_plot + scale_size_discrete(range = c(3, 10)) + geom_point()
    }else{
      pca_plot <- pca_plot + geom_point(size=4)
    }
    
    pca_plot
    
  })
  #Displayable ui for the pcaPlot
  output$principle_component_plots_ui <- renderUI({
    if(is.null(vst_Obj())){
      return_ui <- span("No Data Yet",style="color: red;")
      return_ui
    }else{
        
        plotOutput('pca_plot1', height= "700px" ,width = "100%")
        
    }
  })
  #Change included genes widget
  output$change_n_pca_plot <- renderUI({
    
    m <- 500
    start <- 0
    
    if(is.null(vst_counts())){
      m <- 500
      start <- 500
    }else{
      m <- nrow(vst_counts())
      if(nrow(vst_counts()) >= 500){
        start <- 500
      }else(
        start <- m
      )
    }
    
    div(
      p("This number determines how many genes are included in the PCA ranked by variance."),
      p("If the default 500 is selected the top 500 genes with the most variance will be used."),
      p(paste("Max: ", m)),
      p("Min: 2"),
      numericInput('pca_nrow', '', 
                   value=start, min=2, max=m)
    )
  })
  #Change the metadata columns that should be included in the pcaPlot
  output$change_mColumns_for_pca <- renderUI({
    
    div(
      
      checkboxGroupInput('pca_cond', 'Meta data conditions to view', 
                       choices=colnames(metaData()),
                       selected=colnames(metaData()),
                       inline=FALSE)
    )
    
    
  })
  
  
#______________________________Page 5____________________________#
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-------------heat_map_plot--------------#
  #~~~~~~~~~~~~~~~~plot~~~~~~~~~~~~~~~~~~~~#
  
  #The heat map plot
  output$heatmap_plot1 <- renderPlot({
    if(is.null(vst_counts())){return()}
      
    vst_cor <- cor(vst_counts())
    
    #annotations
    to_annotate <- intersect(colnames(metaData()), input$heatmap_cond)
    
    annotations <- lapply(to_annotate, function(an) {
      
      selected_color <- colorRampPalette(c(input[[paste0("color1_", an)]], 
                                            input[[paste0("color2_", an)]]))(nlevels(factor(metaData()[[an]])))
     
      selected_color <- setNames(selected_color, levels(factor(metaData()[[an]])))
        
      HeatmapAnnotation(df = metaData()[, an, drop = FALSE], 
                        which = "col", 
                        annotation_name_side = "left",
                        annotation_name_gp = gpar(fontsize = 17, color='black'),
                        gap = unit(1, 'cm'),
                        height = unit(1, 'cm'),
                        col = setNames(list(selected_color), an),
                        show_legend = FALSE
                        )
    })
    
    combined_annotations <- do.call(c, annotations)
    
    legends <- lapply(to_annotate, function(an) {
      selected_color <- colorRampPalette(c(input[[paste0("color1_", an)]], 
                                           input[[paste0("color2_", an)]]))(nlevels(factor(metaData()[[an]])))
      
      selected_color <- setNames(selected_color, levels(factor(metaData()[[an]])))
      
      Legend(at = levels(factor(metaData()[[an]])),
             labels = levels(factor(metaData()[[an]])),
             legend_gp = gpar(fill = selected_color),
             title = an,
             title_gp = gpar(fontsize = 17,  col = 'black'),
             labels_gp = gpar(fontsize = 17, col = 'black'),
             ncol = 3,
             direction = "horizontal",
             grid_width = unit(1, "cm"),
             grid_height = unit(0.5,"cm")
      )
      
    })
    
    #Colors Body
    col_fun <- colorRamp2(c(0.92,1), c('white','black'))
    
    l <- NULL
    if(input$heat_body_color == "blue-red"){
      breaks <- c(0.92, 0.95, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 1)
      colors <- c("#003366", "#336699", "#6699cc", "#ffff99", "#ffe066", "#ffcc33", "#ff9933", "#ff6666", "#cc0000")
      col_fun <- colorRamp2(breaks, colors)
      
    }else if(input$heat_body_color == "REV-Rainbow"){
      col_fun <- colorRamp2(c(0.86, 0.90, 0.925, 0.935, 0.945, 0.955, 0.965, 0.975, 0.98, 0.985, 0.99, 0.995, 1), rev(rainbow(13)))
      l <- Legend(col_fun = col_fun, title_gp = gpar(fontsize = 20,  col = 'black'),
                  labels_gp = gpar(fontsize = 20, col = 'black'),legend_height = unit(6, "cm"), grid_width=unit(1,"cm"),
                  at = c(0.90, 0.925, 0.935, 0.945, 0.955, 0.965, 0.975, 0.98, 0.985, 0.99, 0.995, 1),
                  title_gap = unit(0.4, "cm"))
      
    }else if(input$heat_body_color == "white-black"){
      breaks <- c(0.94, 0.96, 0.975, 0.980,0.990,0.995, 1)
      colors <- c('white','lightgrey','grey50','ivory4','lavenderblush4', 'grey25', 'black')
      col_fun <- colorRamp2(breaks, colors)
    }
    
    x <- Heatmap(vst_cor,
            name = "Value",
            col=col_fun,

            top_annotation = combined_annotations,
            
            cluster_rows = TRUE,
            row_dend_side = "right",
            cluster_columns = TRUE,
            
            row_dend_gp = gpar(color='black', lwd = 2.25 ),
            column_dend_gp = gpar(color='black', lwd = 2.25),
            row_dend_width = unit(1.5, 'cm'),
            column_dend_height = unit(1.5, 'cm'),
            
            
            show_row_names = TRUE,
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 20, color='black'),
            
            show_column_names = TRUE, 
            column_names_side = "bottom",
            column_names_rot = 45,
            column_names_gp = gpar(fontsize = 20, color='black'),
            
            width = unit(0.7, "npc"),
            height = unit(18, "cm"),
            
            border_gp = gpar(col = "black", lwd = 1),
            rect_gp = gpar(col = "black", lwd = 1),
            
            show_heatmap_legend = FALSE
            # heatmap_legend_param = list(
            #   title_gp = gpar(fontsize=23, color='black'),
            #   title_position = 'topleft',
            #   labels_gp = gpar(fontsize=18, color='black'),
            #   legend_height = unit(6, "cm"),
            #   grid_width = unit(1, "cm")
            # )
    )
    #draw body
    draw(x, padding = unit(c(0.1, 0.03, 0.15, 0.35), "npc"))
    
    #set and draw body legend
    if(is.null(l)){
      l <- Legend(col_fun = col_fun, title = "Value", title_gp = gpar(fontsize = 17,  col = 'black'),
                  labels_gp = gpar(fontsize = 17, col = 'black'), legend_height = unit(6, "cm"), grid_width= unit(1,"cm"),
                  title_gap = unit(0.40, "cm"))
    }
    draw(l, x= unit(0.85, "npc"), y = unit(0.5, "npc"), just = c("center", "right"))
    
    #pack and draw annotation legends
    legends <- packLegend(legends, direction = "vertical")
    for (i in seq_along(legends)) {
      draw(legends[[i]], x = unit(0.90, "npc"), y = unit(0.85 - i * 0.05, "npc"), just = c("right", "top"))
    }
    
    # Add title
    grid.text(input$title_heatmap_plot, x = unit(0.3, "npc"), y = unit(0.95, "npc"), just = c("center", "top"), gp = gpar(fontsize = 30))
    # Add subtitle
    grid.text(input$subtitle_heatmap_plot, x = unit(0.3, "npc"), y = unit(0.9, "npc"), just = c("center", "top"), gp = gpar(fontsize = 20))
    
    # Add caption
    grid.text(input$caption_heatmap_plot, x = unit(0.5, "npc"), y = unit(0.05, "npc"), just = c("center", "bottom"), gp = gpar(fontsize = 15, fontface = "italic"))
    
  })
  #The ui to display the heatmap plot
  output$heatmap_plots_ui <- renderUI({
    
    if(is.null(vst_counts())){
      return_ui <- span("No Data Yet",style="color: red;")
      return_ui
    }else{
      
      fluidPage(
          plotOutput('heatmap_plot1', height= "900px" ,width = "100%")
        )
    }
  })
  #UI for selecting annotations
  output$heatmap_annotations <- renderUI({
    div(
      checkboxGroupInput('heatmap_cond', 'Meta data conditions to view', 
                         choices=colnames(metaData()),
                         selected=colnames(metaData()),
                         inline=FALSE)
    )
  })
  #displays the correct color editors for the selected annotations
  observe({
    output$color_pickers_heat <- renderUI({
      lapply(input$heatmap_cond, function(cond) {
        div(
          colourInput(inputId = paste0("color1_", cond),
                     label = paste("Choose color1 for", cond),
                     value = "grey"),
          colourInput(inputId = paste0("color2_", cond),
                      label = paste("Choose color2 for", cond),
                      value = "blue")
        )
      })
    })
  })
 
#______________________________Page 6____________________________#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-------------gene_counts_search---------#
  #~~~~~~~~~~~~~~~~ui~~~~~~~~~~~~~~~~~~~~~~#
  #res_cond_data <- merge(as.data.frame(res_cond), as.data.frame(counts(ddsc, normalized=TRUE)), by="row.names", sort=FALSE)  
  #is normalized different than vst
  
  
}#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X

