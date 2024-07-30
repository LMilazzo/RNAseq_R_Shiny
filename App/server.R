#_________________________________Shiny Server__________________________________________
server <- function(input, output) {

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_____________Reactive Values______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----
#----
  # Determines when to run DESeq2
  run_question <- reactiveVal(NULL) 
  
  # Uploaded raw counts data
  raw_counts <- reactiveVal(NULL)
  
  # Extracted gene names from the raw counts
  gene_names <- reactiveVal(NULL)
  
  # Filtered counts based on some criteria
  filtered_counts <- reactiveVal(NULL)
  
  # Uploaded metadata file
  metaData <- reactiveVal(NULL)
  
  # The DESeq2 object
  ddsc <- reactiveVal(NULL)
  
  # Dataframe of DESeq2 results
  results_ddsc <- reactiveVal(NULL)
  
  # Normalized counts from the DESeq2 object
  normalized_counts <- reactiveVal(NULL)
  
  # Variance stabilized counts object from DESeq2
  vst_Obj <- reactiveVal(NULL)
  
  # Variance stabilized counts from the vst_Obj
  vst_counts <- reactiveVal(NULL)
#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_______________Observables________________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----  
  # Reading and Setting data from raw counts upload #----
  #----- Observe UI Event --------# 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Reads the uploaded file and processes the counts matrix.
  # Sets gene_names() and raw_counts() from the processed data.
  # Filters rows with sum < 10 and sets filtered_counts().
  # Errors:
  # - Invalid file type (.csv or .tsv required)
  # - Missing 'gene_name' or 'gene_id' columns
  # - Duplicate 'gene_id' values
  # Notes: 
  # - 'gene_id' must be unique
  # - Displays the number of rows filtered
  observeEvent(input$merged_gene_counts_uploaded_file,{

    req(input$merged_gene_counts_uploaded_file)  
    
    # Inline function to set and clean counts
    setCleanCounts <- function(raw_counts_table) {

      if (!grepl('\\.(csv|tsv)$', raw_counts_table, ignore.case = TRUE)) {
        return("Bad File Type Error")
      }
      # Read as .tsv or .csv
      
      if (grepl('\\.tsv$', raw_counts_table, ignore.case = TRUE)) {
        counts <- read.csv(raw_counts_table, sep = "\t")
      } else {
        counts <- read.csv(raw_counts_table)
      }
      # Check for required columns
      if (!'gene_id' %in% colnames(counts) || !'gene_name' %in% colnames(counts)) {
        return("Missing Column Error")
      }
      
      #Check for duplicate ids
      if (anyDuplicated(counts$gene_id)) {
        return("Duplicate ids")
      }
      
      gene_names <- counts %>% select(gene_name, gene_id)
      
      # Set rownames to gene_id and remove all columns that don't start with '.'
      counts <- counts %>% tibble::column_to_rownames('gene_id')
      sample_cols <- grep("^\\.", colnames(counts), value = TRUE)
      print(sample_cols)
      non_sample_cols <- grep("^[^.]", colnames(counts), value = TRUE)
      
      # Error if there are less than 2 sample columns
      if (length(sample_cols) < 2) {
        return("Insufficient Sample Columns Error")
      }

      # Keep only sample columns
      counts <- counts %>% select(all_of(sample_cols))
      
      return(list(as.matrix(counts), gene_names))

    }
    
    # Inline function to filter counts
    filterCounts <- function(counts){
      
      keep <- rowSums(counts) > 10
      
      counts <- counts[keep,]
      
      counts <- counts
      
      counts <- as.matrix(counts)
      
      return(counts)
    }
    
    #Expect as return c(raw_counts, gene_names)           
    func_return <- setCleanCounts(input$merged_gene_counts_uploaded_file$datapath)
    
    # Successful output should be length 2 
    if (length(func_return) == 1 || is.null(func_return)) {
      showErrorModal("Possible Errors: 
      - Raw Counts upload was not valid file type (.tsv, .csv)
      - Raw Counts upload did not include required 'gene_name' and 'gene_id' columns
      - There were duplicate identifiers in the gene_id column (ids can be arbitrary but must be unique)
      - There were not sufficient sample columns, must be denoted with '.' ")
      return()
    }
    
    # ====== Reactive Values ======
    raw_counts(data.frame(func_return[1]))
    gene_names(data.frame(func_return[2]))
    filtered_counts( filterCounts(raw_counts()) )
    
    #How many rows were filtered
    diff <- nrow( raw_counts() ) - nrow( filtered_counts() )
    showModal(modalDialog(div(p(span(diff, style="color: red;"), 
        " rows were removed from the data set with row sums < 10")), 
      easyClose = TRUE, footer=NULL))

  })
  #----
  # Reading the Metadata file #----
  #----- Observe UI Event --------# 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Reads the uploaded metadata file and processes the data.
  # Sets metaData() to a dataframe with factored conditions from column 2.
  # Row names for metaData() will be the sample names from column 1.
  # Errors:
  # - Invalid file type (.csv or .tsv required)
  # - Insufficient columns (at least 2 required)
  # Notes:
  # - Only the first two columns are considered; more columns are allowed.
    observeEvent(input$meta_data_conditions_uploaded_file, {
  
      req(input$meta_data_conditions_uploaded_file)
      
      metaDataFilePath <- input$meta_data_conditions_uploaded_file$datapath
      
      # Error: Not a readable file
      if (!grepl('\\.(csv|tsv)$', metaDataFilePath, ignore.case = TRUE)) {
        showErrorModal("Error: Meta Data Table upload expected a .csv or .tsv file")
        return()
      }
      
      if (grepl('\\.tsv$', metaDataFilePath, ignore.case = TRUE)) { 
        md <- read.csv(metaDataFilePath, sep = "\t") 
      } else { 
        md <- read.csv(metaDataFilePath) 
      }
      
      # Error: Insufficient columns
      if (ncol(md) < 2 || !"sample" %in% colnames(md)) {
        showErrorModal("Error: Incorrect file format")
        return()
      }
      
      #Enforce sample naming convention
      md$sample <- ifelse(!grepl("^\\.", md$sample), paste0(".", md$sample), md$sample)
      md <- md %>%
              tibble::column_to_rownames('sample')
      
      
      md[] <- lapply(md, as.factor)
      
      # ====== Reactive Assignment ======
      metaData(md) 
      
    })
  #----
  # Running DESeq #----
  #----- Observe UI Event --------# 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Handles the initiation and execution of DESeq2 upon user interaction.
  # Checks for required data availability, prompts user confirmation, and processes the data.
  # Errors:
  # - Displays an error if `filtered_counts()` is NULL or has zero rows.
  # - Displays an error if `metaData()` is NULL.
  # Notes:
  # - Locks out DEG upload option upon running DESeq2.
    interaction_counter <- reactiveVal(0)
    design_concat <- reactiveVal(NULL)
  # - Sets DESeq2 results and related objects as reactive values for further use.
    observeEvent(input$run_DESeq2, {

      # Check for required data
      if (is.null(filtered_counts()) || nrow(filtered_counts()) <= 0) {
        showErrorModal("Merged counts table required")
        return()
      }
      
      if (is.null(metaData())) {
        showErrorModal("Meta Data table required")
        return()
      }
      
      # Prompt user confirmation
      showModal(modalDialog(
        div(
          p("Would you like to generate a DEG dataset with the given raw counts & meta data?"), 
          p("This action will lock you out of uploading an already generated DEG table and will run DESeq2 with default parameters for your data."),
          p("If you have previously uploaded a DEG file, it will be scrapped from the application. All tables will display the new DESeq2 data, and visuals and graphs will also be updated accordingly.")
        ),
        footer = tagList(
          actionButton('cancel', 'Cancel'),
          actionButton('run', 'Run')
        )
      ))

      # Observe cancel and run buttons
      observeEvent(input$cancel, {
        removeModal()
      }, once = TRUE) # `once = TRUE` ensures this observer is only triggered once per modal

      observeEvent(input$run, {

        # ====== Reactive Assignment ======
        kill_deg(TRUE) 
        removeModal()
        
        #Design Creation A lot of parts / quite complex 
        #----
        vars <- colnames(metaData())
        
        output$inter_space <- renderUI({
          if(interaction_counter() <= 0 ){
            return()
          }
          lapply(1:interaction_counter(), function(i) {
            createInteractionUI(i, vars)
          })
        })
        createInteractionsUI <- function(index, var){
          fluidRow(
            column(6, selectInput(paste0("interaction_", index, "_a"), 
                                  paste("Interaction Term", index, "Variable 1:"), 
                                  choices = choices)),
            column(6, selectInput(paste0("interaction_", index, "_b"), 
                                  paste("Interaction Term", index, "Variable 2:"), 
                                  choices = choices))
          )
        }
        
        showModal(modalDialog(div(
          p("Experimental design"),
          renderText({design_concat()}),
          checkboxGroupInput('dvars', 'Variables', choices = vars, inline = TRUE),
          uiOutput("inter_space"),
          actionButton('addinter', "Add interaction"),
          actionButton('mininter', "Remove interaction"),
          actionButton('design_submit', "RUN")
        ),
        easyClose = FALSE))
        
        refresh <- reactive({
          invalidateLater(100)
          Sys.time()
        })
        observeEvent(refresh(),{design_concat(concatDesign(input, output, session, interaction_counter()))})
        observeEvent(input$addinter, {
          interaction_counter(interaction_counter() + 1)
        })
        observeEvent(input$mininter, {
          if(interaction_counter() > 0){
            interaction_counter(interaction_counter() - 1)
          }
        })
        concatDesign <- function(input, output, session, i){
          
          main_effects <- input$dvars
          
          dparts <- c(main_effects)
          
          if(i > 0){
            for (j in 1:i) {
              var1 <- input[[paste0("interaction_", j, "_a")]]
              var2 <- input[[paste0("interaction_", j, "_b")]]
              
              # If both variables for the interaction are selected, add the interaction term
              if (!is.null(var1) && !is.null(var2)) {
                interaction_term <- paste(var1, var2, sep = ":")
                dparts <- c(dparts, interaction_term)
              }
            }
          }
          
          design_formula <- paste(dparts, collapse = " + ")
          
          return(design_formula)
          
        }
        
        #----
        
        observeEvent(input$design_submit,{
        
          # Run DESeq2 process
          metadata <- metaData()
          countdata <- filtered_counts()

          # condition <- metadata[,1]

          tryCatch({
            
            des <- paste0("~ ", design_concat())
            des <- as.formula(des)
            print(des)
            
            ddsc <- DESeqDataSetFromMatrix(countData = round(countdata),
                                          colData = metadata,
                                          design = des)

            showModal(modalDialog("Running DESeq...", footer = NULL))

            deseqdataset <- DESeq(ddsc)

            # ====== Reactive Assignment ======
            ddsc(deseqdataset)
            results_ddsc(data.frame(results(ddsc())))
            vstObject <- vst(ddsc(), blind = TRUE)
            normalized_counts(data.frame(counts(ddsc(), normalized = TRUE)))
            vst_counts(data.frame(assay(vstObject)))
            vst_Obj(vstObject)

            showModal(modalDialog("DESeq complete with no fatal errors", easyClose = TRUE, footer = NULL))

          }, error = function(e) {

            showErrorModal(paste("Error in DESeq2 process:", e$message))

          })
          
        }, once = TRUE)
      }, once = TRUE) # `once = TRUE` ensures this observer is only triggered once per modal
    
    })
  #----
  # Handling a Case of running with DEG data #----
  #----- Observe UI Event --------# 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Handles the uploading and processing of precomputed DEG analysis data.
  # Sets the results_ddsc(), gene_names(), normalized_counts(), vst_counts(), and vst_Obj() 
  # reactive values based on the uploaded data.
  # Errors:
  # - Displays an error if the file is not a .csv or .tsv.
  # - Displays an error if required columns (gene_id, gene_name, log2FoldChange, padj) are missing.
  # - Displays an error if there are duplicate values in the gene_id column.
  # - Displays an error if the number of rows in the counts data does not match the number of rows in gene names.
  # - Displays a warning if sample columns are missing.
  # Notes:
  # - vst_Obj in this case holds a DESeqTransform object instead of a VST object.
  # - determines the status of what has been uploaded and processed PREVENTS reprocessing the same data
    degStatus <- reactiveVal(1)
    # 1 metadata has been uploaded waiting for DEG or no files uploaded
    # 2 DEG has been uploaded wating for metadata
    # 3 Both files are uploaded a file has potentially been changed
  # - determines whether or not this process is available to the user
    kill_deg <- reactiveVal(NULL)
  # - Recommends uploading a proper metadata file for full analysis.
    observeEvent(list(metaData(), input$DEG_analysis_data),{
      
      req(input$DEG_analysis_data)

      if( degStatus() == 1  || degStatus() == 3 ){

        filePath <- input$DEG_analysis_data$datapath

        #Handle for filetype
        if (!grepl('\\.(csv|tsv)$', filePath, ignore.case = TRUE)) {
          showErrorModal("Error: DEG data upload expected a .csv or .tsv file")
          return()
        }

        #Read file
        if (grepl('\\.tsv$', filePath, ignore.case = TRUE)) { 
          data <- read.csv(filePath, sep = "\t") 
        } else { 
          data <- read.csv(filePath) 
        }
        
        #Get data column names
        req_col <- c( "gene_id" , "gene_name" , "log2FoldChange" , "padj" )
        possible_col <- c("gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
        col <- colnames(data)
        
        common_col <- intersect(possible_col, col)
  
        if(! all(req_col %in% common_col) ){
          showErrorModal("Error: data missing either ('gene_id', 'gene_name', 'log2FoldChange', 'padj')")
          return()
        }

        df <- data %>%
                  select( all_of(common_col) )

        # ====== Reactive Assignment ======
        
        gene_names(df %>% select(gene_name, gene_id)) 
        
        #Check for duplicate ids
        if (anyDuplicated(gene_names()$gene_id)) {
          showErrorModal("Error: There appears to be duplicate gene names")
          return()
        }

        # Process DEG data
        df <- df %>%
          tibble::column_to_rownames('gene_id') %>%
          select(-gene_name)
        
        #====== Reactive Assignment ======
        results_ddsc(df) 

        # Process sample columns
        sample_columns <- data[, grep("^\\.", colnames(data), value = TRUE)] # all columns starting with '.'

        if (ncol(sample_columns) == 0) {
          showErrorModal("Warning: Some features not available because of missing sample data")
          return()
        }

        if (!nrow(sample_columns) == nrow(gene_names())) {
          showErrorModal("Error: number of rows in the counts data != rows in gene names list")
          return()
        }

        # Set row names for sample columns
        rownames(sample_columns) <- gene_names()$gene_id

        #====== Reactive Assignment ======
        normalized_counts(sample_columns) 
        degStatus(2)

        # Metadata upload reminder if no metadata return and cancel process
        if (is.null(metaData())) {
          showErrorModal("Warning: For a full analysis, it is recommended to return to page one and submit a proper meta data file")
          return()
        }
        
      }
      
      if( degStatus() == 2 || degStatus() == 3 ){

        if( nrow(metaData()) != ncol(normalized_counts()) ){
          showErrorModal("Error: Number of samples in meta data is not the same as the number of samples in the counts data")
          return()
        }
      
        counts <- round(normalized_counts())    
        
        v <- varianceStabilizingTransformation(as.matrix(counts), blind=TRUE)
        se <- SummarizedExperiment(assays = list(counts = v), colData = metaData())
        vsd <- DESeqTransform(se)

        #====== Reactive Assignment ======
        vst_counts(data.frame(v))
        vst_Obj(vsd)
        degStatus(3)
      }else{
        return()
      }
      
    })
  


#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######____________Output $ Objects______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----
  #______________________________Page 1____________________________#----

    # Raw Counts data table 
    output$raw_counts_PreviewTable <- renderUI({
        
        df <- raw_counts()
        
        table <- renderUI({span("No Raw counts to display",style="color: red;")})

        if(!is.null(df)){
          
          #Some editing to display the gene names as well
          df$gene_id <- row.names(df)
          df <- left_join(df, gene_names()) %>% 
                tibble::column_to_rownames('gene_id')
          
          #Organize columns
          df <- data.frame(df[, c((ncol(df)), 1:(ncol(df)-1))])
          
          table <- renderDT({df}, rownames = TRUE, options = list(pageLength=5))
        }

        div(table)

    })
    
    # render table for metadata
    output$sample_conditions_PreviewTable <- renderUI({
        if(!is.null(metaData())){ renderTable({metaData()}, rownames  = TRUE) }
        else{span("Upload Meta Data",style="color: red;")}
    })
    

  #______________________________Page 2____________________________#----
    
    # Option for uploading DEG data (Only available additionally)
    output$deg_upload <- renderUI({
        
      if(is.null(kill_deg())){
        fileInput("DEG_analysis_data", 
            div(
              p("Differential Expression data .tsv/.csv"),
              p("Required data is: an abitrary but unique gene_id column, a gene_name column, a log2FoldChange column, and an adjusted P-value columns"),
              p("If you upload sample counts you MUST denote sample colums with sample names starting with ' . ' eg: ' .T96.s1 '"),
              p("Sample names in meta data file should match the order of the sample columns in the DEG file, names do not have to match and will be overrided by the ones in the meta data file"),
              p("Sample counts must be a variance stabalized transformation of your raw counts")
            )
          )
        }else{
          return()
        }
    })
    
    # 3 tables feature genes that fit in expression categories
    output$expression_tables <- renderUI({

      if(is.null(results_ddsc()) || is.null(gene_names())){
        return(renderUI({span("No Data Yet",style="color: red;")}))
      }
      
      #Splits the DEG list into groups by fold change 
      func_return <- splitByExpr(results_ddsc(), gene_names(), input$cutOffs )
      up <- data.frame(func_return[1]) #up regulated
      down <- data.frame(func_return[2]) #downregulated
      noR <- data.frame(func_return[3]) #Inbetween (no regulation)

      p <- as.numeric(input$pvaluePg2)
      
      present_columns <- colnames(up)
      display_columns  <- intersect(present_columns, input$display_col)
    
      t1 <- up %>% 
                select('gene_name', display_columns, 'padj') %>% 
                filter(padj < p)

      t2 <- down %>% 
                select('gene_name', display_columns, 'padj') %>% 
                filter(padj < p)

      t3 <- noR %>% 
                select('gene_name', display_columns, 'padj') %>% 
                filter(padj < p)
      div(
        h2("Up Regulated Genes: ",  nrow(t1) ), 
        
        datatable(t1, rownames=FALSE, option=list(pageLength=7)),
        
        h2("Down Regulated Genes: ", nrow(t2) ),
        
        datatable(t2, rownames=FALSE, option=list(pageLength=7)),
        
        h2("Genes Between Selected Cutoffs: ", nrow(t3) ),
        
        datatable(t3, rownames=FALSE, option=list(pageLength=7))
      )
    
    })
    
    # A simple histogram that shows the distribution of fold change
    output$Distribution_Histogram_ui <- renderUI({
        if(is.null(results_ddsc())){
          return()
        }
        
        plot <- renderPlot({
          p <- as.numeric(input$pvaluePg2)
          data <- results_ddsc() %>% 
                              filter(padj < p) %>% 
                              select('log2FoldChange')
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
                  axis.title.x = element_text(size = 15)
          ) 
        })

        div(
          HTML("<h4>Distribution of Fold Change</h4>"),
          plot
        )
        
      })
    
  #______________________________Page 3____________________________#----
    
    # A data table showing normalized counts
    output$normalized_counts_PreviewTable <- renderUI({
      
      if(is.null(normalized_counts())){
        return(span("No Data Yet",style="color: red;"))
      }
      
      p <- as.numeric(input$pvaluePg3)

      normalized <- normalized_counts() %>% 
                                      round(digits = 4) %>% 
                                      mutate(gene_id = rownames(normalized_counts())) 
      
      pvals <- results_ddsc() %>% 
                            mutate(gene_id = rownames(results_ddsc())) %>%
                            left_join(gene_names()) %>% 
                            select(gene_name, padj)
      
      normalized <- left_join(normalized, gene_names()) %>% 
                                                select(-gene_id) 
                                                
      normalized <- left_join(normalized, pvals, 
                    relationship="many-to-many") %>% 
                    filter(padj < p)
      
      normalized <- normalized[order(normalized$padj),c((ncol(normalized) - 1), ncol(normalized), 1:(ncol(normalized) - 2))]
      
      datatable(normalized, rownames=FALSE, options = list(pageLength=15))
      
    })
    
    # A sidebar showing extra info like size factors if they are available
    output$Extra_normalized_count_info <- renderUI({
      
      if(is.null(normalized_counts())){
        return()
      }
      
      if(!is.null(kill_deg())){
        
        size <- as.matrix(sizeFactors(ddsc()) %>% 
                                            round(digits = 3))
        div(
          h4("Size Factors"), 
          renderTable({t(size)})
        )
      }else{
        p("No Size Factors for pre-transformed data")
      }
    })
    
  #______________________________Page 4____________________________#----
    
    # The principle component plot
    output$principle_components <- renderPlot({
      
      if(is.null(vst_Obj())){return()}
      if(length(input$pca_cond) == 0){return()}
      
      pca_plot <- principlePlot(vst_Obj(), 
                                input$pca_cond, 
                                input$pca_nrow, 
                                input$title_pca_plot, 
                                input$subtitle_pca_plot, 
                                input$caption_pca_plot)
      
      pca_plot
      
    })
    # UI holding the plot
    output$principle_component_plots_ui <- renderUI({
      
      if(is.null(vst_Obj()) || length(input$pca_cond) == 0){
        return(span("No Data Yet, Or No Selected Factors",style="color: red;"))
      }

      plotOutput('principle_components', height= "700px", width = "100%")
    
    })
    # Widget for changing genes using in plot
    output$pca_n_counts <- renderUI({
      
      #Sets the number of genes to 500, or the maximum
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
    #Widget for changing the metadata columns
    output$pca_metadata <- renderUI({
      div(
        checkboxGroupInput('pca_cond', 'Meta data conditions to view', 
                        choices=colnames(metaData()),
                        selected=colnames(metaData()),
                        inline=FALSE)
      )
    })
    
  #______________________________Page 5____________________________#----
    
    # The heatmap for correlation analysis
    output$heatmap_plot1 <- renderPlot({
      if(is.null(vst_counts())){return()}
        
      #Correlation Matrix
      vst_cor <- cor(vst_counts())
      
      #Annotations
      to_annotate <- intersect(colnames(metaData()), input$heatmap_cond)
      
      annotations <- lapply(to_annotate, function(an) {
        
        selected_color <- colorRampPalette(c(input[[paste0("color1_", an)]], 
                                              input[[paste0("color2_", an)]])
                                          )(nlevels(factor(metaData()[[an]])))
      
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
      

      #Legends for Annotations
      legends <- lapply(to_annotate, function(an) {
        
        selected_color <- colorRampPalette(c(input[[paste0("color1_", an)]], 
                                            input[[paste0("color2_", an)]])
                                          )(nlevels(factor(metaData()[[an]])))
        
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
      
      #Colors Body default case
      col_fun <- colorRamp2(c(0.92,1), c('white','black'))
    
      if(input$heat_body_color == "blue-red"){
        breaks <- c(0.92, 0.95, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 1)
        colors <- c("#003366", "#336699", "#6699cc", "#ffff99", "#ffe066", "#ffcc33", "#ff9933", "#ff6666", "#cc0000")
        col_fun <- colorRamp2(breaks, colors)
        l <- l <- Legend(col_fun = col_fun, title = "Value", title_gp = gpar(fontsize = 17,  col = 'black'),
                    labels_gp = gpar(fontsize = 17, col = 'black'), legend_height = unit(6, "cm"), grid_width= unit(1,"cm"),
                    title_gap = unit(0.40, "cm"))
      }
      
      else if(input$heat_body_color == "REV-Rainbow"){
        
        col_fun <- colorRamp2(c(0.86, 0.90, 0.925, 0.935, 0.945, 0.955, 0.965, 0.975, 0.98, 0.985, 0.99, 0.995, 1), rev(rainbow(13)))
        l <- Legend(col_fun = col_fun, title_gp = gpar(fontsize = 20,  col = 'black'),
                    labels_gp = gpar(fontsize = 20, col = 'black'),legend_height = unit(6, "cm"), grid_width=unit(1,"cm"),
                    at = c(0.90, 0.925, 0.935, 0.945, 0.955, 0.965, 0.975, 0.98, 0.985, 0.99, 0.995, 1),
                    title_gap = unit(0.4, "cm"))
        
      }
      
      else if(input$heat_body_color == "white-black"){
        breaks <- c(0.94, 0.96, 0.975, 0.980,0.990,0.995, 1)
        colors <- c('white','lightgrey','grey50','ivory4','lavenderblush4', 'grey25', 'black')
        col_fun <- colorRamp2(breaks, colors)
        l <- l <- Legend(col_fun = col_fun, title = "Value", title_gp = gpar(fontsize = 17,  col = 'black'),
                    labels_gp = gpar(fontsize = 17, col = 'black'), legend_height = unit(6, "cm"), grid_width= unit(1,"cm"),
                    title_gap = unit(0.40, "cm"))
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
            
      )

      #draw body
      draw(x, padding = unit(c(0.1, 0.03, 0.15, 0.35), "npc"))
      
      #Draw body legend
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
    # The ui to display the heatmap plot
    output$heatmap_plots_ui <- renderUI({
      
      if(is.null(vst_counts())){
        return(span("No Data Yet",style="color: red;"))
      }  
      fluidPage(
          plotOutput('heatmap_plot1', height= "900px" ,width = "100%")
      )
    })
    # UI for selecting annotations
    output$heatmap_annotations <- renderUI({
      div(
        checkboxGroupInput('heatmap_cond', 'Meta data conditions to view', 
                          choices=colnames(metaData()),
                          selected=colnames(metaData()),
                          inline=FALSE)
      )
    })
    # displays the correct color editors for the selected annotations
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
  
  #______________________________Page 6____________________________#----
    
    # A ui bar on the left displaying some important genes 
    output$leftbar_gene_plots <- renderUI({
      
      if(is.null(results_ddsc()) || is.null(normalized_counts())){
        return()
      }

      p <- as.numeric(input$pvaluePg6)
      
      counts <- data.frame(round(normalized_counts())) %>%
                            mutate(gene_id = rownames(normalized_counts()))
      
      counts <- left_join(counts, gene_names()) 
      counts <- left_join(counts, (results_ddsc() %>%
                                     mutate(gene_id = rownames(results_ddsc())) %>%
                                     select(padj, log2FoldChange, gene_id))
                          ) %>% filter(padj <= p)
      
      info_col <- c('gene_id', 'gene_name', 'padj', 'log2FoldChange')
      counts_col <- counts %>% select(-all_of(info_col))
      counts <- counts[c(info_col, colnames(counts_col))]
      
      tp <- counts[order(counts$padj),] %>% head(1)
      lf <- counts[order(-counts$log2FoldChange),] %>% head(1)
      hf <- counts[order(counts$log2FoldChange),] %>% head(1)
      
      top_pval <- gplop(tp, metaData(), colnames(metaData()))
      
      low_fold <- gplop(lf, metaData(), colnames(metaData()))
      
      high_fold <- gplop(hf, metaData(), colnames(metaData()))
      
      div(
        
        h2("Lowest adjusted pval"),
        
        renderPlot({top_pval}),
        
        h2("Lowest Fold Change"),
        
        renderPlot({high_fold}),
        
        h2("Highest Fold Change"),
        
        renderPlot({low_fold})
        
      )

    })
    # A ui bar on the right displaying a datatable of all the genes
    output$gene_name_list <- renderUI({
      
      if(is.null(gene_names())){

        span("No Data Yet",style="color: red;")
        return()
      }

      names <- datatable(gene_names() %>% 
                        select('gene_name'),
                        rownames=FALSE, 
                        option=list(pageLength=25))
      div(
        h3("Gene List"),
        names
      )

    })
    # A ui for displaying the search method 
    output$gene_count_search <- renderUI({
      
      if(is.null(results_ddsc()) || is.null(normalized_counts())){
        return(span("No Data Yet",style="color: red;"))
      }
      
      counts <- data.frame(round(normalized_counts())) %>%
        mutate(gene_id = rownames(normalized_counts()))
      
      counts <- left_join(counts, gene_names()) 
      counts <- left_join(counts, (results_ddsc() %>%
                                     mutate(gene_id = rownames(results_ddsc())) %>%
                                     select(padj, log2FoldChange, gene_id))
      )
      
      #sort df 
      info_col <- c('gene_id', 'gene_name', 'padj', 'log2FoldChange')
      counts_col <- counts %>% select(-all_of(info_col))
      counts <- counts[c(info_col, colnames(counts_col))]
      
      div(
        
        textInput("gene_count_plot_search",
                  "Search and plot a gene", 
                  NULL),
        
        checkboxGroupInput('gene_count_plot_search_cond', 
                            'Meta data conditions to view, first 4 are considered with mapping order (x, color, shape, size)', 
                            choices=colnames(metaData()),
                            selected=colnames(metaData()),
                            inline=TRUE),
        renderUI({
          
          if(!is.null(input$gene_count_plot_search)){
            gene <- counts %>% filter(gene_name == input$gene_count_plot_search)
            
            if(is.null(gene) || ncol(gene) == 0 || nrow(gene) == 0){
              return()
            }
            if(nrow(gene) > 1 ){
              gene <- gene[1,]
            }
            
            x <- gplop(gene, metaData(), input$gene_count_plot_search_cond)
            renderPlot({x})
            
          }else{
            p("Nothing Here Yet")
          }
        })

      )
    })
    
  #______________________________Page 7____________________________#----
  
    # The plot for the volcano plot page
    output$volcano_plot <- renderPlot({
      
      if(is.null(results_ddsc())){ return() }
      
      if(!is.null(input$volc_search)){
        search <- strsplit(input$volc_search, ",")[[1]]
        search <- trimws(search)
      }else{
        search <- NULL
      }
      
      deg <- results_ddsc() %>% 
        mutate(gene_id = rownames(results_ddsc()))
      
      deg <- left_join(deg, gene_names()) %>% 
        select(gene_name, log2FoldChange, padj)
      
      deg$padj <- as.numeric(deg$padj) 
      
      lowcut <- input$volcano_cutoffs[1]
      highcut <- input$volcano_cutoffs[2] 
      pcut <- input$pvaluePg7 
      pop_score <- input$volcano_pop 
      lab_score <- input$volcano_lab_density 
      title <- input$title_volc_plot 
      subtitle <- input$subtitle_volc_plot 
      caption <- input$caption_volc_plot
      
      plot <- erupt( deg, 
                     lowcut, 
                     highcut, 
                     pcut, 
                     pop_score, 
                     lab_score, 
                     search, 
                     title, 
                     subtitle, 
                     caption)
      
      plot
      
    })
    # A ui for displaying the volcano plot
    output$volcano_plot_ui <- renderUI({
      if(is.null(results_ddsc())){
        return(span("No Data Yet",style="color: red;"))
      }
      
      plotOutput('volcano_plot', height= "900px", width = "100%")
      
    })
     
   
  
#----
}#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X

