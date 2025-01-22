#_________________________________Shiny Server__________________________________________
server <- function(input, output, session) {

  #reload app
  observeEvent(input$reload_app,{
    runjs("location.reload();")
  })
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_____________Reactive Values______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----
  
  # Uploaded raw counts data
  RCV.RAW_COUNTS <- reactiveVal(NULL)
  
  # Extracted gene names from the raw counts
  RCV.GENE_NAMES <- reactiveVal(NULL)
  
  # Filtered counts based on some criteria
  RCV.FILTERED_COUNTS <- reactiveVal(NULL)
  
  # Uploaded metadata file
  RCV.META_DATA <- reactiveVal(NULL)
  
  # The DESeq2 object
  RCV.DDSC <- reactiveVal(NULL)
  
  # Dataframe of DESeq2 results
  RCV.RESULTS_DDSC <- reactiveVal(NULL)
  
  # Normalized counts from the DESeq2 object
  RCV.NORMALIZED_COUNTS <- reactiveVal(NULL)
  
  # Variance stabilized counts object from DESeq2
  RCV.VST_OBJ <- reactiveVal(NULL)
  
  # Variance stabilized counts from the vst_Obj
  RCV.VST_COUNTS <- reactiveVal(NULL)
  
  # The options for selectable contrasts
  RCV.CONTRAST_OPTIONS <- reactiveVal(NULL)
  
  # Pathfinder results
  RCV.PATHWAY_RESULTS <- reactiveVal(NULL)
  
  # Pathfinder old experiment uploads
  RCV.PATHWAY_ABUNDANCE <- reactiveVal(NULL)
  
  #Abunance data gene names
  RCV.PATHWAY_GENES <- reactiveVal(NULL)
  
  #Pathfinder meta data
  RCV.PATHWAY_META_DATA <- reactiveVal(NULL)
  
#---- 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_______________Running Workflow___________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Differential Expression Things______________________________________________

  #-----------NEW EXPERIMENT-------
  
  #UPLOAD MODAL MESSAGE
  observeEvent(input$start_new_experiment,{
    showModal(
      modalDialog(
        div( 
          h5("Counts Matrix"),
          p("Accepted File Types:  .csv   .tsv "),
          p("Format: ( gene_id, gene_name, .samples... ) "),
          p("The gene id column should contain a unique identifier for every row"),
          fileInput("counts_file", "Counts Matrix"),
        
          h5("Sample Conditions Table"),
          p("Accepted File Types:  .csv   .tsv "),
          p("Format: (sample, conditions...) "),
          p("Sample column should include sample names that are found in your experiment data"),
          fileInput("deg_meta_1_file", "Sample Conditions")
        ),
        footer = actionButton("read_cnts_and_meta1", "Finish"), easyClose = TRUE
      )
    )
  })
  
  #READING FILES RAW COUNTS AND META DATA
  observeEvent(input$read_cnts_and_meta1, {
    
    removeModal()
    
    #               REQUIRE BOTH IMPORTS
    #_____________________________________________________
    if( is.null(input$counts_file) || is.null(input$deg_meta_1_file)){
      showErrorModal('A file was not uploaded properly try again')
      return()
    }
  
    #               READ THE COUNTS FILE FIRST
    #     return(gene_names, raw_counts, filtered_counts)
    #____________________________________________________________
    df_list <- readCountsUpload(input$counts_file$datapath)
    if(is.null(df_list)){ return() }
    
    gene_names <- df_list[[1]]
    raw_counts <- df_list[[2]]
    filtered_counts <- df_list[[3]]
    
    #               READ THE META DATA FILE
    #__________________________________________________________________
    meta_data <- readMetaData(input$deg_meta_1_file$datapath)
    if(is.null(meta_data)){ return() } 
  
    #                 SET GLOBAL VARIABLES
    #_____________________________________________________
                    RCV.GENE_NAMES(gene_names)
                    RCV.RAW_COUNTS(raw_counts)
                RCV.FILTERED_COUNTS(filtered_counts)
                     RCV.META_DATA(meta_data)
    
    #             SHOW HOW MANY WERE FILTERED
    #_____________________________________________________
    diff <- nrow(raw_counts) - nrow(filtered_counts)
    
    showModal(modalDialog(div(p(span
      (diff, style="color: red;"), " rows were filtered from the dataset"
    )),easyClose = TRUE, footer=NULL))
    
    #               ADJUST UI
    #_____________________________________________________
    shinyjs::hide('Page1_Upload_Options')
    shinyjs::show('run_DESeq2')
    shinyjs::show('pg1table2')
    shinyjs::show('pg1table1')
    shinyjs::hide('about_DESeq2')
    
  })
  
  #RUNNING DESEQ2
  
  #FORMULA INTERACTION UI TRACKER
  RCV.INT_COUNTER <- reactiveVal(0)
  
  #FORMULA DESIGN
  RCV.DESIGN_FORMULA_STRING <- reactiveVal(NULL)
  
  #DESIGN FORMULA CREATION MODAL
  observeEvent(input$run_DESeq2,{
    
    #     SET LOCAL SAMPLE NAMES
    #__________________________________
    vars <- colnames( RCV.META_DATA() )
    
    #INTERACTION ASSEMBLY===============================================
    
    #         FUNCTION TO CREATE AN INTERACTION SELECTION ELEMENT
    #___________________________________________________________________
    createInteractionUI <- function(index, vars){
      
      fluidRow(
        column(6, selectInput(paste0("interaction_", index, "_a"), 
                              paste("Interaction Term", index, "Variable 1:"), 
                              choices = vars)),
        column(6, selectInput(paste0("interaction_", index, "_b"), 
                              paste("Interaction Term", index, "Variable 2:"), 
                              choices = vars))
      )
      
    }
    
    # UI ADDING A INTERACTION SELECTION ELEMENT OF EACH INTERACTION NEEDED
    #___________________________________________________________________
    output$inter_space <- renderUI({
      
      #IF NO INTERACTIONS SHOW NOTHING
      if( RCV.INT_COUNTER() <= 0 ){ return() }
      
      #CREATE A UI INTERACTION ELEMENT FOR ALL INTERACTIONS NEEDED
      lapply(1: RCV.INT_COUNTER(), function(index) { createInteractionUI(index, vars) })
      
    })
    
    #ADJUSTING INTERACTIVE DESIGN======================================
    
    #           CLOCK USED TO REFRESH FORMULA PERIODICALLY
    #___________________________________________________________________
    refresh <- reactive({
      invalidateLater(100)
      Sys.time()
    })
    
    #     FUNCTION FOR CONCATANATION OF DATA INTO FORMULA STRING
    #___________________________________________________________________
    concatDesign <- function(input, output, session, INT_COUNTER){
      
      #MAIN FORMULA ELEMENTS  X + Y
      main_effects <- input$main_effects
      design_parts <- c(main_effects)
      
      #ADD ANY INTERACTION FORMULAS TO THE DESIGN PART LIST
      if(INT_COUNTER > 0){
        for (j in 1:INT_COUNTER) {
          var1 <- input[[paste0("interaction_", j, "_a")]]
          var2 <- input[[paste0("interaction_", j, "_b")]]
          # If both variables for the interaction are selected, add the interaction term
          if (!is.null(var1) && !is.null(var2)) {
            interaction_term <- paste(var1, var2, sep = ":")
            design_parts <- c(design_parts, interaction_term)
          }
          
        }
      }
      
      #     FORMULATE A STRING OF ALL DESIGN PARTS SEP BY +
      #_______________________________________________________
      design_formula <- paste(design_parts, collapse = " + ")
      return(design_formula)
      
    }
    
    #         UPON CLOCK UPDATE REFRESH THE FORMULA STRING
    #___________________________________________________________________
    observeEvent(refresh(),{
      RCV.DESIGN_FORMULA_STRING(concatDesign(input, output, session, RCV.INT_COUNTER()))
    })
    
    #     ADD AND SUBTRACT INTERACTION SELECTION UI ELEMENTS
    #___________________________________________________________________
    observeEvent(input$plusinter, {
      RCV.INT_COUNTER(RCV.INT_COUNTER() + 1)
    })
    observeEvent(input$mininter, {
      if(RCV.INT_COUNTER() > 0){
        RCV.INT_COUNTER(RCV.INT_COUNTER() - 1)
      }
    })
    
    #MAIN UI MODAL FOR==================================================
    
    showModal(modalDialog(
      div(
        p("Experimental design"),
        renderText({RCV.DESIGN_FORMULA_STRING()}),
        checkboxGroupInput('main_effects', 'Variables', choices = vars, inline = TRUE),
        uiOutput("inter_space"),
        actionButton('plusinter', "Add interaction"),
        actionButton('mininter', "Remove interaction"),
      ),
      footer = tagList(actionButton('design_submit', "Run", 
      style = "background-color: #4CAF50; color: black; border-color: white; background-image: none;")), 
      easyClose = TRUE
    ))
  
  })
  
  #SUBMIT DESIGN AND RUN DESeq2
  observeEvent(input$design_submit,{
    
    #   COPY INFO TO LOCAL VARIABLES
    # ___________________________________
    metadata <- RCV.META_DATA()
    countdata <- RCV.FILTERED_COUNTS() %>% as.matrix()
    
    #   CREATE DESIGN AS A FORMULA
    #___________________________________
    des_formula <- paste0("~ ", RCV.DESIGN_FORMULA_STRING())
    des_formula <- as.formula(des_formula)
    
    #   DEBUGGING PRINTING 
    #_________________________
    print(des_formula)
    print(dim(countdata))
    
    
    tryCatch({
      
      print(metadata)
      print(class(metadata[,1]))
      
      #         CREATE THE DESeq DATA SET   
      #_____________________________________________________________
      ddsc <- DESeqDataSetFromMatrix(countData = round(countdata),
                                     colData = metadata,
                                     design = des_formula)
      print(ddsc)
      
      showModal(modalDialog("Running DESeq...", footer = NULL))
      
      YourDESeqDataSet <- DESeq(ddsc)
      
      #                     SELECT CONTRASTS
      #_____________________________________________________________
      contrast_list <- list()
      #Make options
      for (col_name in intersect(colnames(colData(ddsc)), all.vars(design(ddsc))) ){
        levels <- levels(colData(ddsc)[[col_name]])
        combinations <- combn(levels, 2, simplify = FALSE)
        for (combo in combinations) {
          combo_name1 <- paste(combo, collapse = " vs. ")
          combo_name2 <- paste(rev(combo), collapse = " vs. ")
          contrast_list[[combo_name1]] <- c(col_name, combo)
          contrast_list[[combo_name2]] <- c(col_name, rev(combo))
        }
      }
      
      RCV.CONTRAST_OPTIONS(contrast_list)
    
      def <- contrast_list[[1]]
      
      #         SET GLOBAL VARIABLES FOR THE EXPERIMENT RESULTS
      #___________________________________________________________________
      RCV.DDSC(YourDESeqDataSet)
                                
      RCV.RESULTS_DDSC( as.data.frame(results(YourDESeqDataSet, contrast = def)))
      
      RCV.VST_OBJ( vst(YourDESeqDataSet, blind = TRUE, nsub = 50) )
    
      RCV.VST_COUNTS( as.data.frame(assay(RCV.VST_OBJ())) )
      
      RCV.NORMALIZED_COUNTS(as.data.frame(counts(YourDESeqDataSet, normalized = TRUE)) )
      
      #           EXIT MESSAGE
      #_____________________________
      showModal(modalDialog("DESeq complete with no fatal errors", easyClose = TRUE, footer = NULL))
      
      #           UPDATE UI
      # _______________________________
      shinyjs::hide("run_DESeq2")
      shinyjs::hide('preRun_Data_Preview')
      shinyjs::show("TabSet1_Diffrential_Expression_Analysis")
        
        
    }, error = function(e) {
      
      showErrorModal(paste("Error in DESeq2 process:", e$message))
      
    })
    
  })
  
  observeEvent(input$contrast, {
    
    #LOCALIZE VARIABLES
    options <- RCV.CONTRAST_OPTIONS()
    chosen <- input$contrast
    dataset <- RCV.DDSC()
    
    #SET CONTRAST
    print(options)
    cont <- options[[chosen]]
    
    RCV.RESULTS_DDSC(as.data.frame(results(dataset, contrast = cont)))
    
  })
  
  #-----------OLD EXPERIMENT-------
  
  #UPLOAD MODAL MESSAGE
  observeEvent(input$start_old_experiment,{
    showModal(
      modalDialog(
        div( 
          h5("DEG Experiment"),
          p("Accepted File Types:  .csv   .tsv "),
          p("Format: ( gene_id, gene_name, log2FoldChange, padj, .samples... ) "),
          p("- The gene id column should contain a unique identifier for every row"),
          p("- Sample counts may be included and are recommended in reviewing an experiment for all features"),
          p("- Samples should each have their own column with the sample name proceeding a '.' (standard dot/period)"),
          fileInput("deg_file", "DEG Experiment"),
          
          h5("Sample Conditions Table"),
          p("Accepted File Types:  .csv   .tsv "),
          p("Format: (sample, conditions...) "),
          p("Sample column should include sample names that are found in your experiment data"),
          fileInput("deg_meta_2_file", "Sample Conditions")
        ),footer = actionButton("read_deg_and_meta2", "Finish"), easyClose = TRUE
      )
    )
  })
  
  #READING DEG AND META DATA
  observeEvent(input$read_deg_and_meta2, {
    
    removeModal()
    
    #               REQUIRE BOTH IMPORTS
    #_____________________________________________________
    if( is.null(input$deg_file) || is.null(input$deg_meta_2_file)){
      showErrorModal('A file was not uploaded properly try again')
      return()
    }
    
    
    #               READ META DATA FILE
    #____________________________________________________
    meta_data <- readMetaData(input$deg_meta_2_file$datapath)
    if(is.null(meta_data)){ return() } 
    
    
    #               READ DEG FILE
    # (results, gene names, normalized counts, vst counts, vst object)
    #____________________________________________________
    df_list <- readDEGFile(input$deg_file$datapath, meta_data)
    if(is.null(df_list)){ return() }
    
    #               SETTING GLOBAL VARIABLES
    #_____________________________________________________
    RCV.RESULTS_DDSC(df_list[[1]])
    
    RCV.GENE_NAMES(as.data.frame(df_list[[2]]))
    
    RCV.NORMALIZED_COUNTS(df_list[[3]])
    
    RCV.RAW_COUNTS(df_list[[3]])
    
    RCV.VST_COUNTS(df_list[[4]])
    
    RCV.VST_OBJ(df_list[[5]])
    
    RCV.META_DATA(meta_data)
    
    #             UPDATE UI
    #____________________________________________
    shinyjs::hide('Page1_Upload_Options')
    shinyjs::hide('preRun_Data_Preview')
    shinyjs::show('pg1table2')
    shinyjs::show('pg1table1')
    shinyjs::show("TabSet1_Diffrential_Expression_Analysis")
    shinyjs::hide('about_DESeq2')
    
  })
  
#----
#Pathway Analysis Things_____________________________________________________
  
  #-----------NEW EXPERIMENT-------
  
  #Continue running with data
  observeEvent(input$Run_pathfinder,{
    
    #LOCALIZE VARIABLES
    d.results <- RCV.RESULTS_DDSC()
    norm <- RCV.NORMALIZED_COUNTS()
    meta_data <- RCV.META_DATA()
    
    #       REQUIRE DEG RESULTS
    #_________________________________
    if( is.null( d.results )){
      showErrorModal('No diffrentially expressed gene experiment has been completed')
      return()
    }
    
    #     FORK TO ALLOW COUNTS UPLOAD
    #__________________________________
    if(is.null(norm)){
      showModal(
        modalDialog(
          div(
            p("There is no count data for samples in your expirement yet!"),
            p("If you proceed without them you may miss out on some features"),
            p("(Optional)"),
            HTML('<p>Normalized count data includes (gene_name) and samples where each sample name starts with "."</p>'),
            fileInput('continue_upload_ab_file', 'Abundance Data'),
          ), footer = actionButton("continue_upload_ab", "Continue"),
          easyClose = FALSE
        )
      )
    }
    
  
    else{
      
      #       SET UP VARIABLE OBJECTS
      #___________________________________
      
      #DATA
      data <- d.results
      data$Gene_symbol <- rownames(data)
      data <- data %>% 
        select(Gene_symbol, logFC = log2FoldChange, Padj = padj) %>%
        filter(!is.na(Padj))
      
      #ABUNDANCE DATA / COUNTS
      counts <- norm
      counts$Gene_symbol <- rownames(counts)
      counts <- counts %>% select(Gene_symbol, everything())
      
      #GENE NAMES
      gene_names <- counts %>% select(Gene_symbol)
      
      #META DATA
      md <- meta_data
      
      #         RUNNING PATHFINDER
      #___________________________________
      res <- runPathfindRFunc(data)
      if(is.null(res)){ return() }

      #       SET GLOBAL VARIABLES
      #_____________________________________
      RCV.PATHWAY_RESULTS(res)
      RCV.PATHWAY_ABUNDANCE(counts)
      RCV.PATHWAY_GENES(gene_names)
      RCV.PATHWAY_META_DATA(md)
      
      #ALTER UI
      shinyjs::show("TabSet2_Pathway_Analysis")
      shinyjs::hide('pathfinder_option_buttons')
      shinyjs::hide('about_pathfinder')
    
    }
      
  })
  
  # Extra check for abundance data upload
  observeEvent(input$continue_upload_ab,{
    
    #        REQUIRE AN INPUT TO READ
    #________________________________________
    if( ! is.null(input$continue_upload_ab_file) ){
     
     #READ COUNTS
     objects <- readPathfinderCountsUpload(input$continue_upload_ab_file$datapath)
     if(is.null(objects)){ return() }
     
     #SET GLOABL VARIABLES
     RCV.PATHWAY_ABUNDANCE(objects[1])
     RCV.NORMALIZED_COUNTS(objects[1])
     RCV.PATHWAY_GENES(objects[2])
     
    }
    
    #LOCALIZE VARIABLES
    d.results <- RCV.RESULTS_DDSC()
    meta_data <- RCV.META_DATA()
    
    #       SET UP VARIABLE OBJECTS
    #___________________________________
    #DATA
    data <- results
    data$Gene_symbol <- rownames(data)
    data <- data %>% 
     select(Gene_symbol, logFC = log2FoldChange, Padj = padj) %>%
     filter(!is.na(Padj))
    
    #      RUNNNING PATHWAY ANALYSIS
    #________________________________________
    res <- runPathfindRFunc(data)
    if(is.null(res)){ return() }
    
    #      SET GLOBAL VARIABLES
    #___________________________________
    RCV.PATHWAY_RESULTS(res)
    RCV.PATHWAY_META_DATA(meta_data)
    
    shinyjs::show("TabSet2_Pathway_Analysis")
    shinyjs::hide('pathfinder_option_buttons')
    shinyjs::hide('about_pathfinder')
    
  })
  
  #-----------OLD EXPERIMENT-------
  
  #UPLOAD AN OLD EXPERIMENT
  observeEvent(input$review_pathfinder_new_data,{
    
    #LOCALIZE VARIABLES
    old_md <- RCV.META_DATA()
    
    
    #UI DISPLAY 1 IF NEED NEW METADATA
    if(is.null(old_md)){
      metadata_input <- div(
          h5("Sample Conditions Table"),
          p("Accepted File Types:  .csv   .tsv "),
          p("Format: (sample, conditions...) "),
          p("Sample column should include sample names that are found in your experiment data"),
          fileInput("path_meta_1_file", "Sample Conditions")
      )
    }
    #UI DISPLAY 2 THERES ALREADY METADATA OPTIONAL TO UPLOAD A NEW ONE
    else{
      
      output$tempTableXX <- renderDT({old_md}, rownames  = TRUE)
      
      metadata_input <- div(
          h5("Sample Conditions Table"),
          p('***Optional A valid file has previously been uploaded***'),
          p('Previous upload'),
          tableOutput('tempTableXX'),
          p("File uploaded now will override the previous"),
          p("Accepted File Types:  .csv   .tsv "),
          p("Format: (sample, conditions...) "),
          p("Sample column should include sample names that are found in your experiment data"),
          fileInput("path_meta_1_file", "Sample Conditions")
        )
    }
    
    #UI DISPLAY FOR UPLOADING THE FILES
    showModal(
      modalDialog(
        div( 
        
          p('Pathway Analysis Results'),
          HTML('<p>This is the standard output from a pathfindR clustered experiment.</p>'),
          p("Accepted File Types:  .csv   .tsv "),
          HTML('<p>Must contain all following columns:</p>'),
          HTML('<p>(ID,	Term_Description,	Fold_Enrichment, occurrence, support,	lowest_p,	highest_p,	non_Signif_Snw_Genes,	Up_regulated,	Down_regulated,	all_pathway_genes,	num_genes_in_path,	Cluster, Status)</p>'),
          fileInput('pathfinder_output_file', 'Pathway Analysis Results'),
        
          HTML('<p>Normalized Gene Counts'),
          p("Accepted File Types:  .csv   .tsv "),
          HTML('<p>Format: (gene_name, samples...)</p>'),
          p('All sample column headers should start with a "." character to be recognized'),
          fileInput('pathfinder_output_counts', 'Gene Counts'),
          
          metadata_input
          
        ), footer = actionButton("finish_uploading_old_pathfinder", 'Finish'),
        easyClose = TRUE
      )
    )
    
  })

  #COMPLETE UPLOADS AND READ THE FILES
  observeEvent(input$finish_uploading_old_pathfinder,{
    
    #LOCALIZE VARIABLES
    old_md <- RCV.META_DATA()
    
    #    CHECK IF THE NEEDED FILES WERE UPLOADED
    #________________________________________________
    if(is.null(input$pathfinder_output_file) || (is.null(input$path_meta_1_file) && is.null(old_md))){
      showErrorModal('A file was not uploaded properly try again')
      return()
    }
    
    #
    removeModal()
    #
    
    #         READ META DATA FILE
    #___________________________________
    if(!is.null(input$path_meta_1_file)){
      
      md <- readMetaData(input$path_meta_1_file$datapath)
      if(is.null(md)){ return() }
      
    }else{
      md <- old_md
    }
    
    #           READ DATA FILE
    #__________________________________
    res <- readPathfinderOutput(input$pathfinder_output_file$datapath)
    if(is.null(res)){ return() }
    
    #         CHECK / READ COUNTS
    #_________________________________
    if(!is.null(input$pathfinder_output_counts)){
      
      objects <- readPathfinderCountsUpload(input$pathfinder_output_counts$datapath)
      if(is.null(objects)){ return() }
    
      #     SET GLOBAL VARIABLES
      #______________________________
      RCV.PATHWAY_ABUNDANCE(objects[[1]])
      RCV.PATHWAY_GENES(objects[[2]])
      
    }
    
    #       SET GLOBAL VARIABLES
    #_________________________________
    RCV.PATHWAY_RESULTS(res)
    RCV.PATHWAY_META_DATA(md)
    if(is.null(RCV.PATHWAY_GENES())){
      
      data <- res %>%
        mutate(gene_list = lapply(all_pathway_genes, function(x) trimws(unlist(strsplit(x, ","))))) %>%
        select(gene_list)
      
      dataDF <- c()
      
      for(i in 1:nrow(data) ){
        dataDF <- c(dataDF, unlist(data[i,1]))
      }
      
      data <- data.frame(gene_name = dataDF)
      
      RCV.PATHWAY_GENES(data %>% unique())
      
    }
    
    #ALTER UI
    shinyjs::hide('about_pathfinder')
    shinyjs::show('TabSet2_Pathway_Analysis')
    
  })
#----
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
#~~~~~######____________Save Functionality____________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----

  #Unlink temp directory on session end
  session$onSessionEnded(function() {
    # List all files in the directory
    files <- list.files("SavedFiles", full.names = TRUE)
    
    # Remove all files in the directory
    unlink(files, recursive = TRUE, force = TRUE)
    print("unlinked temp dir")
  })
  
  #list(output, plot, class, h, w )
  #SAVE ASSISTANCE REACTIVE SAVES, Reflects currently displayed visuals
  d.PLOT.PCA <- reactiveVal(NULL) 
  d.PLOT.CORRELATION <- reactiveVal(NULL)
  d.PLOT.GENE <- reactiveVal(NULL)
  d.PLOT.VOLCANO <- reactiveVal(NULL)
  
  p.PLOT.ENRICHMENT <- reactiveVal(NULL)
  p.PLOT.LARGEMAP <- reactiveVal(NULL)
  p.PLOT.CASEENRICHMENT <- reactiveVal(NULL)
  p.PLOT.PATHWAY <- reactiveVal(NULL)
  
  #Iterator for event 
  i.EventID <- reactiveVal(0)
  #Update temp directoru with new plot
  observeEvent(input$save,{
    
    #iterate
    i.EventID(i.EventID() + 1)
    
    #find current page
    getSection <- input$App
    page <- switch(
              getSection, 
              "Diffrentially Expressed Genes" = input$TabSet1,
              "Pathway Analysis" = input$TabSet2
            )
    
    selected <- switch(
            page,
            "Principle Component Plots" =  d.PLOT.PCA(),
            "Correlation Anaylsis" = d.PLOT.CORRELATION(),
            "Gene Counts" = d.PLOT.GENE(),
            "Volcano Plot" = d.PLOT.VOLCANO(),
            
            "Term Enrichment" = p.PLOT.ENRICHMENT(),
            "Term Gene Heatmap" = p.PLOT.LARGEMAP(),
            "Case vs. Control" = p.PLOT.CASEENRICHMENT(),
            "Single Pathway vs. Samples" = p.PLOT.PATHWAY(),
            NULL
          )
    
    # Ensure a plot is selected
    if (is.null(selected)) {
      #Nothing to save
      return()
    }
    
    #Pixels per inch = 150 standard
    h.inches <- as.numeric(gsub("[^0-9]", "", selected$h)) / 72
    w.inches <- as.numeric(gsub("[^0-9]", "", selected$w)) / 50
    
    if(!is.ggplot(selected$plot)){
      toplot <- as.ggplot(selected$plot)
    }else{
      toplot <- selected$plot
    }
    
    #Formulate event ID
    eventID <- paste0("confirm_save_", i.EventID())
    
    showModal(
      modalDialog(
        textInput("save_plot_name", "File Name"),
        footer = actionButton(eventID, "Save"),
        easyClose = TRUE
      )
    )
    
    # Wait for the user to confirm
    observeEvent(input[[eventID]], {
      # Ensure the modal is closed
      removeModal()
      
      # Validate file name and save plot
      if (!is.null(input$save_plot_name) && input$save_plot_name != "") {
        
        ggsave(
          paste0(input$save_plot_name, ".pdf"), 
          path = "SavedFiles",
          plot = toplot, 
          height = h.inches, 
          width = w.inches, 
          units = "cm", 
          dpi = 230,
          limitsize = FALSE
        )
        
        showNotification("Plot saved successfully!", type = "message")
      } else {
        showNotification("Please provide a valid file name", type = "error")
      }
    })
    
  })
  
  
  output$downloadZip <- downloadHandler(
    filename = function() {
      "test.zip"
    },
    content = function(file) {
      
      files <- list.files("SavedFiles", full.names = TRUE)
      print(files)
      
      zip(zipfile = file, files = files)
      
    },
    contentType = "application/zip"
  )
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
#~~~~~######____________Output $ Objects______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----
  
  #Help Page
  #_______________________Interactive Flow Map_____________________#----
  
    output$interactive_flowchart <- renderPlotly({
      
      img_width = 1600 
      img_height = 900 
      scale_factor = 1 
      
      
      # Add invisible scatter trace. 
      # This trace is added to help the autoresize logic work. 
      fig <- plot_ly(width=img_width * scale_factor, 
                     height=img_height * scale_factor 
      ) %>% 
        add_trace( x= c(0, img_width * scale_factor), 
                   y= c(0, img_height * scale_factor), 
                   type = 'scatter',  mode = 'markers', alpha = 0) 
      
      # Configure axes 
      xconfig <- list( 
        title = "", 
        zeroline = FALSE, 
        showline = FALSE, 
        showticklabels = FALSE, 
        showgrid = FALSE, 
        range = c(0, img_width * scale_factor) 
      ) 
      
      yconfig <- list( 
        title = "", 
        zeroline = FALSE, 
        showline = FALSE, 
        showticklabels = FALSE, 
        showgrid = FALSE, 
        range = c(0, img_height * scale_factor), 
        scaleanchor="x" 
      ) 
      
      fig <- fig %>% layout(xaxis = xconfig, yaxis = yconfig) 
      
      # Add image 
      
      fig <- fig %>% layout( 
        images = list(  
          list(  
            source =  "Workflow_Cropped.PNG",  
            x=0, 
            sizex=img_width * scale_factor, 
            y=img_height * scale_factor, 
            sizey=img_height * scale_factor, 
            xref="x", 
            yref="y", 
            opacity=1.0, 
            layer="below", 
            sizing="stretch" 
          )  
        )) 
      
      # Configure other layout 
      
      m = list(r=0, l=0, b=0, t=0) 
      fig <- fig %>% layout(margin = m) %>%
        layout(plot_bgcolor='#e5ecf6',  
               xaxis = list(  
                 zerolinecolor = '#ffff',  
                 zerolinewidth = 2,  
                 gridcolor = 'ffff'),  
               yaxis = list(  
                 zerolinecolor = '#ffff',  
                 zerolinewidth = 2,  
                 gridcolor = 'ffff')  
        )
      fig
    })
  #_______________________File Help Page_____________________#----
    output$file_usage_documentation <- renderUI({
      includeHTML('www/FileUsageDocumentation.HTML')
    })
  
  #----

  #Page1 pretabSET----
    output$about_DESeq2 <- renderUI({
      includeHTML('www/DESeq2_Intro.HTML')
    })
    output$preRun_preview1 <- renderUI({
      #SET LOCAL DF
      df <- RCV.FILTERED_COUNTS()
      
      #SET DEFAULT CASE FOR NO UPLOAD
      if(is.null(df)){
        return()
      }
      
      #transfer gene names
      df$gene_name <- rownames(df)
      
      #Organize columns
      df <- df %>% select(gene_name, everything())
      
      table <- renderDT({df}, rownames = FALSE, options = list(pageLength=7, scrollX = TRUE))
      
      div(table)
    })
    output$preRun_preview2 <- renderUI({
      meta_data <- RCV.META_DATA()
      
      if(!is.null(meta_data)){ 
        
        renderDT({meta_data}, rownames  = TRUE) 
        
      }else{
        
        return()
        
      }
    })
    output$contrast_selection_ui <- renderUI({
      
      #LOCALIZE VARIABLES
      options <- RCV.CONTRAST_OPTIONS()
      
      if(is.null(options)){
        return()
      }
      
      div(selectInput("contrast", "Select Contrast", choices = names(options)))
      
    })
  #----
  
  #TabSet 1
  #______________________________Diff 1________________________#----
  
    # Raw Counts data table 
    output$raw_counts_PreviewTable <- renderUI({
        
        #SET LOCAL DF
        df <- RCV.FILTERED_COUNTS()
        
        #SET DEFAULT CASE FOR NO UPLOAD
        if(is.null(df)){
          return()
        }
  
        #transfer gene names
        df$gene_name <- rownames(df)
        
        #Organize columns
        df <- df %>% select(gene_name, everything())
        
        table <- renderDT({df}, rownames = FALSE, options = list(pageLength=7, scrollX = TRUE))

        div(table)

    })
    
    # render table for metadata
    output$sample_conditions_PreviewTable <- renderUI({
      
      meta_data <- RCV.META_DATA()
      
      if(!is.null(meta_data)){ 
        
        renderDT({meta_data}, rownames  = TRUE) 
        
      }else{
        
       return()
      
      }
    })
    
    
    # A data table showing normalized counts
    output$normalized_counts_PreviewTable <- renderUI({
      
      #LOCALIZE
      norm <- RCV.NORMALIZED_COUNTS()
      
      if(is.null(norm)){
        return(span("No Data Yet",style="color: red;"))
      }
      
      #ROUND AND TRANSFER NAMES
      norm <- norm %>% round(digits = 4)
      norm$gene_name <- rownames(norm)
      
      #PULL PVALUES FROM RESULTS
      res <- RCV.RESULTS_DDSC()
      res$gene_name <- rownames(res)
      pvals <- res %>% select(padj, gene_name)
      
      #MERGE PVALUES TO NORMALIZED COUNTS
      normalized <- left_join(norm, pvals)

      #ORDER BY PVAL
      normalized <- normalized %>% arrange(padj)
      
      #REORDER COLS
      normalized <- normalized %>% select(gene_name, padj, everything())

      renderDT({normalized}, rownames=FALSE, options = list(pageLength=7, scrollX = TRUE))
      
    })
    

  #______________________________Diff 2________________________#----
    
    # 3 tables feature genes that fit in expression categories
    output$expression_tables <- renderUI({

      #LOCALIZE VARIABLES
      results <- RCV.RESULTS_DDSC()
      
      #ASSERT NEEDED OBJECTS
      if( is.null(results) ){
        return(renderUI({span("No Data Yet",style="color: red;")}))
      }
      
      
      #Splits the DEG list into groups by fold change 
      func_return <- splitByExpr(results, input$cutOffs )
      up <- as.data.frame(func_return[1]) #up regulated
      down <- as.data.frame(func_return[2]) #downregulated
      noR <- as.data.frame(func_return[3]) #Inbetween (no regulation)

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
      
      #LOCALIZE VARIABLES
      results <- RCV.RESULTS_DDSC()  
      
      #ASSERT NEEDED THINGS
      if(is.null(results)){
          return()
      }
        
      plot <- renderPlot({
        
        p <- as.numeric(input$pvaluePg2)
        
        data <- results %>% 
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
    
  #_____________PCA______________Diff 4________________________#----
    
    # The principle component plot
    output$principle_component_plots_ui <- renderUI({
      
      #LOCALIZE VARIABLES
      vst <- RCV.VST_OBJ()
      
      if(is.null(vst)){return()}
      if(length(input$pca_cond) == 0){return()}
      
      pca_plot <- principlePlot(vst, 
                                input$pca_cond, 
                                input$pca_nrow, 
                                input$title_pca_plot, 
                                input$subtitle_pca_plot, 
                                input$caption_pca_plot)
      
      #RENDER
      output$principle_components <- renderPlot({pca_plot})
      
      #CHECK BEFORE SHOW
      if(is.null(vst) || length(input$pca_cond) == 0){
        return(span("No Data Yet, Or No Selected Factors",style="color: red;"))
      }
      
      #CREATE OUTPUT
      p.out <- plotOutput('principle_components', height= "700px", width = "100%")
      
      #SAVE
      d.PLOT.PCA(
        list(
          output = p.out, 
          plot = pca_plot, 
          class = class(pca_plot), 
          h = "700px",
          w = "800px"
        )
      )
      
      #PRINT
      p.out
      
    })
    
    # Widget for changing genes used in plot
    output$pca_n_counts <- renderUI({
      
      #LOCALIZE VARIABLES
      counts <- RCV.VST_COUNTS()
      
      #Sets the number of genes to 500, or the maximum
      if(is.null(counts)){
        
        m <- 500
        start <- 500
        
      }else{
        
        m <- nrow(counts)
        
        if(nrow(counts) <= 500){
          start <- m
        }else{
          start <- 500
        }
        
      }
      
      div(
        p("This number determines how many genes are included in the PCA ranked by variance."),
        p("If the default 500 is selected the top 500 genes with the most variance will be used."),
        p(paste("Max: ", m)),
        p("Min: 2"),
        numericInput(
          'pca_nrow', 
          '', 
          value = start, 
          min = 2, 
          max = m
        )
      )
    })
    
    #Widget for changing the metadata columns
    output$pca_metadata <- renderUI({
      
      #LOCALIZE VARIABLES
      meta_data <- RCV.META_DATA()
      
      div(
        checkboxGroupInput('pca_cond', 'Meta data conditions to view', 
                        choices=colnames(meta_data),
                        selected=head(colnames(meta_data),3),
                        inline=FALSE)
      )
    })
    
  #_____________Correlation______Diff 5________________________#----
    
    # The heatmap for correlation analysis
    output$heatmap_plot_ui <- renderUI({
      
      #LOCALIZE VARIABLES
      counts <- RCV.VST_COUNTS()
      meta_data <- RCV.META_DATA()
      corr_anns <- input$corr_annotations
      
      #ASSERT NEEDED VARIABLES
      if(is.null(counts)){return()}
        
      #         CREATE CORRELATION MATRIX
      # ___________________________________________
      vst_cor <- cor(counts)
      
      generate_palette <- function(factor_levels, base_colors = c("cornsilk3", "orchid")) {
        # Check the number of unique factor levels
        unique_levels <- unique(factor_levels)
        num_levels <- length(unique_levels)
        
        # Generate the palette using colorRampPalette
        palette <- colorRampPalette(base_colors)(num_levels)
        
        # Return a named vector mapping levels to colors
        setNames(palette, unique_levels)
      }
      
      ann_colors <- list()
      
      for(i in corr_anns){
        if(nlevels(meta_data[[i]]) > 2){
          add <- setNames(list(generate_palette(meta_data[[i]], base_colors = c("cadetblue", "white", "indianred"))), i)
          ann_colors <- c(ann_colors, add)
        }else{
          add <- setNames(list(generate_palette(meta_data[[i]])), i)
          ann_colors <- c(ann_colors, add)
        }
      }
      
      toAnnotate <- meta_data %>% select(all_of(head(corr_anns, 5))) 
      ann_colors <- head(ann_colors, 5)
      
      p <- pheatmap(
              vst_cor, 
              angle_col = "45", 
              main = input$title_heatmap_plot, 
              annotation_row = if (nrow(toAnnotate) > 0 && ncol(toAnnotate) > 0) 
                toAnnotate else NA,
              annotation_colors = if (nrow(toAnnotate) > 0 && ncol(toAnnotate) > 0) 
                ann_colors else NA,
              fontsize = 20,
              cellwidth = 650 / nrow(vst_cor),
              cellheight = 650 / nrow(vst_cor)
            )
      
      output$heatmap_plot_1 <- renderPlot({(p)})
      
      #RENDER UI
      p.out <- plotOutput('heatmap_plot_1', height = "1200px", width = "100%")
      
      #SAVE
      d.PLOT.CORRELATION(
        list(
          output = p.out, 
          plot = p, 
          class = class(p), 
          h = "1200px", 
          w = "1200px"
        )
      )
      
      #print
      p.out
      
    })
  
    output$choose_annotations_ui <- renderUI({
      
      #LOcalize variables
      meta_columns <- colnames(RCV.META_DATA()) 

      checkboxGroupInput(
        'corr_annotations',
        'Max 5 Shown',
        choices = meta_columns,
        selected = meta_columns[[1]]
      )
    })
    
  #_____________GPlot____________Diff 6________________________#----
    
    # A ui bar on the left displaying some important genes 
    output$leftbar_gene_plots <- renderUI({
      
      #LOCALIZE VARIABLES
      results <- RCV.RESULTS_DDSC()
      norm <- RCV.NORMALIZED_COUNTS()
      meta_data <- RCV.META_DATA()
      
      #ASSERT NEEDED VARIABLES
      if(is.null(results) || is.null(norm)){ return() }

      p <- as.numeric(input$pvaluePg6)
      
      #SETUP COUNTS
      counts <- round(norm) 
      counts$gene_name <- rownames(counts)
      
      #SETUP RESULTS
      results$gene_name <- rownames(results)
      results <- results %>% select(gene_name, padj, log2FoldChange)
      
      #MERGE COUNTS AND RESULTS
      data <- left_join(results, counts)
      
      #FILTER PVALUES
      data <- data %>% filter(padj <= p)
    
      tp <- data %>% arrange(padj) %>% head(1)
      lf <- data %>% arrange(log2FoldChange) %>% head(1)
      hf <- data %>% arrange(desc(log2FoldChange)) %>% head(1)
      
      top_pval <- gplop(tp, meta_data, colnames(meta_data))
      
      low_fold <- gplop(lf, meta_data, colnames(meta_data))
      
      high_fold <- gplop(hf, meta_data, colnames(meta_data))
      
      div(
        
        h2("Lowest adjusted pval"),
        
        renderPlot({top_pval}),
        
        h2("Highest Fold Change"),
        
        renderPlot({high_fold}),
        
        h2("Lowest Fold Change"),
        
        renderPlot({low_fold})
        
      )

    })
    
    # A ui bar on the right displaying a datatable of all the genes
    output$gene_name_list <- renderUI({
      
      #LOCALIZE VARIABLE
      names <- RCV.GENE_NAMES()
      
      if(is.null(names)){
        span("No Data Yet",style="color: red;")
        return()
      }

      names <- datatable(names,
                        rownames=FALSE, 
                        option=list(pageLength=25))
      div(
        h3("Gene List"),
        names
      )

    })
    
    # A ui for displaying the search method 
    output$gene_count_search <- renderUI({
      
      #LOCALIZE VARIABLES
      results <- RCV.RESULTS_DDSC()
      norm <- RCV.NORMALIZED_COUNTS()
      meta_data <- RCV.META_DATA()
      
      if(is.null(results) || is.null(norm)){
        return(span("No Data Yet",style="color: red;"))
      }
      
      #SETUP COUNTS
      counts <- data.frame(round(norm)) 
      counts$gene_name <- rownames(counts) %>% tolower()
      
      #SETUP RESULTS
      results$gene_name <- rownames(results) %>% tolower()
      results <- results %>% select(gene_name, padj, log2FoldChange)
      
      #MERGE COUNTS AND RESULTS
      data <- left_join(results, counts)
      
      div(
        
        textInput("gplop_searched",
                  "Search and plot a gene", 
                  NULL),
        
        checkboxGroupInput('gene_count_plot_search_cond', 
                            'Meta data conditions to view, first 4 are considered with mapping order (x, color, shape, size)', 
                            choices=colnames(meta_data),
                            selected=head(colnames(meta_data),4),
                            inline=TRUE),
        renderUI({
          
          if(!is.null(input$gplop_searched)){
            
            #EDIT SEARCH
            g <- input$gplop_searched %>% tolower()
            
            #PULL SEARCH
            res <- data %>% filter(gene_name == g)
            print(res)
            
            #ASSERT SINGLE ROW
            if(is.null(res) || ncol(res) == 0 || nrow(res) == 0){
              return()
            }
            if(nrow(res) > 1 ){
              res <- res[1,]
            }
            
            plot <- gplop(res, meta_data, input$gene_count_plot_search_cond)
            output$sGene <- renderPlot({plot})
            p.out <- plotOutput('sGene')
            
            #saving
            d.PLOT.GENE(
              list(
                output = p.out, 
                plot = plot, 
                class = class(plot), 
                h = "800px", 
                w = "800px"
              )
            )

            #print
            p.out
            
          }else{
            p("Nothing Here Yet")
          }
          
        })

      )
    })
    
  #_____________Volcano__________Diff 7________________________#----
  
    # The plot for the volcano plot page
    output$volcano_plot_ui <- renderUI({
      
      #LOCALIZE VARIABLES
      results <- RCV.RESULTS_DDSC()
      
      #             REQUIRE RESULTS
      #____________________________________
      if(is.null(results)){ return() }
      
      #         HANDLE GENE SEARCHING
      #____________________________________
      if(!is.null(input$volc_search)){
        search <- strsplit(input$volc_search, ",")[[1]]
        search <- tolower(trimws(search))
      }else{
        search <- NULL
      }
      
      #         MODIFY RESULTS DF
      #________________________________________
      results$gene_name <- rownames(results)
      results <- results %>% select(gene_name,log2FoldChange, padj)
      #Pvalue adjustment
      results$padj <- as.numeric(results$padj) 
      
      #         LOCALIZE PLOT PARAMETERS
      #_______________________________________
      #down regulation cut off
      lowcut <- input$volcano_cutoffs[1]
      #down regulation cut off
      highcut <- input$volcano_cutoffs[2]
      #Pvalue cuttoff
      pcut <- input$pvaluePg7 
      #scatter population
      pop_score <- input$volcano_pop 
      #Label density
      lab_score <- input$volcano_lab_density 
      
      #labs
      title <- input$title_volc_plot 
      subtitle <- input$subtitle_volc_plot 
      caption <- input$caption_volc_plot
      
      #PLOT
      plot <- erupt( results, 
                     lowcut, 
                     highcut, 
                     pcut, 
                     pop_score, 
                     lab_score, 
                     search, 
                     title, 
                     subtitle, 
                     caption)
      
      output$volcano_plot <- renderPlot({plot})
      p.out <- plotOutput('volcano_plot', height= "900px", width = "100%")
      
      #save
      d.PLOT.VOLCANO(
        list(
          output = p.out, 
          plot = plot, 
          class = class(plot), 
          h = "900px", 
          w = "1000px"
        )
      )
      
      #print
      p.out
      
    })

  #----
    
  #Page2 pretabSET----
    output$about_pathfinder <- renderUI({
      includeHTML('www/PathfindR_Intro.HTML')
    })
  #----  
  
  #Tabset 2
  #______________________________Path 1________________________#----
    
    output$pathfinderPreview <- renderUI({
      
      #LOCALIZE VARIABLES
      res <- RCV.PATHWAY_RESULTS()
      
      if( is.null(res) ){ return() }
      
      #ORGANIZE
      data <- res %>% 
        select(-ID, -Up_regulated,	-Down_regulated, -all_pathway_genes)
      
      renderDT({data}, rownames = FALSE)
      
    })
  
  #______________Enrichment______Path 2________________________#----
    
    output$enrichmentUI <- renderUI({
      
      #LOCALIZE VARIABLES
      res <- RCV.PATHWAY_RESULTS()
      
      #         REQUIRE RESULTS
      #_______________________________
      if(is.null(res)){
        return(span("No Data Yet",style="color: red;"))
      }

      #       HANDLE SEARCHES
      #_____________________________
      
      #GENES
      if( input$enrichment_genes != c("")){
        
        search_genes <- strsplit(input$enrichment_genes, ",")[[1]]
        search_genes <- trimws(search_genes)
      
      } else { search_genes <- NULL }
      
      #PATHWAYS
      if( input$enrichment_paths != c("")){
        
        search_paths <- strsplit(input$enrichment_paths, ",")[[1]]
        search_paths <- trimws(search_paths)
        
      } else { search_paths <- NULL }
      
      #CLUSTERS
      if( input$enrichment_clusters != c("")){
        
        search_clusters <- strsplit(input$enrichment_clusters, ",")[[1]]
        search_clusters <- as.numeric(trimws(search_clusters))
        
      }
      #DEFAULT DISPLAY CLUSTERS
      else if( is.null(search_genes) && is.null(search_paths) ){
        
        search_clusters <- seq(1, input$clusters_on_enrichment)
        
      }
      else{ search_clusters <- NULL }
      
      #           PLOT
      #_____________________________
      chart <- enricher(res, 
                       search_clusters,
                       search_genes,
                       search_paths) 
      
      #         EDIT PLOT
      #_____________________________
      chart <- chart +
        theme(axis.text.y = element_text(size = 17),
              axis.text.x = element_text(size = 15),
              legend.text = element_text(size = 15),
              legend.title = element_text(size = 15),
              strip.text = element_text(size = 15))
      
      #         DYNAMIC SIZING
      #________________________________
      hits <- nrow(chart$data)
      h <- paste0(450 + hits * 23, "px")
      w <- "1000px"
      
      
      #           RENDER PLOT
      #________________________________
      output$enrichmentChart <- renderPlot({ chart })
      
      
      #           CREATE UI
      #________________________________
      p.out <- plotOutput('enrichmentChart', height = h, width = w)
      
      #save
      p.PLOT.ENRICHMENT(
        list(
          output = p.out, 
          plot = chart, 
          class = class(chart), 
          h = h, 
          w = w
        )
      )
      
      #print
      p.out
      
    })
    
    output$enrichment_clusters_shown <- renderUI({

      #LOCALIZE VARIABLES
      res <- RCV.PATHWAY_RESULTS()
      
      #           REQUIRE RESULTS
      # ___________________________________
      if( is.null(res) ){ return() }
      
      #MAKE SLIDER UI
      if( input$enrichment_genes == "" && 
          input$enrichment_paths == "" && 
          input$enrichment_clusters == "" ){

        sliderInput('clusters_on_enrichment', 
                    "# Clusters Shown", 
                    min = 1, 
                    max = as.numeric(max(res$Cluster)),
                    value = round(1 + (0.15 * max(res$Cluster))))
      }
      
    })
    
    output$pathwaysDT9 <- renderUI({
      
      #LOCALIZE VARIABLES
      res <- RCV.PATHWAY_RESULTS()
      
      #     REQUIRE RESULTS
      #__________________________
      if( is.null(res) ){ return() }
      
      data <- res %>% 
        select(Term_Description)
      
      renderDT({data}, rownames = FALSE, options = list(pageLength = 5))
      
    })
    
    output$genes_in_paths_DT9 <- renderUI({
      
      #LOCALIZE VARIABLES
      res <- RCV.PATHWAY_RESULTS()
      
      #       REQUIRE RESULTS
      #____________________________
      if( is.null(res)){ return() }
      
      #       CREATE A LIST OF ALL GENEs
      #_____________________________________
      data <- res %>%
        mutate(gene_list = lapply(all_pathway_genes, function(x) trimws(unlist(strsplit(x, ","))))) %>%
        select(gene_list)
      
      dataDF <- c()
      
      for(i in 1:nrow(data) ){
        dataDF <- c(dataDF, unlist(data[i,1]))
      }
    
      data <- data.frame(gene = dataDF)
      
      #MAKE COLUMN FOR HOW MANY PATHS EACH GENE IS IN
      #_______________________________________________
      data <- data %>% 
        group_by(gene) %>%
        mutate(in_n_pathways = n()) %>% 
        distinct(gene, .keep_all = TRUE) %>%
        arrange(desc(in_n_pathways))
      
      #RENDER DATA TABLE
      renderDT({data}, rownames = FALSE, options = list(pageLength = 7))
      
    })
    
  #______________Heatmap_________Path 3________________________#----
  
    plot10 <- reactiveVal(NULL)
    
    #HIDE SEARCH BY PATH
    observeEvent(input$pathway_heatmap_genes,{
      if(nchar(input$pathway_heatmap_genes) > 0){
        shinyjs::hide('heatmap_paths')
      }else{
        shinyjs::show('heatmap_paths')
      }
    })
    #HIDE SEARCH BY GENE
    observeEvent(input$heatmap_paths,{
      if(nchar(input$heatmap_paths) > 0){
        shinyjs::hide('pathway_heatmap_genes')
      }else{
        shinyjs::show('pathway_heatmap_genes')
      }
    })
    
    #OPEN PLOTLY IN NEW TAB
    observeEvent(input$open_in_new_tab10,{
      
      saveWidget(plot10(), 'plot.html')
      browseURL('plot.html')
      
    })
    
    output$open_in_new_tab10 <- renderUI({
      
      #LOCALIZE VARIABLES
      res <- RCV.PATHWAY_RESULTS()
      
      #REQUIRE RESULTS
      if(is.null(res)){ return() }
      
      #BUTTON
      actionButton('open_in_new_tab10', 
                   "Open in new tab")    
    })
    
    #DATA TABLE OF PATHWAYS
    output$pathwaysDT10 <- renderUI({
      
      #LOCALIZE VARIABLES
      res <- RCV.PATHWAY_RESULTS()
      
      #     REQUIRE RESULTS
      #__________________________
      if( is.null(res) ){ return() }
      
      data <- res %>% 
        select(Term_Description)
      
      renderDT({data}, rownames = FALSE, options = list(pageLength = 5))
      
    })
    #DATA TABLE OF GENES
    output$genes_in_paths_DT10 <- renderUI({
      
      #LOCALIZE VARIABLES
      res <- RCV.PATHWAY_RESULTS()
      
      #       REQUIRE RESULTS
      #____________________________
      if( is.null(res)){ return() }
      
      #       CREATE A LIST OF ALL GENEs
      #_____________________________________
      data <- res %>%
        mutate(gene_list = lapply(all_pathway_genes, function(x) trimws(unlist(strsplit(x, ","))))) %>%
        select(gene_list)
      
      dataDF <- c()
      
      for(i in 1:nrow(data) ){
        dataDF <- c(dataDF, unlist(data[i,1]))
      }
      
      data <- data.frame(gene = dataDF)
      
      #MAKE COLUMN FOR HOW MANY PATHS EACH GENE IS IN
      #_______________________________________________
      data <- data %>% 
        group_by(gene) %>%
        mutate(in_n_pathways = n()) %>% 
        distinct(gene, .keep_all = TRUE) %>%
        arrange(desc(in_n_pathways))
      
      #RENDER DATA TABLE
      renderDT({data}, rownames = FALSE, options = list(pageLength = 7))
      
    })
    
    #PLOT UI
    output$pathway_heatmap <- renderUI({
      
      #LOCALIZE VARIABLES
      res <- RCV.PATHWAY_RESULTS()
      
      #     REQUIRE RESULTS
      #______________________________
      if(is.null(res)){ return(span("No Data Yet",style="color: red;")) }
      
      #     SET UP PARAMETERS
      #_______________________________
      #GENES
      genes <- as.character(input$pathway_heatmap_genes)
      if(nchar(genes) > 0){
        genes <- strsplit(genes, ",")[[1]]
        genes <- trimws(genes)
      } else {genes <- NULL}
      
      #PATHWAYS
      pathways <- input$heatmap_paths
      if(nchar(pathways) > 0){
        if(grepl("^[1-9][0-9]*$", pathways)){
          pathways <- as.numeric(pathways)
        }else{
          pathways <- strsplit(pathways, ",")[[1]]
          pathways <- trimws(pathways)
        }
      } else {pathways <- NULL}
      
      #DEAFULT CASE
      plot <- NULL
      
      
      #   CREATE PLOT BASED ON SEARCH
      #__________________________________
      if(!is.null(genes)){
        plot <- geneheatmap(res, genes)
      }
      
      else if(!is.null(pathways)){
        plot <- pathwayheatmap(res, pathways)
      }
      
      else{
        plot <- pathwayheatmap(res, 5)
      }
      
      #   ASSERT PLOT EXISTS
      #___________________________________
      if(is.null(plot)){
        plot10(NULL)
        return(span("No searched genes show diffrential expression",style="color: red;"))
      }
      
      #THE PLOTLY OBJECT THAT CAN BE OPENED IN BROWSER
      plotlie <- plot[[2]]
      #Edit dims
      plotlie$x$layout$width <- plotlie$x$layout$width + 200
      plotlie$x$layout$height <- plotlie$x$layout$height + 200
      plot10(plotlie)
      
      
      edits <- plot[[1]] +
        theme(
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          plot.title = element_text(size = 24, margin = margin(0,0,15,0))
        ) +
        labs(title = input$title_pathwayVgene_heatmap) + 
        scale_x_discrete(position = 'top')
      
      #CREATE PLOT OUTPUT
      output$plot_heatmap_10 <- renderPlot({
        edits
      }) 
      
      #UI OUTPUT
      p.out <- plotOutput('plot_heatmap_10', 
                 width = 'auto', #paste0(plot[[2]]$x$layout$width + 200, "px"), 
                 height = paste0(plot[[2]]$x$layout$height + 275, "px"))
      
      #save
      p.PLOT.LARGEMAP(
        list(
          output = p.out, 
          plot = edits, 
          class = class(edits), 
          h = paste0(plot[[2]]$x$layout$height + 275, "px"), 
          w = paste0(plot[[2]]$x$layout$width + 200, "px")
        )
      )
      
      #print
      p.out
      
    })
    
  #______________Case____________Path 4________________________#----
  
    output$sample_conditions_PreviewTable11 <- renderUI({
      
      #LOCALIZE VARIABLES
      meta_data <- RCV.PATHWAY_META_DATA()
      
      #REQUIRE DATA
      if(is.null(meta_data)){ 
        return(span("Upload Meta Data",style="color: red;"))
      }
      
      #RENDER TABLE
      
      renderDT({meta_data}, rownames  = TRUE) 
      

    })
    
    output$cases_select_box <- renderUI({
      
      #LOCALIZE VARIABLES
      meta_data <- RCV.PATHWAY_META_DATA()
      
      #REQUIRE DATA
      if(is.null(meta_data)){ 
        return(span("Upload Meta Data",style="color: red;"))
      }
      
      caseoptions <- rownames(meta_data)
      
      checkboxGroupInput('score_terms_cases', 'Select Case Samples', 
                         choices = caseoptions, inline = TRUE)
      
    })
    
    output$case_plot_ui <- renderUI({

      #LOCALIZE VARIABLES
      res <- RCV.PATHWAY_RESULTS()
      md <- RCV.PATHWAY_META_DATA()
      norm <- RCV.NORMALIZED_COUNTS()
      p.ab <- RCV.PATHWAY_ABUNDANCE()
      
      #     REQUIRE DATA
      #_______________________
      
      if(is.null(res) || is.null(md)){ return() }
      if(is.null(p.ab) && is.null(norm)){
        return(span("Abundance or counts data is needed for this visual",style="color: red;"))
      }
      
      #   SET UP CASE VS CONTROL PARAM
      #____________________________________
      if(is.null(input$score_terms_cases)){
        caseSamples <- NULL
      } else {
        caseSamples <- input$score_terms_cases
      }
      
      #     FIND COUNTS SOURCE
      #__________________________________
      #PRIORITIZE ABUNDANCE UPLOAD
      if(!is.null(p.ab)){
        ab <- p.ab
      }
      #RESORT TO NORMALIZED COUNTS
      else if (!is.null(norm)){
        ab <- norm
        ab$Gene_symbol <- rownames(ab)
        ab <- ab %>% select(Gene_symbol, everything())
      }
      #FINALLY GIVE UP
      else{ return() }
      
      print(head(res, n=1))
      print(head(ab, n=1))
      #       CREATE THE PLOT
      #____________________________
      plotdata <- score_pathway_terms(
        data = res,
        abundance = ab,
        cases = caseSamples,
        repOnly = input$repOnly, 
        pathways = NULL)
      
      #GATHER DIMENSIONS
      h <- plotdata[[2]]
      w <- plotdata[[3]]
      
      #RENDER PLOT
      output$case_map_plot <- renderPlot({
        plotdata[[1]] + labs(title = input$case_plot_title) + 
          theme(plot.title = element_text(size =20, margin = margin(0,0,10,0)))
      })
      
      p.out <- plotOutput('case_map_plot', height = h, width = w)
      
      #save
      p.PLOT.CASEENRICHMENT(
        list(
          output = p.out, 
          plot = plotdata[[1]], 
          class = class(plotdata[[1]]), 
          h = h, 
          w = w
        )
      )
      
      #print
      p.out
      
      
    })

    
  #______________Pathway_________Path 5________________________#----

    output$Single_pathway_plot_ui <- renderUI({
      
      #LOCALIZE VARIABLES
      res <- RCV.PATHWAY_RESULTS()
      p.ab <- RCV.PATHWAY_ABUNDANCE()
      norm <- RCV.NORMALIZED_COUNTS()
      deg.res <- RCV.RESULTS_DDSC()
      annotations <- input$path_heatmap_annotations
      meta_data <- RCV.PATHWAY_META_DATA()
      
      #     REQUIRE DATA
      #__________________________
      
      if(is.null(res)){ return() }
      if(is.null(p.ab) && is.null(norm)){
        return(span("Abundance or counts data is needed for this visual",style="color: red;"))
      }
      if(!tolower(input$Single_pathway_plot_search) %in% tolower(res$Term_Description)){
        return(span("The searched pathway was not found",style="color: red;"))
      }

      #   FIND COUNTS SOURCE
      #___________________________
      #PRIORITIZE ABUNDANCE UPLOAD
      if(!is.null(p.ab)){
        
        ab <- p.ab
        ab$gene_name <- tolower(ab$Gene_symbol)
        ab <- ab %>% select(-Gene_symbol)
        
      }
      #RESORT TO NORMALIZED COUNTS
      else if(!is.null(norm)){
        
       ab <- norm
       ab$gene_name <- tolower(rownames(ab))
                
      }
      #FINALLY GIVE UP
      else{ return() }
      
      #   SET UP DEG SOURCE
      #__________________________
      #PRIORITIZE DEG EXPERIMENT
      if(!is.null(deg.res)){
        
        d <- deg.res
        d$gene_name <- tolower(rownames(d))
        d_arg <- d %>% left_join(ab)
        
      }
      #RESORT TO NULL Padj
      else{
        
        ab$padj <- 0
        d_arg <- ab
        
      }
      
      #   CREATE THE PLOT DATA
      #_________________________________
      data <- single_pathway_heatmap(
        input$Single_pathway_plot_search, 
        d_arg, 
        res, 
        genes_listed = input$Single_pathway_plot_num_points)
      
      #     ASSERT THE FUNCTION WORKED
      #________________________________
      if(!is.matrix(data)){ return() }
      
      #   SET UP PLOT DIMENSIONS
      #_________________________________
      w <- max(ncol(data) * 100, 900)
      h <- max(nrow(data) * 20, 700)
      w <- paste0(w, "px")
      h <- paste0(h, "px")
      
      #   SET UP ANNOTATIONS
      #_____________________________
      generate_palette <- function(factor_levels, base_colors = c("cornsilk3", "orchid")) {
        # Check the number of unique factor levels
        unique_levels <- unique(factor_levels)
        num_levels <- length(unique_levels)
        
        # Generate the palette using colorRampPalette
        palette <- colorRampPalette(base_colors)(num_levels)
        
        # Return a named vector mapping levels to colors
        setNames(palette, unique_levels)
      }
      
      ann_colors <- list()
      for(i in annotations){
        if(nlevels(meta_data[[i]]) > 2){
          add <- setNames(list(generate_palette(meta_data[[i]], base_colors = c("cadetblue", "white", "indianred"))), paste0(i, '   '))
          ann_colors <- c(ann_colors, add)
        }else{
          add <- setNames(list(generate_palette(meta_data[[i]])), paste0(i, '   '))
          ann_colors <- c(ann_colors, add)
        }
      }
      
      toAnnotate <- meta_data %>% select(all_of(head(annotations, 5)))
      if(ncol(toAnnotate) > 0){
        colnames(toAnnotate) <- sapply(colnames(toAnnotate), function(x) paste0(x, '   '))
      }
      ann_colors <- head(ann_colors, 5)
      
      #   PLOT THE DATA
      #______________________________
      plot <- pheatmap(data, 
               scale = 'row', 
               fontsize = 15, 
               angle_col = "45", 
               legend_breaks = c(-2, -1, 0, 1, 2), 
               main = input$Single_pathway_plot_search,
               annotation_col = if (nrow(toAnnotate) > 0 && ncol(toAnnotate) > 0) 
                 toAnnotate else NA,
               annotation_colors = if (nrow(toAnnotate) > 0 && ncol(toAnnotate) > 0) 
                 ann_colors else NA)
      
      #RENDER PLOT
      output$Single_pathway_plot <- renderPlot({plot})
      
      p.out <- plotOutput('Single_pathway_plot', height = h, width = 'auto')
      
      #save
      p.PLOT.PATHWAY(
        list(
          output = p.out, 
          plot = plot, 
          class = class(plot), 
          h = h, 
          w = w
        )
      )
      
      #print
      p.out
      
    })
    
    output$choose_annotations_ui_p <- renderUI({
      
      #LOcalize variables
      meta_columns <- colnames(RCV.PATHWAY_META_DATA()) 
      
      checkboxGroupInput(
        'path_heatmap_annotations',
        'Max 5 Shown',
        choices = meta_columns,
        selected = meta_columns[[1]]
      )
    })
    
#----
} #X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X

