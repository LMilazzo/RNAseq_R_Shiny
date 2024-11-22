#_________________________________Shiny Server__________________________________________
server <- function(input, output) {

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
#~~~~~######_______________WORKFLOW________________######~~~~~#
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
    hide('Page1_Upload_Options')
    show('run_DESeq2')
    show('pg1table2')
    show('pg1table1')
    hide('about_DESeq2')
    
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
    
    # UI ADDING A INTERACTION SELECTION ELEMNT OF EACH INTERACTION NEEDED
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
            design_part <- c(design_parts, interaction_term)
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
      if(RCV.INT_COUNTER > 0){
        RCV.INT_COUNTER(RCV.INT_COUNTER - 1)
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
      style = "background-color: #4CAF50; color: #4CAF50; border-color: #4CAF50;")), 
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
      
      #         CREATE THE DESeq DATA SET   
      #_____________________________________________________________
      ddsc <- DESeqDataSetFromMatrix(countData = round(countdata),
                                     colData = metadata,
                                     design = des_formula)
      
      showModal(modalDialog("Running DESeq...", footer = NULL))
      
      YourDESeqDataSet <- DESeq(ddsc)
      
      #                     SELECT CONTRASTS
      #_____________________________________________________________
      contrast_list <- list()
      #Make options
      for (col_name in colnames(colData(ddsc))) {
        levels <- levels(colData(ddsc)[[col_name]])
        combinations <- combn(levels, 2, simplify = FALSE)
        for (combo in combinations) {
          combo_name1 <- paste(combo, collapse = " vs. ")
          combo_name2 <- paste(rev(combo), collapse = " vs. ")
          contrast_list[[combo_name1]] <- c(col_name, combo)
          contrast_list[[combo_name2]] <- c(col_name, rev(combo))
        }
      }
      
      showModal(modalDialog(
        div(selectInput("contrast", "Select Contrast", choices = names(contrast_list))),
        footer = actionButton('finish', "Complete")
      ))
      
      #         FINISH THE PROCESS
      #___________________________________________________________
      observeEvent(input$finish,{
        
        cont <- contrast_list[[input$contrast]]
        
        #         SET GLOBAL VARIABLES FOR THE EXPERIMENT RESULTS
        #___________________________________________________________________
        RCV.DDSC(YourDESeqDataSet)
                                  
        RCV.RESULTS_DDSC( as.data.frame(results(YourDESeqDataSet, contrast = cont)))
        
        RCV.VST_OBJ( vst(YourDESeqDataSet, blind = TRUE, nsub = 50) )
      
        RCV.VST_COUNTS( as.data.frame(assay(RCV.VST_OBJ())) )
        
        RCV.NORMALIZED_COUNTS( as.data.frame(counts(YourDESeqDataSet, normalized = TRUE)) )
        
        #           EXIT MESSAGE
        #_____________________________
        showModal(modalDialog("DESeq complete with no fatal errors", easyClose = TRUE, footer = NULL))
        
        #           UPDATE UI
        # _______________________________
        hide("run_DESeq2")
        hide('preRun_Data_Preview')
        show("TabSet1_Diffrential_Expression_Analysis")
        
      })
        
    }, error = function(e) {
      
      showErrorModal(paste("Error in DESeq2 process:", e$message))
      
    })
    
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
    hide('Page1_Upload_Options')
    hide('preRun_Data_Preview')
    show('pg1table2')
    show('pg1table1')
    show("TabSet1_Diffrential_Expression_Analysis")
    hide('about_DESeq2')
    
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
      show("TabSet2_Pathway_Analysis")
      hide('pathfinder_option_buttons')
      hide('about_pathfinder')
    
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
    
    show("TabSet2_Pathway_Analysis")
    hide('pathfinder_option_buttons')
    hide('about_pathfinder')
    
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
      
      output$tempTableXX <- renderTable({old_md}, rownames  = TRUE)
      
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
        
          HTML('<p>(csv or tsv) containing (ID,	Term_Description,	Fold_Enrichment, occurrence, support,	lowest_p,	highest_p,	non_Signif_Snw_Genes,	Up_regulated,	Down_regulated,	all_pathway_genes,	num_genes_in_path,	Cluster, Status,)</p>'),
          HTML('<br>This is the standard output from a pathfindR clustered experiment.'),
          fileInput('pathfinder_output_file', 'Pathway Information'),
        
          HTML('<p>Normalized count data includes (gene_name) and samples where each sample name starts with "."</p>'),
          fileInput('pathfinder_output_counts', 'Abundance Data'),
          
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
    hide('about_pathfinder')
    show('TabSet2_Pathway_Analysis')
    
  })
#----
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######____________Output $ Objects______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----
  
  #Help Tab
  #_______________________Interactive Flow Map_____________________#----
  
    output$interactive_image <- renderPlotly({
      
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
  
  
  
  
  #Page1 pretab----
    output$about_DESeq2 <- renderUI({
      includeHTML('www/DESeq2_Intro.HTML')
    })
    output$preRun_preview1 <- renderUI({
      #SET LOCAL DF
      df <- RCV.FILTERED_COUNTS()
      
      #SET DEFAULT CASE FOR NO UPLOAD
      if(is.null(df)){
        return(span("No Data Yet",style="color: red;"))
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
        
        renderTable({meta_data}, rownames  = TRUE) 
        
      }else{
        
        span("Upload a file",style="color: red;")
        
      }
    })
  #----
  
  #TabSet 1
  #______________________________Diff 1_______DONE_________________#----
  
    # Raw Counts data table 
    output$raw_counts_PreviewTable <- renderUI({
        
        #SET LOCAL DF
        df <- RCV.FILTERED_COUNTS()
        
        #SET DEFAULT CASE FOR NO UPLOAD
        if(is.null(df)){
          return(span("No Data Yet",style="color: red;"))
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
        
        renderTable({meta_data}, rownames  = TRUE) 
        
      }else{
        
        span("Upload a file",style="color: red;")
      
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
    

  #______________________________Diff 2_______DONE_________________#----
    
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
    
  #______________________________Diff 4_______DONE_________________#----
    
    # The principle component plot
    output$principle_components <- renderPlot({
      
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
      
      pca_plot
      
    })
    
    # UI holding the plot
    output$principle_component_plots_ui <- renderUI({
      
      #LOCALIZE VARiABLES
      vst <- RCV.VST_OBJ()
      
      if(is.null(vst) || length(input$pca_cond) == 0){
        return(span("No Data Yet, Or No Selected Factors",style="color: red;"))
      }

      plotOutput('principle_components', height= "700px", width = "100%")
    
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
        numericInput('pca_nrow', '', 
                    value=start, min=2, max=m)
      )
    })
    
    #Widget for changing the metadata columns
    output$pca_metadata <- renderUI({
      
      #LOCALIZE VARIABLES
      meta_data <- RCV.META_DATA()
      
      div(
        checkboxGroupInput('pca_cond', 'Meta data conditions to view', 
                        choices=colnames(meta_data),
                        selected=colnames(meta_data),
                        inline=FALSE)
      )
    })
    
  #______________________________Diff 5_______DONE_________________#----
    
    # The heatmap for correlation analysis
    output$heatmap_plot1 <- renderPlot({
      
      #LOCALIZE VARIABLES
      counts <- RCV.VST_COUNTS()
      meta_data <- RCV.META_DATA()
      
      #ASSERT NEEDED VARIABLES
      if(is.null(counts)){return()}
        
      #         CREATE CORRELATION MATRIX
      # ___________________________________________
      vst_cor <- cor(counts)
      
      #   CREATE LIST OF CONDITIONS TO ANNOTATE
      #____________________________________________
      to_annotate <- intersect(colnames(meta_data), input$heatmap_cond)
      
      #CREATE ANNOTATION ELEMENTS
      annotations <- lapply(to_annotate, function(an) {
        
        #CREATE RAMP PALLETS
        selected_color <- colorRampPalette(c(input[[paste0("color1_", an)]], 
                                              input[[paste0("color2_", an)]])
                                          )(nlevels(factor(meta_data[[an]])))
        
        selected_color <- setNames(selected_color, levels(factor(meta_data[[an]])))
        
        #CREATE HEATMAP ANNOTATION ELEMENT
        HeatmapAnnotation(df = meta_data[, an, drop = FALSE], 
                          which = "col", 
                          annotation_name_side = "left",
                          annotation_name_gp = gpar(fontsize = 17, color='black'),
                          gap = unit(1, 'cm'),
                          height = unit(1, 'cm'),
                          col = setNames(list(selected_color), an),
                          show_legend = FALSE
                          )
      })
      
      #COMBINE ELEMENTS
      combined_annotations <- do.call(c, annotations)
  
      #CREATE LEGENDS FOR ANNOTATIONS
      legends <- lapply(to_annotate, function(an) {
        
        selected_color <- colorRampPalette(c(input[[paste0("color1_", an)]], 
                                            input[[paste0("color2_", an)]])
                                          )(nlevels(factor(meta_data[[an]])))
        
        selected_color <- setNames(selected_color, levels(factor(meta_data[[an]])))
        
        Legend(at = levels(factor(meta_data[[an]])),
              labels = levels(factor(meta_data[[an]])),
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
      
      #                   COLOR SCHEMES
      #___________________________________________________
      
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
      
      
      x <- Heatmap(vst_cor,   name = "Value",  col=col_fun,

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

      #           DRAW PLOT
      #______________________________________________________
      draw(x, padding = unit(c(0.1, 0.03, 0.15, 0.35), "npc"))
      
      #           DRAW BODY LEGEND
      #_______________________________________________________
      draw(l, x= unit(0.85, "npc"), y = unit(0.5, "npc"), 
           just = c("center", "right"))
      
      #         PACK AND DRAW ANNOTATION LEGENDS
      #______________________________________________________
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
    
     # UI for selecting annotations
    output$heatmap_annotations <- renderUI({
      
      #LOCALIZE VARIABLES
      meta_data <- RCV.META_DATA()
      
      if(is.null(meta_data)){
        return()
      }
      
      div(
        checkboxGroupInput('heatmap_cond', 'Meta data conditions to view', 
                          choices=colnames(meta_data),
                          selected=colnames(meta_data),
                          inline=FALSE)
      )
    })
    
    # The ui to display the heatmap plot
    output$heatmap_plots_ui <- renderUI({
      
      #LOCALIZE VARIABLES
      counts <- RCV.VST_COUNTS()
      
      if(is.null(counts)){
        return(span("No Data Yet",style="color: red;"))
      }  
      
      if( (ncol(counts) * 45) > 800){
        h <- ncol(counts) * 45 + 200
      }else{
        h <- 800
      }
      
      h <- paste0(h, "px")
     
      fluidPage(
        plotOutput('heatmap_plot1', height= h,width = "100%")
      )
      
    })
    
    # displays the correct color editors for the selected annotations
    observe({
      
      output$color_pickers_heat <- renderUI({
        
        lapply(input$heatmap_cond, function(cond) {
          
          div(
            colourpicker::colourInput(inputId = paste0("color1_", cond),
                      label = paste("Choose color1 for", cond),
                      value = "grey"),
            colourpicker::colourInput(inputId = paste0("color2_", cond),
                        label = paste("Choose color2 for", cond),
                        value = "blue")
          )
        })
      })
    })
  
  #______________________________Diff 6_______DONE_________________#----
    
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
        
        h2("Lowest Fold Change"),
        
        renderPlot({high_fold}),
        
        h2("Highest Fold Change"),
        
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
                            selected=colnames(meta_data),
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
            
            x <- gplop(res, meta_data, input$gene_count_plot_search_cond)
            renderPlot({x})
            
          }else{
            p("Nothing Here Yet")
          }
          
        })

      )
    })
    
  #______________________________Diff 7_______DONE_________________#----
  
    # The plot for the volcano plot page
    output$volcano_plot <- renderPlot({
      
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
      
      plot
      
    })
    # A ui for displaying the volcano plot
    output$volcano_plot_ui <- renderUI({
      
      #LOCALIZE VARIABLES
      results <- RCV.RESULTS_DDSC()
      
      if(is.null(results)){
        return(span("No Data Yet",style="color: red;"))
      }
      
      plotOutput('volcano_plot', height= "900px", width = "100%")
      
    })

  #----
    
  #Page2 pretab----
    output$about_pathfinder <- renderUI({
      includeHTML('www/PathfindR_Intro.HTML')
    })
  #----  
  
  #Tabset 2
  #______________________________Path 1_______DONE_________________#----
    
    output$pathfinderPreview <- renderUI({
      
      #LOCALIZE VARIABLES
      res <- RCV.PATHWAY_RESULTS()
      
      if( is.null(res) ){ return() }
      
      #ORGANIZE
      data <- res %>% 
        select(-ID, -Up_regulated,	-Down_regulated, -all_pathway_genes)
      
      renderDT({data}, rownames = FALSE)
      
    })
  
  #______________________________Path 2_______DONE_________________#----
    
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
      div( plotOutput('enrichmentChart', height = h, width = w) )
      
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
    
  #______________________________Path 3_______DONE_________________#----
  
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
      
      #CREATE PLOT OUTPUT
      output$plot_heatmap_10 <- renderPlot({
        plot[[1]] +
          theme(
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10)
          ) +
          scale_x_discrete(position = 'top')
      }) 
      
      #UI OUTPUT
      plotOutput('plot_heatmap_10', 
                 width = paste0(plot[[2]]$x$layout$width + 200, "px"), 
                 height = paste0(plot[[2]]$x$layout$height + 275, "px"))
      
      
    })
    
  #______________________________Path 4_______DONE_________________#----
  
    output$sample_conditions_PreviewTable11 <- renderUI({
      
      #LOCALIZE VARIABLES
      meta_data <- RCV.PATHWAY_META_DATA()
      
      #REQUIRE DATA
      if(is.null(meta_data)){ 
        return(span("Upload Meta Data",style="color: red;"))
      }
      
      #RENDER TABLE
      
      renderTable({meta_data}, rownames  = TRUE) 
      

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
      output$case_map_plot <- renderPlot({plotdata[[1]]})
      
      #RENDER UI
      div(
        plotOutput('case_map_plot', height = h, width = w)
      )
      
    })

    
  #______________________________Path 5_______DONE_________________#----

    output$Single_pathway_plot_ui <- renderUI({
      
      #LOCALIZE VARIABLES
      res <- RCV.PATHWAY_RESULTS()
      p.ab <- RCV.PATHWAY_ABUNDANCE()
      norm <- RCV.NORMALIZED_COUNTS()
      deg.res <- RCV.RESULTS_DDSC()
      
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
      
      #   PLOT THE DATA
      #______________________________
      p <- pheatmap(data, 
               scale = 'row', 
               fontsize = 15, 
               angle_col = "45", 
               legend_breaks = c(-2, -1, 0, 1, 2), 
               main = input$Single_pathway_plot_search)
      
      #RENDER PLOT
      output$Single_pathway_plot <- renderPlot({p})
      
      #RENDER UI
      div(
        plotOutput('Single_pathway_plot', height = h, width = w)
      )
      
    })
    
    
    
#----
} #X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X

