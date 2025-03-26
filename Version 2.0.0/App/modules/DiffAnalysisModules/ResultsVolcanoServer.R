ResultsVolcanoServer <- function(id, Data, Uploads){
  moduleServer(id, function(input, output, session){
    
    ns <- session$ns
    
    #Data ~~~~~~~~~
    # DiffResults <- list(
    #    "contrasts" = list, can be null
    #    "ddsc" = Reactive, can be null
    #    "Results" = Reactive
    #    "vstObject" = Reactive
    #    "vstCounts" = Reactive
    #    "NormCounts" = Reactive
    #    "status" = Reactive 0, 1
    # )
    #~~~~~~~~~~~~~~~
    #Uploads ~~~~~~~~~
    # DiffUploads <- list(
    #   "RawCounts" = Reactive,
    #   "GeneNames" = Reactive,
    #   "MetaData" = Reactive,
    #   "PrevResults" = Reactive list,
    #   "status" = Reactive 0, 1, 2
    # )
    # #~~~~~~~~~~~~~~~
    
    clean_data <- 
      function(data, pval, lab_density, population, genes, up_reg, down_reg){
      
      if(is.null(data)){return()}
      
      #______________Handle pvals of 0_____________
      #----
      
      #They are set to 1e-300
      data <- data %>%
        mutate( padj = ifelse( padj == 0 , 1e-300 , padj ) )
      
      #----
      
      data$gene_name <- rownames(data)
      
      #_____Create a row with expression direction with a character representation____
      #----
      
      #down regulated genes are "DOWN"
      #up regulated genes are "UP"
      #genes found in the search vector are "FOUND"
      #genes between expression cutoffs are "NO"
      
      data <- data %>%
        
        mutate(ex = case_when(
          
          tolower(gene_name) %in% genes ~ "FOUND",
          
          log2FoldChange > up_reg & padj < pval ~ "UP",
          
          log2FoldChange < down_reg & padj < pval ~ "DOWN",
          
          .default = as.character("NO"))
          
        )
      
      #----
      
      #___Exclude genes found in the search vector from downward filtering___
      #----
      
      found_genes <- NULL
      
      if(! is.null(genes) ){
        
        found_genes <- data %>%
          filter(ex == "FOUND") %>%
          mutate(glabel = gene_name)
        
        data <- data %>% filter(ex != "FOUND")
        
      }
      
      #----
      
      #___Filtering and sorting the data_____
      #----
      
      data <- data %>%
        arrange(
          ex == "NO", #Genes with no expression are given lowest priority (put at bottom)
          padj,       #Genes with lowest Pvalues have highest priority (put at top)
          desc( abs( log2FoldChange ) ) #fold change is used in case where there is similar pvalue
        )
      #----
      
      #____Extract a percentage of the data to plot given by the population argument___
      #----
      
      #In case of 0%, 20 genes are plotted anyways
      
      data <- head( data , ( ( round(nrow(data) * population) ) + 20 ) )
      
      #----
      
      #___Extract a percentage of the genes to name based on the lab density argument___
      #----
      
      # 100% will always plot a maximum of 230 labels
      # any  genes found withing the top 200 genes by pvalue
      # and any genes found within top 30 genes by foldchange
      # (with overlap)
      # This value can be increased past 100% to label more than 230 genes but plot may
      # get cluttered
      #
      # Genes selected by P value MUTS also have some expression
      #
      # NOTE: these dataframes only consist of gene names
      p_dens <- round( 10 + (200 - 10) * lab_density )
      l_dens <- round( 5 + (30 - 5) * lab_density )
      
      p_top <- data %>%
        filter( ex != "NO" ) %>%
        arrange(padj) %>%
        pull(gene_name) %>%
        head(n = p_dens)
      
      l_top <- data %>%
        arrange( desc(abs(log2FoldChange)) ) %>%
        pull(gene_name) %>%
        head(n = l_dens)
      
      #----
      
      #_____Create a label column in the dataframe___
      #----
      
      #Genes found in either p_top or l_top are labeled by their gene_name
      #Genes not found in either are NA
      
      data$glabel <- ifelse(
        
        data$gene_name %in% c(p_top, l_top),
        
        data$gene_name, #value for if
        
        NA #value for else
      )
      
      return(list(data, found_genes))
      
      
    }
    
    #look up table
    output$Lookup <- renderDT({
      
      #df <- data.frame(Uploads$GeneNames())
      df <- data.frame(1:1000)
      colnames(df) <- "Genes"
      
      datatable(
        df, rownames = FALSE, 
        options = list(
          pagelength = 10, dom = "ftp", pagingType = "simple",
          language = list(
            search = ""
          )
        )
      )
      
    })
    
    Title <- debounce(reactive(input$Title), 750)
    Subtitle <- debounce(reactive(input$Subtitle), 750)
    Caption <- debounce(reactive(input$Caption), 750)
    Search <- 
      debounce(reactive(
        strsplit(input$search,split=",",fixed=TRUE)[[1]] %>% 
          trimws("both") %>% tolower()), 
      750)
    
    Points <- debounce(reactive(input$points), 750)
    Labels <- debounce(reactive(input$labels), 750)
    min <- debounce(reactive(input$min), 750)
    max <- debounce(reactive(input$max), 750)
    pline <- debounce(reactive(input$PLine), 750)
    vData <- reactiveVal("foo")
    

    observeEvent(list(Points(), Labels(), min(), max(), Search(), Data$Results()), {
      vData(
        clean_data(Data$Results(), pline(), Labels(), Points(), Search(), max(), min())
      )
    })
    
    output$main <- renderPlot({
      
      if(is.null(vData())){return()}
      
      plot <- ggplot(data = vData()[[1]], 
             aes(
               x = log2FoldChange, 
               y = -log10(padj), 
               color = ex)
             ) +
        
        geom_point(size = 2) +
        geom_vline(xintercept = min(), linetype = "dotted") + 
        geom_vline(xintercept = max(), linetype = "dotted") +
        geom_hline(yintercept = -log10(pline()), linetype = "dotted") + 
        
        geom_label_repel(aes(x = log2FoldChange,y = -log10(padj), label = glabel), 
                         alpha = 0.8, max.overlaps = Inf)+
        
        scale_color_manual(
          values = c("DOWN" = "#2171b5", "NO" = "grey", "UP" = "#bb0c00"),
          labels = c("DOWN" = "Downregulated",
                     "NO" = "Not significant",
                     "UP" = "Upregulated")) +
        
        scale_x_continuous() +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        
        labs(
          x = expression("logFC"),
          y = expression("-log"[10]*"p-value"),
          color = "Genes",
          title = Title(),
          subtitle = Subtitle(),
          caption = Caption()
        ) + 
        
        theme(
          plot.margin = margin(10, 10, 10, 10, "pt"),
         
          text = element_text(size = input$Font, color = "black"),
          plot.title = element_text(size = input$Font + 5, color = "black", margin=margin(20,20,5,10,"pt")),
          plot.subtitle = element_text(size = input$Font + 2, color = "black", margin=margin(5,5,15,10,"pt")),
          plot.caption = element_text(size = input$Font - 2, color = "black", margin=margin(10, 10, 10, 10, "pt")),
          
          legend.margin = margin(15, 15, 15, 15, "pt"),

          panel.background = element_rect(fill = "grey98"),
          panel.grid  = element_line(color = "grey80"),
          axis.line = element_line(color = "grey40"),
        ) + 
        
        eval(parse(text=as.character(input$Grid))) +
        #Add custom Aspect ratio
        theme(aspect.ratio = as.numeric(input$Aspect))
      
      #Process Searched stuff
      if(length(Search) > 0){
        plot <- plot + 
          geom_point(data = vData()[[2]],
            aes(x = log2FoldChange,
                y = -log10(padj)
            ),
            color='mediumseagreen',
            size=3,
            show.legend = FALSE) +

          geom_label_repel(data = vData()[[2]],
                           aes(x = log2FoldChange,y = -log10(padj), label = gene_name), 
                           alpha = 0.8, max.overlaps = Inf,
                           color = 'mediumseagreen')
      }

      
      plot
      
    },
      width = function() {ifelse(is.na(input$W), 1000, input$W)},
      height = function() {ifelse(is.na(input$H), 750, input$H)}
    )
    
  })
  
}