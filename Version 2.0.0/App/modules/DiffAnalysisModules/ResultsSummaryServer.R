ResultsSummaryServer <- function(id, Data, Uploads){
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
   
    #Pvalue Table stats
    output$PvalStats <- renderDT({
      
      if(is.null(Data$Results())  || nrow(Data$Results()) < 1){return()}
      
      stats <- Data$Results() %>% na.exclude() %>%
        filter(padj <= input$Pval) %>%
        summarize(
          "Max" = max(padj),
          "Mean" = mean(padj),
          "Min" = min(padj)
        )
      
      datatable(
        t(stats), 
        colnames = c("Adj. P-val"),
        options = list(
          dom = "t"
        )
      )
      
    })
    
    #Pval Histogram
    output$PvalHisto <- renderPlot({
     
      if(is.null(Data$Results()) || nrow(Data$Results()) < 1){return()}
      
      ggplot(Data$Results() %>% filter(padj <= input$Pval), aes(pmin(padj, 0.5))) +
        geom_histogram(bins = 10, fill = '#677780') +
        theme_minimal() +
        theme(
          axis.text.x = element_text(color = 'black', size = 16),
          axis.text.y = element_blank(),
          axis.title = element_text(color = 'black', size = 17),
          panel.grid.major = element_line(color = 'grey25'),
          panel.grid.minor = element_blank(),
          plot.margin = margin(5,5,5,5,"pt")
        ) +
        labs(x = "Adjusted P-value", y = "Count") +
        scale_x_continuous(breaks = c(0, 0.1, 0.5),
                           labels = c("0", "0.1", "> 0.5"))

    })
    
    #Fold Change Stats
    output$FoldStats <- renderDT({
      
      if(is.null(Data$Results()) || nrow(Data$Results()) < 1){return()}
      
      stats <- Data$Results() %>% na.exclude() %>%
        filter(padj <= input$Pval) %>%
        summarize(
          "Max" = max(log2FoldChange),
          "Mean" = mean(log2FoldChange),
          "Min" = min(log2FoldChange)
        )
      
      datatable(
        t(stats), 
        colnames = c("Fold Change"),
        options = list(
          dom = "t"
        )
      )
      
    })
  
    #Fold Change Histogram
    output$FoldHisto <- renderPlot({
     
      if(is.null(Data$Results()) || nrow(Data$Results()) < 1){return()}
      ggplot(Data$Results() %>% filter(padj <= input$Pval), aes(pmax(pmin(log2FoldChange, 5), -5))) + 
        geom_histogram(bins = 12, fill = '#677780') + 
        theme_minimal() + 
        theme(
          axis.text.x = element_text(color = 'black', size = 16),
          axis.text.y = element_blank(),
          axis.title = element_text(color = 'black', size = 17),
          panel.grid.major = element_line(color = 'grey25'),
          panel.grid.minor = element_blank(),
          plot.margin = margin(5,5,5,5,"pt")
        ) +
        labs(x = "Fold Change", y = "Count") +
        scale_x_continuous(breaks = c(-5, -1, 0, 1, 5), labels = c( "< -5", "-1", "0", "1", "> 5"))
      
    })
    
    #Full Results DT
    output$Res <- renderDT({
      
      if(is.null(Data$Results()) || nrow(Data$Results()) < 1){return()}
      
      datatable(
        Data$Results() %>% 
          filter(log2FoldChange >= input$FoldRange[1] & log2FoldChange <= input$FoldRange[2]) %>%
          filter(padj <= input$Pval), 
        rownames = TRUE, 
        options = list(pageLength = 15)
      )
    })
    
    #Show Extra Rows 
    
  })
}