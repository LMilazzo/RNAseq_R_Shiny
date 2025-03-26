ResultsGeneCountsServer <- function(id, Data, Uploads){
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
  
    #Debounce Text Inputs
    Title <- debounce(reactive(input$Title), 750)
    Subtitle <- debounce(reactive(input$Subtitle), 750)
    Caption <- debounce(reactive(input$Caption), 750)
    Main <- debounce(reactive(input$main), 750)
    
    #look up table
    output$Lookup <- renderDT({
      
      df <- data.frame(Uploads$GeneNames())
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
    
    #Annotation selection
    output$x <- renderUI({
      
      selectInput(
        ns("x"),
        NULL,
        choices = colnames(Uploads$MetaData()),
        selected = 
          ifelse(
            is.null(Data$ddsc()),
            colnames(Uploads$MetaData())[1],   # NULL case
            all.vars(design(Data$ddsc()))   #Not NULL case
          ),
        width = "100%"
      )
      
    })
    output$color <- renderUI({
      
      selectInput(
        ns("color"),
        NULL,
        choices = c("None", colnames(Uploads$MetaData())),
        selected = c("None"),
        # ifelse(
        #   is.null(Data$ddsc()),
        #   colnames(Uploads$MetaData())[1],   # NULL case
        #   all.vars(design(Data$ddsc()))   #Not NULL case
        # )
        width = "100%"
      )
      
    })
    output$shape <- renderUI({
      
      selectInput(
        ns("shape"),
        NULL,
        choices = c("None", colnames(Uploads$MetaData())),
        selected = c("None"),
        # ifelse(
        #   is.null(Data$ddsc()),
        #   colnames(Uploads$MetaData())[1],   # NULL case
        #   all.vars(design(Data$ddsc()))   #Not NULL case
        # )
        width = "100%"
      )
      
    })
    output$size <- renderUI({
      
      selectInput(
        ns("size"),
        NULL,
        choices = c("None", colnames(Uploads$MetaData())),
        selected = c("None"),
        # ifelse(
        #   is.null(Data$ddsc()),
        #   colnames(Uploads$MetaData())[1],   # NULL case
        #   all.vars(design(Data$ddsc()))   #Not NULL case
        # )
        width = "100%"
      )
      
    })

    output$main <- renderPlot({
      #print(Main())
      #print(Uploads$GeneNames()[1])
      if(Main() == ""){
        main <- Uploads$GeneNames()[1]
      }else{
        main <- Main()
      }
      
      data <- validateSearch(Data, Uploads, main)
      
      a <- aes(y = gene)
      
      # Dynamically add aesthetics based on the variables
      if (!is.na(input$x) && input$x != "None") { a$x <- as.name(input$x) }
      if (!is.na(input$color) && input$color != "None") { a$colour <- as.name(input$color) }
      if (!is.na(input$shape) && input$shape != "None") { a$shape <- as.name(input$shape) }
      if (!is.na(input$size) && input$size != "None") { a$size <- as.name(input$size) }
      
      plot <- ggplot(data, a) +
        labs(
          title = ifelse(is.null(Title()) || Title() == "", main, Title()),
          subtitle = Subtitle(), 
          caption = Caption(),
          x = input$x,
          y = "Count",
          colour = input$color,
          shape = input$shape
        ) + 
        #Add custom Aspect ratio
        theme(aspect.ratio = as.numeric(input$Aspect)) +
        #Add the custom font size
        theme(
          text = element_text(size = input$Font, color = "black"),
          plot.title = element_text(size = input$Font + 5, color = "black", margin=margin(20,20,5,10,"pt")),
          plot.subtitle = element_text(size = input$Font + 2, color = "black", margin=margin(5,5,15,10,"pt")),
          plot.caption = element_text(size = input$Font - 2, color = "black", margin=margin(10, 10, 10, 10, "pt"))
        ) +
        #EXtra theme
        theme(
          plot.margin = margin(10,10,10,10,"pt"),
          panel.background = element_rect(fill = "grey98"),
          panel.grid  = element_line(color = "grey80"),
          axis.line = element_line(color = "grey40"),
          legend.margin = margin(15, 15, 15, 15, "pt"),
          axis.title = element_text(margin = margin(15, 15, 15, 15, "pt"))
        ) +
        
        #Add The custum grid lines
        eval(parse(text=as.character(input$Grid)))
        
        
      #ADD Final size scale
      if(input$size != "None"){
        plot <- plot + scale_size_discrete(range = c(3, 10)) + geom_jitter(width = 0.2, height = 0)
      }else{
        plot <- plot + geom_jitter(size=4, width = 0.2, height = 0)
      }
      
      plot
      
    },
      width = function() {ifelse(is.na(input$W), 750, input$W)},
      height = function() {ifelse(is.na(input$H), 450, input$H)}
    ) 
    
    #Validates a search by merging nessecary things and converting to lower case
    validateSearch <- function(Data, Uploads, search){
      
      search <- tolower(search)
      
      print(Data$NormCounts() %>% head())
      
      y <- Uploads$MetaData() %>% 
        mutate(gene_name = rownames(Uploads$MetaData())) %>%
        mutate(gene = 
          (Data$NormCounts() %>% 
            mutate(gene_name = tolower(rownames(Data$NormCounts()))) %>%
            filter(gene_name == search) %>% select(-gene_name) %>%t())[,1]
        )
      
      return(y)
      
    }
    
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  })
  
}