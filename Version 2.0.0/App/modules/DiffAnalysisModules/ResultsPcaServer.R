ResultsPcaServer <- function(id, Data, Uploads){
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
    
    #Number Of rows to use               #REMOVE RESTAINT
    output$N <- renderUI({
      numericInput(
        ns("n"), NULL, min = 500, 
        max = nrow(Data$vstCounts()),
        value = 500,
        width = "80%"
      )
      
    })
    
    #Annotation selection
    output$color <- renderUI({
      
      selectInput(
        ns("color"),
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
    
    #Plot output
    output$PCA <- renderPlot(
      {
        
        getPCA(
          Data$vstObject(), 
          input
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
          #Add The custum grid lines
          eval(parse(text=as.character(input$Grid))) + 
          #Labels
          labs(
            title = Title(), 
            subtitle = Subtitle(), 
            caption = Caption()
          )
      
      },
      width = function() {ifelse(is.na(input$W), 750, input$W)},
      height = function() {ifelse(is.na(input$H), 450, input$H)}
    )
    
    #Helper Function
    #vst = vstObject
    #metadata list of conditions to apply
    #n number of elements
    getPCA <- function(vst, input){
      
      intg <- setdiff( c(input$color, input$shape, input$size), c("None"))
      
      print(intg)
      
      pca.data <- plotPCA(
        vst, intgroup = intg, ntop = ifelse(is.na(input$n), 500, input$n), returnData = TRUE
      )
      
      percentVar <- round(100 * attr(pca.data, "percentVar"))
      
      #X Y aes
      a <- aes(x = PC1, y = PC2)
      
      #Set additional aes  Color shape size
      if (!is.na(input$color) && input$color != "None") { a$colour <- as.name(input$color) }
      if (!is.na(input$shape) && input$shape != "None") { a$shape <- as.name(input$shape) }
      if (!is.na(input$size) && input$size != "None") { a$size <- as.name(input$size) }
      
      #Make Plot Barebones
      plot <- ggplot(pca.data, a)+
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        labs(
          colour = input$color,
          shape = input$shape,
          size = input$shape
        ) + 
        theme(
          plot.margin = margin(10,10,10,10,"pt"),
          panel.background = element_rect(fill = "grey98"),
          panel.grid  = element_line(color = "grey80"),
          axis.line = element_line(color = "grey40"),
          legend.margin = margin(15, 15, 15, 15, "pt"),
          axis.title = element_text(margin = margin(15, 15, 15, 15, "pt"))
        ) 
  
      #Change size of points according to present number of variables
      if(input$size != "None"){
        plot <- plot + scale_size_discrete(range = c(3, 10)) + geom_point()
      }else{
        plot <- plot + geom_point(size=4)
      }

      return(plot)
    }
    
  })
}