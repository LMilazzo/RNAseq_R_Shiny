
pcaPlot <- function(input, vst, advSettings){
  
  #InputHandling
  annotations <- input$pca_cond
  print(annotations)
  n <- input$pca_nrow
  aspect <- ifelse(advSettings$Aspect == "auto", 1, advSettings$Aspect)

  #Handle missing data
  if(is.null(vst)){return()}
  if(length(annotations) == 0){return()}
  
  #Plot the data
  plot_object <- principlePlot(
                   vst,
                   annotations,
                   n,
                   advSettings$Title,
                   advSettings$Subtitle, 
                   advSettings$Caption,
                   advSettings$TSA
                 )
  
  #Add aspect ratio
  plot_object <- plot_object + theme(aspect.ratio = aspect) 
  
  #Set up plotting dimensions
  
  if(advSettings$Height == "auto"){ 
    height <- paste0(advSettings$DefaultHeight, "px")
  }else{
    height <- paste0(advSettings$Height, "px")
    
  }
  
  if(advSettings$Width == "auto"){
    width <-  paste0(advSettings$DefaultWidth, "px")
  }else{
    width <- paste0(advSettings$Width, "px")
  }
  

  #RENDER
  renderedPlot <- renderPlot({plot_object})
  
  return(
    list(
      "plot" = plot_object,
      "rendered" = renderedPlot,
      "w" = width,
      "h" = height
    )
  )
  
}


#Creates the ggplot for principle component analysis
principlePlot <- function(data, cond, genes, title, sub, cap, TSA){
  
  pca_data <- plotPCA(data, intgroup=cond, ntop=genes, returnData=TRUE)
  
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  aesthetics <- aes(x = pca_data$PC1, y = pca_data$PC2 )
  
  x <- length(cond)
  
  if(x >= 1){
    aesthetics$colour = as.name(cond[1])
  }
  if(x >= 2){
    aesthetics$shape = as.name(cond[2])
  }
  if(x >= 3){
    aesthetics$size = as.name(cond[3])
  }
  
  pca_plot <- ggplot(pca_data,aesthetics)+
    
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    
    theme_minimal() +
    
    labs(title = title,
         subtitle = sub,
         caption = cap,
         color = if (x >= 1) cond[1] else NULL,
         shape = if (x >= 2) cond[2] else NULL,
         size = if (x >= 3) cond[3] else NULL
    ) +
    
    theme(plot.margin = margin(10, 10, 10, 10, "pt"),
          axis.title.x = element_text(color='black', size = 20 + TSA, margin = margin(15, 15, 15, 15, "pt")),
          axis.title.y = element_text(color='black', size= 20 + TSA, margin = margin(15, 15, 15, 15, "pt")),
          axis.text.y = element_text(color='black', size= 15 + TSA),
          axis.text.x = element_text(color='black', size= 15 + TSA),
          panel.grid.minor.y = element_blank(), 
          panel.grid.minor.x = element_blank(),
          axis.line.x = element_line(color='grey'),
          axis.line.y = element_line(color='grey'),
          
          legend.text = element_text(color='black', size=20 + TSA),
          legend.title = element_text(color='black', size=20 + TSA),
          legend.margin = margin(15, 15, 15, 15, "pt"),
          
          plot.title = element_text(color='black', size=30 + TSA, margin=margin(20,20,5,10,"pt")),
          plot.subtitle = element_text(color='black', size=20 + TSA, margin=margin(5,5,15,10,"pt")),
          plot.caption = element_text(color='black', size=15 + TSA, margin=margin(10, 10, 10, 10, "pt")),
          
    ) +
    coord_fixed(1)
  
  #Change size of points according to present number of variables
  if(x >= 3){
    pca_plot <- pca_plot + scale_size_discrete(range = c(3, 10)) + geom_point()
  }else{
    pca_plot <- pca_plot + geom_point(size=4)
  }
  
  return(pca_plot)
}