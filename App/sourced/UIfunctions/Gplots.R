singleGenePlot <- function(data, counts, md, cond, searched, settings){

  #Check data
  if(is.null(data) || is.null(counts) || 
     is.null(searched) || length(searched) < 1 ||
     is.null(cond) || length(cond) < 1){return()}
  
  #Set up data
  counts <- data.frame(round(counts))
  counts$gene_name <- rownames(counts) %>% tolower()

  #Set up results
  data$gene_name <- rownames(data) %>% tolower()
  data <- data %>% select(gene_name, padj, log2FoldChange)
  
  #Merge
  combined <- left_join(data, counts)
  
  #find searched
  g <- combined %>% filter(gene_name == tolower(searched))
  
  #ASSERT SINGLE ROW
  if(is.null(g) || ncol(g) == 0 || nrow(g) == 0){
    return()
  }
  if(nrow(g) > 1 ){
    g <- g[1,]
  }
  
  #Plotting
  plot <- gplop(g, md, cond, settings$TSA)
  if(!is.null(settings$Title)){plot <- plot + labs(title = settings$Title)}
  if(!is.null(settings$Subtitle)){plot <- plot + labs(subtitle = settings$Subtitle)}
  if(!is.null(settings$Caption)){plot <- plot + labs(caption = settings$Caption)}
  
  if(!settings$Aspect == "auto") {
    plot <- plot + theme(aspect.ratio = settings$Aspect)
  }
  
  #Render
  rendered <- renderPlot({plot})
  
  #Work out sizing
  if(settings$Width == "auto"){w <- settings$DefaultWidth}
  else{w <- paste0(settings$Width, "px")}
  if(settings$Height == "auto"){h <- settings$DefaultHeight}
  else{h <- paste0(settings$Height, "px")}
  
  return(
    list(
      "plot" = plot,
      "rendered" = rendered,
      "w" = w,
      "h" = h
    )
  )
}