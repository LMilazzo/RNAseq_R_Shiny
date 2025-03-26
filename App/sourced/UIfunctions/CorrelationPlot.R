correlate <- function(counts, md, annotations, advSettings){

  #ASSERT NEEDED VARIABLES
  if(is.null(counts)){return()}
  
  #Create correlation matrix
  vst_cor <- cor(counts)

  #get annotation colors
  ann_colors <- setUpAnnotations(annotations, md)
  toAnnotate <- md %>% select(all_of(annotations))
  
 
  #Cell sizing
  if(advSettings$CellWidth == "auto"){
    cell_w <- 700 / nrow(vst_cor)
  }#auto
  else{cell_w <- advSettings$CellWidth} 
  
  if(advSettings$CellHeight == "auto"){
    cell_h <- 700 / nrow(vst_cor)
  } #auto
  else{cell_h <- advSettings$CellHeight}
  
  if(advSettings$Height == "auto"){
    height <- paste0(advSettings$DefaultHeight, "px")
  }else{
    height <- paste0(advSettings$Height, "px")
  }
  
  if(advSettings$Width == "auto"){
    width <- paste0(advSettings$DefaultWidth, "px")
  }else{
    width <- paste0(advSettings$Width, "px")
  }
  
  
  #create plot
  p <- pheatmap(
    vst_cor, 
    angle_col = "45", 
    main = advSettings$Title, 
    annotation_row = if (nrow(toAnnotate) > 0 && ncol(toAnnotate) > 0) 
      toAnnotate else NA,
    annotation_colors = if (nrow(toAnnotate) > 0 && ncol(toAnnotate) > 0) 
      ann_colors else NA,
    fontsize = 20 + advSettings$TSA,
    cellwidth = cell_w,
    cellheight = cell_h
  )
  print(input$H)
  
  p <- as.ggplot(p) 
  renderedPlot <- renderPlot({(p)})
  
  print("++++++++++++++++++++++++++++++++++++++++++++")
  print(advSettings)
 
  return(
    list(
      "plot" = p,
      "rendered" = renderedPlot,
      "w" = width,
      "h" = height
    )
  )
  
  
  
}


generate_palette <- function(factor_levels, base_colors = c("cornsilk3", "orchid")) {
  # Check the number of unique factor levels
  unique_levels <- unique(factor_levels)
  num_levels <- length(unique_levels)
  
  # Generate the palette using colorRampPalette
  palette <- colorRampPalette(base_colors)(num_levels)
  
  # Return a named vector mapping levels to colors
  setNames(palette, unique_levels)
}


setUpAnnotations <- function(annotations, md){
  ann_colors <- list()
  for(i in annotations){
    if(nlevels(md[[i]]) > 2){
      add <- setNames(list(generate_palette(md[[i]], base_colors = c("cadetblue", "white", "indianred"))), i)
      ann_colors <- c(ann_colors, add)
    }else{
      add <- setNames(list(generate_palette(md[[i]])), i)
      ann_colors <- c(ann_colors, add)
    }
  }
  
  return(ann_colors)
}