
library(readr)
library(dplyr)

# Primitive function for displaying a modal that displays and error
showErrorModal <- function(message) {
  showModal(modalDialog(tags$p(style = "color: red;", message), easyClose = TRUE, footer = NULL))
}


# Splits a DESeq result dataset into 3 dataframes base on fold change
splitByExpr <- function(data, gene_names, cut){
  
  lowerBound <- cut[1]
  upperBound <- cut[2]
  
  data$gene_id <- row.names(data)
  
  data <- left_join(data, gene_names) %>% select(-gene_id)
  
  
  res <- data.frame(data)

  up <- res %>% filter(log2FoldChange > upperBound)
  up <- up[order(up$padj) , ]
  
  down <- res %>% filter(log2FoldChange < lowerBound)
  down <- down[order(down$padj) , ]
  
  noR <- res %>% 
         filter(log2FoldChange > lowerBound) %>% 
         filter(log2FoldChange < upperBound)
  
  noR <- noR[order(noR$padj) , ]
  
  return( list( up, down, noR ) )
  
}


# Returns a ggplot of the principle components
principlePlot <- function(data, cond, genes, title, sub, cap){
  
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
              axis.title.x = element_text(color='black', size = 20, margin = margin(15, 15, 15, 15, "pt")),
              axis.title.y = element_text(color='black', size=20, margin = margin(15, 15, 15, 15, "pt")),
              axis.text.y = element_text(color='black', size=15),
              axis.text.x = element_text(color='black', size=15),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.line.x = element_line(color='grey'),
              axis.line.y = element_line(color='grey'),

              legend.text = element_text(color='black', size=20),
              legend.title = element_text(color='black', size=20),
              legend.margin = margin(15, 15, 15, 15, "pt"),

              plot.title = element_text(color='black', size=30, margin=margin(20,20,5,10,"pt")),
              plot.subtitle = element_text(color='black', size=20, margin=margin(5,5,15,10,"pt")),
              plot.caption = element_text(color='black', size=15, margin=margin(10, 10, 10, 10, "pt")),

              aspect.ratio = 1) +
        coord_fixed()

      #Change size of points according to present number of variables
      if(x >= 3){
        pca_plot <- pca_plot + scale_size_discrete(range = c(3, 10)) + geom_point()
      }else{
        pca_plot <- pca_plot + geom_point(size=4)
      }

      return(pca_plot)
}

