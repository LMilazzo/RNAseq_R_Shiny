
library(readr)
library(dplyr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_______________UI Functions_______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_____________Server Functions_____________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-----Set Clean Counts-----#
#~~~~~~~~~~~Func~~~~~~~~~~~#
#Reads from the File upload containing the counts matrix. Counts matrix is 
#handled to return to the global data frame raw_counts in the correct and 
#following format:
  #axis.names| Sample1 | Sample 2 | Sample 3
  #Gene_ID_1 |    x    |    x     |    x
  #Gene_ID_2 |    x    |    x     |    x
#Input format is expected to be:
#gene_id  | gene_name | sample1	| sample2	
#ENSG...	|   XYXY    |    x	  |   x	
#ENSG...  |   XYXY    |	   x	  |   x	
# 
# @param the input list containing all ui inputs
# 
# @return raw_counts matrix
# @return gene_names data frame
Set_Clean_Counts <- function(raw_counts_table){
  
  if( !grepl('.tsv', raw_counts_table, fixed=TRUE) && !grepl('.csv', raw_counts_table, fixed=TRUE)){
    return("Bad File Type Error")
  }
  
  #READ as a .tsv or .csv
  if(grepl('.tsv', raw_counts_table, fixed=TRUE)){counts <- read.csv(raw_counts_table, sep="\t")} #read as tsv
  else{counts <- read.csv(raw_counts_table)} #read as csv

  #gene names and ids are recommended 
  if(!'gene_id' %in% colnames(counts) || !'gene_name' %in% colnames(counts)){
    return("Missing Column Error")
  }
  
  #Check for duplicated Gene IDS
  g_ids <- counts %>% select(gene_id)
  g_ids_unique <- g_ids %>% unique()
  if(nrow(g_ids) != nrow(g_ids_unique)){
    return("Duplicate Gene IDS")
  }
  
  gene_names <- counts %>% 
                select(gene_name , gene_id)
  
  counts <- counts %>% 
           select(-gene_name) %>% 
           tibble::column_to_rownames('gene_id')
  
  return( list( as.matrix(counts) , gene_names ) )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------Filter Counts------#
#~~~~~~~~~~~Func~~~~~~~~~~~#
#Takes raw_counts as input and returns a filtered version that removes low 
#counts genes
#@param raw_counts matrix
#@return the matrix of filtered counts
#@return the number of rows deleted
filterCounts <- function(counts){
  
  keep <- rowSums(counts) > 10
  
  counts <- counts[keep,]
  
  counts <- counts
  
  counts <- as.matrix(counts)
  
  return(counts)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------Split by Expr------#
#~~~~~~~~~~~Func~~~~~~~~~~~#
#Splits the results by expression level orders by P val
#@param the data from the results
#@param the list of geneids and gene names
#@return a list of the three new dataframes
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------plot gene count----#
#~~~~~~~~~~~Func~~~~~~~~~~~#
#Takes a single row with a gene name, id, foldchange, and padj
#metadata file, and conditions to view
#returns a plot with the counts as y value and other conditions suplementing
plotGeneCounts <- function(gene, metadata, conditions){
  
  name <- (gene %>% select(gene_name))[1]
  id <- (gene %>% select(gene_id))[1]
  p <- format((gene %>% select(padj))[1], scientific = TRUE, digits = 8)
  fold <- round((gene %>% select(log2FoldChange))[1], digits=8)
  counts <- (gene[,5:ncol(gene)])
  mean <- round(mean(as.numeric(counts)))
  
  if(is.null(metadata)){
    return()
  }
  
  toplot <- metadata %>% mutate(counts = data.frame(t(counts))[,1])
  
  print(toplot)
  
  print("_____________________")
  
  cond <- intersect(conditions, colnames(metadata))
  
  if(is.null(cond) || length(cond) <= 0){
    print("select a factor")
    return()
  }
  
  aesthetics <- aes(x = as.factor(toplot[[cond[1]]]) , y = toplot$counts)
  
  a <- length(cond)
  if(a >= 2){
    aesthetics$colour = as.factor(toplot[[cond[2]]])
  }
  if(a >= 3){
    aesthetics$shape = as.factor(toplot[[cond[3]]])
  }
  if(a >= 4){
    aesthetics$size = as.factor(toplot[[cond[4]]])
  }
  
  x <- ggplot(toplot, aesthetics) +
    xlab(cond[1]) +
    ylab("Count") +
  
    theme_minimal() +
    
    labs(title = name,
         subtitle = paste("Fold Change: ", fold),
         caption = paste("Mean Count: ",mean, "   Padj: ", p ),
         color = if (a >= 2) cond[2] else NULL,
         shape = if (a >= 3) cond[3] else NULL,
         size = if (a >= 4) cond[4] else NULL
    ) + 
    
    theme(plot.margin = margin(3, 2, 2, 2, "pt"),
          axis.title.x = element_text(color='black', size = 15, margin = margin(2, 2,2, 2, "pt")),
          axis.title.y = element_text(color='black', size=15, margin = margin(2, 2, 2, 2, "pt")),
          axis.text.y = element_text(color='black', size=15),
          axis.text.x = element_text(color='black', size=15),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.line.x = element_line(color='grey'),
          axis.line.y = element_line(color='grey'),
          
          legend.text = element_text(color='black', size=20),
          legend.title = element_text(color='black', size=20),
          legend.margin = margin(2, 2, 2, 2, "pt"),
          
          plot.title = element_text(color='black', size=20, margin=margin(2,2,2,2,"pt")),
          plot.subtitle = element_text(color='black', size=15, margin=margin(2,2,5,2,"pt")),
          plot.caption = element_text(color='black', size=10, margin=margin(2, 80, 2, 1, "pt")),
          
    )
  
  #Change size of points according to present number of variables
  if(a >= 4){
    x <- x + scale_size_discrete(range = c(3, 10)) + geom_beeswarm(cex = 3)
  }else{
    x <- x + geom_beeswarm(size=3, cex = 3)
  }
  
  return(renderPlot({x}))
}


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

