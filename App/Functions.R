
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
  if(grepl('.tsv', raw_counts_table, fixed=TRUE)){
    counts <- read.csv(raw_counts_table, sep="\t")
  }else{
    counts <- read.csv(raw_counts_table)
  }

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
