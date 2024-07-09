
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
  
  if( !grepl('.tsv', raw_counts_table, fixed=TRUE)){
      return("Not .tsv")
  }
  
  count <- read.csv(raw_counts_table, sep="\t")

  if(!'gene_id' %in% colnames(count) || !'gene_name' %in% colnames(count)){
    return("Bad .tsv")
  }

  gene_names <- count %>% 
                select(gene_name , gene_id) %>%
                mutate(gene = gene_id) %>%
                select(-gene_id)
  
  
  count <- count %>% select(-gene_name)
  
  count <- count %>% tibble::column_to_rownames('gene_id')
  
  count <- as.matrix(count)
  
  return( list(count, gene_names) )
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

  fcounts <- counts
  
  keep <- rowSums(fcounts) > 10
  
  fcounts <- fcounts[keep,]
  
  fcounts <- round(fcounts)
  
  fcounts <- as.matrix(fcounts)
  
  return(fcounts)
}

