
library(readr)
library(dplyr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_______________UI Functions_______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------File Upload--------#
#~~~~~~~~~~Widget~~~~~~~~~~#
# Creates a widget for a file upload
# @param Identifing name for the uploaded file
# @param statement to be displayed near the widget
FileUpload <- function(widget_Identfier, widget_Statement){
  
  fileInput(widget_Identfier, widget_Statement)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------File Preview-------#
#~~~~~~~~~~Widget~~~~~~~~~~#
# Creates a widget for a selecting what preview style should appear 
# @param Identifing name for the preview type describing which object it belongs to
# @param statement to be displayed near the widget
FilePreviewSize <- function(widget_Identfier, widget_Statement){
  
  #Options <- A vector of character choices defining the display types
  options <- c("None", "Head", "Full")
  
  selectInput(widget_Identfier, widget_Statement, options)
  
}


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
    
    return("Error: raw counts matrix upload expected a .tsv file")
  
  }
  
  count <- read.csv(raw_counts_table, sep = "\t")

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
#----Make Preview Table----#
#~~~~~~~~~~~Func~~~~~~~~~~~#
#based on the table and view settings returns the right rows/columns to turn
#into a data table
#@param the data frame or matrix to view
#@param the view setting
makePreviewTable <- function(subject, setting){
  if(setting == "Full"){
    subject
  }
  else if(setting == "Head"){
    head(subject)
  }
  else if(setting == "Column Names"){
    colnames(subject)
  }
  else if(setting == "Row Names"){
    rownames(subject)
  }
  else{
  }
}

