#Import Shiny and needed packages
library(shiny)
library(shinythemes)
library(DT)
library(plyr)
library(dplyr)
library(readr)
source("Functions.R")
library(crayon)

options(shiny.maxRequestSize=100*1024^2)  # Limits file upload size to 100 MB

ui <- navbarPage("DESeq2",  
   theme = shinytheme("sandstone"),
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 1----------------------------------#
#-----------------------------File Upload-------------------------------#
#---------------------------Widget Count: 2-----------------------------#
#-----------------------Expected output Count: 3------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   tabPanel("File Upload",
      fluidPage(
        sidebarLayout(
#______________________________Side Bar_________________________________#
          
          sidebarPanel(
            
            HTML("<h1>Upload Data</h1>")
            
            ,
            
            #-----------------input----------------#
            #---merged_gene_counts_uploaded_file---#
            #--------------------------------------#
            # An upload file widget for the merged gene counts tsv,csv
            #          Object Name    
            FileUpload("merged_gene_counts_uploaded_file", 
            #          Widget Message 
                       "Salmon Merged Counts .tsv file")
            
            ,

            
            HTML("<h1>Upload Meta Data</h1>")
            
            ,
            
            #-----------------input----------------#
            #--meta_data_conditions_uploaded_file--#
            #--------------------------------------#
            # An upload file widget for the metadata .csv file
            #          Object Name    
            FileUpload("meta_data_conditions_uploaded_file", 
                       #          Widget Message 
                       "Conditions")
            
            ,
            
          ), ##XX##~~~Side Panel End~~~##XX##
          
#_____________________________Main Panel________________________________#
          
          mainPanel(
            
            HTML("<h1>Raw Counts Table</h1>"),
            HTML("<p1>Upload your salmon merged counts matrix on the left. Your 
                 uploaded file must be a .tsv file. Below you can explore and search
                 through your unfiltered data once uploaded.")
            
            ,
            
            
            #Text output for error 1.0
            #Error: Raw Counts upload must be .tsv file
            uiOutput('errorMessagesPG1.0')
            
            ,
            
            
            #----------------output----------------#
            #--------raw_counts_PreviewTable-------#
            #--------------------------------------#
            #A Data Table preview of what is inside the uploaded counts table 
            #after formating correctly
              #Expected Format:
                #axis.names| geneName| Sample 2 | Sample 3
                #Gene_ID_1 |    x    |    x     |    x
                #Gene_ID_2 |    x    |    x     |    x
            DTOutput('raw_counts_PreviewTable')
            
            ,
            
            HTML("<h1>Conditions (meta data) Table</h1>"),
            HTML("<p1>Upload a .csv file with the samples in column 1 and their 
                 conditions in column 2.")
            
            ,
            
            #Text output for error 1.1
            #Error: Meta Data Table upload expected a .csv file
            uiOutput('errorMessagesPG1.1')
            
            ,
            
            #Text output for error 1.2
            #"Error: Incorrect .csv format
            uiOutput('errorMessagesPG1.2')
            
            ,
            
            #----------------output----------------#
            #----sample_conditions_PreviewTable----#
            #--------------------------------------#
            #A Data Table Preview of the samples and their respective conditions
            tableOutput('sample_conditions_PreviewTable')
            
          )##XX##~~~Main Panel End~~~##XX##
          
        )##XX##~~~Side Bar Layout Closing Bracket~~~##XX##
      ) ##XX##~~~Fluid Page Closing Bracket~~~##XX##
    ) ##XX##~~~Tab Panel Closing Bracket~~~##XX##
    
    ,

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 2----------------------------------#
#-----------------------------Filter Info-------------------------------#
#---------------------------Widget Count: 0-----------------------------#
#-----------------------Expected output Count: 2------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    tabPanel("DESeq2 and Filtering Info",
      fluidPage(
            
            HTML("<h2>Filter Notes</h2>")
        
            ,
            #----------------output----------------#
            #---------------nFiltered--------------#
            #--------------------------------------#
            #A message stating how many rows were filtered
            textOutput('nFiltered')
            
            ,
            
            #----------------output----------------#
            #-----filtered_counts_PreviewTable-----#
            #--------------------------------------#
            #A Data Table preview of what is inside the filtered counts matrix 
            #after filtering correctly
              #Expected Format:
              #axis.names| geneName| Sample 2 | Sample 3
              #Gene_ID_1 |    x    |    x     |    x
              #Gene_ID_2 |    x    |    x     |    x
            #DTOutput('filtered_counts_PreviewTable')
            
      ) ##XX##~~~Fluid Page Closing Bracket~~~##XX##
    ) ##XX##~~~Tab Panel Closing Bracket~~~##XX##
            
                 
)#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X
   
   
     
            
