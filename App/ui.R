#Import Shiny and needed packages
library(shiny)
library(shinythemes)
library(DT)
library(plyr)
library(dplyr)
library(readr)
source("Functions.R")

options(shiny.maxRequestSize=100*1024^2)  # Limits file upload size to 100 MB

ui <- navbarPage("DESeq2",  
   theme = shinytheme("sandstone"),
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 1----------------------------------#
#-----------------------------File Upload-------------------------------#
#---------------------------Widget Count: 3-----------------------------#
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
            
            HTML("<h1>View Data</h1>")
            
            ,
            

            #-----------------input----------------#
            #-------raw_counts_matrix_Preview------#
            #--------------------------------------#
            # A preview settings select option to change size of the displayed
            # raw_counts matrix returned from the Set Clean Counts function.
            #                       Object Name      
            FilePreviewSize("raw_counts_matrix_Preview", 
            #                       Widget Message                
                            "Counts Matrix Preview")
            
            ,

            #-----------------input----------------#
            #--------geneID_geneName_Preview-------#
            #--------------------------------------#
            # A preview settings select option to change the size of the displayed
            # dataframe containing gene ids and matching gene names
            #                       Object Name      
            FilePreviewSize("geneID_geneName_Preview", 
            #                       Widget Message                
                            "Gene Data Preview")
          
          ), ##XX##~~~Side Panel End~~~##XX##
          
#_____________________________Main Panel________________________________#
          
          mainPanel(
            
            #Text output for any error will appear at the top
            textOutput('errorMessagesPG1.0')
            
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

            #----------------output----------------#
            #-----geneID_geneName_PreviewTable-----#
            #--------------------------------------#
            #A Data Table preview of all gene ids and their gene names
              #Expected Format:
                #axis.names|    gene    | gene_name
                #   1      | Gene_ID_1  |   Actin    
                #   2      | Gene_ID_2  |   COOR4_7    
            DTOutput('geneID_geneName_PreviewTable')

            
            
            #----------------output----------------#
            #----sample_conditions_PreviewTable----#
            #--------------------------------------#
            #A Data Table Preview of the samples and their respective conditions
              #Expected Format:
              # #axis.names|  Condition    
              #   T0.s1    |  T0  
              #   T96.s2   |  T96 
            #Condition column stored as a factored category.
            #TO-DO
            
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
   
   
     
            
