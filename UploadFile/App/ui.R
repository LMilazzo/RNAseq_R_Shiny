#Import Shiny and needed packages
library(shiny)
library(shinythemes)
library(DT)
library(plyr)
library(dplyr)
library(readr)
source("Functions.R")


options(shiny.maxRequestSize=100*1024^2)  # Limits file upload size to 100 MB

#_________________________Create Navigation Bar System__________________________

ui <- navbarPage("Application",
                 
#_________________________Set Theme To sandstone________________________________
   
   theme = shinytheme("sandstone"),
   
##############################Page <- File Upload###############################
   
   tabPanel("File Upload",
      
      fluidPage(
        
        #Sidebar layout displays widgets on a left sidebar and output object in 
        #the main right panel
        sidebarLayout(
        
##############################Side Bar##########################################
      
          
          #The sidebar panel contains widgets for this page
          
          #Expected Widgets for this sidebar
          
          # merged_gene_counts_uploaded_file
            # An upload file widget for the merged gene counts tsv,csv
          # raw_counts_matrix_Preview
            # A preview settings select option to change size of the displayed
            # raw_counts matrix returned from the Set Clean Counts function.

          sidebarPanel(
            
            # merged_gene_counts_uploaded_file
            #          Object Name    
            FileUpload("merged_gene_counts_uploaded_file", 
            #          Widget Message 
                       "Upload Counts")
            
            ,
            
            # raw_counts_matrix_Preview
            #                       Object Name      
            FilePreviewSize("raw_counts_matrix_Preview", 
            #                       Widget Message                
                            "Counts Matrix Preview")
            
            ,
            
            # geneID_geneName_Preview
            #                       Object Name      
            FilePreviewSize("geneID_geneName_Preview", 
            #                       Widget Message                
                            "Gene Data Preview")
          
            ),
          
##########################Main Panel############################################
          
          #The main panel contains output objects for this page
          
          #Expected possible output objects:
          
          #output$raw_counts_PreviewTable
            #A Data Table preview of what is inside the uploaded counts table 
            #after formating correctly
              #Expected Format:
                #axis.names| Sample1 | Sample 2 | Sample 3
                #Gene_ID_1 |    x    |    x     |    x
                #Gene_ID_2 |    x    |    x     |    x
            
          #output$geneID_geneName_PreviewTable
            #A Data Table preview of all gene ids and their gene names
              #Expected Format:
                #axis.names|    gene    | gene_name
                #   1      | Gene_ID_1  |   Actin    
                #   2      | Gene_ID_2  |   COOR4_7    
            
          #output$sample_conditions_PreviewTable
            #A Data Table Preview of the samples and their respective conditions
              #Expected Format:
              # #axis.names|  Condition    
              #   T0.s1    |  T0  
              #   T96.s2   |  T96 
            #Condition column stored as a factored category.
          
          mainPanel(
            
            DTOutput('raw_counts_PreviewTable')
            
            ,
            
            DTOutput('geneID_geneName_PreviewTable')
            
          )
          
        )
      )
    )
)
   
   
     
            