#Import Shiny and needed packages
library(shiny)
library(shinythemes)
library(DT)
library(plyr)
library(dplyr)
library(readr)
source("Functions.R")
library(crayon)
library(DESeq2)

options(shiny.maxRequestSize=100*1024^2)  # Limits file upload size to 100 MB

ui <- navbarPage("DESeq2",  
   theme = shinytheme("cyborg"),
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 1----------------------------------#
#-----------------------------File Upload-------------------------------#
#---------------------------Widget Count: 2-----------------------------#
#-----------------------Expected output Count: 2------------------------#
#-----------------------Possible Errors Count: 3------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   tabPanel("File Upload",
      fluidPage(
        sidebarLayout(
#______________________________Side Bar_________________________________#
          
          sidebarPanel(
            
            HTML("<p>(Optional If you still need to run the DEG process for your data)</p>"),
            HTML("<p>If you already have this done and a file with data skip this page</p>")
            
            ,
            
            HTML("<h2>Upload Data</h1>")
            
            ,
            
            #-----------------input----------------#
            #---merged_gene_counts_uploaded_file---#
            #--------------------------------------#
            # An upload file widget for the merged gene counts tsv,csv
            #          Object Name    
            fileInput("merged_gene_counts_uploaded_file", 
            #          Widget Message 
                       "Salmon Merged Counts .tsv/.csv file")
            
            ,
           
            HTML("<p>(Not Optional)</p>"),
            HTML("<p>Meta Data and conditions is needed for analysis</p>")
            
            ,
            
            HTML("<h2>Upload Meta Data</h1>")
            
            ,
            
            #-----------------input----------------#
            #--meta_data_conditions_uploaded_file--#
            #--------------------------------------#
            # An upload file widget for the metadata .csv file
            #          Object Name    
            fileInput("meta_data_conditions_uploaded_file", 
                       #  Widget Message 
                       "Conditions")
            
            ,
            
            #-----------------input----------------#
            #---------------run_DESeq2-------------#
            #--------------------------------------#
            # Action button that starts DSeq2 disabling the upload of
            # a more advanced dataset
            actionButton("run_DESeq2", "Run Differential Expression")
            
            
          ), ##XX##~~~Side Panel End~~~##XX##
          
#_____________________________Main Panel________________________________#
          
          mainPanel(
            
            HTML("<h1>Raw Counts Table</h1>"),
            HTML("<p1>Upload your salmon merged counts matrix on the left. Your 
                 uploaded file must be a .tsv file. Below you can explore and search
                 through your unfiltered data once uploaded.")
            
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
#-----------------------------DESeq2 Info-------------------------------#
#---------------------------Widget Count: 1-----------------------------#
#-----------------------Expected output Count: 1------------------------#
#-----------------------Possible Errors Count: 0------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    tabPanel("DEG Analysis",
      fluidPage(
        sidebarLayout(
#______________________________Side Bar_________________________________#
        
          sidebarPanel(
            
            #-----------------input----------------#
            #----------DEG_analysis_data-----------#
            #--------------------------------------#
            # An upload file widget for the DEG data
            # If the raw counts and metadata were uploaded to perform DESeq2
            # This will not be visable
            uiOutput('pg2_display_upload_DEG_data')
            
            ,
            
            #-----------------input----------------#
            #----------------pvalue----------------#
            #--------------------------------------#
            # A select input giving the choice to change the displayed genes off pvalue
            # not just resort lists but actually display more or less genes
            selectInput("pvalue", "Select adjusted P value cutoff for following tables", 
                        choices=c(1.0 ,0.5, 0.05, 0.01, 0.001))
            
            
          ), ##XX##~~~Side Panel End~~~##XX##
          
#_____________________________Main Panel________________________________#

          mainPanel(      
            
            #-----------------input----------------#
            #----------------cutOffs----------------#
            #--------------------------------------#
            # A slider to adjust what cutoffs will be shown for the tables
            sliderInput('cutOffs', 'Cut-off values for the genes being displayed as a certain direction of diffrential expression', min=-12, max=12, step=0.10, value = c(-0.5, 0.5), width='100%')
            
            ,
            
            #----------------output----------------#
            #---DESeq_Expression_Analysis_Tables---#
            #--------------------------------------#
            #A single Div containing three tables one for up-regulated, down-reg, 
            #and no regulated genes
            uiOutput('DESeq_Expression_Analysis_Tables')
            
          )##XX##~~~Main Panel End~~~##XX##

        ) ##XX##~~~Side Bar Layout Closing Bracket~~~##XX##
      ) ##XX##~~~Fluid Page Closing Bracket~~~##XX##
    ) ##XX##~~~Tab Panel Closing Bracket~~~##XX##
                 
)#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X
   
   
     
            
