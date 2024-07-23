#Import Shiny and needed packages
library(shiny)
library(shinythemes)
library(DT)
library(plyr)
library(dplyr)
library(readr)
source("Functions.R")
library(DESeq2)
library(ggplot2)
library(ggplotify)
library(patchwork)
library(matrixStats)
library(S4Vectors)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(grid)
library(colourpicker)
library(RColorBrewer)

options(shiny.maxRequestSize=100*1024^2)  # Limits file upload size to 100 MB

ui <- navbarPage("DESeq2",  
   theme = shinytheme("cyborg"),
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 1----------------------------------#
#-----------------------------File Upload-------------------------------#
#---------------------------Widget Count: 2-----------------------------#
#-----------------------Expected output Count: 2------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   tabPanel("File Upload",
      fluidPage(
        sidebarLayout(
#______________________________Side Bar_________________________________#
          
          sidebarPanel(
            
            HTML("<p>(Optional If you still need to run the DEG process for your data)</p>"),
            HTML("<p>If you already have this done and a file with data skip this page</p>"),
            HTML("<h2>Upload Data</h1>"),
            
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
            HTML("<p>Meta Data and conditions is needed for analysis</p>"),
            HTML("<h2>Upload Meta Data</h1>"),
            
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
                 through your unfiltered data once uploaded."),
            
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
                 conditions in column 2."),
            
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
#---------------------------Widget Count: 3-----------------------------#
#-----------------------Expected output Count: 2------------------------#
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
            
            ,
            
            checkboxGroupInput('display_col', 'Values of Interest (may take longer to load if more selected)', 
                               c('Base Mean' = 'baseMean', 
                               'Log2 Fold Change' ='log2FoldChange',
                               'lfcSE' ='lfcSE',
                               'Stat' ='stat',	
                               'P-Value' ='pvalue'),
                               selected=c('log2FoldChange'),
                               inline=FALSE)
            ,
            
            uiOutput("Distribution_Histogram_ui")
            
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
    
    ,

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 3----------------------------------#
#---------------------------Normalized Counts---------------------------#
#----------------------------Widget Count: 1----------------------------#
#-----------------------Expected output Count: 2------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
tabPanel("Normalized Counts Transformation",
         fluidPage(
#_____________________________Main Panel________________________________#
             
             mainPanel( 
               
               #-----------------input----------------#
               #----------------pvalue----------------#
               #--------------------------------------#
               # A select input giving the choice to change the displayed genes off pvalue
               # not just resort lists but actually display more or less genes
               selectInput("pvaluePg3", "Select adjusted P value cutoff for following tables", 
                           choices=c(1.0 ,0.5, 0.05, 0.01, 0.001))
               
               ,
               
               HTML("<h3>Normalized Counts Data</h3>")
               
               ,
               
               #----------------output----------------#
               #--------vst_counts_PreviewTable-------#
               #--------------------------------------#
               #A ui element containing a datatable with the padj 
               #gene names and vst counts
               uiOutput('vst_counts_PreviewTable')
             
              )##XX##~~~Main Panel End~~~##XX##
             
             ,
             
             sidebarPanel(
               
               uiOutput('Extra_vst_count_info')
               
             )
         ) ##XX##~~~Fluid Page Closing Bracket~~~##XX##
       ) ##XX##~~~Tab Panel Closing Bracket~~~##XX##

      ,

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 4----------------------------------#
#-----------------------------PCA Plots---------------------------------#
#---------------------------Widget Count: 5-----------------------------#
#-----------------------Expected output Count: 1------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
tabPanel("Principle Component Plots",
         fluidPage(
#______________________________Side Bar_________________________________#
             
             sidebarPanel( width=2,
               
               HTML("<h3>Labels</h3>")    
               
               ,
               
               #-----------------input----------------#
               #----------------Title-----------------#
               #--------------------------------------#           
               textInput('title_pca_plot', 'Title', value = 'Title')
               
               ,
               
               #-----------------input----------------#
               #----------------subTitle--------------#
               #--------------------------------------#
               textInput('subtitle_pca_plot', 'Sub Title', value = 'Sub Title')
               
               
               ,
               
               #-----------------input----------------#
               #----------------caption---------------#
               #--------------------------------------#
               textAreaInput('caption_pca_plot', 'Caption', value = 'caption', width=200, rows=3)
               
             ), ##XX##~~~Side Panel End~~~##XX##
             
#_____________________________Main Panel________________________________#
             
             mainPanel( width=7,     
               
                #-----------------output---------------#
                #----------------pca_plot--------------#
                #--------------------------------------#
                ##Principle component Analysis plot
                uiOutput('principle_component_plots_ui')
               
               
             ),##XX##~~~Main Panel End~~~##XX##
             
#______________________________Side Bar_________________________________#
             
             sidebarPanel( width=3,
              
               HTML("<h3>Data</h3>")              
               
               ,
               
               #-----------------input----------------#
               #----# of things included in pca-------#
               #--------------------------------------#
               #The number of things with top variance included in pca
               uiOutput('change_n_pca_plot')
               
               ,
               
               #-----------------input----------------#
               #-------columns shown in pca-----------#
               #--------------------------------------#
               # a Check box group of each column in the metaData()
               # Selected on makes it included in the pca up to three
               uiOutput('change_mColumns_for_pca')
               
             ), ##XX##~~~Side Panel End~~~##XX##
             
         ) ##XX##~~~Fluid Page Closing Bracket~~~##XX##
       ) ##XX##~~~Tab Panel Closing Bracket~~~##XX##

      ,

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 5----------------------------------#
#-------------------------Correlation Analysis--------------------------#
#---------------------------Widget Count: 0-----------------------------#
#-----------------------Expected output Count: 0------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
tabPanel("Correlation Anaylsis",
         fluidPage(
#______________________________Side Bar_________________________________#
           
           sidebarPanel( width=2,
             
             HTML("<h3>Labels</h3>")    
             
             ,
             
             textInput('title_heatmap_plot', 'Title', value = 'Title')
             
             ,
             
             textInput('subtitle_heatmap_plot', 'Sub Title', value = 'Sub Title')
             
             
             ,
             
             textAreaInput('caption_heatmap_plot', 'Caption', value = 'caption', width=200, rows=3)
             
           ), ##XX##~~~Side Panel End~~~##XX##
           
#_____________________________Main Panel________________________________#
           
           mainPanel( width=7,     
                      
              uiOutput('heatmap_plots_ui')
                      
           ),##XX##~~~Main Panel End~~~##XX##
           
#______________________________Side Bar_________________________________#
           
           sidebarPanel( width=3,
              
              uiOutput('heatmap_annotations')
              #minimalize feature           
              
           ), ##XX##~~~Side Panel End~~~##XX##
           
         ) ##XX##~~~Fluid Page Closing Bracket~~~##XX## 
      ) ##XX##~~~Tab Panel Closing Bracket~~~##XX##


)#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X
   
   
     
            
