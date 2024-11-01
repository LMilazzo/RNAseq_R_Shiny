#Import Shiny and needed packages
#----
#----
library(shiny)
library(shinythemes)
library(DT)
library(dplyr)
library(readr)
source("Functions.R")
library(DESeq2)
library(ggplot2)
library(ggplotify)
library(patchwork)
library(matrixStats)
library(SummarizedExperiment)
library(circlize)
library(ComplexHeatmap)
library(colourpicker)
library(ggbeeswarm)
library(ggrepel)
library(BioVis)
library(pathfindR)
library(tidyr)
library(shinyjs)
library(plotly)
library(htmlwidgets)
#----

options(shiny.maxRequestSize=100*1024^2)  # Limits file upload size to 100 MB

ui <- navbarPage("DESeq2",  
   theme = shinytheme("cyborg"),
   useShinyjs(),
#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 1----------------------------------#
#-----------------------------File Upload-------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
   tabPanel("File Upload",
      fluidPage(
        sidebarLayout(
#______________________________Side Bar_________________________________#
          
          sidebarPanel(
            
            HTML("<p>Sample counts matrix:</p>"),
            HTML("<p>Supported file types:  .csv   .tsv"),
            HTML("<p>Format: ( gene_id, gene_name, .samples... ) </p>"),
            HTML("<p>gene_id must be unique for all genes</p>"),
            
            #-----------------input----------------#
            #---merged_gene_counts_uploaded_file---#
            #--------------------------------------#
            # An upload file widget for the merged gene counts tsv,csv
            fileInput("merged_gene_counts_uploaded_file", 
                      "Counts Matrix"),

            HTML("<p>Sample conditions / meta data : </p>"),
            HTML("<p>Supported file types:  .csv   .tsv"),
            HTML("<p>Format: ( sample, condions... ) </p>"),
            HTML("<p>The sample column should match the list of samples found in the counts matrix</p>"),
            
            #-----------------input----------------#
            #--meta_data_conditions_uploaded_file--#
            #--------------------------------------#
            # An upload file widget for the metadata .csv file 
            fileInput("meta_data_conditions_uploaded_file", 
                      "Sample Conditions"),
            
            #-----------------input----------------#
            #---------------run_DESeq2-------------#
            #--------------------------------------#
            # Action button that starts DSeq2 disabling the upload of
            # a post experiment dataset
            actionButton("run_DESeq2", "Run Differential Expression")
            
          ), ##XX##~~~Side Panel End~~~##XX##
          
#_____________________________Main Panel________________________________#
          
          mainPanel(
            
            HTML("<h5>Counts Matrix</h5>"),
            # HTML("<p1>Upload your salmon merged counts matrix on the left. Your 
            #      uploaded file must be a .tsv file. Below you can explore and search
            #      through your unfiltered data once uploaded."),
            
            uiOutput('raw_counts_PreviewTable'),
            
            HTML("<h5>Sample Conditions</h5>"),
            # HTML("<p1>Upload a .csv file with the samples in column 1 and their 
            #      conditions in column 2."),
            
            uiOutput('sample_conditions_PreviewTable')
            
          )##XX##~~~Main Panel End~~~##XX##
          
        )##XX##~~~Side Bar Layout Closing Bracket~~~##XX##
      ) ##XX##~~~Fluid Page Closing Bracket~~~##XX##
    ), ##XX##~~~Tab Panel Closing Bracket~~~##XX##
#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 2----------------------------------#
#-----------------------------DESeq Info--------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
    tabPanel("DEG Analysis",
      fluidPage(
        sidebarLayout(
#______________________________Side Bar_________________________________#
        
          sidebarPanel(
            
            uiOutput('deg_upload'),
            
            #-----------------input----------------#
            #----------------pvalue----------------#
            #--------------------------------------#
            # A select input giving the choice to change the displayed genes off pvalue
            # not just resort lists but actually display more or less genes
            selectInput("pvaluePg2", "Select adjusted P value cutoff for following tables", 
                        choices=c(1.0 ,0.5, 0.05, 0.01, 0.001)),
            
            #-----------------input----------------#
            #------------displayed vals------------#
            #--------------------------------------#
            #Determines what dataset values are displayed in the tables
            checkboxGroupInput('display_col', 
                               'Values of Interest (may take longer to load if more selected)', 
                               c('Base Mean' = 'baseMean', 
                               'Log2 Fold Change' ='log2FoldChange',
                               'lfcSE' ='lfcSE',
                               'Stat' ='stat',	
                               'P-Value' ='pvalue'),
                               selected=c('log2FoldChange'),
                               inline=FALSE
                              ),
            
            uiOutput("Distribution_Histogram_ui")
            
          ), ##XX##~~~Side Panel End~~~##XX##
          
#_____________________________Main Panel________________________________#

          mainPanel(
            
            #-----------------input----------------#
            #----------------cutOffs----------------#
            #--------------------------------------#
            # A slider to adjust what cutoffs will be shown for the tables
            sliderInput('cutOffs', 
            'Cut-off values for the genes being displayed as a certain direction of diffrential expression', 
            min=-12, max=12, step=0.10, value = c(-0.5, 0.5), width='100%'),
            
            uiOutput('expression_tables')
            
          )##XX##~~~Main Panel End~~~##XX##

        ) ##XX##~~~Side Bar Layout Closing Bracket~~~##XX##
      ) ##XX##~~~Fluid Page Closing Bracket~~~##XX##
    ), ##XX##~~~Tab Panel Closing Bracket~~~##XX##
#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 3----------------------------------#
#---------------------------Normalized Counts---------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
  tabPanel("Normalized Counts Preview",
    fluidPage(
#_____________________________Main Panel________________________________#
             
        mainPanel(
          
          #-----------------input----------------#
          #----------------pvalue----------------#
          #--------------------------------------#
          # A select input giving the choice to change the displayed genes off pvalue
          # not just resort lists but actually display more or less genes
          selectInput("pvaluePg3", "Select adjusted P value cutoff for following tables", 
                      choices=c(1.0 ,0.5, 0.05, 0.01, 0.001)),
          
          HTML("<h3>Normalized Counts Data</h3>"),
          
          uiOutput('normalized_counts_PreviewTable')
          
        ), ##XX##~~~Main Panel End~~~##XX##
             
#______________________________Side Bar_________________________________#

        sidebarPanel(
          
          uiOutput('Extra_normalized_count_info')
          
        ) ##XX##~~~Side Panel End~~~##XX##
      ) ##XX##~~~Fluid Page Closing Bracket~~~##XX##
    ), ##XX##~~~Tab Panel Closing Bracket~~~##XX##
#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 4----------------------------------#
#-----------------------------PCA Plots---------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
  tabPanel("Principle Component Plots",
    fluidPage(
#______________________________Side Bar_________________________________#
             
      sidebarPanel( width=2,
                    
        HTML("<h3>Labels</h3>"),
        
        #-----------------input----------------#
        #---------title, subtitle, caption-----#
        #--------------------------------------#           
        textInput('title_pca_plot', 'Title', value = 'Title'),
        textInput('subtitle_pca_plot', 'Sub Title', value = 'Sub Title'),
        textAreaInput('caption_pca_plot', 'Caption', value = 'caption', width=200, rows=3)
        
      ), ##XX##~~~Side Panel End~~~##XX##
             
#_____________________________Main Panel________________________________#
             
      mainPanel( width=7,     
        
        uiOutput('principle_component_plots_ui')
        
      ), ##XX##~~~Main Panel End~~~##XX##
             
#______________________________Side Bar_________________________________#
             
      sidebarPanel( width=3,
                    
        HTML("<h3>Data</h3>"),
    
        uiOutput('pca_n_counts'),
        
        uiOutput('pca_metadata')
        
      ), ##XX##~~~Side Panel End~~~##XX##
             
    ) ##XX##~~~Fluid Page Closing Bracket~~~##XX##
  ), ##XX##~~~Tab Panel Closing Bracket~~~##XX##
#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 5----------------------------------#
#-------------------------Correlation Analysis--------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
  tabPanel("Correlation Anaylsis",
    fluidPage(
#______________________________Side Bar_________________________________#
           
      sidebarPanel( width=2,
                    
        HTML("<h3>Labels</h3>"),
        
        #-----------------input----------------#
        #---------title, subtitle, caption-----#
        #--------------------------------------#                 
        textInput('title_heatmap_plot', 'Title', value = 'Title'),
        textInput('subtitle_heatmap_plot', 'Sub Title', value = 'Sub Title'),
        textAreaInput('caption_heatmap_plot', 'Caption', value = 'caption', width=200, rows=3)
          
      ), ##XX##~~~Side Panel End~~~##XX##
           
#_____________________________Main Panel________________________________#
           
      mainPanel( width=7,     
                
        uiOutput('heatmap_plots_ui')
                
      ),##XX##~~~Main Panel End~~~##XX##
           
#______________________________Side Bar_________________________________#
           
      sidebarPanel( width=3,
                    
        uiOutput('heatmap_annotations'),
        
        #-----------------input----------------#
        #----------------color-----------------#
        #--------------------------------------#
        # A select input giving the choice to change the displayed genes off pvalue
        # not just resort lists but actually display more or less genes
        selectInput("heat_body_color", "Select color palette", 
                    choices=c("blue-red", "REV-Rainbow", "white-black")),
      
        uiOutput('color_pickers_heat')
        
      ), ##XX##~~~Side Panel End~~~##XX##
           
    ) ##XX##~~~Fluid Page Closing Bracket~~~##XX## 
  ), ##XX##~~~Tab Panel Closing Bracket~~~##XX##
#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 6----------------------------------#
#-----------------------------Gene Search-------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
  tabPanel("Gene Counts",
    fluidPage(
#______________________________Side Bar_________________________________#
      sidebarPanel(width=4,
                   
        #-----------------input----------------#
        #----------------pvalue----------------#
        #--------------------------------------#
        # A select input giving the choice to change the displayed genes off pvalue
        # not just resort lists but actually display more or less genes
        selectInput("pvaluePg6", "Select adjusted P value cutoff for following plots", 
                    choices=c(1.0 ,0.5, 0.05, 0.01, 0.001)),
        
        uiOutput('leftbar_gene_plots')
        
      ), ##XX##~~~Side Panel End~~~##XX##
#_____________________________Main Panel________________________________#
      mainPanel(width=6,
                
          uiOutput('gene_count_search')
          
      ), ##XX##~~~Main Panel End~~~##XX##

#______________________________Side Bar_________________________________#
      sidebarPanel(width=2,
                   
        uiOutput("gene_name_list")

      ) ##XX##~~~Side Panel End~~~##XX##
           
    ) ##XX##~~~Fluid Page Closing Bracket~~~##XX## 
  ) ,  ##XX##~~~Tab Panel Closing Bracket~~~##XX##
#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 7----------------------------------#
#-----------------------------Volcano plot------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
  tabPanel("Volcano Plot",
    fluidPage(
      
#_____________________________Side Bar__________________________________#
      
      sidebarPanel(width=2,
                   
         #-----------------input----------------#
         #---------title, subtitle, caption-----#
         #--------------------------------------#                
         textInput('title_volc_plot', 'Title', value = 'Title'),
         textInput('subtitle_volc_plot', 'Sub Title', value = 'Sub Title'),
         textAreaInput('caption_volc_plot', 'Caption', value = 'caption', width=200, rows=3)         
        
      ),##XX##~~~Side Panel End~~~##XX##
#_____________________________Main Panel________________________________#
      
      mainPanel(width=7,

        #-----------------input----------------#
        #----------------cutOffs---------------#
        #--------------------------------------#
        # A slider to adjust what cutoffs for the vert lines will be
        sliderInput('volcano_cutoffs', 
        'Cut of values for diffrential expression (vertivle lines on Volcano Plot)', 
        min=-12, max=12,step=0.01, value = c(-1, 1), width='100%'),
        
        uiOutput('volcano_plot_ui')
        
        
      ), ##XX##~~~Main Panel End~~~##XX##
      
#_____________________________Side Bar__________________________________#
      
      sidebarPanel(width=3,
        
        #-----------------input----------------#
        #----------------pvalue----------------#
        #--------------------------------------#
        # A select input giving the choice to change the displayed genes off pvalue
        # not just resort lists but actually display more or less genes
        selectInput("pvaluePg7", "Select adjusted P value cutoff for following plot", 
                    choices=c(1.0 ,0.5, 0.05, 0.01, 0.001), selected=0.05),
        
        #-----------------input----------------#
        #--------------population--------------#
        #--------------------------------------#
        # A slider to adjust the population of the graph
        sliderInput('volcano_pop', 
                    'Change the population of the Volcano Plot', 
                    min=0, max=1, value = 0.3, width='100%', ticks=FALSE),
        
        #-----------------input----------------#
        #-------------label density------------#
        #--------------------------------------#
        # A slider to adjust the label population of the graph
        sliderInput('volcano_lab_density', 
                    'Change the label density', 
                    min=0, max=1, value = 0.3, width='100%', ticks=FALSE),
        
        #-----------------input----------------#
        #--------------gene_to_search----------#
        #--------------------------------------#
        # A text input later converted to vector of genes to highlight
        textInput('volc_search',
                  'Genes to search/highlight (can be comma seperated list)',
                  value=NULL)

      )##XX##~~~Side Panel End~~~##XX##

    )##XX##~~~Fluid Page Closing Bracket~~~##XX## 
  ), ##XX##~~~Tab Panel Closing Bracket~~~##XX##
#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 8----------------------------------#
#--------------------------Pathway Analysis-----------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
  tabPanel("Pathway Analysis",
   fluidPage(
     
#_____________________________Side Bar__________________________________#
      sidebarPanel(width = 3,
        
        HTML('<h3>Pathway Analysis</h3><br><br><p>Continue Experiment</p>'),
        actionButton('Run_pathfinder', 'Run'),
        
        HTML('<br><br><p>Upload a past experiment<p>'),
        actionButton('Run_pathfinder_new_data', 'Review old experiment')
        
        
      ),##XX##~~~Side Panel End~~~##XX##
#_____________________________Main Panel________________________________#
      mainPanel(width = 9,
        
        
        uiOutput('pathfinderPreview')
                
                
      )##XX##~~~Main Panel End~~~##XX##

   )##XX##~~~Fluid Page Closing Bracket~~~##XX## 
  ),##XX##~~~Tab Panel Closing Bracket~~~##XX##
#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 9----------------------------------#
#--------------------------Pathway Enrichment---------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
  tabPanel("Term Enrichment",
    fluidPage(
      
      uiOutput('enrichment_clusters_shown'),
#_____________________________Side Bar__________________________________#
      sidebarPanel(width = 3,
      
        #-----------------input----------------#
        #--------------gene_to_search----------#
        #--------------------------------------#
        # A text input later converted to vector of genes to highlight
        textInput('enrichment_genes',
                  'Genes to include in the enrichment chart',
                  value=NULL),
        
        #-----------------input----------------#
        #--------------gene_to_search----------#
        #--------------------------------------#
        # A text input later converted to vector of genes to highlight
        textInput('enrichment_paths',
                  'Pathway terms to be included in the enrichment chart',
                  value=NULL),
        
        #-----------------input----------------#
        #--------------gene_to_search----------#
        #--------------------------------------#
        # A text input later converted to vector of genes to highlight
        textInput('enrichment_clusters',
                  'Specific clusters to view in the enrichment chart',
                  value=NULL),
        
        uiOutput('pathwaysDT9'),
        
        uiOutput('genes_in_paths_DT9')
        
      
      ), ##XX##~~~Side Panel End~~~##XX##

#_____________________________Main Panel________________________________#
      mainPanel(width = 9,
                
        uiOutput('enrichmentUI')
        
      )##XX##~~~Main Panel End~~~##XX##
    )##XX##~~~Fluid Page Closing Bracket~~~##XX## 
  ), ##XX##~~~Tab Panel Closing Bracket~~~##XX##
#----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 10---------------------------------#
#--------------------------Pathway Heatmaps-----------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
  tabPanel("Term Gene Heatmap",
   fluidPage(
#_____________________________Side Bar__________________________________#
     sidebarPanel(width = 3,
      
        uiOutput('open_in_new_tab10'), 
        
        #-----------------input----------------#
        #--------------gene_to_search----------#
        #--------------------------------------#
        # A text input later converted to vector of genes to highlight
        textInput('pathway_heatmap_genes',
                  'Genes to include in the heatmap chart',
                  value=NULL),
        
        #-----------------input----------------#
        #--------------gene_to_search----------#
        #--------------------------------------#
        # A text input later converted to vector of genes to highlight
        textInput('heatmap_paths',
                  'Pathway terms to be included in the heatmap chart',
                  value=NULL),
    
        uiOutput('pathwaysDT10'),
        
        uiOutput('genes_in_paths_DT10')
                  
                  
     ), ##XX##~~~Side Panel End~~~##XX##
     
#_____________________________Main Panel________________________________#
     mainPanel(width = 9,
               
        uiOutput('pathway_heatmap')
               
     )
   )##XX##~~~Fluid Page Closing Bracket~~~##XX## 
  ),##XX##~~~Tab Panel Closing Bracket~~~##XX##


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------PAGE 10---------------------------------#
#--------------------------Pathway Heatmaps-----------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#----
tabPanel("example",
   fluidPage(
     
     uiOutput('cases_select_box'),
    
    #_____________________________Side Bar__________________________________#
    sidebarPanel(width = 4,
        
        uiOutput('sample_conditions_PreviewTable11')
  
     ), ##XX##~~~Side Panel End~~~##XX##
    
    #_____________________________Main Panel________________________________#
    mainPanel(width = 8,
    
       uiOutput('bsexample')
  
    )##XX##~~~Main Panel End~~~##XX##
   )##XX##~~~Fluid Page Closing Bracket~~~##XX## 
  )##XX##~~~Tab Panel Closing Bracket~~~##XX##
)#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X




