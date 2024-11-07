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
library(grid)
#----

options(shiny.maxRequestSize=100*1024^2)  # Limits file upload size to 100 MB

ui <- tagList(tags$style(HTML("
    /* Custom CSS for DataTable */
    .dataTables_wrapper .dataTables_length, 
    .dataTables_wrapper .dataTables_filter, 
    .dataTables_wrapper .dataTables_info, 
    .dataTables_wrapper .dataTables_paginate,
    table.dataTable {
      color: #ffffff;  /* Change text color */
    }
    table.dataTable thead {
      background-color: #444444; /* Change header background */
      color: #ffffff; /* Change header text color */
    }
    table.dataTable tbody {
      background-color: #333333; /* Change body background */
    }
  ")),
  useShinyjs(),
  
  navbarPage("", selected = "Diffrentially Expressed Genes", 
    theme = shinytheme("slate"),
#----
  
# Page 1__________________________________________________________________#----

  tabPanel("Help"),

# Page 2__________________________________________________________________#----

   tabPanel("Diffrentially Expressed Genes",
      fluidPage(
        
        fluidRow(
          column(1,
            div(
              actionButton('reload_app', 
                           'Reload Application', 
                           style = "background-color: #FF6F61; color: #FF6F61; border-color: #FF6F61;")
            )
          ),
        
          column(6, offset = 1,
        
            div(id = 'Page1_Upload_Options',
              fluidRow(
        
                column(2,
                  actionButton('start_new_experiment', 'Start New Experiment',
                               style = "background-color: #4CAF50; color: #4CAF50; border-color: #4CAF50;"
                               )
                ),
        
                column(2, offset = 2,
                  actionButton('start_old_experiment', 'Review A Previous Experiment')
                )
              )
            ),
        
            hidden(
              actionButton("run_DESeq2", "Run Differential Expression",
                           style = "background-color: #4CAF50; color: #4CAF50; border-color: #4CAF50;"
                           )
            )
          )
        ),
        
            
# Tabset 1 Diffrential Expression Analysis________________________________#
          br(),
          mainPanel(width = 12,
            
            uiOutput('about_DESeq2'),
                    
            div(id = 'preRun_Data_Preview',
                uiOutput('preRun_preview1'),
                uiOutput('preRun_preview2')
            ),
                    
            hidden(
              div(id="TabSet1_Diffrential_Expression_Analysis",
                tabsetPanel(

#TABS 
# Preview_________________________________________________________________#----
tabPanel("Data Preview", 
         
  sidebarPanel(width = 6,
         
   HTML("<h3>Counts Matrix</h3>"),
   
   uiOutput('raw_counts_PreviewTable'),
  
  ),
  
  sidebarPanel(width = 6,
  
    HTML("<h3>Normalized Counts</h3>"),
    
    uiOutput('normalized_counts_PreviewTable')
               
  ),

 HTML("<h3>Sample Conditions</h3>"),
 
 uiOutput('sample_conditions_PreviewTable'),
   
   
),
# DEG analysis____________________________________________________________#----
tabPanel("Diffrential Analysis",
  fluidPage(
    sidebarLayout(
    
      sidebarPanel(
      
        #Pvalue box
        selectInput("pvaluePg2", 
        "Select adjusted P value cutoff for following tables", 
        choices=c(1.0 ,0.5, 0.05, 0.01, 0.001)),
        
        #Stat check box
        checkboxGroupInput('display_col', 
              'Values of Interest (may take longer to load if more selected)', 
              c('Base Mean' = 'baseMean', 
                'Log2 Fold Change' ='log2FoldChange',
                'lfcSE' ='lfcSE',
                'Stat' ='stat',	
                'P-Value' ='pvalue'),
              selected=c('log2FoldChange'),
              inline=FALSE),
        
        #Distribution
        uiOutput("Distribution_Histogram_ui")
      
      ),
      
      mainPanel(
      
        #Cutoff Slider
        sliderInput('cutOffs', 
        'Cut-off values for the genes being displayed as a certain direction of diffrential expression', 
        min=-12, max=12, step=0.10, value = c(-0.5, 0.5), width='100%'),
        
        #Data tables
        uiOutput('expression_tables')
        
      )
    )
  )
),
# Principle Component Analysis____________________________________________#----
tabPanel("Principle Component Plots",
  fluidPage(
  
    sidebarPanel(width=3,
                 
      HTML("<h3>Data</h3>"),
     
      uiOutput('pca_n_counts'),
     
      uiOutput('pca_metadata'),
     
      HTML("<h3>Labels</h3>"),
      
      textInput('title_pca_plot', 'Title', value = 'Principle Component Analysis of Sample Varience'),
      textInput('subtitle_pca_plot', 'Sub Title', value = 'Colored by sample traits'),
      textAreaInput('caption_pca_plot', '', value = '', width=200, rows=3)
             
    ),
       
    mainPanel(width=9,     
          
      uiOutput('principle_component_plots_ui')
          
    )
    
  )
),
# Correlation Analysis____________________________________________________#----
tabPanel("Correlation Anaylsis",
  fluidPage(
  
    sidebarPanel( width=3,
                  
      uiOutput('heatmap_annotations'),
      
      selectInput("heat_body_color", "Select color palette", 
                  choices=c("blue-red", "REV-Rainbow", "white-black")),
      
      uiOutput('color_pickers_heat'),
      
      HTML("<h3>Labels</h3>"),
    
      textInput('title_heatmap_plot', 'Title', value = 'Sample Correlation'),
      textInput('subtitle_heatmap_plot', 'Sub Title', value = ''),
      textAreaInput('caption_heatmap_plot', 'Caption', value = '', width=200, rows=3)
      
    ),
    
    mainPanel( width=9,     
    
      uiOutput('heatmap_plots_ui')
    
    ),

  ) 
),
# Gene Counts_____________________________________________________________#----
tabPanel("Gene Counts",
  fluidPage(
  
    sidebarPanel(width=4,
                
      selectInput("pvaluePg6", "Select adjusted P value cutoff for following plots", 
                  choices=c(1.0 ,0.5, 0.05, 0.01, 0.001)),
      
      uiOutput('leftbar_gene_plots')
                
    ), 
    
    mainPanel(width=6,
             
      uiOutput('gene_count_search')
             
    ), 
    
    sidebarPanel(width=2,
                
      uiOutput("gene_name_list")
                
    )
  )
), 
# Volcano Plot____________________________________________________________#----
tabPanel("Volcano Plot",
  fluidPage(
  
    sidebarPanel(width=3,
    
     sliderInput('volcano_cutoffs', 
                 'Cut of values for diffrential expression (verticle lines on Volcano Plot)', 
                 min=-12, max=12,step=0.01, value = c(-1, 1), width='100%'),
     
     selectInput("pvaluePg7", "Select adjusted P value cutoff for following plot", 
                 choices=c(1.0 ,0.5, 0.05, 0.01, 0.001), selected=0.05),
     
     sliderInput('volcano_pop', 
                 'Change the population of the Volcano Plot', 
                 min=0, max=1, value = 0.3, width='100%', ticks=FALSE),
     
     sliderInput('volcano_lab_density', 
                 'Change the label density', 
                 min=0, max=1, value = 0.3, width='100%', ticks=FALSE),
     
     textInput('volc_search',
               'Genes to search/highlight (can be comma seperated list)',
               value=NULL),
                 
                 
      textInput('title_volc_plot', 'Title', value = 'Gene Expression'),
      textInput('subtitle_volc_plot', 'Sub Title', value = ''),
      textAreaInput('caption_volc_plot', 'Caption', value = '', width=200, rows=3)         
    
    ),
    
    mainPanel(width=9,
      
      uiOutput('volcano_plot_ui')
    
    ), 
    
  ) 
)
#----
              ) #Tabset 1
            )
          )
        )
      )
    ),

# Page 3__________________________________________________________________#----
  tabPanel("Pathway Analysis",
    fluidPage(
        
      fluidRow(
        column(1,
           div(
             actionButton('reload_app', 
                          'Reload Application',
                          style = "background-color: #FF6F61; color: #FF6F61; border-color: #FF6F61;")
           )
        ),
        div(id="pathfinder_option_buttons",
          column(9, offset = 1,
            fluidRow(
              column(4,
                actionButton('Run_pathfinder', 'Continue Experiment With Pathway Analysis',
                             style = "background-color: #4CAF50; color: #4CAF50; border-color: #4CAF50;"
                             ),
              ),
              column(3, offset = 1,
                actionButton('review_pathfinder_new_data', 'Review A Pathway Analysis Experiment')
              )
            )
          )
        )
      ),
# Tabset 2 Pathway Analysis ________________________________#   
      br(),
      mainPanel(width = 12,
        
        uiOutput('about_pathfinder'),        
                
        hidden(
          div(id="TabSet2_Pathway_Analysis",
            tabsetPanel(

#TABS
# Preview_________________________________________________________________#----
tabPanel("Pathway Analysis Table",
  fluidPage(
    mainPanel(
      
      uiOutput('pathfinderPreview')
      
    )
  )
),
# Pathway Enrichment Plot_________________________________________________#----
tabPanel("Term Enrichment",
  fluidPage(
     
    uiOutput('enrichment_clusters_shown'),
    
    sidebarPanel(width = 3,
                    
      textInput('enrichment_genes',
                'Genes to include in the enrichment chart',
                value=NULL),
      
      textInput('enrichment_paths',
                'Pathway terms to be included in the enrichment chart',
                value=NULL),
      
      textInput('enrichment_clusters',
                'Specific clusters to view in the enrichment chart',
                value=NULL),
      
      uiOutput('pathwaysDT9'),
      
      uiOutput('genes_in_paths_DT9')
                    
    ), 
       
    mainPanel(width = 9,
             
      uiOutput('enrichmentUI')
             
    )
  )
), 
# Heatmaps________________________________________________________________#----
tabPanel("Term Gene Heatmap",
  fluidPage(
  
    sidebarPanel(width = 3,
    
      uiOutput('open_in_new_tab10'), 
      
      textInput('pathway_heatmap_genes',
      'Genes to include in the heatmap chart',
      value=NULL),
      
      textInput('heatmap_paths',
      'Pathway terms to be included in the heatmap chart',
      value=NULL),
      
      uiOutput('pathwaysDT10'),
      
      uiOutput('genes_in_paths_DT10')
    
    ),
    
    mainPanel(width = 9,
    
      uiOutput('pathway_heatmap')
    
    )
  )
),
# Case Map________________________________________________________________#----
tabPanel("Case vs. Control",
  fluidPage(
    
    uiOutput('cases_select_box'),
    
    sidebarPanel(width = 4,
              
      uiOutput('sample_conditions_PreviewTable11')
              
    ),
    
    mainPanel(width = 8,
         
      uiOutput('case_plot_ui')
    
    )
  )
),
# Single Pathway heatmap__________________________________________________#----
tabPanel("Single Pathway vs. Samples",
  fluidPage(
   
    sidebarPanel(width = 4,

      textInput('Single_pathway_plot_search', 
                'Pathway to visualize', 
                value = ""),
      
      numericInput('Single_pathway_plot_num_points', 
                   'Number of genes to visualize (priority to lowest pvals',
                   value = 40,
                   min = 10)
      
    ),
    
    mainPanel(width = 8,
              
      uiOutput('Single_pathway_plot_ui')
      
    )
  )
)
#----       
            ) #Tabset 2
          )
        )
      )
    )
  )


  )#End of tabPanel list
)#End of tag list ui


