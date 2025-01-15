
#Sourced Files
source("sourced/Functions.R")
source("sourced/LoadPackages.R")
source("sourced/Buttons.R")

#Load packages
runLibs()

#App options and settings ----

options(shiny.maxRequestSize=100*1024^2)  # Limits file upload size to 100 MB

#CSS
tags.css <- tags$style(HTML("             
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
    
    body {
    -moz-transform: scale(0.8, 0.8); /* Moz-browsers */
    zoom: 0.8; /* Other non-webkit browsers */
    zoom: 80%; /* Webkit browsers */
    }
  "))

#App UI tag List ----

ui <- 
  
  tagList(
    
    tags.css, # add css
  
    useShinyjs(), #initialize shiny javascript
    
#Top level navigation bar setup ---- 

navbarPage(
  id = "App",
  title = "", 
  selected = "Diffrentially Expressed Genes", #Initial child page
  theme = shinytheme("slate"),
  
#______Navigation page 1 "Help"______ ----

tabPanel(
  title = "Help",
  
  fluidPage(
    
    plotlyOutput('interactive_image')
  
  )
     
),
    

#______Navigation page 2 "Differential Expressed Genes"_____ ----

tabPanel("Diffrentially Expressed Genes",

  fluidPage(

#Page 2 functional button row ----

#Top row of buttons [ reload, new exp, old exp, run, contrasts]
fluidRow(
  
  #Reload application
  column(
    width = 1, 
    div(reloadApplication())
  ),
    
  #Experiment starting point options
  div(
    id = 'Page1_Upload_Options',
    column(
      width = 1,
      offset = 1,
      actionButton(
        'start_new_experiment', 
        'Start New Experiment',
        style = "background-color: #4CAF50; color: #4CAF50; border-color: #4CAF50;"
      )
    ),
    column(
      width = 2,
      offset = 1,
      actionButton(
        'start_old_experiment', 
        'Review A Previous Experiment'
      )
    )
  ),
  
  #Run button ( start hidden )
  column(
    width = 1,
    offset = 1,
    hidden(
      actionButton(
        "run_DESeq2", 
        "Run Differential Expression",
        style = "background-color: #4CAF50; color: #4CAF50; border-color: #4CAF50;"
      )
    )
  ),
  
  #Contrast selection drop down
  column(
    width = 2,
    offset = 1,
    uiOutput('contrast_selection_ui')
  )
),
        
#Page 2 Main panel set up ----

br(),

mainPanel(
  width = 12,
  
#Page 2 HTML Info Page ( pre data only) ----
  
uiOutput('about_DESeq2'), #HTML Info Page

#Page 2 upload previews ----

div(
  id = 'preRun_Data_Preview',
  uiOutput('preRun_preview1'),
  uiOutput('preRun_preview2')
),

#Page 2 Bulk TabSet setup ----  

hidden( #Start hidden
  
  div(
    id = "TabSet1_Diffrential_Expression_Analysis", 
    
    tabsetPanel(
      id = "TabSet1",

#Preview Experiment tab 1 ----

tabPanel(
  title = "Data Preview", 

  #Side Bar Raw Counts preview
  sidebarPanel(
    width = 6,
    HTML("<h3>Counts Matrix</h3>"),
    uiOutput('raw_counts_PreviewTable')
  ),
  
  #Side Bar Normalized counts preview
  sidebarPanel(
    width = 6,
    HTML("<h3>Normalized Counts</h3>"),
    uiOutput('normalized_counts_PreviewTable')
  ),

  #Bottom page Condition Preview
  HTML("<h3>Sample Conditions</h3>"),
  uiOutput('sample_conditions_PreviewTable'),
  
),

#DEG analysis tab 2 ----

tabPanel(
  title = "Diffrential Analysis",
  
  #Page setup
  fluidPage(
    sidebarLayout(
    
      #Side bar 
      sidebarPanel(
        #Pvalue select box
        selectInput(
          "pvaluePg2", 
          "Select adjusted P value cutoff for following tables", 
          choices=c(1.0 ,0.5, 0.05, 0.01, 0.001)
        ),
        
        #Stat check box
        checkboxGroupInput(
          'display_col', 
          'Values of Interest (may take longer to load if more selected)', 
          choices = c(
            'Base Mean' = 'baseMean', 
            'Log2 Fold Change' ='log2FoldChange',
            'lfcSE' ='lfcSE',
            'Stat' ='stat',	
            'P-Value' ='pvalue'
          ),
          selected=c('log2FoldChange'),
          inline=FALSE
        ),
        
        #Distribution of fold change
        uiOutput("Distribution_Histogram_ui")
      ),
      
      #Main Panel
      mainPanel(
        #Cutoff Slider
        sliderInput(
          'cutOffs', 
          'Cut-off values for the genes being displayed as a certain direction of diffrential expression', 
          min=-12, 
          max=12, 
          step=0.10, 
          value = c(-0.5, 0.5), 
          width='100%'
        ),
        
        #Data tables
        uiOutput('expression_tables')
      )
      
    )#Sidebar
  )#fluid page
),

#PCA tab 3 ----
tabPanel(
  "Principle Component Plots",
  
  fluidPage(
  
    #side bar
    sidebarPanel(
      width=3,
      
      HTML("<h3>Data</h3>"),
      uiOutput('pca_n_counts'),
      uiOutput('pca_metadata'),
      
      HTML("<h3>Labels</h3>"),
      
      textInput(
        'title_pca_plot', 
        'Title', 
        value = 'Principle Component Analysis of Sample Varience'
      ),
      textInput(
        'subtitle_pca_plot', 
        'Sub Title', 
        value = 'Colored by sample traits'
      ),
      textAreaInput(
        'caption_pca_plot', 
        '', 
        value = '', 
        width=200, 
        rows=3
      ),
      
      savePlotButton()
      
    ),
    
    #main panel
    mainPanel(
      width = 9,     
      uiOutput('principle_component_plots_ui')
    )
    
  )#fluidpage
),

#Correlation Analysis tab 4 ----
tabPanel(
  "Correlation Anaylsis",
  
  fluidPage(
    
    #side bar
    sidebarPanel(
      width=2,
      HTML("<h3>Labels</h3>"),
      textInput(
        'title_heatmap_plot', 
        'Title', 
        value = 'Sample Correlation'
      ),
      HTML("<h4>Annotations to Include</h4>"),
      uiOutput('choose_annotations_ui'),
      savePlotButton()
    ),
    
    #main panel
    mainPanel(
      width = 10,     
      uiOutput('heatmap_plot_ui')
    )

  )#fluid page
),
#Gene Counts tab 5 ----
tabPanel(
  "Gene Counts",
  
  fluidPage(
  
    #Side bar
    sidebarPanel(
      width=4,
      
      #Pvalue cut off 
      selectInput(
        "pvaluePg6", 
        "Select adjusted P value cutoff for following plots", 
        choices=c(1.0 ,0.5, 0.05, 0.01, 0.001)
      ),
      
      uiOutput('leftbar_gene_plots')
    ), 
    
    #main panel
    mainPanel(
      width=6,
      uiOutput('gene_count_search'),
      br(),
      savePlotButton()
    ), 
    
    #side panel
    sidebarPanel(
      width=2,
      uiOutput("gene_name_list")
    )
  
  )#fluid page
), 
#Volcano Plot tab 6 ----
tabPanel(
  "Volcano Plot",
  
  fluidPage(
  
    #side bar 
    sidebarPanel(
      width=3,
      
      #diff expression cutoffs
      sliderInput(
        'volcano_cutoffs', 
        'Cut of values for diffrential expression (verticle lines on Volcano Plot)', 
        min = -12, 
        max = 12,
        step = 0.01, 
        value = c(-1, 1), 
        width = '100%'
      ),
     
     #P value cut off
     selectInput(
        "pvaluePg7", 
        "Select adjusted P value cutoff for following plot", 
        choices = c(1.0 ,0.5, 0.05, 0.01, 0.001), 
        selected = 0.05
      ),
     
     #Volcano Population
      sliderInput(
        'volcano_pop', 
        'Change the # of data points', 
        min = 0, 
        max = 1, 
        value = 0.9, 
        width = '100%', 
        ticks = FALSE
      ),
     
      #Label density
      sliderInput(
        'volcano_lab_density', 
        'Change # of points labeled', 
        min = 0, 
        max = 1, 
        value = 0.3, 
        width = '100%', 
        ticks = FALSE
      ),
     
     #Searched genes
      textInput(
        'volc_search',
        'Genes to search/highlight (can be comma seperated list)',
        value = NULL
      ),
     
      #Plot labels
      textInput('title_volc_plot', 'Title', value = 'Gene Expression'),
      textInput('subtitle_volc_plot', 'Sub Title', value = ''),
      textAreaInput('caption_volc_plot', 'Caption', value = '', width = 200, rows = 3),
      savePlotButton(),
    ),
    
    #main panel
    mainPanel(
      width = 9,
      uiOutput('volcano_plot_ui')
    ), 
    
  )#fluid page
)
#Page 2 Closing brackets ----
          ) #Tabset 1
        ) #Div
      ) #Hidden
    ) #Main panel below button
  ) #High level fluid page
), #High level tab panel

#______Navigation page 3 "Pathway Analysis"______ ----
tabPanel("Pathway Analysis",

  fluidPage(

#Page 3 Functional button row ----

#Top row of buttons [ reload, run, old exp]
fluidRow(
  
  #Reload application
  column(
    width = 1, 
    div(reloadApplication())
  ),
  
  div(
    id = "pathfinder_option_buttons",
    column(
      width = 2,
      offset = 1,
      actionButton(
        'Run_pathfinder', 
        'Continue Experiment With Pathway Analysis',
        style = "background-color: #4CAF50; color: #4CAF50; border-color: #4CAF50;"
      ),
    ),
    
    column(
      width = 2, 
      offset = 1,
      actionButton(
        'review_pathfinder_new_data', 
        'Review A Pathway Analysis Experiment'
      )
    )
  )
),

#Page 3 Main panel set up ----  

br(),

mainPanel(
  width = 12,

#Page 3 HTML Info Page ( pre data only ) ----

uiOutput('about_pathfinder'), #HTML Info Page 

#Page 3 Bulk TabSet setup ----

hidden( #Start hidden
  
  div(
    id = "TabSet2_Pathway_Analysis",
    
    tabsetPanel(
      id = "TabSet2",

#Preview Experiment tab 1 ----
tabPanel(
  "Pathway Analysis Table",

  fluidPage(
    mainPanel(
      uiOutput('pathfinderPreview')
    )
  )
  
),
#Pathway Enrichment Plot tab 2 ----
tabPanel(
  "Term Enrichment",
  
  fluidPage(
    
    #Slider for number of clusters included 
    fluidRow(
      column(
        width = 4,
        uiOutput('enrichment_clusters_shown')
      ),
      column(
        width = 1, 
        offset = 2,
        savePlotButton()
      )
    ),
    
    #sidebar
    sidebarPanel(
      width = 3,
      
      #genes to include
      textInput(
        'enrichment_genes',
        'Genes to include in the enrichment chart',
        value = NULL
      ),
      
      #Pathways to include
      textInput(
        'enrichment_paths',
        'Pathway terms to be included in the enrichment chart',
        value = NULL
      ),
      
      #search specific clusters
      textInput(
        'enrichment_clusters',
        'Specific clusters to view in the enrichment chart',
        value = NULL
      ),
      
      #A list of all the pathways
      uiOutput('pathwaysDT9'),
      
      #A list of all the genes
      uiOutput('genes_in_paths_DT9')
    ), 
    
    #main panel
    mainPanel(
      width = 8,
      uiOutput('enrichmentUI')
    )
  )#fluid page
), 
#Heatmaps tab 3 ----
tabPanel(
  "Term Gene Heatmap",
  
  fluidPage(
    
    #side bar
    sidebarPanel(
      width = 3,
      
      #link to open interactive plotly or save plot
      fluidRow(column(2,uiOutput('open_in_new_tab10')), column(2, offset = 4,savePlotButton())), 
      br(),
      
      #search for genes text box
      textInput(
        'pathway_heatmap_genes',
        'Genes to include in the heatmap chart',
        value = NULL
      ),
      
      #search for pathways text box
      textInput(
        'heatmap_paths',
        'Pathway terms to be included in the heatmap chart',
        value = NULL
      ),
      
      #Title for the plot
      textInput(
        'title_pathwayVgene_heatmap',
        'Plot Title',
        value = ''
      ),
      
      #A table of all pathways
      uiOutput('pathwaysDT10'),
      
      #A table of all genes
      uiOutput('genes_in_paths_DT10')
    ),
    
    #main panel
    mainPanel(
      width = 9,
      uiOutput('pathway_heatmap')
    )
    
  )#fluid page
),
#Case Map tab 4 ----
tabPanel(
  "Case vs. Control",
  
  fluidPage(
    
    #a list of all samples to choose whch ones are part of the case group
    uiOutput('cases_select_box'),
    
    #side bar
    sidebarPanel(
      width = 4,
      
      #A check box to decide if only representetive pathways are used
      checkboxInput(
        'repOnly', 
        "Representative pathways only", 
        value = TRUE
      ),
      
      #plot title
      textInput(
        'case_plot_title',
        'Plot Title',
        ''
      ),
      
      # a table with the meta data file for reference
      uiOutput('sample_conditions_PreviewTable11'),
      
      savePlotButton()
    ),
    
    #main panel
    mainPanel(
      width = 8,
      uiOutput('case_plot_ui')
    )
    
  )#fluid page
),
#Single Pathway heatmap tab 5 ----
tabPanel(
  "Single Pathway vs. Samples",
  
  fluidPage(
    
    #side bar
    sidebarPanel(
      width = 4,
      
      #The pathway to visualze text input
      textInput(
        'Single_pathway_plot_search', 
        'Pathway to visualize', 
        value = ""
      ),
      
      #The number of genes for that pathway to show
      numericInput(
        'Single_pathway_plot_num_points', 
        'Number of genes to visualize (priority to lowest pvals',
        value = 40,
        min = 10
      ),
      
      HTML("<h4>Annotations to Include</h4>"),
      uiOutput('choose_annotations_ui_p'),
      
      savePlotButton()
    ),
    
    #main panel
    mainPanel(
      width = 8,
      uiOutput('Single_pathway_plot_ui')
    )
    
  ) #fluid page
)
#Page 3 Closing brackets ----       
          ) #Tabset 2
        ) #Div
      ) #Hidden
    ) #Main Panel
  ) #High level fluid page
), #High level tab panel

#High level tab panel closing ----
  #download button
  header = tags$ul(
    class = "nav navbar-nav navbar-right", # Align to the right
    tags$li(
      style = "padding: 0px 45px 0 0; margin-bottom: 0px;", # Top, right, bottom, left,
      downloadZipButton()
    )
  )
)#End of tabPanel list

#UI tag list closing ----
)#End of tag list ui


