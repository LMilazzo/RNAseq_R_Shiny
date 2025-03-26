ResultsVolcanoUI <- function(id){
  
  ns <- NS(id)
  
  fluidPage(
    br(),
    
    sidebarPanel(
      width = 3,
      tabsetPanel(
        
        #----
        #   DATA    #TO-DO
        #----    
        tabPanel(
          "Data",                   
          # div(id = ns("Data-Stuff"),
          #   div(
          #     p("P-value: ", style = "margin-bottom: 18px; font-size: 18px;"),
          #     selectInput(ns("Pval"), label = NULL, 
          #                 choices = c("0.5", "0.1", "0.05", "0.01", "0.005", "0.001", "Ignore"),
          #                 selected = "Ignore",
          #                 width = "50%"
          #     ),
          #     style = "display: flex; align-items: center; gap: 5%; width: 100%; margin-top: 10px; "
          #   ),
          #   style = "border: 2px solid #272b30; border-radius: 15px; padding-left: 5px; margin-top: 5px;"
          # ),
          div(id = ns("Search-Stuff"),
            p("Highlight Searched Genes", style = "margin-bottom: 10px; font-size: 18px;"),
            p("Input as comma seperated list", style = "margin-bottom: 3px; font-size: 12px;"),
            textAreaInput(ns("search"), NULL, width = "90%"),
            style = "border: 2px solid #272b30; border-radius: 15px; padding-left: 5px; margin-top: 5px; overflow: auto;"
          ),
          div(id = ns("Visual-Stuff"),
            p("# of Labels", style = "margin-bottom: 8px; font-size: 18px;"),
            p("(as a % of the maximum # of labels, 230. 200 Are selected by pval and 30 selelcted by fold change.)", style = "margin-bottom: 3px; font-size: 12px;"),
            div(
              shinyWidgets::noUiSliderInput(ns("labels"), NULL, max = 1, min = 0, value = 0.25, width = "90%", color = "#ff6f61"),
              style = "display: flex;"
            ),
            p("# of Points", style = "margin-bottom: 8px; font-size: 18px;"),
            p("(as a % of the total number of genes)", style = "margin-bottom: 3px; font-size: 12px;"),
            div(
              shinyWidgets::noUiSliderInput(ns("points"), NULL, max = 1, min = 0, value = 0.90, width = "90%", color = "#0275d8"),
              style = "display: flex;"
            ),
            div(
              p("P-Value Line", style = "margin-bottom: 18px; font-size: 18px;"),
              numericInput(ns("PLine"), NULL, min = 0, max = 1, value = 0.5, step = 0.05, width = "40%"),
              style = "display: inline-flex; align-items: center; gap: 5%; width: 95%;"
            ),
            p("Fold Change Lines", style = "margin-bottom: 18px; font-size: 18px;"),
            div(
              numericInput(ns("min"), NULL, min = -100, max = 10, value = -2, step = 0.10, width = "40%"),
              p("~"),
              numericInput(ns("max"), NULL, min = -10, max = 100, value = 2, step = 0.10, width = "40%"),
              style = "display: flex; align-items: center; gap: 5%; width: 95%; justify-content: center;"
            ),
            style = "border: 2px solid #272b30; border-radius: 15px; padding-left: 5px; margin-top: 5px; overflow: auto;"
          ),
          style = "overflow: auto;"
        ),
        #----
        #   SETTINGS ----
        tabPanel(
          "Settings", 
          div(id = ns("Text-Stuff"),
            #Plot Title----
            div(
              p("Plot Title:", style = "margin-right: 5%; padding-bottom: 2%"),
              textInput(ns("Title"), label =  NULL, value = "", width = "50%"),
              style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center;"
            ),
            #Subtitle----
            div(
              p("Sub Title:", style = "margin-right: 5%; padding-bottom: 2%"),
              textInput(ns("Subtitle"), label = NULL, value = "", width = "50%"),
              style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center;"
            ),
            #Caption----                   
            div(
              p("Caption:", style = "margin-right: 5%; padding-bottom: 2%"),
              textInput(ns("Caption"), label = NULL, value = "", width = "50%"),
              style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center;"
            ),
            #text size----           
            div(
              p("Font:", style = "margin-right: 5%; padding-bottom: 2%"),
              numericInput(ns("Font"), max = 32, min = 4, value = 18, label = NULL, width = "30%"),
              style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center;
          width: 100%;"
            ),
            style = "border: 2px solid #272b30; border-radius: 15px; padding-left: 5px; margin-top: 5px;"
          ),
          div(id=ns("Theme-Stuff"),
            div(
              p("Grid Lines:", style = "margin-right: 5%; padding-bottom: 2%"),
              selectInput(
                ns("Grid"), label = NULL, width = "100%",
                choices = c(
                  "None" = "theme(panel.grid = element_blank())", 
                  "Major" = "theme(panel.grid.major = element_line(), panel.grid.minor = element_blank())", 
                  "Minor" = "theme(panel.grid.major = element_blank(), panel.grid.minor = element_line())", 
                  "All" = "theme(panel.grid.major = element_line(), panel.grid.minor = element_line())"
                ), 
                selected = c("theme(panel.grid.major = element_line(), panel.grid.minor = element_blank())")
              ),
              style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center;"
            ),
            style = "border: 2px solid #272b30; border-radius: 15px; padding-left: 5px; margin-top: 5px;"
          ),
          div(id=ns("Size-Stuff"),
            #Aspect Ratio (DETERMINED BY CELL SIZINGS----
            div(
              p("Aspect Ratio:", style = "margin-right: 5%; padding-bottom: 2%"),
              selectInput(
                ns("Aspect"),label = NULL,
                choices = c("4:3" = 4/3, "3:4" = 3/4, "16:9" = 16/9, "9:16" = 9/16, "1:1" = 1),
                selected = c(4/3), width = "100%"
              ),
              style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center;"
            ),
            #width----
            div(
              p("Width:", style = "margin-right: 5%; padding-bottom: 2%"),
              numericInput(ns("W"), max = 4280, min = 0, value = 900, label = NULL, width = "50%"),
              style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center; width: 100%;"
            ),
            #height----
            div(
              p("Height:", style = "margin-right: 5%; padding-bottom: 2%"),
              numericInput(ns("H"), max = 4280, min = 0, value = 700, label = NULL, width = "50%"),
              style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center; width: 100%;"
            ),
            style = "border: 2px solid #272b30; border-radius: 15px; padding-left: 5px; margin-top: 5px;"
          ),
          style = "overflow: auto;"
        ),
        tabPanel(
          "Search",
          #Lookup table ----
          div(  
            DTOutput(ns("Lookup"), width = "100%"),
            style = "display: flex; justify-content: center; align-items: center; text-align: center; overflow: auto; margin-top: 10px;",
            tags$style(HTML(".dataTables_filter input {width: 70% !important;}"))
          )
          #----
        )
      )
    ),
    
    mainPanel(
      width = 9,
      
      div(
        withSpinner(
          plotOutput(ns("main")),
          type = 6
        )
      )
      
    )
  )
  
  
  
}