ResultsGeneCountsUI <- function(id){
  
  ns <- NS(id)
  
  fluidPage(
    br(),
    
    sidebarPanel(
      width = 3,
      tabsetPanel(
        tabPanel(
          "Data",
          #Gene to plot----
          div(
            p("Search: ", style = "margin-right: 5%; padding-bottom: 2%"),
            textInput(ns("main"), label = NULL, value = "", width = "50%"),
            style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center;"
          ),
          #Annotations to show x----
          div(class = "collapsible-container",
            tags$details(
              tags$summary(
                p("X axis", style = "margin-right: 5%; padding-bottom: 2%")
              ),
              tags$div(class = "content", 
                       uiOutput(ns("x"), width = "100%"), 
                       style = ""
              )
            ),
            style = "margin-top: 5px;"
          ),
          #----
          #Annotations to show Color----
          div(class = "collapsible-container",
              tags$details(
                tags$summary(
                  p("Color", style = "margin-right: 5%; padding-bottom: 2%")
                ),
                tags$div(class = "content", 
                         uiOutput(ns("color"), width = "100%"), 
                         style = ""
                )
              ),
              style = "margin-top: 5px;"
          ),
          #----
          #Annotations to show Shape----
          div(class = "collapsible-container",
              tags$details(
                tags$summary(
                  p("Shape", style = "margin-right: 5%; padding-bottom: 2%")
                ),
                tags$div(class = "content", 
                         uiOutput(ns("shape"), width = "100%"), 
                         style = ""
                )
              ),
              style = "margin-top: 5px;"
          ),
          #----
          #Annotations to show Size----
          div(class = "collapsible-container",
              tags$details(
                tags$summary(
                  p("Size", style = "margin-right: 5%; padding-bottom: 2%")
                ),
                tags$div(class = "content", 
                         uiOutput(ns("size"), width = "100%"), 
                         style = ""
                )
              ),
              style = "margin-top: 5px;"
          )
          #----
        ),
        tabPanel(
          "Settings",
          #Settings----
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
            #grid lines----
            div(
              p("Grid Lines:", style = "margin-right: 5%; padding-bottom: 2%"),
              selectInput(
                ns("Grid"), label = NULL, width = "100%",
                choices = c(
                  "None" = "theme(panel.grid = element_blank(), )", 
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
            #Aspect Ratio----
            div(
              p("Aspect Ratio:", style = "margin-right: 5%; padding-bottom: 2%"),
              selectInput(
                ns("Aspect"),label = NULL,
                choices = c("4:3" = 4/3, "3:4" = 3/4, "16:9" = 16/9, "9:16" = 9/16, "1:1" = 1), 
                selected = c(3/4), width = "100%"
              ),
              style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center;"
            ),
            #width----
            div(
              p("Width:", style = "margin-right: 5%; padding-bottom: 2%"),
              numericInput(ns("W"), max = 4280, min = 0, value = 750, label = NULL, width = "50%"),
              style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center; width: 100%;"
            ),
            #height----
            div(
              p("Height:", style = "margin-right: 5%; padding-bottom: 2%"),
              numericInput(ns("H"), max = 4280, min = 0, value = 450, label = NULL, width = "50%"),
              style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center; width: 100%;"
            ),
            style = "border: 2px solid #272b30; border-radius: 15px; padding-left: 5px; margin-top: 5px;"
        ),
          #----
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
          plotOutput(ns("main"), width = "100%", height = "450px"),
          type = 6
        )
      )
      
      
      
    )
    
    
    
    
    
  )
  
  
  
}