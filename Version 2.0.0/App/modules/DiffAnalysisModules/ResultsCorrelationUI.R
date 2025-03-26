ResultsCorrelationUI <- function(id){
  
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
          div(id = ns("Data-Stuff"),
            #Display Numbers
            div(id="NumDisplay",
              p("Display Values: ", style = "margin-right: 15px; padding-bottom: 12px"),
              checkboxInput(ns("NumDisplay"), "", value = FALSE),
              style = "font-size: 17px; display: inline-flex; align-text: center; justify-content: center; align-items: center; overflow: auto;"
            ),  
            #Meta Data to show----
            div(class = "collapsible-container",
                tags$details(
                  tags$summary(
                    p("Annotations", style = "margin-right: 5%; padding-bottom: 2%"),
                  ),
                  tags$div(class = "content", uiOutput(ns("Meta")))
                ),
                style = "padding-right: 5px; margin-top: 10px;"
            ),
            style = "border: 2px solid #272b30; border-radius: 15px; padding-left: 5px; margin-top: 5px; overflow: auto;"
          )
          
          
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
          #NOT USED YET
          # div(id=ns("Theme-Stuff"),
          #   style = "border: 2px solid #272b30; border-radius: 15px; padding-left: 5px; margin-top: 5px;"
          # ),
          div(id=ns("Size-Stuff"),
            # #Aspect Ratio (DETERMINED BY CELL SIZINGS---- 
            # div(
            #   p("Aspect Ratio:", style = "margin-right: 5%; padding-bottom: 2%"),
            #   selectInput(
            #     ns("Aspect"),label = NULL,
            #     choices = c("4:3" = 4/3, "3:4" = 3/4, "16:9" = 16/9, "9:16" = 9/16, "1:1" = 1), 
            #     selected = c(3/4), width = "100%"
            #   ),
            #   style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center;"
            # ),
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
            #Cell-width----
            div(
              p("Cell Width:", style = "margin-right: 5%; padding-bottom: 2%"),
              numericInput(ns("cW"), max = 4280, min = 0, value = 40, label = NULL, width = "50%"),
              style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center; width: 100%;"
            ),
            #Cell-height----
            div(
              p("Cell Height:", style = "margin-right: 5%; padding-bottom: 2%"),
              numericInput(ns("cH"), max = 4280, min = 0, value = 40, label = NULL, width = "50%"),
              style = "margin-top: 10px; margin-bottom: 5px; display: inline-flex; font-size: 17px; align-items: center; width: 100%;"
            ),
            style = "border: 2px solid #272b30; border-radius: 15px; padding-left: 5px; margin-top: 5px;"
          ),
          style = "overflow: auto;"
        )
        #----
      )
    ),
    
    mainPanel(
      width = 9,
      
      div(
        withSpinner(
          plotOutput(ns("test")),
          type = 6
        )
      )
      
    )
  )
  
  
  
}