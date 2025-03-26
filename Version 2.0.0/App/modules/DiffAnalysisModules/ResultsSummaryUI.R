ResultsSummaryUI <- function(id){

  ns <- NS(id)

  fluidPage(
    
    sidebarPanel(
      width = 3,
      
      #CHANGE PVAL
      div(
        p("P-value: ", style = "margin-bottom: 18px; font-size: 18px;"),
        selectInput(ns("Pval"), label = NULL, 
          choices = c("0.5", "0.1", "0.05", "0.01", "0.005", "0.001", "Ignore"),
          selected = "Ignore",
          width = "50%"
        ),
        style = "display: flex; align-items: center; gap: 5%; width: 100%; "
      ),
      
      #PVAL STATS
      div(
        DTOutput(ns("PvalStats")), #DT
        style = "width: 90%; text-align: center; font-size: 16px; overflow: auto;"
      ),
      
      div(
        plotOutput(ns("PvalHisto"), height = "100%", width = "90%"), #HISTOGRAM
        style = "max-height: 300px; height: 300px; margin-top: 15px; overflow: auto;"
      ),
      
      #FOLD CHANGE STATS
      div(
        p("Fold Change: ", style = "margin-bottom: 18px; margin-top: 16px; font-size: 18px; text-align: left;"),
        DTOutput(ns("FoldStats")), #DT
        style = "width: 90%; text-align: center; font-size: 16px; overflow: auto;"
      ),
      
      div(
        plotOutput(ns("FoldHisto"), height = "100%", width = "90%"), #HISTOGRAM
        style = "max-height: 300px; height: 300px; margin-top: 15px; overflow: auto;"
      ),
      
    ),
    
    mainPanel(
      width = 9,
      
      #FILTER RANGE
      div(
        h4("Filter Results By Fold Change Range:", style = "padding-bottom: 20px;"),
        shinyWidgets::noUiSliderInput(
          ns("FoldRange"),
          color = "#FF6F61",
          value = c(-30,30),
          width = "100%",
          range = list("min" = -100, "8%" = -30, "14%" = -15, "20%" = -10, "26%" = -7,
                       "32%" = -5, "38%" = -3, "44%" = -1, "50%" = 0, "56%" = 1, "62%" = 3,
                       "68%" = 5, "74%" = 7, "80%" = 10, "86%" = 15, "92%" = 30, "max" = 100),
          pips = list(
            mode = "values",
            values = list(-30,-15, -10, -7, -5, -3, -1, 0, 1, 3, 5, 7, 10, 15, 30),
            density = 40
          ),
          update_on = "end"
        )
      ),
      
      div(
        DTOutput(ns("Res"))
      )
    )
    
    
    # sidebarLayout(
    #   sidebarPanel(
    #     width = 2,
    #     
    #     div(
    #       id = "GeneNamesTableContainer",
    #       DTOutput(ns("GeneNames")),
    #       style = "overflow: auto; text-align: left;"
    #     ),
    #     
    #     tags$style(HTML(
    #       "GeneNamesTableContainer {
    #         width: 100%; display: flex; justify-content: left; align-items: left; test-align: left;
    #       }
    #       
    #       .dataTables_filter input {
    #         width: 85% !important;
    #         box-sizing: border-box;
    #         display: flex; justify-content: left; align-items: left; text-align: left;
    #       }"
    #     ))
    #     
    #     
    #   ),
    #   mainPanel(
    #     
    #   )
    # )
    
    
  )
    
}