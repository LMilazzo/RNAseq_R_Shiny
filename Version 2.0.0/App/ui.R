library(shiny)
library(shinycssloaders)
source("Pages-Themes-Assets/CSS/Styles.R")

files <- list.files("modules", "*.R")  # locate all .R files
sapply(paste0("modules/",files), source)
files <- list.files("modules/DiffAnalysisModules", "*.R")
sapply(paste0("modules/DiffAnalysisModules/", files), source)
files <- list.files("sourcedFunctions/UI Elements", "*.R")
sapply(paste0("sourcedFunctions/UI Elements/", files), source)
files <- list.files("sourcedFunctions/Scripts", "*.R")
sapply(paste0("sourcedFunctions/Scripts/", files), source)


runLibs()

ui <- 
  
  tagList(
    
    getStyles(),
    
    useShinyjs(), 

    fluidPage(
      
      navbarPage(id = "Landing", title = "🧬",
      theme = shinytheme("slate"),
      selected = "Differential Gene Expression",
      collapsible = TRUE,
      
      
      header = tagList(
        #Reload and save buttons
        div(
          actionButton("ReloadApp", "Restart Session",
            style = 
            "background-color: #FF6F61; color: black; border-color: black; background-image: none; margin-right: 5%"
          ),
          shiny::downloadButton(
            outputId = 'downloadZip', 
            label = 'Download Saved files',
            style = " border-color: #4CAF50; background-image: none; margin-right: 5%"
          ),
          style = "text-align: right"
        )
      ),
      
      helpTab(),
    
      tabPanel(
        title = "Differential Gene Expression",
        
        #hidden(
          div(id = "Diff-Uploads", diffUploadingDataTab("Diff-Uploads")),#######
        #),
        
        hidden(
          div(id = "Diff-Intermediate", diffIntermediateTab("Diff-Intermediate"))
        ),
        
        hidden(
          div(id = "Diff-Results", diffResultsTab("Diff-Results"))
        )
        
      ),
      
      tabPanel(
        title = "Pathway Analysis",
        
        div(id = "Path-Uploads", pathUploadingDataTabUI("Path-Uploads"))
      )
      
      )
    )
  )
