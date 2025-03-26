
pathUploadingDataTabUI <- function(id){

  ns <- NS(id)
  
  fluidPage(
    
    fluidRow(
      actionButton(
        ns('ContinueSessionWithPathway'), 
        'Continue Session',
        style = "background-color: #4CAF50; color: black; 
        border-color: white; background-image: none;
        margin-top: 0px ; margin-bottom: 0px; margin-right: 5%"
      ),
      actionButton(
        ns('UploadPrevPathData'), 
        'Review An Experiment',
        style = "background-color: #007bff ; color: black; 
        border-color: white; background-image: none;
        margin-top: 0px ; margin-bottom: 0px; margin-left: 5%"
      ),
      style = "text-align: center; margin-top: 5%; margin-bottom: 0px;
      z-index: 2; position: relative;"
    ),
    
    div(
      includeHTML("Pages-Themes-Assets/HTML Pages/PathFindR_Intro.HTML"),
      style = "margin-top: 10px; padding: 0px; z-index: 1;"
    )
    
  )
  
  
}