
#A button that reloads the application when clicked
reloadApplication <- function(){
  return(
    actionButton(
      'reload_app', 
      'Reload Application', 
       style = "background-color: #FF6F61; color: #FF6F61; border-color: #FF6F61;"
    )  
  )
}

#A button that allows for the download of the temp zip file of saved plots
downloadZipButton <- function(){
  return(
    shiny::downloadButton(
      outputId = 'downloadTemp', 
      label = 'Download Saved files'
    )
  )
}