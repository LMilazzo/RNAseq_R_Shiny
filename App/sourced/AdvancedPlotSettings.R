advancedSettingsModal <- function(
    eventID, 
    settings,
    input,
    wantedSettings
){
  
  orgT <- settings$Title
  orgS <- settings$Subtitle
  orgC <- settings$Caption
  orgAR <- settings$Aspect
  orgTSA <- settings$TSA
  orgWidth <- settings$Width
  orgHeight <- settings$Height
  orgCellWidth <- settings$CellWidth
  orgCellHeight <- settings$CellHeight

  
  #TEXT ADDITIONS
  title.input <- if("tit" %in% wantedSettings) returnTitle(orgT, eventID) else NULL
  caption.input <- if("cap" %in% wantedSettings) returnCaption(orgC, eventID) else NULL
  subtitle.input <- if("sub" %in% wantedSettings) returnSubtitle(orgS, eventID) else NULL
  if(any(!sapply(list(
    title.input, 
    caption.input, 
    subtitle.input), is.null))){
    textAdditions <- div(
      h4("Text Additions"),
      title.input,
      subtitle.input,
      caption.input
    )
  }else{textAdditions <- NULL}
  
  #SIZE ADJUSTMENTS
  textSizeAdjustment.input <- if("tsa" %in% wantedSettings) returnTextAdjustment(orgTSA, eventID) else NULL
  
  #Find previous dims
  if(orgWidth == "auto"){orgWidth <- settings$DefaultWidth}
  if(orgHeight == "auto"){orgHeight <- settings$DefaultHeight}
  plotWidth.input <- if("w" %in% wantedSettings) returnWidth(orgWidth, eventID) else NULL
  plotHeight.input <- if("h" %in% wantedSettings) returnHeight(orgHeight, eventID) else NULL
  
  
  plotRatio.input <- if("ar" %in% wantedSettings) returnAspectRatio(orgAR, eventID) else NULL
  #Heatmap Cell Sizing
  cellWidth.input <- if("cw" %in% wantedSettings) returnCellWidth(eventID, settings) else NULL
  observeEvent(input[[paste0(eventID,"_autoCellWidth")]],{
    
    if(input[[paste0(eventID,"_autoCellWidth")]] == FALSE){
      shinyjs::show(paste0(eventID, "_cellWidth"))
    }
    if(input[[paste0(eventID,"_autoCellWidth")]] == TRUE){
      shinyjs::hide(paste0(eventID, "_cellWidth"))
    }
    
  })
  cellHeight.input <- if("ch" %in% wantedSettings) returnCellHeight(eventID, settings) else NULL
  observeEvent(input[[paste0(eventID,"_autoCellHeight")]],{
    
    if(input[[paste0(eventID,"_autoCellHeight")]] == FALSE){
      shinyjs::show(paste0(eventID, "_cellHeight"))
    }
    if(input[[paste0(eventID,"_autoCellHeight")]] == TRUE){
      shinyjs::hide(paste0(eventID, "_cellHeight"))
    }
    
  })
  
  if(any(!sapply(list(
    textSizeAdjustment.input, 
    plotWidth.input, 
    plotHeight.input, 
    plotRatio.input,
    cellWidth.input,
    cellHeight.input), is.null))){
    sizeAdjustments <- div(
      h4("Size Adjustments"),
      textSizeAdjustment.input,
      plotRatio.input,
      plotWidth.input,
      plotHeight.input,
      cellWidth.input,
      cellHeight.input
    )
  }else{sizeAdjustments <- NULL}
  
  
  modalStuff <- div(
    
    textAdditions,
    
    sizeAdjustments
    
  )
  
  
  showModal(
    modalDialog(
      modalStuff,
      footer = tagList(
        actionButton(inputId = paste0(eventID,"_scrapSettings"), label ="Scrap"),
        actionButton(inputId = paste0(eventID,"_applySettings"), label = "Apply")
      ),
      easyClose = FALSE
    )
  )
  
}

returnTitle <- function(orgT, eventID){
  title.input <- textInput(paste0(eventID, "_title"),
                           "Title",
                           value = orgT)
  return(title.input)
}
returnCaption <- function(orgC, eventID){
  caption.input <- textInput(paste0(eventID, "_caption"),
            "Caption",
            value = orgC)
  return(caption.input)
}
returnSubtitle <- function(orgS, eventID){
  subtitle.input <- textInput(paste0(eventID, "_subtitle"),
            "Subtitle",
            value = orgS)
  return(subtitle.input)
}
returnTextAdjustment <- function(orgTSA, eventID){
  textSizeAdjustment.input <-numericInput(paste0(eventID, "_textSizeAdjustment"),
               "Text Size Adjustment",
               value = orgTSA,
               min = -5,
               max = 15)
  return(textSizeAdjustment.input)
}
returnWidth <- function(orgWidth, eventID){
  plotWidth.input <- shinyWidgets::noUiSliderInput(paste0(eventID, "_width"),
                                                   "Plot Width",
                                                   value = orgWidth,
                                                   min = 100,
                                                   step = 10,
                                                   max = orgWidth * 2,
                                                   width = "100%",
                                                   color = "indianred"
  )
  return(plotWidth.input)
}
returnHeight <- function(orgHeight, eventID){
  plotHeight.input <- shinyWidgets::noUiSliderInput(paste0(eventID, "_height"),
                                                   "Plot Height",
                                                   value = orgHeight,
                                                   min = 100,
                                                   step = 10,
                                                   max = orgHeight * 2,
                                                   width = "100%",
                                                   color = "indianred"
  )
  return(plotHeight.input)
}
returnAspectRatio <- function(orgAR, eventID){
  orgAR <- switch(as.character(orgAR),
                  "1" = "1:1",
                  "0.75" = "3:4",
                  "1.33333333333333" = "4:3",
                  "1.77777777777778" = "16:9",
                  "0.5625" = "9:16",
                  "auto" = "auto")
  plotRatio.input <- selectInput(paste0(eventID, "_aspect"),
                                 "Aspect Ratio",
                                 choices = c(
                                   "auto",
                                   "1:1", 
                                   "4:3", 
                                   "3:4",
                                   "16:9",
                                   "9:16"
                                 ),
                                 selected = orgAR
  )
  return(plotRatio.input)
}
returnCellWidth <- function(eventID, settings){
  
  if(settings$CellWidth != "auto"){
    start <- settings$CellWidth
  }else{
    start <- 75
  }
  
  d <- div(
    checkbox <- checkboxInput(paste0(eventID,"_autoCellWidth"), "Auto Cell Width", value = settings$CellWidth == "auto"),
    cellWidth.input <- shinyWidgets::noUiSliderInput(paste0(eventID, "_cellWidth"),
                                                   "Cell Width",
                                                   value = start,
                                                   min = 10,
                                                   max = 200,
                                                   step = 1,
                                                   width = "100%",
                                                   color = "indianred")
  )
  
  return(d)
  
}
returnCellHeight <- function(eventID, settings){
  
  if(settings$CellHeight != "auto"){
    start <- settings$CellHeight
  }else{
    start <- 75
  }
  
  d <- div(
    checkbox <- checkboxInput(paste0(eventID,"_autoCellHeight"), "Auto Cell Height", value = settings$CellHeight == "auto"),
    cellHeight.input <- shinyWidgets::noUiSliderInput(paste0(eventID, "_cellHeight"),
                                                     "Cell Height",
                                                     value = start,
                                                     min = 10,
                                                     max = 200,
                                                     step = 1,
                                                     width = "100%",
                                                     color = "indianred")
  )
  
  return(d)
}