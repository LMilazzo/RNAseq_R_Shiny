
# Basic function that returns and displays a modal with the given message in red text
showErrorModal <- function(message) {
  showModal(modalDialog(tags$p(style = "color: red;", message), easyClose = TRUE, footer = NULL))
}
