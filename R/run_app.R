#' Run the ShinyPLS Application
#'
#' @description Launches the graphical user interface.
#' @export
#' @import shiny
#' @return A Shiny application object
run_app <- function() {

  library(shiny)
  library(bslib)
  library(DT)
  library(dplyr)
  library(readxl)
  library(openxlsx)
  library(digest)
  library(cSEM)
  library(genpathmox)
  library(DiagrammeR)
  library(ggplot2)
  library(ggforce)
  library(Rtsne)

  options(bitmapType = "cairo")

  shiny::shinyApp(ui = app_ui(), server = app_server, options = list(launch.browser = TRUE))
}
