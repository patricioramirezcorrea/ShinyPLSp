#' Run the ShinyPLS Application
#'
#' @description Launches the graphical user interface.
#' @export
#' @import shiny
#' @return A Shiny application object
run_app <- function() {

  # Al poner las librerías aquí, nos aseguramos de que toda la aplicación
  # tenga acceso a ellas (dplyr, readxl, etc) sin tener que cambiar tu código.
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

  # Lanza la app usando las funciones que acabamos de crear
  shiny::shinyApp(ui = app_ui(), server = app_server)
}
