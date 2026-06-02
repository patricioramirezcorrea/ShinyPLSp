# ==============================================================================
# LIBRERIAS Y CONFIGURACION GLOBAL
# ==============================================================================

# Core de Shiny y UI
library(shiny)
library(bslib)
library(DT)

# Manipulación y lectura de datos
library(dplyr)
library(readxl)
library(openxlsx)
library(digest)

# Estimación SEM y Segmentación
library(cSEM)
library(genpathmox)

# Visualización y Machine Learning (WPI / Árboles)
library(DiagrammeR)
library(ggplot2)
library(ggforce)
library(Rtsne)

# Útil en servidores Linux para evitar errores de renderizado de gráficos
options(bitmapType = "cairo")
