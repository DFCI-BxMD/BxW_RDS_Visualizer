#increase filesize limit
options(shiny.maxRequestSize=10000*1024^2) #max 10gb upload

# Install required packages
packages <- c("shiny", "shinydashboard", "shinydashboardPlus", "purrr", "stringr", "Seurat", "DT", "shinyalert", "ggplot2")

# Check if each package is already installed, if not, install it
for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

# Load the installed packages
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(purrr)
library(stringr)
library(Seurat)
library(DT)
library(shinyalert)
library(ggplot2)
