#increase filesize limit
base::options(shiny.maxRequestSize=1000000*1024^2) #max 10gb upload
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(purrr)
library(stringr)
library(Seurat)
library(DT)
library(shinyalert)
library(ggplot2)


function(input, output, session) {

reactivevalue=reactiveValues(RDS_directory=NULL,
                             SeuratObject=NULL,
                             Loaded=F,
                             metadata=NULL,
                             genes_name=NULL)

rds_file <- list.files(pattern = "\\.rds$", full.names = TRUE)

reactivevalue$RDS_directory=rds_file
source('operator.R',local = T)

}



