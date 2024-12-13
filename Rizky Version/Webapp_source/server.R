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

shinyFileChoose(input, "files", roots=c(wd="/home/dnanexus/project/"), filetypes = c("", "rds"))

observeEvent(input$files, {
    # Parse the selected file path
    fileinfo <- parseFilePaths(c(wd="/home/dnanexus/project/"), input$files)
    selected_file <- as.character(fileinfo$datapath)
    
    # Check if a file was selected
    if (length(selected_file) > 0 && file.exists(selected_file)) {
      # Load the selected RDS file
      seurat_obj <- readRDS(selected_file)
      
      # Display the selected file path
      output$loadedFile <- renderText({
        paste("Loaded file:", selected_file)
      })
    }
  })

reactivevalue$RDS_directory=seurat_obj

source('operator.R',local = T)

}
