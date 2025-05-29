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
library(waiter)
library(tiledbsoma)
library(tiledb)
library(ggplotify)
library(patchwork)




function(input, output, session) {

reactivevalue=reactiveValues(RDS_directory=NULL,
                             SeuratObject=NULL,
                             Loaded=F,
                             metadata=NULL,
                             genes_name=NULL
                             )

waitress <- Waitress$new(theme = "overlay-percent")

shinyFileChoose(input, "files", 
                roots=c(wd="/home/dnanexus/project/"), 
                #roots=c(wd="/Users/zhaorong/Maynard/"), 
                
                filetypes = c("", "rds"))

observeEvent(input$files, {

  # Parse the selected file path
  fileinfo <- parseFilePaths(
    c(wd = "/home/dnanexus/project/"), 
    #c(wd = "/Users/zhaorong/Maynard/"), 
    
    input$files)
  selected_file <- as.character(fileinfo$datapath)
  
  # Check if a file was selected
  if (length(selected_file) > 0 && file.exists(selected_file)) {

    waitress$start()
    
    # Store the file path in reactive values
    reactivevalue$RDS_directory <- selected_file
    
    # Load the selected RDS file
    seurat_obj <- readRDS(selected_file)
    reactivevalue$SeuratObject <- seurat_obj

    for(i in 1:10){
      Sys.sleep(.5)
      waitress$inc(10) 
    }
    
    output$loadedFile <- renderText({
      paste("Loaded file:", selected_file)
    })
    
    reactivevalue$Loaded <- FALSE

    waitress$close()

   
  }
})


source('operator.R',local = T)

}
