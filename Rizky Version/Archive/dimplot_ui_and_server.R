#ui

tab_dim <- tabItem(tabName = "dim",
                   h2("Dim Plot"),
                   br(),
                   p("Select a variable to split the Dimplot By:"),
                   # selectInput("variableInput", label = NULL,
                   #             choices = c("Clusters", "Treatments", "Clonotype", "Cytotoxic" )),
                   selectInput("variableInput", label = "Select",
                               choices = "Input File For Dropdown Options"),
                   
                   br(),
                   p("Here are the available categories that the Dimplot can be split by:"),
                   br(),
                   div(verbatimTextOutput("sampleCols"), class = "full-width-cols"),
                   
                   # p("Metadata split.by params:"),
                   # br(),
                   # verbatimTextOutput("samplesplit"),
                   br(),
                   actionButton("plotButton_Dim","Generate DimPlot"),
                   br(),
                   plotOutput("featurePlotDim", height = 600, width = 900)
                   
)

#server

###update dimplot dropdown###
observeEvent(input$seuratFile, {
  req(seuratData()) # wait for file input
  seurat_data <- seuratData()
  
  #fetch metadata colnames
  obj_meta <- seurat_data@meta.data
  
  #categorical colnames only
  categorical_cols <- names(obj_meta)[sapply(obj_meta, is.factor)]
  
  #reassign dropdown options
  updateSelectInput(session, "variableInput",label = "Select a Field to Split Dimplot By", choices=categorical_cols)
  
})

###dimplot logic### 
observeEvent(input$plotButton_Dim, {
  req(seuratData()) # wait for file input
  seurat_data <- seuratData()
  if (is.null(seurat_data)) {
    output$errormessage <- renderText("Null file")
  } else {
    # variable <- switch(input$variableInput,
    #                    "Clusters" = "seurat_clusters",
    #                    "Treatments" = "sample",
    #                    "Clonotype" = "exact_subclonotype_id",
    #                    "Cytotoxic" = "cytotoxic")
    
    # Generate dim plot
    dim_plot <- DimPlot(object = seurat_data, reduction = "umap", split.by = input$variableInput) 
    
    # make sure dimplot is not null before rendering plot
    if (is.null(dim_plot)) {
      output$nodimplot <- renderText("Bad argument")
    } else {
      output$featurePlotDim <- renderPlot({ dim_plot })
    }
  }
})
