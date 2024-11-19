#increase filesize limit
base::options(shiny.maxRequestSize=1000000*1024^2) #max 10gb upload

#########           SERVER START          ##########


#define server
function(input, output, session) {

  ### input logic ###
  
  # Read Seurat object
  seuratData <- reactive({
    rds_file <- list.files(pattern = "\\.rds$", full.names = TRUE)[1]
    if (!is.na(rds_file)) {
        seurat_obj <- readRDS(rds_file)
    } else {
      stop("No .rds file found in the current directory.")
    }
    return(seurat_obj)
  })
  
  # Render table-categorical only
  output$table_cat <- DT::renderDataTable({
    #retrieve seurat object from reactive expression in above code block
    seurat_data <- seuratData() 
    if (is.null(seurat_data))
      return(NULL)
    
    #get metadata
    meta_cat <- as.data.frame(seurat_data@meta.data)
    
    # include these collumns if they dont exist
    columns_to_check <- c("sample", "cytotoxic", "exact_subclonotype_id")
    
    # Find and convert columns that match the specified names
    for (col in columns_to_check) {
      if (col %in% colnames(meta_cat) && !is.factor(meta_cat[[col]])) {
        meta_cat[[col]] <- as.factor(meta_cat[[col]])
      }
    }
    
    #filter metadata by categorical
    categorical_metadata <- meta_cat[, sapply(meta_cat, is.factor)]
    
    #table of cat metadata
    datatable(categorical_metadata, options=list(scrollX=TRUE))
  })
  
  # Render table-metadata
  output$table_meta <- DT::renderDataTable({
    #retrieve seurat object from reactive expression in above code block
    seurat_data <- seuratData() 
    if (is.null(seurat_data))
      return(NULL)
    # Metadata
    seurat_df <- as.data.frame(seurat_data@meta.data)
    datatable(seurat_df, options=list(scrollX=TRUE))
  })
  
  ### featureplot logic ###
  #umap featureplot logic
  observeEvent(input$plotButton_UMAP, {
    req(seuratData(),input$featuresInput) #wait for file input and feature input
    seurat_data <- seuratData()
    if (is.null(seurat_data))
      output$errormessage <- renderText("Null file")
    
    # Split the input string @ commas into individual gene names using regex
    #[[1]] extracts the first element in list of gene names, which gives us a chr vector 
    #containing the individaul gene names as separate elements
    gene_names <- strsplit(input$featuresInput, ",\\s*")[[1]] 
    gene_names <- trimws(gene_names)  # Trim leading and trailing spaces from gene names
    
    # Check if any of the requested genes are missing
    missing_genes <- setdiff(gene_names, rownames(seurat_data@assays$RNA@data))
    
    if (length(missing_genes) > 0) {
      # output$nofeaturefound <- renderText(paste("The following genes were not found:", paste(missing_genes, collapse = ", "))) #this worked... trying to phase out in favor of warning message
      #warning message
      shinyalert::shinyalert(
        title = "Warning: Invalid Input",
        text = paste("The following genes were not found:", paste(missing_genes, collapse = ", ")),
        type = "warning")
    } else {
      # Generate feature plot
      feature_plot <- FeaturePlot(object = seurat_data, features = gene_names,reduction = "umap")
      
      #assign featureplot to global env
      assign("feature_plot", "new", envir = .GlobalEnv) #this doesnt seem to work... try declaring reactive vals at top of server script
      
      #make sure featureplot is not null before rendering plot
      if (is.null(feature_plot)) {
        output$nogenesfound <- renderText("No genes found!")
      } else {
        output$featurePlotUMAP <- renderPlot({ feature_plot })
      }
    }
  })
  
  #pca featureplot logic
  observeEvent(input$plotButton_PCA, {
    req(seuratData(),input$featuresInput) #wait for file input and feature input
    seurat_data <- seuratData()
    if (is.null(seurat_data))
      output$errormessage <- renderText("Null file")
    
    # Split the input string @ commas into individual gene names using regex
    #[[1]] extracts the first element in list of gene names, which gives us a chr vector containing the individaul gene names as separate elements
    gene_names <- strsplit(input$featuresInput, ",\\s*")[[1]] 
    gene_names <- trimws(gene_names)  # Trim leading and trailing spaces from gene names
    
    # Check if any of the requested genes are missing
    missing_genes <- setdiff(gene_names, rownames(seurat_data@assays$RNA@data))
    
    if (length(missing_genes) > 0) {
      # output$nofeaturefound <- renderText(paste("The following genes were not found:", paste(missing_genes, collapse = ", "))) #this worked... trying to phase out in favor of warning message
      #warning message
      shinyalert::shinyalert(
        title = "Warning: Invalid Input",
        text = paste("The following genes were not found:", paste(missing_genes, collapse = ", ")),
        type = "warning")
    } else {
      # Generate feature plot
      feature_plot <- FeaturePlot(object = seurat_data, features = gene_names, reduction = "pca")
      
      #assign featureplot to global env
      # assign("feature_plot", "new", envir = .GlobalEnv) #this doesnt seem to work... try declaring reactive vals at top of server script
      
      #make sure featureplot is not null before rendering plot
      if (is.null(feature_plot)) {
        output$nogenesfound <- renderText("No genes found!")
      } else {
        output$featurePlotPCA <- renderPlot({ feature_plot })
      }
    }
  })
  
  #sampleGene logic
  output$sampleGenes <- renderText({
    seurat_data <- seuratData()
    if (is.null(seurat_data))
      return(NULL)
    # Access the RNA assay
    rna_assay <- seurat_data@assays$RNA
    
    # Get the gene names
    gene_names <- rownames(rna_assay@data)
    
    # Set the random seed for reproducibility (optional)
    set.seed(123)
    
    # Select 10 random genes
    random_genes <- sample(gene_names, 10)
    
    # Return the random genes as a character string
    paste(random_genes, collapse = ", ")
  })
  
  
  ###update dimplot/violinplot dropdown###
  observeEvent(input$seuratFile, {
    req(seuratData()) # wait for file input
    seurat_data <- seuratData()
    
    #fetch metadata colnames
    obj_meta <- seurat_data@meta.data
    
    # Specify the metadata columns you want to check and convert
    columns_to_check <- c("sample", "cytotoxic", "exact_subclonotype_id")
    
    # Find and convert columns that match the specified names
    for (col in columns_to_check) {
      if (col %in% colnames(obj_meta) && !is.factor(obj_meta[[col]])) {
        obj_meta[[col]] <- as.factor(obj_meta[[col]])
      }
    }
    
    #assign column factors to categorical_cols var
    # categorical_cols <- names(obj_meta)[sapply(obj_meta, is.factor)]
    categorical_cols <- names(obj_meta)[sapply(obj_meta, function(col) is.factor(col) && nlevels(col) < 5)]
    
    #reassign dropdown options
    updateSelectInput(session, "variableInput",label = "Select a Field to Split Dim Plot By", choices=categorical_cols)
    
    updateSelectInput(session, "violinInput",label = "Select a Field to Split Violin Plot By", choices = categorical_cols)
    
  })
  
  #new dimplot logic
  observeEvent(input$plotButton_Dim, {
    req(seuratData(), input$variableInput)
    seurat_data <- seuratData()
    
    if (is.null(seurat_data)) {
      output$errormessage <- renderText("Null file")
    } else {
      if (!input$splitToggle) {
        # Unsplit dim plot colored by group
        dim_plot <- DimPlot(object = seurat_data, reduction = "umap", group.by = input$variableInput)
        
        if (is.null(dim_plot)) {
          output$nodimplot <- renderText("Bad argument")
        } else {
          output$featurePlotDim <- renderPlot({ dim_plot })
          output$splitToggle <- renderUI({ NULL })
        }
      } else {
        # Split dim plot colored by group
        dim_plot <- DimPlot(object = seurat_data, reduction = "umap", split.by = input$variableInput, group.by = input$variableInput)
        
        if (is.null(dim_plot)) {
          output$nodimplot <- renderText("Bad argument")
        } else {
          output$featurePlotDim <- renderPlot({ dim_plot })
          output$splitToggle <- renderUI({ checkboxInput("splitToggle", "Split Dim Plot", value = TRUE) })
        }
      }
    }
  })
  
  
  ###
  
  #update violinplot dropdown#
  
  
  ### violinplot logic ###
  observeEvent(input$plotButton_violin, {
    req(seuratData(),input$featuresInput_vln,input$violinInput) #wait for all inputs
    seurat_data <- seuratData()
    
    # Split the input string @ commas into individual gene names using regex
    gene_names <- strsplit(input$featuresInput_vln, ",\\s*")[[1]] 
    gene_names <- trimws(gene_names)  # Trim leading and trailing spaces from gene names
    
    # Check if any of the requested genes are missing
    missing_genes_vln <- setdiff(gene_names, rownames(seurat_data@assays$RNA@data))
    
    if (length(missing_genes_vln) > 0) {
      #warning message
      shinyalert::shinyalert(
        title = "Warning: Invalid Input",
        text = paste("The following genes were not found:", paste(missing_genes_vln, collapse = ", ")),
        type = "warning")
    } else {
      # Generate feature plot
      violin_plot <- VlnPlot(object = seurat_data, features = gene_names, split.by =input$violinInput)
      
      #make sure featureplot is not null before rendering plot
      if (is.null(violin_plot)) {
        output$badviolinplot <- renderText("Bad argument")
      } else {
        output$violinplot <- renderPlot({ violin_plot })
      }
    }
  })
} #end of server

#########           SERVER END          ##########
