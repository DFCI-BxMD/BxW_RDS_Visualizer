#increase filesize limit
base::options(shiny.maxRequestSize=1000000*1024^2) #max 10gb upload

#########           SERVER START          ##########


#define server
function(input, output, session) {

  ### input logic ###
  
  seuratData <- reactiveVal(NULL)

  shinyFileChoose(input, "files", roots=c(wd="../project"), filetypes = c("", "rds"))

  observeEvent(input$files, {
      # Parse the selected file path
      fileinfo <- parseFilePaths(c(wd="../project"), input$files)
      selected_file <- as.character(fileinfo$datapath)
      
      # Check if a file was selected
      if (length(selected_file) > 0 && file.exists(selected_file)) {
        # Load the selected RDS file
        seurat_obj <- readRDS(selected_file)
        
        # Store the loaded Seurat object in a reactive value
        seuratData(seurat_obj)
        
        # Display the selected file path
        output$loadedFile <- renderText({
          paste("Loaded file:", selected_file)
        })
      }
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
    gene_names <- strsplit(tolower(input$featuresInput), ",\\s*")[[1]] 
    gene_names <- trimws(gene_names)  # Trim leading and trailing spaces from gene names
    
    # Get original Seurat gene names
    original_seurat_gene_names <- rownames(seurat_data@assays$RNA@data)
    
    # Create a mapping of lowercase to original gene names
    seurat_gene_names_lower <- tolower(original_seurat_gene_names)
    gene_name_map <- setNames(original_seurat_gene_names, seurat_gene_names_lower)
    
    # Check for missing genes using lowercase comparison
    missing_genes <- setdiff(gene_names, seurat_gene_names_lower)
    
    if (length(missing_genes) > 0) {
      # output$nofeaturefound <- renderText(paste("The following genes were not found:", paste(missing_genes, collapse = ", "))) #this worked... trying to phase out in favor of warning message
      #warning message
      shinyalert::shinyalert(
        title = "Warning: Invalid Input",
        text = paste("The following genes were not found:", paste(missing_genes, collapse = ", ")),
        type = "warning")
    } else {

      matched_gene_names <- gene_name_map[gene_names]

      # Generate feature plot
      feature_plot <- FeaturePlot(object = seurat_data, features = matched_gene_names,reduction = "umap")
      
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
    gene_names <- strsplit(tolower(input$featuresInput), ",\\s*")[[1]] 
    gene_names <- trimws(gene_names)  # Trim leading and trailing spaces from gene names
    
    # Get original Seurat gene names
    original_seurat_gene_names <- rownames(seurat_data@assays$RNA@data)
    
    # Create a mapping of lowercase to original gene names
    seurat_gene_names_lower <- tolower(original_seurat_gene_names)
    gene_name_map <- setNames(original_seurat_gene_names, seurat_gene_names_lower)
    
    # Check for missing genes using lowercase comparison
    missing_genes <- setdiff(gene_names, seurat_gene_names_lower)
    
    if (length(missing_genes) > 0) {
      # output$nofeaturefound <- renderText(paste("The following genes were not found:", paste(missing_genes, collapse = ", "))) #this worked... trying to phase out in favor of warning message
      #warning message
      shinyalert::shinyalert(
        title = "Warning: Invalid Input",
        text = paste("The following genes were not found:", paste(missing_genes, collapse = ", ")),
        type = "warning")
    } else {

      matched_gene_names <- gene_name_map[gene_names]

      # Generate feature plot
      feature_plot <- FeaturePlot(object = seurat_data, features = matched_gene_names, reduction = "pca")
      
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
    
    # Get the gene names
    gene_names <- rownames(seurat_data@assays$RNA@data)
    
    # Set the random seed for reproducibility (optional)
    set.seed(123)
    
    # Select 10 random genes
    random_genes <- sample(gene_names, 10)
    
    # Return the random genes as a character string
    paste(random_genes, collapse = ", ")
  })
  
  
  ###update dimplot/violinplot dropdown###

  # Update dropdowns for DimPlot and ViolinPlot when Seurat object is loaded
  observe({
    req(seuratData()) # Ensure Seurat object is loaded
    seurat_data <- seuratData()
    
    # Fetch metadata column names
    obj_meta <- seurat_data@meta.data
    
    # Convert specified columns to factors if needed
    columns_to_check <- c("sample", "cytotoxic", "exact_subclonotype_id")
    for (col in columns_to_check) {
      if (col %in% colnames(obj_meta) && !is.factor(obj_meta[[col]])) {
        obj_meta[[col]] <- as.factor(obj_meta[[col]])
      }
    }
    
    # Assign column factors to categorical_cols var
    categorical_cols <- names(obj_meta)[sapply(obj_meta, function(col) is.factor(col) && nlevels(col) < 5)]
    
    # Reassign dropdown options
    updateSelectInput(session, "variableInput", label = "Select a Field to Split Dim Plot By", choices = categorical_cols)
    updateSelectInput(session, "violinInput", label = "Select a Field to Split Violin Plot By", choices = categorical_cols)
  })

  # DimPlot logic
  observeEvent(input$plotButton_Dim, {
    req(seuratData(), input$variableInput)
    seurat_data <- seuratData()
    
    if (is.null(seurat_data)) {
      output$errormessage <- renderText("Null file")
    } else {
      if (!input$splitToggle) {
        # Unsplit dim plot colored by group
        dim_plot <- DimPlot(object = seurat_data, reduction = "umap", group.by = input$variableInput)
      } else {
        # Split dim plot colored by group
        dim_plot <- DimPlot(object = seurat_data, reduction = "umap", split.by = input$variableInput, group.by = input$variableInput)
      }
      
      if (is.null(dim_plot)) {
        output$nodimplot <- renderText("Bad argument")
      } else {
        output$featurePlotDim <- renderPlot({ dim_plot })
      }
    }
  })

  # ViolinPlot logic
  observeEvent(input$plotButton_violin, {
    req(seuratData(), input$featuresInput_vln, input$violinInput)
    seurat_data <- seuratData()
    
    # Split the input string @ commas into individual gene names using regex
    gene_names <- strsplit(input$featuresInput_vln, ",\\s*")[[1]] 
    gene_names <- trimws(gene_names)
    
    # Check if any of the requested genes are missing
    missing_genes_vln <- setdiff(gene_names, rownames(seurat_data@assays$RNA@data))
    
    if (length(missing_genes_vln) > 0) {
      shinyalert::shinyalert(
        title = "Warning: Invalid Input",
        text = paste("The following genes were not found:", paste(missing_genes_vln, collapse = ", ")),
        type = "warning"
      )
    } else {
      # Generate violin plot
      violin_plot <- VlnPlot(object = seurat_data, features = gene_names, split.by = input$violinInput)
      
      if (is.null(violin_plot)) {
        output$badviolinplot <- renderText("Bad argument")
      } else {
        output$violinplot <- renderPlot({ violin_plot })
      }
    }
  })
}

#########           SERVER END          ##########
