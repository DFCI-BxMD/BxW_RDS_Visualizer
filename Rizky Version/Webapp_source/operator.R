#### Check if Seurat Object exist and whether it is loaded: if yes and no, load, if yes and no, stop, if no and no, stop.####


waitress <- Waitress$new(theme = "overlay-percent")


observe(if (length(reactivevalue$RDS_directory)!=0&(!reactivevalue$Loaded)) {
  print(reactivevalue$Loaded)
  reactivevalue$SeuratObject=readRDS(reactivevalue$RDS_directory)
  reactivevalue$SeuratObject=UpdateSeuratObject(reactivevalue$SeuratObject)
  reactivevalue$metadata=reactivevalue$SeuratObject@meta.data
  DefaultAssay(reactivevalue$SeuratObject)='RNA'
  reactivevalue$SeuratObject=NormalizeData(reactivevalue$SeuratObject)
  reactivevalue$genes_name=rownames(reactivevalue$SeuratObject)
  
  reactivevalue$reduction=names((reactivevalue$SeuratObject@reductions))


  # Convert Seurat object to soma experiment

  # EXPERIMENT_URI <- sprintf("dx://project-xxxx:/path/to/file")
  # write_soma(reactivevalue$SeuratObject, uri = EXPERIMENT_URI)

  #reactivevalue$reduction=reactivevalue$reduction[!grepl('pca',reactivevalue$reduction,ignore.case = T)]
  #reactivevalue$reduction=reactivevalue$reduction[!grepl('harmony',reactivevalue$reduction,ignore.case = T)]
  
  output$MainFigure=renderPlot(DimPlot(reactivevalue$SeuratObject)+coord_fixed()+ theme(plot.title = element_text(size = 20), 
                                                                                                    axis.title = element_text(size = 18, face = "bold"), 
                                                                                                    axis.text = element_text(size = 17),
                                                                                                    legend.title = element_text(size = 18),
                                                                                                    legend.text = element_text(size = 17, face = "bold")))
  
  
  output$Metadata=DT::renderDataTable(DT::datatable(reactivevalue$metadata,
                                                            options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),filter = list(position = "top")),server = T)
  
  
  Bar_Graph_Columns=c()
  for (i in colnames(reactivevalue$metadata)) {
    if (typeof((reactivevalue$metadata[,i]))!='double') {
      if (typeof(reactivevalue$metadata[,i])=='character') {
        Bar_Graph_Columns=c(Bar_Graph_Columns,i)
        
      } else {
        if (is.factor(reactivevalue$metadata[,i])) {
          Bar_Graph_Columns=c(Bar_Graph_Columns,i)
          
        }
      }
    }
  }
  Numerical_features=c()
  for (i in colnames(reactivevalue$metadata)) {
    if (typeof((reactivevalue$metadata[,i]))=='double') {
      Numerical_features=c(Numerical_features,i)
    } else {
      if (typeof((reactivevalue$metadata[,i]))=='integer') {
        
        if (!is.factor(reactivevalue$metadata[,i])) {
          
          Numerical_features=c(Numerical_features,i)
        }
        
      }
    }
  }
  
  updateSelectizeInput(session = session,inputId = 'Bar_Graph_y',choices =Bar_Graph_Columns,selected = NULL,server=T)
  updateSelectizeInput(session = session,inputId = 'Bar_Graph_fill',choices =Bar_Graph_Columns,selected = NULL,server=T)

  # Feature Plots
  default_reduction <- ifelse("umap" %in% reactivevalue$reduction, "umap", 
                                ifelse("pca" %in% reactivevalue$reduction, "pca", 
                                       reactivevalue$reduction[1]))
    
  updateSelectizeInput(session = session,inputId = 'FeaturePlot_GeneInput',choices=reactivevalue$genes_name,selected = NULL,server = T)
  updateSelectizeInput(session = session,inputId = 'FeaturePlot_MetaInput',choices=Numerical_features,selected = NULL,server = T)
  updateSelectizeInput(session = session,inputId = 'GeneFeaturePlot_reduction',choices=reactivevalue$reduction,selected = default_reduction,server = T)
  updateSelectizeInput(session = session,inputId = 'MetaFeaturePlot_reduction',choices=reactivevalue$reduction,selected = default_reduction,server = T)
  
  updateSelectizeInput(session = session,inputId = 'DimPlot_group_by',choices =Bar_Graph_Columns,selected = NULL,server=T)
  updateSelectizeInput(session = session,inputId = 'DimPlot_split_by',choices =c(Bar_Graph_Columns,''),selected = '',server=T)
  updateSelectizeInput(session = session,inputId = 'DimPlot_reduction',choices =reactivevalue$reduction,selected = default_reduction,server=T)
  
  updateSelectizeInput(session = session,inputId = 'VlnPlot_GeneInput',choices=reactivevalue$genes_name,selected = NULL,server = T)
  updateSelectizeInput(session = session,inputId = 'VlnPlot_MetaInput',choices=Numerical_features,selected = NULL,server = T)
  
  updateSelectizeInput(session = session,inputId = 'GeneVlnPlot_group_by',choices=Bar_Graph_Columns,selected = NULL,server = T)
  updateSelectizeInput(session = session,inputId = 'MetaVlnPlot_group_by',choices=Bar_Graph_Columns,selected = NULL,server = T)
    
  
  reactivevalue$Loaded=T
}
)


BarGraphListener <- reactive({
  list(input$Bar_Graph_y,input$Bar_Graph_fill,input$Bar_Graph_Percentage)
})


observeEvent(BarGraphListener(),{
  
  if (reactivevalue$Loaded) {
  if (!is.null(reactivevalue$metadata)&input$Bar_Graph_y!=''&input$Bar_Graph_fill!='') {
    if (input$Bar_Graph_y!=input$Bar_Graph_fill){
      Reference=data.frame(table(reactivevalue$metadata[,input$Bar_Graph_fill]))
      colnames(Reference)=c("Variable2",'CellNumber')
      Reference$Variable1='Reference'
      Reference=Reference[,c("Variable1","Variable2",'CellNumber')]
      temp=data.frame(table(reactivevalue$metadata[,input$Bar_Graph_y],reactivevalue$metadata[,input$Bar_Graph_fill]))
      colnames(temp)=c("Variable1","Variable2",'CellNumber')
      temp=rbind(Reference,temp)
      temp$Variable1=factor(temp$Variable1)
      temp$Variable1=relevel(temp$Variable1,ref = 'Reference')
      if (input$Bar_Graph_Percentage){
        output$BarGraph=renderPlot(ggplot(temp,aes(
          x=CellNumber,y=Variable1,fill=Variable2
        ))+geom_bar(stat = 'identity',position = 'fill')+xlab('Cell Number')+ylab(input$BarGraph1)+labs(fill=input$BarGraph2))} else {
          output$BarGraph=renderPlot(ggplot(temp,aes(
            x=CellNumber,y=Variable1,fill=Variable2
          ))+geom_bar(stat = 'identity',position = 'stack')+xlab('Cell Number')+ylab(input$BarGraph1)+labs(fill=input$BarGraph2) )
        }
      output$BarGraph_Table=DT::renderDataTable(DT::datatable(temp,
                                                              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),filter = list(position = "top")),server = T)
    } else {
      temp=data.frame(table(reactivevalue$metadata[,input$Bar_Graph_y]))
      colnames(temp)=c("Variable1",'CellNumber')
      if (input$Bar_Graph_Percentage){
        output$BarGraph=renderPlot(ggplot(temp,aes(
          x=CellNumber,y=Variable1
        ))+geom_bar(stat = 'identity',position = 'fill')+xlab('Cell Number')+ylab(input$BarGraph1)+labs(fill=input$BarGraph2))} else {
          output$BarGraph=renderPlot(ggplot(temp,aes(
            x=CellNumber,y=Variable1
          ))+geom_bar(stat = 'identity',position = 'stack')+xlab('Cell Number')+ylab(input$BarGraph1)+labs(fill=input$BarGraph2) )
        }
      output$BarGraph_Table=DT::renderDataTable(DT::datatable(temp,
                                                        options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),filter = list(position = "top")),server = T)
      
    }
  }
  }
})


## Dim Plot 
plotDimplot=eventReactive(input$plotDimPlot_Button, {
    if (reactivevalue$Loaded) {
        waitress$start() 
        if (input$DimPlot_split_by!=''){
            if (length(unique(reactivevalue$metadata[,input$DimPlot_split_by]))>2) {
                number_col=round(sqrt(length(unique(reactivevalue$metadata[,input$DimPlot_split_by]))))
            } else {
                number_col=length(unique(reactivevalue$metadata[,input$DimPlot_split_by]))
            }
            plot=(DimPlot(reactivevalue$SeuratObject,group.by = input$DimPlot_group_by,split.by = input$DimPlot_split_by,
                                    ncol = number_col,reduction = input$DimPlot_reduction) + coord_fixed()+ theme(plot.title = element_text(size = 20), 
                                                                                                                    axis.title = element_text(size = 18, face = "bold"), 
                                                                                                                    axis.text = element_text(size = 17),
                                                                                                                    legend.title = element_text(size = 18),
                                                                                                                    legend.text = element_text(size = 17, face = "bold")))
            } else {
            plot=(DimPlot(reactivevalue$SeuratObject,group.by = input$DimPlot_group_by,reduction = input$DimPlot_reduction) + coord_fixed()+ theme(plot.title = element_text(size = 20), 
                                                                                                                                                    axis.title = element_text(size = 18, face = "bold"), 
                                                                                                                                                    axis.text = element_text(size = 17),
                                                                                                                                                    legend.title = element_text(size = 18),
                                                                                                                                                    legend.text = element_text(size = 17, face = "bold")))
        }

        waitress$close()
    
  }

  output$DimPlot=renderPlot(plot)
  reactivevalue$dimPlot = plot
})

observe(plotDimplot())

# Save the plot in mounted project
shinyFileSave(input, "saveDimPlot", roots =c(wd="/home/dnanexus/project/"), filetypes = c("pdf"))
  observeEvent(input$saveDimPlot, {
    fileinfo <- parseSavePath(c(wd="/home/dnanexus/project/"), input$saveDimPlot)
    if (nrow(fileinfo) > 0) {
      pdf(fileinfo$datapath) 
        print(reactivevalue$dimPlot)
        dev.off()
    }
    })
    


## Feature Plot
plotGeneFeaturePlot=eventReactive(input$plotGeneFeaturePlot_Button, {
  if (reactivevalue$Loaded) {
    waitress$start()

   if (length(unique(reactivevalue$metadata[,input$FeaturePlot_GeneInput]))>1) {
                number_col=round(sqrt(length(unique(reactivevalue$metadata[,input$FeaturePlot_GeneInput]))))
            } else {
                number_col=length(unique(reactivevalue$metadata[,input$FeaturePlot_GeneInput]))
            }

    
    plots <- FeaturePlot(reactivevalue$SeuratObject,
                           features = c(input$FeaturePlot_GeneInput),
                           reduction = input$GeneFeaturePlot_reduction,
                           order = T, ncol = number_col) &
        theme(plot.title = element_text(size = 20),
              axis.title = element_text(size = 18, face = "bold"),
              axis.text = element_text(size = 17),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 17, face = "bold"))
      
      output$GeneFeaturePlot <- renderPlot(plots)
    reactivevalue$featurePlot = plots
    
    waitress$close()
  }
})


observe(plotGeneFeaturePlot())

plotMetaFeaturePlot=eventReactive(input$plotMetaFeaturePlot_Button, {
   if (reactivevalue$Loaded) {
    waitress$start()

    plots <- FeaturePlot(reactivevalue$SeuratObject,
                           features = c(input$FeaturePlot_MetaInput),
                           reduction = input$MetaFeaturePlot_reduction,
                           order = T, ncol = 3) &
        theme(plot.title = element_text(size = 20),
              axis.title = element_text(size = 18, face = "bold"),
              axis.text = element_text(size = 17),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 17, face = "bold"))
      
      output$MetaFeaturePlot <- renderPlot({
        print(plots) 
      }, width = 1300, height = 1200) 
    reactivevalue$featurePlot = plots
    
    waitress$close()
  }
})


observe(plotMetaFeaturePlot())

shinyFileSave(input, "saveFeaturePlot", roots =c(wd="/home/dnanexus/project/"), filetypes = c("pdf"))
  observeEvent(input$saveFeaturePlot, {
    fileinfo <- parseSavePath(c(wd="/home/dnanexus/project/"), input$saveFeaturePlot)
    if (nrow(fileinfo) > 0) {
      pdf(fileinfo$datapath) 
        print(reactivevalue$featurePlot)
        dev.off()
    }
    })



## Violin Plot
plotGeneVlnplot=eventReactive(input$plotGeneVlnPlot_Button, {
  if (reactivevalue$Loaded) {
  if (length(input$VlnPlot_GeneInput) > 2) {
    number_col = round(sqrt(length(input$VlnPlot_GeneInput)))
  } else {
    number_col = length(input$VlnPlot_GeneInput)
  }

  group_by <- ifelse(input$GeneVlnPlot_group_by == "", NULL, input$GeneVlnPlot_group_by)

  plot = VlnPlot(
    reactivevalue$SeuratObject,
    features = c(input$VlnPlot_GeneInput),
    group.by = group_by,
    ncol = number_col,
    same.y.lims = TRUE,
    raster = TRUE
  ) & theme(plot.title = element_text(size = 20), 
                axis.title = element_text(size = 18, face = "bold"), 
                axis.text = element_text(size = 17),
                legend.title = element_text(size = 18),
                legend.text = element_text(size = 17, face = "bold"))
}

output$GeneVlnPlot = renderPlot(plot)
reactivevalue$vlnPlot = plot
 
})

observe(plotGeneVlnplot())


plotMetaVlnplot=eventReactive(input$plotMetaVlnPlot_Button, {
  if (reactivevalue$Loaded) {
  if (length(input$VlnPlot_MetaInput) > 2) {
    number_col = round(sqrt(length(input$VlnPlot_MetaInput)))
  } else {
    number_col = length(input$VlnPlot_MetaInput)
  }

  group_by <- ifelse(input$MetaVlnPlot_group_by == "", NULL, input$MetaVlnPlot_group_by)

  plot = VlnPlot(
    reactivevalue$SeuratObject,
    features = c(input$VlnPlot_MetaInput),
    group.by = group_by,
    ncol = number_col,
    same.y.lims = TRUE,
    raster = TRUE
  ) & theme(plot.title = element_text(size = 20), 
                axis.title = element_text(size = 18, face = "bold"), 
                axis.text = element_text(size = 17),
                legend.title = element_text(size = 18),
                legend.text = element_text(size = 17, face = "bold"))
}

output$MetaVlnPlot = renderPlot(plot)
reactivevalue$vlnPlot = plot
 
})

observe(plotMetaVlnplot())

shinyFileSave(input, "saveViolinPlot", roots =c(wd="/home/dnanexus/project/"), filetypes = c("pdf"))
  observeEvent(input$saveViolinPlot, {
    fileinfo <- parseSavePath(c(wd="/home/dnanexus/project/"), input$saveViolinPlot)
    if (nrow(fileinfo) > 0) {
      pdf(fileinfo$datapath) 
        print(reactivevalue$VlnPlot)
        dev.off()
    }
    })