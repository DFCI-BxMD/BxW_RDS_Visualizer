#### Check if Seurat Object exist and whether it is loaded: if yes and no, load, if yes and no, stop, if no and no, stop.####





observe(if (length(reactivevalue$RDS_directory)!=0&(!reactivevalue$Loaded)) {
  print(reactivevalue$Loaded)
  reactivevalue$SeuratObject=readRDS(reactivevalue$RDS_directory)
  reactivevalue$SeuratObject=UpdateSeuratObject(reactivevalue$SeuratObject)
  reactivevalue$metadata=reactivevalue$SeuratObject@meta.data
  DefaultAssay(reactivevalue$SeuratObject)='RNA'
  reactivevalue$SeuratObject=NormalizeData(reactivevalue$SeuratObject)
  reactivevalue$genes_name=rownames(reactivevalue$SeuratObject)
  
  reactivevalue$reduction=names((reactivevalue$SeuratObject@reductions))
  #reactivevalue$reduction=reactivevalue$reduction[!grepl('pca',reactivevalue$reduction,ignore.case = T)]
  #reactivevalue$reduction=reactivevalue$reduction[!grepl('harmony',reactivevalue$reduction,ignore.case = T)]
  
  output$MainFigure=renderPlot(DimPlot(reactivevalue$SeuratObject))
  
  
  output$Metadata=DT::renderDataTable(DT::datatable(reactivevalue$metadata,
                                                            options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),filter = list(position = "top")),server = T)
  
  
  Bar_Graph_Columns=c()
  for (i in colnames(reactivevalue$metadata)) {
    if (length(unique(reactivevalue$metadata[,i]))<=10) {
      Bar_Graph_Columns=c(Bar_Graph_Columns,i)
    }
  }
  
  updateSelectizeInput(session = session,inputId = 'Bar_Graph_y',choices =Bar_Graph_Columns,selected = NULL,server=T)
  updateSelectizeInput(session = session,inputId = 'Bar_Graph_fill',choices =Bar_Graph_Columns,selected = NULL,server=T)
  updateSelectizeInput(session = session,inputId = 'FeaturePlot_GeneInput',choices=reactivevalue$genes_name,selected = NULL,server = T)
  updateSelectizeInput(session = session,inputId = 'FeaturePlot_reduction',choices=reactivevalue$reduction,selected = NULL,server = T)
  
  updateSelectizeInput(session = session,inputId = 'DimPlot_group_by',choices =Bar_Graph_Columns,selected = NULL,server=T)
  updateSelectizeInput(session = session,inputId = 'DimPlot_split_by',choices =c(Bar_Graph_Columns,''),selected = '',server=T)
  updateSelectizeInput(session = session,inputId = 'DimPlot_reduction',choices =reactivevalue$reduction,selected = NULL,server=T)
  
  updateSelectizeInput(session = session,inputId = 'VlnPlot_GeneInput',choices=reactivevalue$genes_name,selected = NULL,server = T)
  updateSelectizeInput(session = session,inputId = 'VlnPlot_group_by',choices=Bar_Graph_Columns,selected = NULL,server = T)
  
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


plotDimplot=eventReactive(input$plotDimPlot_Button, {
    if (reactivevalue$Loaded) {
    if (input$DimPlot_split_by!=''){
      if (length(unique(reactivevalue$metadata[,input$DimPlot_split_by]))>2) {
        number_col=round(sqrt(length(unique(reactivevalue$metadata[,input$DimPlot_split_by]))))
      } else {
        number_col=length(unique(reactivevalue$metadata[,input$DimPlot_split_by]))
      }
      plot=(DimPlot(reactivevalue$SeuratObject,group.by = input$DimPlot_group_by,split.by = input$DimPlot_split_by,
                             ncol = number_col,reduction = input$DimPlot_reduction))
    } else {
      plot=(DimPlot(reactivevalue$SeuratObject,group.by = input$DimPlot_group_by,reduction = input$DimPlot_reduction))
    }
    
  }

  output$DimPlot=renderPlot(plot)
  reactivevalue$dimPlot = plot
})

observe(plotDimplot())

# Download handler for the plot
output$downloadDimPlot <- downloadHandler(
  filename = function() {
    paste("Dim_plot", Sys.Date(), ".pdf", sep = "")
  },
    content = function(file) {
        pdf(file)
        print(reactivevalue$dimPlot)
        dev.off()
    })






plotFeaturePlot=eventReactive(input$plotFeaturePlot_Button, {
  if (reactivevalue$Loaded) {
    
  plot=FeaturePlot(reactivevalue$SeuratObject,features = input$FeaturePlot_GeneInput,reduction = input$FeaturePlot_reduction,order = T)
  
  output$FeaturePlot=renderPlot(plot)

  reactivevalue$featurePlot = plot
  }
})

observe(plotFeaturePlot())

# Download handler for the plot
output$downloadFeaturePlot <- downloadHandler(
  filename = function() {
    paste("feature_plot", Sys.Date(), ".pdf", sep = "")
  },
  content = function(file) {
        pdf(file)
        print(reactivevalue$featurePlot)
        dev.off()
    })




plotVlnPlot=eventReactive(input$plotVlnPlot_Button, {
  if (reactivevalue$Loaded) {
    
  if (length(input$VlnPlot_GeneInput)>2) {
    number_of_cols=round(sqrt(length(input$VlnPlot_GeneInput)))
  } else {
    number_of_cols=length(input$VlnPlot_GeneInput)
    
  }
  plot=VlnPlot(reactivevalue$SeuratObject,features = input$VlnPlot_GeneInput,group.by = input$VlnPlot_group_by,ncol = number_of_cols,same.y.lims = T,raster = T)
  
  output$VlnPlot=renderPlot(plot)
  reactivevalue$VlnPlot = plot
  }
})

observe(plotVlnPlot())

# Download handler for the plot
output$downloadVlnPlot <- downloadHandler(
  filename = function() {
    paste("Violin_plot", Sys.Date(), ".pdf", sep = "")
  },
  content = function(file) {
        pdf(file)
        print(reactivevalue$VlnPlot)
        dev.off()
    })