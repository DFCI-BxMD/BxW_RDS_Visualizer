#### Check if Seurat Object exist and whether it is loaded: if yes and no, load, if yes and no, stop, if no and no, stop.####





observe(if (length(reactivevalue$RDS_directory)!=0&(!reactivevalue$Loaded)) {
  print(reactivevalue$Loaded)
  reactivevalue$Loaded=T
  reactivevalue$SeuratObject=readRDS(reactivevalue$RDS_directory[1])
  reactivevalue$SeuratObject=UpdateSeuratObject(reactivevalue$SeuratObject)
  reactivevalue$metadata=reactivevalue$SeuratObject@meta.data
  DefaultAssay(reactivevalue$SeuratObject)='RNA'
  reactivevalue$SeuratObject=NormalizeData(reactivevalue$SeuratObject)
  reactivevalue$genes_name=rownames(reactivevalue$SeuratObject)
  
  reactivevalue$reduction=names((reactivevalue$SeuratObject@reductions))
  reactivevalue$reduction=reactivevalue$reduction[!grepl('pca',reactivevalue$reduction,ignore.case = T)]
  reactivevalue$reduction=reactivevalue$reduction[!grepl('harmony',reactivevalue$reduction,ignore.case = T)]
  
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
  updateSelectizeInput(session = session,inputId = 'GeneInput',choices=reactivevalue$genes_name,selected = NULL,server = T)
  
}
)


BarGraphListener <- reactive({
  list(input$Bar_Graph_y,input$Bar_Graph_fill,input$Bar_Graph_Percentage)
})


observeEvent(BarGraphListener(),{
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
    }
  }
})





