  #### Check if Seurat Object exist and whether it is loaded: if yes and no, load, if yes and no, stop, if no and no, stop.####
  volcano.plot.function=function(input=NULL,
                                 logfc.axis.name=NULL,
                                 pvalue.axis.name=NULL,
                                 name.axis.name=NULL,
                                 x.threshold=0.5,
                                 y.threshold=0.05,
                                 genes.to.label=NULL
  ) {
    temp=data.frame(xaxis=input[,logfc.axis.name],
                    yaxis=input[,pvalue.axis.name],
                    rawyaxis=input[,pvalue.axis.name],
                    name=input[,name.axis.name])
    temp$yaxis=ifelse(temp$yaxis==0,yes=.Machine$double.xmin,
                      no=temp$yaxis)
    temp$yaxis=-log10(temp$yaxis)
    temp$score=temp$yaxis*temp$xaxis
    temp$score=ifelse(abs(temp$xaxis)>=1&temp$rawyaxis<y.threshold,yes=temp$score,no = 0)
    score=temp$score
    names(score)=temp$name
    score=score[order(score,decreasing = T)]
    label=genes.to.label
    temp$label=ifelse(temp$name%in%label,yes=temp$name,no='')
    
    require(ggplot2)
    require(ggrepel)
    temp$group='Insignificant'
    temp$group=ifelse(temp$xaxis>x.threshold&temp$rawyaxis<y.threshold,yes='UpRegulated',no=temp$group)
    temp$group=ifelse(temp$xaxis<(-1)*x.threshold&temp$rawyaxis<y.threshold&temp$group=='Insignificant',yes='DownRegulated',no=temp$group)
    if (!is.null(genes.to.label)) {
      temp$group=ifelse(temp$name%in%genes.to.label,yes='Testing',no=temp$group)
      for (i in 1:nrow(temp)) {
        if (temp$group[i]=='Testing') {
          if (temp$rawyaxis[i]<y.threshold) {
            if (temp$xaxis[i]>x.threshold) {
              temp$group[i]='UpRegulated'
            } else if (temp$xaxis[i]<((-1)*x.threshold)) {
              temp$group[i]='DownRegulated'
            }
            else {
              temp$group[i]='Insignificant'
              
            }
          } else {
            temp$group[i]='Insignificant'
            
          }
        }
      }
      
      
    }
    colorscheme=c('darkred','grey','darkblue')
    names(colorscheme)=c('UpRegulated','Insignificant','DownRegulated')
    
    #temp$label=ifelse(as.character(temp$group)=='Insignificant',yes='',no=temp$label)
    plot=ggplot(temp,aes(x=xaxis,y=yaxis,colour=group,label=label))+geom_point()+geom_text_repel(box.padding = 0.5, max.overlaps = Inf)+xlim(
      max(abs(temp$xaxis))*(-1),
      max(abs(temp$xaxis))*(1)
    )+ylim(0,max(temp$yaxis)*1.1)+xlab(logfc.axis.name)+ylab(paste0('-log10(',pvalue.axis.name,')'))+geom_vline(xintercept = x.threshold, linetype="dashed")+geom_vline(xintercept = (-1)*x.threshold, linetype="dashed")+
      geom_hline(yintercept = -log10(y.threshold))+theme_bw()
    plot=plot+scale_colour_manual(values = colorscheme)
    
    
    
    return(plot)
    
  }
  
  waitress <- Waitress$new(theme = "overlay-percent")
  
  
  observe(if (length(reactivevalue$RDS_directory)!=0&(!reactivevalue$Loaded)) {
    print(reactivevalue$Loaded)
    reactivevalue$SeuratObject=readRDS(reactivevalue$RDS_directory)
    reactivevalue$SeuratObject=UpdateSeuratObject(reactivevalue$SeuratObject)
    DefaultAssay(reactivevalue$SeuratObject)='RNA'
    reactivevalue$SeuratObject=NormalizeData(reactivevalue$SeuratObject)
    reactivevalue$metadata=reactivevalue$SeuratObject@meta.data
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
    
    # Differential Expression Analysis
    
    spliter='_thisisastringthatisspecificallydesignedforthesingleviewershinyappbyzhaorongifthisstringappearinthemetadataandcauseproblemthereisreallynotmuchzhaorongcando_'
    
    DGE_List=list()
    for (column in Bar_Graph_Columns) {
      for (value in unique(reactivevalue$metadata[,column])) {
        DGE_List[[paste0(column,' - ',value)]]=paste0(column,spliter,value)
        
      }
    }
    reactivevalue$spliter=spliter
    reactivevalue$DGE_List=DGE_List

    updateSelectizeInput(session = session,inputId = 'DGE_Group_1',choices=names(DGE_List),selected = NULL,server = T)
    updateSelectizeInput(session = session,inputId = 'DGE_Group_2',choices=names(DGE_List),selected = NULL,server = T)
    
    # Feature Plots
    default_reduction <- ifelse("umap" %in% reactivevalue$reduction, "umap", 
                                ifelse("tsne" %in% reactivevalue$reduction, "tsne", 
                                      ifelse("pca" %in% reactivevalue$reduction, "pca", 
                                       reactivevalue$reduction[1])))
    
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
      
      num_genes=length(input$FeaturePlot_GeneInput)
      
      plots=list()
      number_col=1
      if (length(input$FeaturePlot_GeneInput)>2) {
        number_col=ceiling(sqrt(length(input$FeaturePlot_GeneInput)))
      } else {
        number_col=length(unique(input$FeaturePlot_GeneInput))
      }
      number_col_contestant=ceiling(length(input$FeaturePlot_GeneInput)/number_col)
      if (number_col_contestant>number_col) {
        number_col=number_col_contestant
      }
      
      for (i in input$FeaturePlot_GeneInput) {
        if (all(as.numeric(reactivevalue$SeuratObject@assays$RNA$data[i,])==0)) {
          temp=FeaturePlot(reactivevalue$SeuratObject,features = i,reduction = input$GeneFeaturePlot_reduction,order = T,cols = c('grey','grey'))&theme(plot.title = element_text(size = 20),
                                                                                                                                                                axis.title = element_text(size = 18, face = "bold"),
                                                                                                                                                                axis.text = element_text(size = 17),
                                                                                                                                                                legend.title = element_text(size = 18),
                                                                                                                                                                legend.text = element_text(size = 17, face = "bold"))
          
        } else {
          temp=FeaturePlot(reactivevalue$SeuratObject,features = i,reduction = input$GeneFeaturePlot_reduction,order = T,cols = c('grey','blue'))&theme(plot.title = element_text(size = 20),
                                                                                                                                                                axis.title = element_text(size = 18, face = "bold"),
                                                                                                                                                                axis.text = element_text(size = 17),
                                                                                                                                                                legend.title = element_text(size = 18),
                                                                                                                                                                legend.text = element_text(size = 17, face = "bold"))
          
        }
        as_grob(temp)
        
        plots[[i]]=temp
      }
      
      #plots=FeaturePlot(reactivevalue$SeuratObject,features = c(input$FeaturePlot_GeneInput),reduction = input$GeneFeaturePlot_reduction,order = T)
      
      #as_grob(plots)
      plots=patchwork::wrap_plots(plots,ncol=number_col,nrow=number_col)
      
      

      
      output$GeneFeaturePlot <- renderPlot({
        print(plots) 
      }, width = reactive({
        
        if (num_genes == 1) { 
          500
        } else if (num_genes == 2) { 
          1000
        } else {
          1300
        }
      }), height = reactive({
        
        if (num_genes == 1) { 
          500
        } else if (num_genes == 2) { 
          500
        } else {
          1200
        }
      })) 
      reactivevalue$featurePlot = plots
      
      waitress$close()
    }
  })
  
  
  observe(plotGeneFeaturePlot())
  
  plotMetaFeaturePlot=eventReactive(input$plotMetaFeaturePlot_Button, {
    if (reactivevalue$Loaded) {
      waitress$start()
      
      num_genes <- length(input$FeaturePlot_MetaInput)
      
      
      
      plots=list()
      number_col=1
      if (length(input$FeaturePlot_MetaInput)>2) {
        number_col=ceiling(sqrt(length(input$FeaturePlot_MetaInput)))
      } else {
        number_col=length(unique(input$FeaturePlot_MetaInput))
      }
      number_col_contestant=ceiling(length(input$FeaturePlot_MetaInput)/number_col)
      if (number_col_contestant>number_col) {
        number_col=number_col_contestant
      }
      
      for (i in input$FeaturePlot_MetaInput) {
        if (all(as.numeric(reactivevalue$SeuratObject@meta.data[,i])==0)) {
          temp=FeaturePlot(reactivevalue$SeuratObject,features = i,reduction = input$MetaFeaturePlot_reduction,order = T,cols = c('grey','grey'))&theme(plot.title = element_text(size = 20),
                                                                                                                                                        axis.title = element_text(size = 18, face = "bold"),
                                                                                                                                                        axis.text = element_text(size = 17),
                                                                                                                                                        legend.title = element_text(size = 18),
                                                                                                                                                        legend.text = element_text(size = 17, face = "bold"))
          
        } else {
          temp=FeaturePlot(reactivevalue$SeuratObject,features = i,reduction = input$MetaFeaturePlot_reduction,order = T,cols = c('grey','blue'))&theme(plot.title = element_text(size = 20),
                                                                                                                                                        axis.title = element_text(size = 18, face = "bold"),
                                                                                                                                                        axis.text = element_text(size = 17),
                                                                                                                                                        legend.title = element_text(size = 18),
                                                                                                                                                        legend.text = element_text(size = 17, face = "bold"))
          
        }
        as_grob(temp)
        
        plots[[i]]=temp
      }
      
      plots=patchwork::wrap_plots(plots,ncol=number_col,nrow=number_col)
      output$MetaFeaturePlot <- renderPlot({
        print(plots) 
      }, width = reactive({
        
        if (num_genes == 1) { 
          500
        } else if (num_genes == 2) { 
          1000
        } else {
          1300
        }
      }), height = reactive({
        
        if (num_genes == 1) { 
          500
        } else if (num_genes == 2) { 
          500
        } else {
          1200
        }
      })) 
      
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
      waitress$start()
      
      num_genes <- length(input$VlnPlot_GeneInput)
      
      
      if (length(unique(input$VlnPlot_GeneInput)) > 2) {
        number_col = ceiling(sqrt(length(unique(input$VlnPlot_GeneInput))))
      } else {
        number_col = length(unique(input$VlnPlot_GeneInput))
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
    
    output$GeneVlnPlot = renderPlot({
      print(plot) 
    }, width = reactive({
      
      if (num_genes == 1) { 
        500
      } else if (num_genes == 2) { 
        1000
      } else {
        1300
      }
    }), height = reactive({
      if (num_genes == 1) { 
        500
      } else if (num_genes == 2) { 
        500
      } else {
        1200
      }
    })) 
    reactivevalue$vlnPlot = plot
    
    waitress$close()
    
  })
  
  observe(plotGeneVlnplot())
  
  
  plotMetaVlnplot=eventReactive(input$plotMetaVlnPlot_Button, {
    if (reactivevalue$Loaded) {
      waitress$start()
      
      num_genes <- length(input$VlnPlot_MetaInput)
      
      
      if (length(unique(input$VlnPlot_MetaInput)) > 2) {
        number_col = ceiling(sqrt(length(unique(input$VlnPlot_MetaInput))))
      } else {
        number_col = length(unique(input$VlnPlot_MetaInput))
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
    
    output$MetaVlnPlot = renderPlot({
      print(plot) 
    }, width = reactive({
      
      if (num_genes == 1) { 
        500
      } else if (num_genes == 2) { 
        1000
      } else {
        1300
      }
    }), height = reactive({
      
      if (num_genes == 1) { 
        500
      } else if (num_genes == 2) { 
        500
      } else {
        1200
      }
    })) 
    reactivevalue$vlnPlot = plot
    
    waitress$close()
    
  })
  
  observe(plotMetaVlnplot())
  
  shinyFileSave(input, "saveViolinPlot", roots =c(wd="/home/dnanexus/project/"), filetypes = c("pdf"))
  observeEvent(input$saveViolinPlot, {
    fileinfo <- parseSavePath(c(wd="/home/dnanexus/project/"), input$saveViolinPlot)
    if (nrow(fileinfo) > 0) {
      pdf(fileinfo$datapath) 
      print(reactivevalue$vlnPlot)
      dev.off()
    }
  })
  
  
  
  performDGE=eventReactive(input$DGEAction, {
    if (reactivevalue$Loaded) {
      waitress$start()
      temp_DGE=reactivevalue$SeuratObject
      spliter=reactivevalue$spliter
      DGE_List=reactivevalue$DGE_List
      Group1=list()
      for (i in input$DGE_Group_1) {
        Group1[[i]]=unlist(strsplit(DGE_List[[i]],split = spliter))
      }
      Group2=list()
      for (i in input$DGE_Group_2) {
        Group2[[i]]=unlist(strsplit(DGE_List[[i]],split = spliter))
      }
      Group1_cells=c()
      for (i in names(Group1)) {
        Group1_cells=c(Group1_cells,
                       rownames(temp_DGE@meta.data)[temp_DGE@meta.data[,Group1[[i]][1]]==Group1[[i]][2]]
                       )
      }
      
      Group2_cells=c()
      for (i in names(Group2)) {
        Group2_cells=c(Group2_cells,
                       rownames(temp_DGE@meta.data)[temp_DGE@meta.data[,Group2[[i]][1]]==Group2[[i]][2]]
        )
      }

      intersect_cells=intersect(Group1_cells,Group2_cells)
      Group1_cells=Group1_cells[!Group1_cells%in%intersect_cells]
      Group2_cells=Group2_cells[!Group2_cells%in%intersect_cells]
      temp_DGE$group='Others'
      for (i in 1:nrow(temp_DGE@meta.data)) {
        if (rownames(temp_DGE@meta.data)[i]%in%Group1_cells) {
          temp_DGE@meta.data$group[i]='Group1'
        }
        if (rownames(temp_DGE@meta.data)[i]%in%Group2_cells) {
          temp_DGE@meta.data$group[i]='Group2'
        }
      }
      if (length(Group1_cells)>2&length(Group2_cells)>2) {
        DGE_results=FindMarkers(temp_DGE,group.by='group',ident.1='Group1',ident.2='Group2')
        
      } else {
        showNotification('Too little cells to perform differential test.')
      }
      
      
      output$DGE_table=DT::renderDataTable(cbind(data.frame(row.names = rownames(DGE_results),genes=rownames(DGE_results)),
                                                               DGE_results),options = list(dom = 'Bfrtip'),
                                                         rownames= FALSE,
                                                         filter = list(position = "top"), selection='multiple',server = T)
      reactivevalue$DGE_table=cbind(data.frame(row.names = rownames(DGE_results),genes=rownames(DGE_results)),
                                    DGE_results)
      output$volcano=renderPlot(volcano.plot.function(input = reactivevalue$DGE_table,
                                           logfc.axis.name = 'avg_log2FC',
                                           pvalue.axis.name = 'p_val_adj',name.axis.name = 'genes'))
      
      output$downloadDGE <- downloadHandler(
        filename = function() {
          # Use the selected dataset as the suggested file name
          paste0(input$DGE_Group_1,'.vs.',input$DGE_Group_2, ".Differential.test.tsv")
        },
        content = function(file) {
          # Write the dataset to the file that will be downloaded
          write.table(reactivevalue$DGE_table, file,sep = '\t',quote=F,row.names = F)
        }
      )
      
      waitress$close()
    }

    
  })
  
  observe(performDGE())
  
  observeEvent(input$DGE_table_rows_selected, {
    print(input$DGE_table_rows_selected)
    output$volcano=renderPlot(volcano.plot.function(input = reactivevalue$DGE_table,
                                                    logfc.axis.name = 'avg_log2FC',
                                                    pvalue.axis.name = 'p_val_adj',name.axis.name = 'genes',genes.to.label = reactivevalue$DGE_table[input$DGE_table_rows_selected,'genes']))
  })
  
  
  

  