#increase filesize limit
#options(shiny.maxRequestSize=1000000*1024^2) #max 10gb upload

#load libraries
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyFiles)
library(shinycssloaders)
library(purrr)
library(stringr)
library(Seurat)
library(DT)
library(shinyalert)
library(ggplot2)
library(waiter)
library(cowplot)
library(patchwork)
library(qs) 
library(future)


#########           UI START          ##########

# fully defined ui
dashboardPage(
  # webapp layout #
  dashboardHeader(title = span("BxW Single-Cell Analyzer v2.1", style = "font-size: 14px;"), 
                  tags$li(
                    class = "dropdown", 
                    style = "padding-top: 0px;",
                    tags$a(target="_blank", href = "https://github.com/rkafrawi/RDS_Vis_v1_1/tree/main", "Github Repo", .noWS = "outside")
                  )
  ),
  menu_bar <- dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Data Summary", tabName="input_summ", icon=icon("table")),
      menuItem("Dimensionality Reduction Plot", tabName = "DimPlot", icon = icon("cubes")),
      menuItem("Feature Plot", tabName = "Feature_plot_page", icon = icon("chart-bar")),
      menuItem("Violin Plot", tabName = "VlnPlot", icon = icon("record-vinyl")),
      menuItem("Help Page", tabName = "help", icon = icon("question"))
    )
  ),
  dashboardBody(
    tags$head(  
      tags$style(HTML("
        #GeneFeaturePlot-container {
          width: 1300px;
          height: 1200px;
          overflow-x: auto;
          overflow-y: auto;
        }
        #GeneFeaturePlot {
          width: auto;
          height: auto;
        }
        #MetaFeaturePlot-container {
          width: 1300px;
          height: 1200px;
          overflow-x: auto;
          overflow-y: auto;
        }
        #MetaFeaturePlot {
          width: auto;
          height: auto;
        }
        #GeneVlnPlot-container {
           width: 1300px;
          height: 1200px;
          overflow-x: auto;
          overflow-y: auto;
        }
        #GeneVlnPlot {
          width: auto;
          height: auto;
        }
        #MetaVlnPlot-container {
          width: 1300px;
          height: 1200px;
          overflow-x: auto;
          overflow-y: auto;
        }
        #MetaVlnPlot {
          width: auto;
          height: auto;
        }
      "))
    ),
    useWaitress(color = "#0047AB"),
    
    tabItems(
      tab_home <- tabItem(tabName = "home",
                          h2("Select a File"),
                          shinyFilesButton("files", label="Browse", title="Please select a file", multiple=FALSE),
                          verbatimTextOutput("loadedFile"),
                          h2("Home Page"),
                          br(),
                          br(),
                          p("Welcome to the BxW Single-Cell Analyzer (v2.1)!"),
                          br(),
                          p("This web app is designed for the visualization and exploration of single-cell RNA-seq data contained in Seurat objects. 
                             It will provide various plots and features to help you analyze and gain insights from your data. 
                             This web app has been optimized/structured a DNAnexus build, which circumvents the instance size restrictions of shinyapps.io's base plan."),
                          br(),
                          p("This web app currently consists of four pages of note: a Data Summary page, a Feature Plot page, a Dimensionality Reduction Plot page, and a Violin Plot page. 
                            As a general rule of thumb, these pages will not visualize any data if no .rds file has been provided by the user in the Home page.
                            Note also that the drop downs on the Dim Plot and Violin Plot pages will not display updated categories until a file has been provided by the user."),
                          br(),
                          
                          # Add any additional content or UI elements here
      ),
      tab_input <- tabItem(
        tabName = "input_summ",
        h2("Data Summary Page"),
        br(),
        # Update the text to remove file upload instructions
        p("This app generates visualizations for a pre-loaded .rds file."),
        br(),
        # Display panels with data output
        tabsetPanel(
          tabPanel(
            "Main Figure",
            br(),
            imageOutput('MainFigure', height="700px")
          ),
          tabPanel(
            "Metadata",
            br(),
            dataTableOutput("Metadata")
          ),
          tabPanel(
            "Data Distribution",
            br(),
            selectizeInput('Bar_Graph_y','Group in Y-Axis',choices=NULL,selected=NULL),
            selectizeInput('Bar_Graph_fill','Groups that will fill the colors of bar',choices=NULL,selected=NULL),
            br(),
            checkboxInput('Bar_Graph_Percentage','Generate Percentage Bar Graph'),
            br(),
            plotOutput('BarGraph'),
            br(),
            dataTableOutput("BarGraph_Table")
            
          )
        )
      ),
      tab_inputfeatures <- tabItem(tabName = "Feature_plot_page",
                                   h2("Feature Plot"),
                                   tabsetPanel(
                                     tabPanel(
                                       "Gene Expression Feature Plot",
                                       
                                       br(),
                                       p("Enter the gene of interest and color the selected dimension reduction with its expression"),
                                       br(),
                                       selectizeInput("FeaturePlot_GeneInput", "Select Gene:",choices = NULL,selected = NULL,multiple = T,options = list(maxItems = 10)),
                                       #selectizeInput("FeaturePlot_MetaInput", "Select Numerical Metadata Column:",choices = NULL,selected = NULL,multiple = T,options = list(maxItems = 10)),
                                       
                                       selectizeInput("GeneFeaturePlot_reduction", "Select Reduction Name:",choices = NULL,selected = NULL),
                                       actionButton('plotGeneFeaturePlot_Button','Plot FeaturePlot'),
                                       br(),
                                       div(id = "GeneFeaturePlot-container",
                                           plotOutput("GeneFeaturePlot")
                                       ),
                                       br(),
                                       shinySaveButton("saveFeaturePlot", "Download Plot (PDF)", title = "Save Feature Plot", filetype = list(PDF = "pdf"))
                                     ),
                                     
                                     tabPanel(
                                       "Metadata Feature Plot",
                                       
                                       br(),
                                       p("Enter the numerical metadata of interest and color the selected dimension reduction with its values"),
                                       br(),
                                       #selectizeInput("FeaturePlot_GeneInput", "Select Gene:",choices = NULL,selected = NULL,multiple = T,options = list(maxItems = 10)),
                                       selectizeInput("FeaturePlot_MetaInput", "Select Numerical Metadata Column:",choices = NULL,selected = NULL,multiple = T,options = list(maxItems = 10)),
                                       
                                       selectizeInput("MetaFeaturePlot_reduction", "Select Reduction Name:",choices = NULL,selected = NULL),
                                       actionButton('plotMetaFeaturePlot_Button','Plot FeaturePlot'),
                                       br(),
                                       div(id = "MetaFeaturePlot-container",
                                           plotOutput("MetaFeaturePlot")
                                       ),
                                       br(),
                                       shinySaveButton("saveFeaturePlot", "Download Plot (PDF)", title = "Save Feature Plot", filetype = list(PDF = "pdf"))
                                     ))
                                   
                                   
      ),
      tab_dim <- tabItem(tabName = "DimPlot",
                         h2("Dimensionality Reduction Plot"),
                         br(),
                         p("This page generates Dimensionality Reduction Plot that can be split by different groups 
                           found in a dropdown that updates based on the metadata columns in 
                           the uploaded Seurat object. By default, the Dimplot will not be split by a group. 
                           To split the visualization by a group, toggle the checkbox and select a group from the dropdown."),
                         br(),
                         #checkboxInput("splitToggle", "Split Dim Plot", value = FALSE),
                         selectizeInput('DimPlot_group_by','Color the plots by: ',choices=NULL),
                         br(),
                         # selectInput("variableInput", label = "Select a Group:",
                         #             choices = "Input File For Dropdown Options"),
                         selectizeInput('DimPlot_split_by','Split the plots by: ',choices=NULL),
                         selectizeInput('DimPlot_reduction','Plot the plots by : ',choices=NULL),
                         
                         br(),
                         actionButton("plotDimPlot_Button","Generate Plot"),
                         br(),
                         plotOutput("DimPlot", height="700px"),
                         br(),
                         shinySaveButton("saveDimPlot", "Download Plot (PDF)", title = "Save Dim Plot", filetype = list(PDF = "pdf"))
                         
      ),
      tab_violin <- tabItem(tabName = "VlnPlot",
                            h2("Violin Plot"),
                            tabsetPanel(
                              tabPanel(
                                "Gene Expression Violin Plot",
                                br(),
                                p("Generate Violin Plots based on selected genes and groups"),
                                selectizeInput("VlnPlot_GeneInput", "Select Gene:",choices = NULL,selected = NULL,multiple = T,options = list(maxItems = 8)),
                                selectizeInput('GeneVlnPlot_group_by','Group the Violin Plot by: ',choices=NULL),
                                actionButton("plotGeneVlnPlot_Button","Generate Violin Plot"),
                                br(),
                                div(id = "GeneVlnPlot-container",
                                    plotOutput("GeneVlnPlot")
                                ),
                                br(),
                                shinySaveButton("saveViolinPlot", "Download Violin Plot (PDF)", title = "Save Violin Plot", filetype = list(PDF = "pdf"))
                              ),
                              tabPanel(
                                "Metadata Violin Plot",
                                br(),
                                p("Generate Violin Plots based on selected metadata values and groups"),
                                selectizeInput("VlnPlot_MetaInput", "Select Numerical Metadata Column:",choices = NULL,selected = NULL,multiple = T,options = list(maxItems = 8)),
                                selectizeInput('MetaVlnPlot_group_by','Group the Violin Plot by: ',choices=NULL),
                                actionButton("plotMetaVlnPlot_Button","Generate Violin Plot"),
                                br(),
                                div(id = "MetaVlnPlot-container",
                                    plotOutput("MetaVlnPlot")
                                ),
                                br(),
                                shinySaveButton("saveViolinPlot", "Download Violin Plot (PDF)", title = "Save Violin Plot", filetype = list(PDF = "pdf"))
                              ),
                            )
                            
      ),
      
      tab_help <- tabItem(tabName = "help",
                          h2("Help Page"),
                          br(),
                          p("Below are a list of Frequently Asked Questions (FAQs)."),
                          br(),
                          accordion(
                            id = "helpAccordion",
                            
                            accordionItem(
                              title = "How do I download these figures?",
                              "To download the figures, you can click on the 'Download PDF' button and select/create a folder within your project and click save. You can also right-click on each figure and choose the 'Save Image As' option. Then, select your desired location on your computer to save the image."
                            ),
                            accordionItem(
                              title = "Is there any way to view tabular representations of the visualizations?",
                              "Currently, the server does not support table generation to supplement each visual. The only way to inspect tabular representations of the uploaded data is in the Data Summary page."
                            )
                          ),
                          br(),
                          br(),
                          br(),
                          p("Click ", a("here", target="_blank", href="https://satijalab.org/seurat/articles/pbmc3k_tutorial.html", .noWS = "outside"), " to reference the clustering guide used to build this web app.", .noWS = c("after-begin", "before-end")),
                          
                          
      )
    ), 
    div(class = "hr-container",
        br(),
        hr(class = "hr-line")),
    # Footer
    div(class = "foot",
        tags$footer(
          class = "footer",
          p("Developed by the BxMD Group at DFCI"),
          br(),
        ),
        tags$img(
          src ="https://raw.githubusercontent.com/tpathakdfci/Logo/main/Home%20_%20Informatics%20and%20Analytics.png",
          style="display: block; margin-left: auto; margin-right: auto; width:20%; height:15%",
        )
    ),
    
    
    
    # CSS Styling
    tags$style(HTML("
      .hr-container {
        text-align: center;
      }
      
      .hr-line {
        border-top: 1px solid #ccc;
        margin: 20px 0;
      }
      .img {
    
      }
      .footer {
        padding: 10px;
        text-align: center;
      }
      .placeholder {
      text-align: center;
      }
      .full-width-cols {
      overflow:auto;
      }
    "))
  )
)

#########           UI END          ##########

