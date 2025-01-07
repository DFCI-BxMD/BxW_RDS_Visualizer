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


#########           UI START          ##########

# fully defined ui
dashboardPage(
      # webapp layout #
      dashboardHeader(title = "RDS Visualizer v1.1", 
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
      menuItem("Dim Plot", tabName = "DimPlot", icon = icon("cubes")),
      menuItem("Feature Plot", tabName = "Feature_plot_page", icon = icon("chart-bar")),
      menuItem("Violin Plot", tabName = "VlnPlot", icon = icon("record-vinyl")),
      menuItem("Help Page", tabName = "help", icon = icon("question"))
    )
  ),
  dashboardBody(
    useWaitress(color = "#0047AB"),
    tabItems(
      tab_home <- tabItem(tabName = "home",
                          h2("Select a File"),
                          shinyFilesButton("files", label="Browse", title="Please select a file", multiple=FALSE),
                          verbatimTextOutput("loadedFile"),
                          h2("Home Page"),
                          br(),
                          br(),
                          p("Welcome to the RDS Visualizer Web App (v1.1)!"),
                          br(),
                          p("This web app is designed for the visualization and exploration of single-cell RNA-seq data contained in Seurat objects. 
                             It will provide various plots and features to help you analyze and gain insights from your data. 
                             This web app has been optimized/structured a DNAnexus build, which circumvents the instance size restrictions of shinyapps.io's base plan."),
                          br(),
                          p("This web app currently consists of four pages of note: a Data Summary page, a Feature Plot page, a Dim Plot page, and a Violin Plot page. 
                            As a general rule of thumb, these pages will not visualize any data if no .rds file has been provided by the user in the data summary page.
                            Note also that the drop downs on the Dim Plot and Violin Plot pages will not display updated categories until a file has been provided by the user."),
                          br(),
                          #embed href in paragraph sentence
                          p("Click ", a("here", target="_blank", href="https://raw.githubusercontent.com/rkafrawi/RDS_Vis_v1_1/main/docs/RDS_Visualizer_user_guide.pdf", .noWS = "outside"), " for an in depth user guide.", .noWS = c("after-begin", "before-end")),
                          
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
                          imageOutput('MainFigure')
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

                                   br(),
                                   p("Enter the gene of interest and color the selected dimension reduction with its expression"),
                                   br(),
                                   selectInput("FeaturePlot_GeneInput", "Select Gene:",choices = NULL,selected = NULL,multiple = T),
                                   selectizeInput("FeaturePlot_reduction", "Select Reduction Name:",choices = NULL,selected = NULL),
                                   actionButton('plotFeaturePlot_Button','Plot FeaturePlot'),
                                   br(),
                                   plotOutput('FeaturePlot'),
                                   br(),
                                   shinySaveButton("saveFeaturePlot", "Download Feature Plot (PDF)", title = "Save Feature Plot", filetype = list(PDF = "pdf"))

                                   
      ),
      tab_dim <- tabItem(tabName = "DimPlot",
                         h2("Dim Plot"),
                         br(),
                         p("This page generates DimPlots that can be split by different groups 
                           found in a dropdown that updates based on the metadata columns in 
                           the uploaded Seurat object. By default, the Dimplot will not be split by a group. 
                           To split the visualization by a group, toggle the checkbox and select a group from the dropdown."),
                         br(),
                         #checkboxInput("splitToggle", "Split Dim Plot", value = FALSE),
                         selectizeInput('DimPlot_group_by','Color the DimPlots by: ',choices=NULL),
                         br(),
                        # selectInput("variableInput", label = "Select a Group:",
                        #             choices = "Input File For Dropdown Options"),
                        selectizeInput('DimPlot_split_by','Split the DimPlots by: ',choices=NULL),
                        selectizeInput('DimPlot_reduction','Plot the DimPlots by : ',choices=NULL),
                        
                        br(),
                        actionButton("plotDimPlot_Button","Generate DimPlot"),
                        br(),
                        plotOutput("DimPlot"),
                        br(),
                        shinySaveButton("saveDimPlot", "Download Dim Plot (PDF)", title = "Save Dim Plot", filetype = list(PDF = "pdf"))
                         
      ),
      tab_violin <- tabItem(tabName = "VlnPlot",
                            h2("Violin Plot"),
                            br(),
                            selectInput("VlnPlot_GeneInput", "Select Gene:",choices = NULL,selected = NULL,multiple = T),
                            selectizeInput('VlnPlot_group_by','Group the Violin Plot by: ',choices=NULL),
                            actionButton("plotVlnPlot_Button","Generate Violin Plot"),
                            plotOutput("VlnPlot"),
                            br(),
                            shinySaveButton("saveViolinPlot", "Download Violin Plot (PDF)", title = "Save Violin Plot", filetype = list(PDF = "pdf"))
      ),
      
      tab_help <- tabItem(tabName = "help",
                          h2("Help Page"),
                          br(),
                          p("Below are a list of Frequently Asked Questions (FAQs)."),
                          br(),
                          accordion(
                            id = "helpAccordion",
                            accordionItem(
                              title = "Why is my input file not loading?",
                              "Check if the file path is correct and if the file format is compatible with the application. This web app currently only supports .rds files with an upper size limitation of 10 GB."
                            ),
                            accordionItem(
                              title = "How do I download these figures?",
                              "To download the figures, you can right-click on each figure and choose the 'Save Image As' option. Then, select your desired location on your computer to save the image."
                            ),
                            accordionItem(
                              title = "Is there any way to view tabular representations of the visualizations?",
                              "Currently, the server does not support table generation to supplement each visual. The only way to inspect tabular representations of the uploaded data is in the Data Summary page."
                            ),
                            accordionItem(
                              title = "Why are only some of my categorical variables showing in my dropdown options?",
                              "For the sake of visual parity, only categorical variables containing less than 5 levels are included from the uploaded .rds file."
                            )
                          ),
                          br(),
                          br(),
                          br(),
                          p("Click ", a("here", target="_blank", href="https://satijalab.org/seurat/articles/pbmc3k_tutorial.html", .noWS = "outside"), " to reference the clustering guide used to build this web app.", .noWS = c("after-begin", "before-end")),
                          br(),
                          p("Sample seurat objects can be found ", a("here", target="_blank",  href="https://www.dropbox.com/home/Rizky/RShiny/rds/BMS_DGKi", .noWS = "outside"),  ". Feel free to experiment with these example datasets to familiarize yourself with the web app workflow!", .noWS = c("after-begin", "before-end")),
                          
      )
    ), 
    div(class = "hr-container",
        br(),
        hr(class = "hr-line")),
    # Footer
    div(class = "foot",
        tags$footer(
          class = "footer",
          p("Developed for use of the Belfer Center Research Team by the BxMD Group at DFCI."),
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
