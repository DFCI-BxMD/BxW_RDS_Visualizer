#load tab source files


ui <- dashboardPage(
  dashboardHeader(title = "[Tool x]"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Feature Plot", tabName = "feature", icon = icon("chart-bar")),
      menuItem("Volcano Plot", tabName = "volcano", icon = icon("fire")),
      menuItem("Violin Plot", tabName = "violin", icon = icon("record-vinyl")),
      menuItem("Dim Plot", tabName = "dim", icon = icon("cubes")),
      menuItem("Help Page", tabName = "help", icon = icon("question"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "home",
              h2("Home Page"),
              br(),
              br(),
              p("Welcome to the 'x' Tool!"),
              br(),
              p("This web app is designed for the visualization and exploration of
                single-cell RNA-seq data contained in Seurat objects. It provides various interactive plots
                and features to help you analyze and gain insights from your data."),
              br(),
              p("Developed for the Belfer Center Research Team by the BxMD Group at DFCI."),
              br(),
              p("To get started, upload your rds object and then proceed to one of the visualization tabs on the lefthand side of the web app."),
              #file upload
              fileInput("RDS_element", "Upload RDS File")
              # Add any additional content or UI elements here
      ),
      tabItem(tabName = "feature",
              h2("Feature Plot Page"),
              br(),      
              p("Enter gene names:"),
              textInput("geneInput", "Gene Names"),
              br()
              
              # Add feature plot page content here
      ),
      tabItem(tabName = "volcano",
              h2("Volcano Plot Page"),
              br(),
              p("Enter gene names:"),
              textInput("geneInput", "Gene Names"),
              br(),
              p("Select a variable:"),
              selectInput("variableInput", label = NULL,
                          choices = c("Clusters", "Treatments", "Responder States")),
              br()
              # Add volcano plot page content here
      ),
      tabItem(tabName = "violin",
              h2("Violin Plot Page"),
              # Add violin plot page content here
              br(),
              p("Select a variable:"),
              selectInput("variableInput", label = NULL,
                          choices = c("Clusters", "Treatments", "Responder States"))
      ),
      tabItem(tabName = "dim",
              h2("Dim Plot"),
              br(),
              p("Select a variable:"),
              selectInput("variableInput", label = NULL,
                          choices = c("Clusters", "Treatments", "Responder States"))
              # Add any additional content or UI elements here
      ),
      tabItem(tabName = "help",
              h2("Help Page"),
              br(),
              accordion(
                id = "helpAccordion",
                accordionItem(
                  title = "Why is my input file not loading?",
                  "The input file may not be in the correct format or location. Make sure the file exists and is accessible. Check if the file path is correct and if the file format is compatible with the application."
                ),
                accordionItem(
                  title = "How do I download these figures?",
                  "To download the figures, you can right-click on each figure and choose the 'Save Image As' option. Select a location on your computer to save the image."
                ),
                accordionItem(
                  title = "Are there tables corresponding to the plots generated in each tab?",
                  "Yes, there are corresponding tables for the plots generated in each tab. You can find these tables in the 'Data' section of each tab. The tables contain the underlying data used to generate the plots."
                )
              )
      )
    )
  )
)
