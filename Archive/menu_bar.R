menu_bar <- dashboardSidebar(
                sidebarMenu(
                  menuItem("Home", tabName = "home", icon = icon("home")),
                  menuItem("Feature Plot", tabName = "feature", icon = icon("chart-bar")),
                  menuItem("Volcano Plot", tabName = "volcano", icon = icon("fire")),
                  menuItem("Violin Plot", tabName = "violin", icon = icon("record-vinyl")),
                  menuItem("Dim Plot", tabName = "dim", icon = icon("cubes")),
                  menuItem("Help Page", tabName = "help", icon = icon("question"))
                )
              )
