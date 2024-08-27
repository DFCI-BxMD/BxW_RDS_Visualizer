#increase filesize limit
options(shiny.maxRequestSize=10000*1024^2) #max 10gb upload

#load libraries
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(purrr)
library(stringr)
library(Seurat)
library(SeuratObject)
library(DT)
library(shinyalert)
library(grDevices)
library(ggplot2)

#load ui/server source
source("myUI_mod.R")
source("myServer.R")

#launch shinyapp
shinyApp(ui = ui, 
         server = server
         )

# 
# usethis::create_from_github(
#   "https://github.com/rkafrawi/RDS_Vis_v1_1.git",
#   destdir = "/Users/rkafrawi/Desktop/RDS_Vis_v1_1",
#   fork=FALSE
# )

# usethis::create_from_github(
#   "https://github.com/rkafrawi/RDS_Vis_v1_1.git",
#   destdir = "/Users/rkafrawi/Downloads",
#   fork=FALSE
# )
#optimizing app:
# restructure ui to pull for tabs instead of hosting all in one place
# 
# restructure server logic using modules since a lot of the code is reused
#   plot module
#   gene module
#   add a sample gene textbox for violinplot

# git remote add origin git@github.com:rkafrawi/RDS_Vis_v1_1.git
# git branch -M main
# git push -u origin main
