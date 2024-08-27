# base image: rocker/verse (with a specific version of R)
#   has R, RStudio, tidyverse, devtools, tex, and publishing-related packages
FROM rocker/shiny

MAINTAINER Rizky Kafrawi <rkafrawi@dfci.harvard.edu>

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libgsl-dev \
    libxml2-dev \
    libglpk-dev \
    libcurl4-openssl-dev

RUN R -e "install.packages(c('shiny','shinydashboard','shinydashboardPlus', 'purrr', 'stringr', 'DT', 'shinyalert', 'ggplot2'))" 

RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"
RUN R --no-echo --no-restore --no-save -e "install.packages('Seurat')"

# Install SeuratDisk
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('mojaveazure/seurat-disk')"
