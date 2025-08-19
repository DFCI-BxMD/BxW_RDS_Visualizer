# base image: rocker/verse (with a specific version of R)
#   has R, RStudio, tidyverse, devtools, tex, and publishing-related packages
FROM rocker/shiny


# Install system dependencies
RUN apt-get clean && \
     apt-get update && \
    apt-get install -y --no-install-recommends \
    libgsl-dev \
    libxml2-dev \
    libglpk-dev \
    libcurl4-openssl-dev \
    cmake \
    git

RUN R -e "install.packages(c('shiny','shinydashboard','shinydashboardPlus', 'purrr', 'stringr', 'DT', 'shinyalert', 'ggplot2'))" 

RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"
RUN R --no-echo --no-restore --no-save -e "install.packages('Seurat')"

# Install SeuratDisk
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('mojaveazure/seurat-disk')"

RUN R -e "install.packages('shinyFiles')"
RUN R -e "install.packages('waiter')"
RUN R -e "install.packages('shinycssloaders')"

RUN R -e "install.packages('tiledb')"
RUN R -e "Sys.setenv('NOT_CRAN' = 'true')"
RUN R -e "install.packages('arrow')"
RUN R -e "install.packages('tiledbsoma', repos = c('https://tiledb-inc.r-universe.dev','https://cloud.r-project.org'))"
RUN R -e "install.packages('ggrastr')"
RUN R -e "install.packages('ggplotify')"
RUN R -e "install.packages('qs')"
RUN R -e "install.packages('future')"

RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('immunogenomics/presto')"

COPY ./Webapp_source/ /srv/shiny-server/

CMD ["/usr/bin/shiny-server"]
