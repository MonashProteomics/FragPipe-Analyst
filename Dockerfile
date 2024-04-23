# syntax=docker/dockerfile:1
FROM rocker/shiny-verse:4.2.1
RUN apt-get update && apt-get install -yq \
    libbz2-dev \
    libhdf5-dev \
    libnetcdf-dev \
    build-essential \
    libgd-dev \
    libudunits2-dev \
    libproj-dev \
    libgdal-dev  \
    texlive-latex-base \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-extra

RUN Rscript -e 'install.packages(c("devtools", "BiocManager", "tidyverse", "ggrepel", "httr", "rjson", "mvtnorm", "tmvtnorm", \
"imputeLCMD", "plotly", "DT", "testthat", "RColorBrewer", "shiny","shinyalert","shinydashboard", \
"shinyjs", "svglite", "rhandsontable", "shinyBS", "shinyWidgets", "ggVennDiagram", "conflicted", "png", "vegan", \
"shinycssloaders","shiny.info","UpSetR"), dependencies=TRUE)'
### #FROM bioconductor/bioconductor_docker:RELEASE_3_15
RUN Rscript -e 'BiocManager::install(pkgs=c("ensembldb", "EnsDb.Hsapiens.v86", "DEP", "SummarizedExperiment", "limma", "ComplexHeatmap", "impute", "pcaMethods"), ask=F, dependencies=TRUE)'

### RUN Rscript -e 'install.packages("renv")'
COPY ./ /srv/shiny-server/fragpipe-analyst
COPY shiny-server.conf.prod /etc/shiny-server/shiny-server.conf
### WORKDIR /srv/shiny-server/fragpipe-analyst
### RUN Rscript -e 'renv::init()' # This is already run and pushed to github.
### RUN Rscript -e 'renv::restore()'
### RUN Rscript -e 'renv::isolate()'
RUN chmod -R +r /srv/shiny-server/fragpipe-analyst
