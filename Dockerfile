FROM rocker/shiny-verse:4.2.1
RUN apt-get update && apt-get install -yq \
     libhdf5-dev \
     libnetcdf-dev \
     build-essential \
     libgd-dev \
     libbz2-dev \
     libudunits2-dev \
     libproj-dev \
     libgdal-dev

RUN Rscript -e 'install.packages(c("devtools", "mvtnorm", "tmvtnorm","impute", "pcaMethods", "imputeLCMD", "plotly", \
    "DT", "BiocManager","testthat", "RColorBrewer", "shiny","shinyalert","shinydashboard", "shinyjs", "svglite", \
    "rhandsontable", "shinyBS", "shinyWidgets", "ggVennDiagram", "conflicted"), dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install(pkgs=c("DEP", "SummarizedExperiment", "limma", "ComplexHeatmap","ensembldb","EnsDb.Hsapiens.v86", ask=F))'
COPY . /srv/shiny-server/fragpipe-analyst
RUN rm -f /srv/shiny-server/fragpipe-analyst/.Rprofile
RUN chmod -R +r /srv/shiny-server/fragpipe-analyst