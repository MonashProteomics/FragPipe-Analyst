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
    libgdal-dev

### RUN Rscript -e 'install.packages(c("devtools", "BiocManager", "tidyverse", "ggrepel", "httr", "rjson", "mvtnorm", "tmvtnorm", "impute", \
### "pcaMethods", "imputeLCMD", "plotly", "DT", "testthat", "RColorBrewer", "shiny","shinyalert","shinydashboard", \
### "shinyjs", "svglite", "rhandsontable", "shinyBS", "shinyWidgets", "ggVennDiagram", "conflicted", "png"), dependencies=TRUE)'
### #FROM bioconductor/bioconductor_docker:RELEASE_3_15
### RUN Rscript -e 'BiocManager::install(pkgs=c("BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", \
### "GenomicFeatures", "AnnotationDbi", "Biobase", "ensembldb", "EnsDb.Hsapiens.v86", "DEP", "SummarizedExperiment", "limma", "ComplexHeatmap", ask=F), type = "source"))'

RUN Rscript -e 'install.packages("renv")'
COPY ./ /srv/shiny-server/fragpipe-analyst
WORKDIR /srv/shiny-server/fragpipe-analyst
<<<<<<< HEAD
### RUN Rscript -e 'renv::init()' # This is already run and pushed to github.
RUN Rscript -e 'renv::restore()'
RUN chmod -R +r /srv/shiny-server/fragpipe-analyst
=======
### RUN Rscript -e 'renv::init()'
RUN Rscript -e 'renv::restore()'
RUN chmod -R +r /srv/shiny-server/fragpipe-analyst
>>>>>>> 94d3c5fbbcd9d9ca73a7d2e17eb4d06c262d3794
