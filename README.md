[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![R](https://img.shields.io/badge/R-%3E4.4-brightgreen)

# FragPipe-Analyst

A tool for analyzing quantitative proteomics datasets for [FragPipe](https://fragpipe.nesvilab.org/).

## Features

- Differential expression analysis
- Enrichment analysis (GO/Pathways)
- Imputation (optional)
- Data visualization
  1. PCA
  2. Sample correlation
  3. Heatmaps
  4. Missing value inspection
  5. Sample coverage
  6. Protein intensity plots for slected protein(s)
  7. Imputation effect evaluation

- Report and multiple levels of exportable tables for further analysis
  - Table options
    - DE results
    - Unimputed data matrix: Original protein intensities before imputation
    - Imputed data matrix: Protein intensities after performing selected imputation method

## Public servers

There are two server instances
- Production server is hosted at [https://fragpipe-analyst.org/](https://fragpipe-analyst.org/).
- Dev server is also hosted at [http://fragpipe-analyst.nesvilab.org/](http://fragpipe-analyst.nesvilab.org/).

## Documentation

Now documentation is all moved to [here](https://fragpipe-analyst-doc.nesvilab.org/).

## Install on your own machine

Two options are available right now for FragPipe-Analyst: local installation or run through Docker.

### Local Installation

#### Prerequisite
- R >= 4.4
- PDFlatex

Once all the prerequisites are installed, follow steps below to build and run the server locally.

``` sh
# Clone the repository
git clone https://github.com/MonashProteomics/FragPipe-Analys.git

# Move to the folder
cd FragPipe-Analyst

# Inside R console or R studio
> install.packages("renv")
> renv::init(bioconductor = T)

# Install shiny.info, it's removed from CRAN since 2025-03-21 https://cran.r-project.org/web/packages/shiny.info/index.html, so we need to install their github version
> renv::install("Appsilon/shiny.info")

# Execute
> shiny::runApp()
```

### Installation through Docker:

``` sh
# Clone the repository
git clone https://github.com/MonashProteomics/FragPipe-Analyst.git

# Move to the folder
cd FragPipe-Analyst

# Build FragPipe-Analyst (Any name after -t)
docker buildx build -f Dockerfile.local -t fragpipe-analyst  --output=type=docker --platform=linux/amd64 .

# Run FragPipe-Analyst
docker run -it --platform=linux/amd64 -d -p 3838:3838 fragpipe-analyst

# Open local interface
http://localhost:3838/fragpipe-analyst/
```
