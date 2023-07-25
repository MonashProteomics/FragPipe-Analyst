# Define UI for data upload app ----
VERSION <- "v0.25"

ENABLE_PEPTIDE_ANALYSIS <- T
if (ENABLE_PEPTIDE_ANALYSIS) {
  analysis_options <- c("LFQ"="LFQ", "TMT"="TMT", "DIA"="DIA", "TMT (peptide)"="TMT-peptide")
} else {
  analysis_options <- c("LFQ"="LFQ", "TMT"="TMT", "DIA"="DIA")
}

ui <- function(request){shinyUI(
  dashboardPage(
    skin = "blue",
    dashboardHeader(title = "FragPipe-Analyst"),
    # disable = TRUE),# Disable title bar
    dashboardSidebar(
      sidebarMenu(
        id="tabs_selected",
        convertMenuItem(
          tabName = "home",
          menuItem('Home', icon=icon("home"), selected = TRUE, tabName = "home")
        ),
        convertMenuItem(
          tabName = 'analysis',
          menuItem("Analysis", tabName="analysis", icon=icon("flask"),
                   selectInput("exp", "Experiment type:", analysis_options, selected = "LFQ"),
                   conditionalPanel(
                     condition = "input.exp == 'LFQ'",
                     fileInput('lfq_expr', 'Upload FragPipe combined_protein.tsv',
                               accept=c('text/tsv',
                                        'text/tab-separated-values,text/plain',
                                        '.tsv')),
                     fileInput('lfq_manifest', 'Upload sample annotation',
                               accept=c('text/tsv',
                                        'text/tab-separated-values,text/plain',
                                        '.tsv')),
                     radioButtons("lfq_type",
                                  "Intensity Type",
                                  choices = c("Intensity"="Intensity",
                                              "MaxLFQ Intensity"="MaxLFQ",
                                              "Spectral Count"="Spectral Count"),
                                              selected = "Intensity"),
                     tags$hr(),
                     downloadLink("lfq_example", label="Example LFQ data"),
                     br(),
                     downloadLink("lfq_annotation", label="Example annotation")),
                   conditionalPanel(
                     condition = "input.exp == 'TMT'",
                     fileInput('tmt_expr', 'Upload gene-level TMT-I report *.tsv',
                               accept=c('text/tsv',
                                        'text/tab-separated-values,text/plain',
                                        '.tsv')),
                     fileInput('tmt_annot', 'Upload sample annotation',
                               accept=c('text/tsv',
                                        'text/tab-separated-values,text/plain',
                                        '.tsv'))
                     ),
                   conditionalPanel(
                     condition = "input.exp == 'DIA'",
                     fileInput('dia_expr', 'Upload protein group (PG) matrix *.tsv',
                               accept=c('text/tsv',
                                        'text/tab-separated-values,text/plain',
                                        '.tsv')),
                     fileInput('dia_manifest', 'Upload sample annotation',
                               accept=c('text/tsv',
                                        'text/tab-separated-values,text/plain',
                                        '.tsv'))),
                   conditionalPanel(
                     condition = "input.exp == 'TMT-peptide'",
                     fileInput('tmt_pept_expr', 'Upload peptide-level TMT-I report *.tsv',
                               accept=c('text/tsv',
                                        'text/tab-separated-values,text/plain',
                                        '.tsv')),
                     fileInput('tmt_pept_annot', 'Upload sample annotation',
                               accept=c('text/tsv',
                                        'text/tab-separated-values,text/plain',
                                        '.tsv'))
                   ),
                 tags$hr(),
                 menuItem("Advanced Options",tabName="advanced", icon = icon("cogs"),
                          numericInput("min_global_appearance",
                                       "Min percentage of non-missing values globally",
                                       min = 0, max = 100, value = 0),
                          numericInput("min_appearance_each_condition",
                                       "Min percentage of non-missing values in at least one condition",
                                       min = 0, max = 100, value = 0),
                          numericInput("p", 
                                       "Adjusted p-value cutoff",
                                       min = 0.0001, max = 0.1, value = 0.01),
                          numericInput("lfc",
                                       "Log2 fold change cutoff",
                                       min = 0, max = 10, value = 1),
                          # checkboxInput("paired",
                          #               "Paired test", FALSE),
                          radioButtons("normalization",
                                       "Normalization type",
                                       choices = c("No normalization"="none",
                                                   "Variance stabilizing normalization (LFQ and DIA only)"="vsn"
                                                   #, "Median centered"="MD"
                                                   ), selected = "none"),
                          radioButtons("imputation",
                                       "Imputation type",
                                       choices = c("No imputation"="none", "Perseus-type"="man",
                                                   "MLE"="MLE", "knn"="knn", "bpca"="bpca"),
                                       selected = "man"),
                          radioButtons("fdr_correction",
                                       "Type of FDR correction",
                                       choices =  c("Benjamini Hochberg"="BH",
                                                    "Local and tail area-based"="fdrtool"
                                       ), selected= "BH")
                          # checkboxInput("s",
                          #               "Paired test", FALSE),
                          # numericInput("k_number",
                          #              "Number of clusters in heatmap",
                          #              min = 1, max = 10, value = 3)
                 ),
               tags$hr(),
               actionButton("analyze", "Run"),
               #actionButton("load_data", "Load example data")
               tags$script(HTML("
                 $(document).ready(function() {
                    $('#analyze').on('click', function(){$(this).blur()});
                  });
                "))
        )),
        # convertMenuItem(menuItem('Demo', icon=icon("eye"), tabName = "demo"), tabName = "demo"),
        convertMenuItem(
          tabName = "info",
          menuItem('Documentation', icon=icon("question"),
                    tabName = "info"))
        )
    ), # sidebar close

 ################################################################ 
    ## DASHBOARD BODY
 ################################################################ 
    
    dashboardBody(
      useShinyjs(), #imp to use shinyjs functions

      # comment out when testing on local
      tags$head(includeHTML(("google_analytics.html"))),
     
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "./css/custom.css")
      ),

      #  Add logo to the body
      #  tags$img(src="mbpf_logo.jpg",height=50, align="right"),
      
      ## Add tabItems
     # id="body",
      tabItems(
        tabItem(tabName = "home",
                fluidRow(
                  box(
                    title = "Overview",
                    width = 12,
                    solidHeader = TRUE,
                    status = "primary",
                    h3("FragPipe-Analyst"),
                    p(HTML(paste0("FragPipe-Analyst is an easy-to-use, interactive web application developed to perform 
                      differential expression analysis with “one click” and to visualize quantitative proteomic datasets analyzed using ",
                                  a(href="https://fragpipe.nesvilab.org/", target="_blank", "FragPipe"), " computational platform.",
                      " It is compatible with the LFQ-MBR, TMT, and DIA quantification workflows in FragPipe.",
                      " FragPipe-Analyst is based on the original ",
                      a(href="https://bioinformatics.erc.monash.edu/apps/LFQ-Analyst/", target="_blank", "LFQ-Analyst"),
                      " code."))),
                    br(),
                    fluidRow(
                      column(width = 4,
                        h4("Features"),
                        tags$ul(
                          tags$li("Differential expression analysis"),
                          tags$li("Enrichment analysis (GO/Pathways)"),
                          tags$li("Imputation (optional)"),
                          tags$li("Data visualization"),
                          tags$ul(
                            tags$li("PCA"),
                            tags$li("Sample correlation"),
                            tags$li("Heatmaps"),
                            tags$li("Missing values inspection"),
                            tags$li("Sample coverage"),
                            tags$li("Protein intensity plots for selected protein(s)"),
                            tags$li("Imputation effect evaluation")
                          ),
                        )
                      ),
                      column(width = 8,
                        h4("Example Results"),
                        fluidRow(
                          style='margin: 0px;',
                          column(width = 6, box(width="100%", img(src="PCA_plot.png", style="width: 100%;"))),
                          column(width = 6, box(width="100%", img(src="heatmap.png", style="width: 100%;")))
                        )
                      )
                    ),
                    br(),
                    fluidRow(
                      column(width = 4,
                        h4("Sidebar tabs"),
                        tags$ul(
                          tags$li(tags$b("Analysis: "),"perform your own analysis"),
                          tags$li(tags$b("Documentation: "), "Learn more about how to use FragPipe-Analyst")
                          # tags$li(tags$b("Demo: "),"familiarise yourself with FragPipe-Analyst by browsing through pre-analysed results"),
                        )
                      ),
                      column(width = 8,
                        h4("Support"),
                        tags$ul(
                          tags$li(tags$b("Questions/Suggestions/Bug reports: "), "Ask us in ", tags$a(href="https://github.com/Nesvilab/FragPipe-Analyst", target="_blank", "our GitHub forum"), "."),
                          tags$li(tags$b("Documentation/Tutorials: "), "Learn more ", tags$a(href="https://github.com/MonashProteomics/FragPipe-Analyst/tree/main/docs", target="_blank", "here"), "."),
                          tags$li(tags$b("Servers:"), "Our production (stable) server is at ", tags$a(href="https://fragpipe-analyst.org/", target="_blank", "https://fragpipe-analyst.org/"), "but we also provide our latest dev server ", tags$a(href="http://fragpipe-analyst.nesvilab.org/", target="_blank", "http://fragpipe-analyst.nesvilab.org/"), "with most recent updates and bug fixes."),
                        )
                      )
                    )
                  ) # box 1 closed
                 ) # fluidRow close
                ), # home tab close
      tabItem(tabName = "analysis",
      div(id="quickstart_info",
                    fluidPage(
                      box(
                        title = "Getting Started",
                        h3(tags$b(span("Quick Start", style="text-decoration:underline"))),
                        tags$ul(
                          tags$li("Choose the type of experiment you performed. Currently, DDA-based LFQ (MS1-based or spectral count), TMT, and DIA are supported."),
                          tags$li("For DDA LFQ:",
                                  tags$ul(
                                    tags$li("Upload ", tags$b("combined_protein.tsv "), "generated by ",
                                            tags$a(href="https://github.com/Nesvilab/IonQuant", target="_blank", "IonQuant"),
                                            " in ", tags$a(href="https://fragpipe.nesvilab.org/", target="_blank", "FragPipe")),
                                    tags$li("Upload ", tags$b("experiment_annotation.tsv"), "file. Edit the template file generated by FragPipe;",
                                            " Check ", tags$a(href="https://github.com/MonashProteomics/FragPipe-Analyst/tree/main/docs", target="_blank", "here"),
                                            " for details."),
                                    tags$li("Select quantification method (MS1-based Intensity or MaxLFQ Intensity, or spectral counts).")
                                  )
                          ),
                          tags$li("For TMT:",
                                  tags$ul(
                                    tags$li("Upload gene-level report", tags$b("[abundance/ratio]_gene_[normalization].tsv"), "generated by ",
                                            tags$a(href="https://tmt-integrator.nesvilab.org/", target="_blank", "TMT-Integrator (TMT-I)"), "in ",
                                            tags$a(href="https://fragpipe.nesvilab.org/", target="_blank", "FragPipe"),
                                            "(we recommend ", tags$b("abundance_gene_MD.tsv"), "file). If peptide option is enabled in the server, TMT-I peptide report could also be used alternatively."),
                                    tags$li("Upload", tags$b("experiment_annotation.tsv"), "file. Edit the template file generated by FragPipe;",
                                            " Check ", tags$a(href="https://github.com/MonashProteomics/FragPipe-Analyst/tree/main/docs", target="_blank", "here"),
                                            " for details.")
                                  )
                          ),
                          tags$li("For DIA:",
                                  tags$ul(
                                    tags$li("Upload protein group (PG) matrix (", tags$b("diann-output.pg_matrix.tsv"), ") generated by ",
                                            tags$a(href="https://github.com/vdemichev/DiaNN", target="_blank", "DIA-NN"),
                                            " in ", tags$a(href="https://fragpipe.nesvilab.org/", target="_blank", "FragPipe")),
                                    tags$li("Upload ", tags$b("experiment_annotation.tsv"), "file.",
                                            "Check ", tags$a(href="https://github.com/MonashProteomics/FragPipe-Analyst/tree/main/docs", target="_blank", "here"),
                                            " for details.")
                                  )
                          ),
                          tags$li(tags$b("Optional: "),
                                  "Adjust the p-value cut-off, the log2 fold change cut-off, missing value imputation, FDR correction method ",
                                  "in the", tags$b("Advanced Options"),
                                  '. Note that the missing value imputation method is set by default to “Perseus-like” for DDA LFQ and DIA, and to “No imputation” for TMT.'),
                          tags$li("Press ", tags$b("'Run' ")),
                          tags$li(tags$b("Hint: "), " Check the ", tags$b("Documentation ")," tab for a detailed explanation of inputs, 
                                advanced options and outputs"),
                        ),
                        br(),
                        # HTML('<center><img src="./LFQ_analyst.svg" width="500px"></center>'),
                        width = 12,
                        solidHeader = TRUE,
                        status = "danger"
                      )
                    )
        ), # QUICKSTART INFO CLOSE
      shinyjs::hidden(
        div(id="panel_list",
            tabsetPanel(id = "tab_panels",
                        type = "tabs",
                        selected = "LFQ-Analyst",
                        tabPanel("Quantification",
                                 value = "quantification_panel",
                                 br(),
                          fluidRow(
                              box(
                                column(6,uiOutput("downloadTable"),offset = 1), 
                                column(4,uiOutput("downloadButton")), # make the button on same line
                                width = 4),
                              box(
                                column(12, uiOutput("significantBox")),  width = 4
                              ),
                              box(
                                column(5,uiOutput("downloadreport")), # offset for dist between buttons
                                #tags$br(),
                                #column(5,uiOutput('downloadPlots')),
                                width = 4),
                          ), #close first fluidrow
                          # align save button
                          tags$style(type='text/css', "#downloadButton { width:100%; margin-top: 25px;}"), 
                          tags$style(type='text/css', "#downloadreport { width:100%; vertical-align- middle; margin-top: 25px; 
                                     margin-bottom: 25px;}"),
                          #tags$style(type='text/css', "#downloadPlots { width:100%; margin-top: 25px;}"),
                          tags$br(), # Blank lines
                          
                          ## Data table and result plots box
                          fluidRow(id="results_tab",
                              box(
                                title = "Results Table",
                                shinycssloaders::withSpinner(DT::dataTableOutput("contents"),
                                                             color = "#3c8dbc"),
                                #  actionButton("clear", "Deselect Rows"),
                                actionButton("original", "Refresh Table"),
                                width = 6,
                                status = "success",
                                #color=""
                                solidHeader = TRUE
                              ),
                              # column(
                              box(
                                width = 6,
                                collapsible = TRUE,
                                #status="primary",
                                #solidHeader=TRUE,
                                tabBox(
                                  title = "Result Plots",
                                  width = 12,
                                  tabPanel(title = "Volcano plot",
                                           fluidRow(
                                             box(uiOutput("volcano_cntrst"), width = 5),
                                             box(numericInput("fontsize",
                                                              "Font size",
                                                              min = 0, max = 8, value = 4),
                                                 width = 3),
                                             box(checkboxInput("check_names",
                                                               "Display names",
                                                               value = T),
                                                 checkboxInput("p_adj",
                                                               "Adjusted p values",
                                                               value = T),
                                                 width = 4),
                                           ),
                                           fluidRow(
                                             box(
                                               tags$p("Select features from Results Table to highlight them on the plot OR 
                                                      drag the mouse on plot to show expression of features in Table"),
                                               width = 12
                                             )
                                           ),
                                           fluidRow(
                                             shinycssloaders::withSpinner(
                                               plotOutput("volcano",
                                                          height = 600,
                                                          brush = "protein_brush" # click = "protein_click"
                                                        ), color = "#3c8dbc"),
                                             downloadButton('downloadVolcano', 'Save Highlighted Plot'),
                                             actionButton("resetPlot", "Clear Selection")
                                           )),
                                  tabPanel(title= "Heatmap",
                                           fluidRow(
                                             box(checkboxInput("show_row_names",
                                                               "Show rownames",
                                                               value = F),
                                                 width = 6
                                             )
                                           ),
                                           fluidRow(
                                             shinycssloaders::withSpinner(plotOutput("heatmap", height = 600), color = "#3c8dbc")
                                           ),
                                           fluidRow(
                                             # box(numericInput("cluster_number",
                                             #                  "Cluster to download",
                                             #                  min=1, max=6, value = 3), width = 6),
                                             box(
                                               # downloadButton('downloadCluster',"Save Cluster"),
                                                downloadButton('download_hm_svg', "Save svg"),
                                                width = 5)
                                           ),
                                      # align save button
                                      tags$style(type='text/css', "#downloadCluster {margin-top: 25px;}"),
                                      tags$style(type='text/css', "#download_hm_svg {margin-top: 25px;}")
                                  ),
                                  tabPanel(title = "Feature Plot",
                                           fluidRow(
                                             box(radioButtons("type",
                                                              "Plot type",
                                                              choices = c("Box Plot"= "boxplot",
                                                                          "Violin Plot"="violin"),
                                                              selected = "boxplot", 
                                                              inline = TRUE),
                                                 width = 6
                                             ),
                                             box(checkboxInput("check_impute",
                                                               "Show imputed values",
                                                               value = F),
                                                 width = 6
                                             )
                                           ),
                                           fluidRow(
                                             box(
                                                 tags$p("Select one or more rows from Results Table to plot individual 
                                                    protein intesities across conditions and replicates"),
                                                 width = 12
                                             )
                                           ),
                                           fluidRow(
                                             shinycssloaders::withSpinner(plotlyOutput("protein_plot"), color = "#3c8dbc")
                                             # downloadButton('downloadProtein', 'Download Plot')
                                             )
                                           )
                              ) # tabBox end
                            ) # box or column end
         ), # result fluidRow close
        
        ## QC Box
        fluidRow(
          id="qc_tab",
          box(width = 6,
              tabBox(title = "QC Plots", width = 12, id="qc_tabBox", height=800,
                     tabPanel(title = "PCA Plot",
                              fluidRow(
                                column(6, checkboxInput("pca_imputed",
                                                "Show imputed version",
                                                value = F)),
                                column(6, checkboxInput("pca_scale",
                                                        "Show scaled version",
                                                        value = T))
                              ),
                              fluidRow(shinycssloaders::withSpinner(plotlyOutput("pca_plot", height = 600), color = "#3c8dbc"))
                              ),
                     tabPanel(title="Sample Correlation",
                              fluidRow(box(checkboxInput("cor_imputed",
                                                         "Show imputed version",
                                                         value = F),
                                           width = 6
                              )),
                              fluidRow(shinycssloaders::withSpinner(plotOutput("sample_corr", height = 600), color = "#3c8dbc")),
                              fluidRow(downloadButton('download_corr_svg', "Save svg"))
                              ),
                     tabPanel(title= "Sample CVs",
                              fluidRow(
                                box(checkboxInput("cvs_full_range",
                                                  "Show full range",
                                                  value = F),
                                    width = 6
                                )
                              ),
                              fluidRow(
                                shinycssloaders::withSpinner(plotOutput("sample_cvs", height = 600), color = "#3c8dbc")
                              ),
                              fluidRow(
                                downloadButton('download_cvs_svg', "Save svg")
                              )),
                     tabPanel(title = "Feature Numbers",
                              shinycssloaders::withSpinner(plotOutput("numbers", height = 600), color = "#3c8dbc"),
                              downloadButton('download_num_svg', "Save svg")),
                     tabPanel(title = "Sample coverage", value="sample_coverage_tab",
                              shinycssloaders::withSpinner(plotOutput("coverage", height = 600), color = "#3c8dbc"),
                              downloadButton('download_cov_svg', "Save svg")),
                     tabPanel(title = "Missing values - Heatmap",
                              shinycssloaders::withSpinner(plotOutput("missval", height = 600), color = "#3c8dbc"),
                              downloadButton('download_missval_svg', "Save svg")
                              ),
                     tabPanel(title = "Density plot", value="density_tab",
                              shinycssloaders::withSpinner(plotOutput("density", height = 600), color = "#3c8dbc"),
                              downloadButton('download_density_svg', "Save svg")
                              )
                     # ,
                     # tabPanel(title = "p-value Histogram",
                     #          plotOutput("p_hist", height = 600)
                     # )
              ) # Tab box close
          ),
          box(
            width=6,
              tabBox(title = "Enrichment", width = 12, height=800,
                     tabPanel(title= "Pathway enrichment",
                              fluidRow(
                                column(4,uiOutput("contrast_1")),
                                column(4, selectInput("pathway_database", "Pathway database:",
                                                      c("Hallmark"="MSigDB_Hallmark_2020",
                                                        "KEGG"="KEGG_2021_Human",
                                                        "Reactome"="Reactome_2022"))),
                                column(3, radioButtons("pathway_direction",
                                                       "Direction",
                                                       choices = c("Up"="UP", "Down"="DOWN"),
                                                       inline = T,
                                                       selected = "UP")),
                                column(1, shinyWidgets::dropdownButton(
                                  circle = TRUE, status = "default", right = T,
                                  icon = icon("gear"), width = "300px",
                                  numericInput("p_path",
                                               "DE adjusted p-value cutoff",
                                               min = 0, max = 1, value = 0.01),
                                  numericInput("lfc_path",
                                               "DE log2 fold change cutoff",
                                               min = 0, max = 10, value = 1),
                                  checkboxInput("path_adjust",
                                                "Use adjusted p-value",
                                                value = F),
                                  checkboxInput("pathway_whole_proteome",
                                                "Whole proteome as background",
                                                value = T),
                                  tooltip = tooltipOptions(placement="left",
                                                           title = "customize settings"))
                                )
                              ),
                              fluidRow(
                                column(4, actionButton("pathway_analysis", "Run Enrichment")),
                                column(6, tags$p("Note: Currently, only human data is supported."))
                              ),
                              fluidRow(box(width = 12, uiOutput("spinner_pa"), height = 500)),
                              fluidRow(column(12, downloadButton('downloadPA', 'Download Table')))
                     ),
                     tabPanel(title="Gene Ontology",
                              fluidRow(
                                column(4, uiOutput("contrast")),
                                column(4, selectInput("go_database", "GO database:",
                                                      c("Molecular Function"="GO_Molecular_Function_2021",
                                                        "Cellular Component"="GO_Cellular_Component_2021",
                                                        "Biological Process"="GO_Biological_Process_2021"))),
                                column(3, radioButtons("go_direction",
                                                       "Direction",
                                                       choices = c("Up"="UP", "Down"="DOWN"),
                                                       inline=T,
                                                       selected = "UP")),
                                column(1, shinyWidgets::dropdownButton(
                                  circle = TRUE, status = "default", right = T,
                                  icon = icon("gear"), width = "300px",
                                  numericInput("p_go",
                                               "DE adjusted p-value cutoff",
                                               min = 0, max = 1, value = 0.01),
                                  numericInput("lfc_go",
                                               "DE log2 fold change cutoff",
                                               min = 0, max = 10, value = 1),
                                  checkboxInput("go_adjust",
                                                "Use adjusted p-value",
                                                value = F),
                                  checkboxInput("go_whole_proteome",
                                                "Whole proteome as background",
                                                value = T),
                                  tooltip = tooltipOptions(placement="left",
                                                           title = "customize settings")
                                ))
                                ),
                              fluidRow(
                                column(4, actionButton("go_analysis", "Run Enrichment")),
                                column(6, tags$p("Note: Currently, only human data is supported"))
                              ),
                              fluidRow(box(width = 12, uiOutput("spinner_go"), height = 500)),
                              fluidRow(column(12, downloadButton('downloadGO', 'Download Table')))
                              )
                     ) # Tab box close
            ) # box end
          ) # fluidrow qc close
        ),# lfq-analyst panel close
        tabPanel('Absence/Presence',
                 value = "occ_panel",
                 br(),
                 fluidRow(
                   tags$style(
                     ".box {
                                 border-top: none;
                                 box-shadow: 0 0px 0px rgb(0 0 0 / 10%);
                                 }"
                   ),
                   column(3,
                          box(width =NULL,
                              title = "Options",
                              tags$p("Pre-filtering Results table and Venn plot based on preference Filter Condition if has, 
                                                    and/or changing the Sliders of each condition/group below"),
                              br(),
                              div(id = "venn_filter",
                                  tags$h4("Subset Results Table"),
                                  shinyWidgets::prettyCheckboxGroup("filtered_condition_fragpipe",
                                                                    "Filtered Condition",
                                                                    choices = c('Proteins with more than two peptides'),
                                                                    shape = "round",
                                                                    selected = NULL)
                                  ),
                              tags$hr(),
                              tags$h4("Number of samples present"),
                              tags$style(type = "text/css", ".irs-grid-pol.small {height: 0px;}"), # hide minor ticks of a sliderInput
                              uiOutput('sidebar'),
                              status = "success",
                              solidHeader = TRUE)
                   ), # slider bar column closed
                   column(9,
                          box(width = NULL,
                              title = "Results Table",
                              shinycssloaders::withSpinner(DT::dataTableOutput("contents_occ"), color = "#3c8dbc"),
                              downloadButton('download_attendance', 'Download Table'),
                              status = "success",
                              solidHeader = TRUE),
                          box(width = NULL,
                              title = "Venn Plot",
                              tags$p('Select conditions/groups to generate the Venn plot. By default, more than three conditions/groups generates a 3D Venn plot,
                                            set Condition 3 as "NONE" to generate a 2D Venn plot'),
                              column(12,
                                     box(width = 4,id = "con_1",uiOutput("condition_1")),
                                     box(width = 4,id = "con_2", uiOutput("condition_2")),
                                     box(width = 4,id = "con_3", uiOutput("condition_3"))),
                              column(12,
                                     shinycssloaders::withSpinner(plotOutput("venn_plot"), color = "#3c8dbc")),
                              column(12, downloadButton('download_venn_svg', "Save svg")),
                              status = "success",
                              solidHeader = TRUE)
                   ) # Venn plot column closed
                 ) # fuildRow closed
              )  # occurrence panel closed
            ) # panel_list close
          ) # div close
      #bookmarkButton()
        )), #analysis tab close
      
      tabItem(tabName = "info",
              fluidRow( 
               box(
                   title = "Documentation",
                   h3("Need help?"),
                   tags$ul(
                     tags$li("Read our documentation and tutorial ", a(href = 'https://github.com/MonashProteomics/FragPipe-Analyst/tree/main/docs', target='_blank', tags$b('here')), "."), 
                     tags$li("Report issues and ask questions ", a(href = 'https://github.com/Nesvilab/FragPipe-Analyst', target='_blank', tags$b('here')), "."), 
                     tags$li("FragPipe-Analyst is open-source! You are more than welcome to ",a(href = 'https://github.com/MonashProteomics/FragPipe-Analyst', target='_blank', tags$b('contribute')), "."),
                     tags$li('Learn more about our FragPipe', a(href = 'https://fragpipe.nesvilab.org/', target='_blank', tags$b('here')), "."),
                     tags$li("The user manual of original LFQ-Analyst can be accessed",
                             a(href = 'https://bioinformatics.erc.monash.edu/apps/LFQ-Analyst/LFQ-Analyst_manual.pdf', target='_blank', tags$b("here")), ".")
                   ),

                   h4("Contact Us"),
                   p("For any feedback or question regarding FragPipe-Analyst, please contact the 
                     Proteomics & Integrative Bioinformatics Lab (P.I. Alexey Nesvizhskii; University of Michigan):"),
                   tags$ul(
                     tags$li("Professor Alexey Nesvizhskii: ", a(href="mailto: nesvi@med.umich.edu", target='_blank', "nesvi@med.umich.edu"))),

                   h4("News and Updates"),
                   tags$ul(
                     tags$li("12-02-2022: FragPipe-Analyst is first released for beta testing."),
                     tags$li("07-13-2022: FragPipe-Analyst is first created.")
                    ),
                   width = 12,
                   solidHeader = TRUE,
                   status = "primary"
                )
              )
          )# info tab close
     #     , tabItem(tabName = "demo",
     #        div(id="downloadbox_dm",
     #                      fluidRow(
     #                        box(
     #                          column(6,uiOutput("downloadTable_dm"),offset = 1), 
     #                          column(4,uiOutput("downloadButton_dm")), # make the button on same line
     #                          width = 4),
     #                        
     #                        infoBoxOutput("significantBox_dm",width = 4),
     #                        box(
     #                          column(5,uiOutput("downloadreport_dm")), # offset for dist between buttons
     #                          #tags$br(),
     #                         # column(5,uiOutput('downloadPlots_dm')),
     #                          width = 4
     #                        )
     #                      )), #close div and first row 
     #  
     #  # align save button
     #  tags$style(type='text/css', "#downloadButton_dm { width:100%; margin-top: 25px;}"), 
     #  tags$style(type='text/css', "#downloadreport_dm { width:100%; margin-top: 25px; margin-bottom: 25px;}"),
     # # tags$style(type='text/css', "#downloadPlots_dm { width:100%; margin-top: 25px;}"),
     #  
     #  tags$br(), # Blank lines
     #  tags$br(),
     #  
     #  ## Data table and result plots box
     #  fluidRow(
     #    div(id="results_tab_dm",
     #                        box(
     #                          title = "LFQ Results Table",
     #                          DT::dataTableOutput("contents_dm"),
     #                          #  actionButton("clear", "Deselect Rows"),
     #                          actionButton("original_dm", "Refresh Table"),
     #                          width = 6,
     #                          status = "success",
     #                          #color=""
     #                          solidHeader = TRUE
     #                        ),
     #                        # column(
     #                        box(
     #                          width= 6,
     #                          collapsible = TRUE,
     #                          #status="primary",
     #                          #solidHeader=TRUE,
     #                          tabBox(
     #                            title = "Result Plots",
     #                            width = 12,
     #                            tabPanel(title = "Volcano plot",
     #                                     fluidRow(
     #                                       box(uiOutput("volcano_cntrst_dm"), width = 5),
     #                                       box(numericInput("fontsize_dm",
     #                                                        "Font size",
     #                                                        min = 0, max = 8, value = 4),
     #                                           width = 3),
     #                                       box(checkboxInput("check_names_dm",
     #                                                         "Display names",
     #                                                         value = FALSE),
     #                                           checkboxInput("p_adj_dm",
     #                                                         "Adjusted p values",
     #                                                         value = FALSE),
     #                                           width = 4),
     #                                       tags$p("Select protein from LFQ Results Table to highlight on the plot OR 
     #                                              drag the mouse on plot to show expression of proteins in Table")
     #                                       #Add text line
     #                                       # tags$p("OR"),
     #                                       #  tags$p("Drag the mouse on plot to show expression of proteins in Table") 
     #                                       ),
     #                                     
     #                                     fluidRow(
     #                                       plotOutput("volcano_dm", height = 600,
     #                                                  # hover = "protein_hover"),
     #                                                  #),
     #                                                  # click = "protein_click"),
     #                                                  brush = "protein_brush_dm",
     #                                                  click = "protein_click_dm"),
     #                                       downloadButton('downloadVolcano_dm', 'Save Highlighted Plot'),
     #                                       actionButton("resetPlot_dm", "Clear Selection")
     #                                       #)),
     #                                     )),
     #                            tabPanel(title= "Heatmap",
     #                                     fluidRow(
     #                                       plotOutput("heatmap_dm", height = 600)
     #                                     ),
     #                                     fluidRow(
     #                                       box(numericInput("cluster_number_dm",
     #                                                        "Cluster to download",
     #                                                        min=1, max=6, value = 1), width = 6),
     #                                       box(downloadButton('downloadCluster_dm',"Save Cluster"),width = 3)
     #                                     )
     #                            ),
     #                            tabPanel(title = "Protein Plot",
     #                                     fluidRow(
     #                                       box(radioButtons("type_dm",
     #                                                        "Plot type",
     #                                                        choices = c("Box Plot"= "boxplot",
     #                                                                    "Violin Plot"="violin", 
     #                                                                    "Interaction Plot"= "interaction",
     #                                                                    "Intensity Plot"="dot"
     #                                                        ),
     #                                                        selected = "boxplot", 
     #                                                        inline = TRUE),
     #                                           width = 12
     #                                       ),
     #                                       tags$p("Select one or more rows from LFQ Results Table to plot individual 
     #                                              protein intesities across conditions and replicates")
     #                                       ),
     #                                     fluidRow(
     #                                       plotOutput("protein_plot_dm"),
     #                                       downloadButton('downloadProtein_dm', 'Download Plot')
     #                                     )
     #                                     )
     #                            # verbatimTextOutput("protein_info"))
     #                        )
     #                        ) # box or column end
     #    )),
     #  
     #  ## QC Box
     #  fluidRow(
     #    div(id="qc_tab_dm",
     #                        column(
     #                          width=6,
     #                          tabBox(title = "QC Plots", width = 12,
     #                            tabPanel(title = "PCA Plot",
     #                                     plotOutput("pca_plot_dm"), height=600),
     #                            tabPanel(title="Sample Correlation",
     #                                     plotOutput("sample_corr_dm", height = 600)),
     #                            tabPanel(title= "Sample CVs",
     #                                     plotOutput("sample_cvs_dm", height = 600)),
     #                            conditionalPanel(condition="input.exp != 'TMT'",
     #                                             tabPanel(title = "Protein Numbers",
     #                                                      plotOutput("numbers_dm", height = 600))),
     #                            conditionalPanel(condition="input.exp != 'TMT'",
     #                                             tabPanel(title = "Sample coverage",
     #                                                      plotOutput("coverage_dm", height = 600))),
     #                            tabPanel(title = "Normalization",
     #                                     plotOutput("norm_dm", height = 600)),
     #                            # tabPanel(title = "Missing values - Quant",
     #                            #          plotOutput("detect_dm", height = 600)
     #                            # ),
     #                            tabPanel(title = "Missing values - Heatmap",
     #                                     plotOutput("missval_dm", height = 600)),
     #                            tabPanel(title = "Imputation",
     #                                     plotOutput("imputation_dm", height = 600))
     #                            #,
     #                            # tabPanel(title = "p-value Histogram",
     #                            #          plotOutput("p_hist_dm", height = 600)
     #                            # )
     #                          ) # Tab box close
     #                        ),
     #                        column(
     #                          width=6,
     #                          tabBox(title = "Enrichment", width = 12,
     #                                 tabPanel(title="Gene Ontology",
     #                                          box(uiOutput("contrast_dm"), width = 5),
     #                                        box(
     #                                          selectInput("go_database_dm", "GO database:",
     #                                                    c("Molecular Function"="GO_Molecular_Function_2017b",
     #                                                      "Cellular Component"="GO_Cellular_Component_2017b",
     #                                                      "Biological Process"="GO_Biological_Process_2017b")),
     #                                          width= 5),
     #                                       actionButton("go_analysis_dm", "Run Enrichment"),
     #                                          plotOutput("go_enrichment_dm", height=600),
     #                                          downloadButton('downloadGO_dm', 'Download Table'),
     #                                       height=600
     #                                 ),
     #                                 tabPanel(title= "Pathway enrichment",
     #                                          box(uiOutput("contrast_dm_1"), width = 5),
     #                                          box(
     #                                            selectInput("pathway_database_dm", "Pathway database:",
     #                                                        c("KEGG"="KEGG_2016",
     #                                                          "Reactome"="Reactome_2016")),
     #                                            width= 5),
     #                                          actionButton("pathway_analysis_dm", "Run Enrichment"),
     #                                          plotOutput("pathway_enrichment_dm", height=600),
     #                                          downloadButton('downloadPA_dm', 'Download Table'),
     #                                          height=600
     #                                 ) #### Tab demo closed
     #                                 
     #                          ) # Tab box close
     #    ))) # fluidrow qc close
     #  ) # Tab items close
    ),
        fluidRow(
          tags$div(
            tags$footer(
              tags$p("Proteomics & Integrative Bioinformatics Lab at the University of Michigan (P.I. Alexey Nesvizhskii) and the Monash Proteomics & Metabolomics Facility, Monash University (P.I. Ralf Schittenhelm)."),
              align = "left",
              style = "margin-left: 20px;")
              # style = "position:absolute;
              #         bottom:0;
              #         width:100%;
              #         height:50px;   /* Height of the footer */")
            ),
          shiny.info::version(position = "bottom left")
          )
      ) # Dasbboardbody close
    ) #Dashboard page close
  )#Shiny U Close
}
