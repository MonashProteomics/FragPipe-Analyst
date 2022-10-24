# Define UI for data upload app ----
library(dplyr)

ui <- function(request){shinyUI(
  dashboardPage(
    skin = "blue",
    dashboardHeader(title = "Fragpipe-Analyst"),
    # disable = TRUE),# Disable title bar
    dashboardSidebar(
      sidebarMenu(
        id="tabs_selected",
        convertMenuItem(menuItem('Home', icon=icon("home"), selected = TRUE, tabName = "home"), tabName = "home"),
        convertMenuItem(menuItem("Analysis",  tabName="analysis", icon=icon("flask"),
                                 selectInput("exp", "Experiment type:", c("LFQ"="LFQ", "TMT"="TMT", "DIA"="DIA"), selected = "LFQ"),
                                 conditionalPanel(
                                   condition = "input.exp == 'LFQ'",
                                   fileInput('lfq_expr',
                                             'Upload FragPipe combined_protein.tsv',
                                      accept=c('text/tsv',
                                               'text/tab-separated-values,text/plain',
                                               '.tsv')),
                                   fileInput('lfq_manifest',
                                      'Upload Manifest',
                                      accept=c('text/tsv',
                                               'text/tab-separated-values,text/plain',
                                               '.tsv',
                                               '.fp-manifest')),
                                   tags$hr(),
                                   downloadLink("lfq_example", label="Example LFQ data"),
                                   br(),
                                   downloadLink("lfq_manifest", label="Example FragPipe Manifest")
                                   # p(a("Example LFQ data", target= "_blank",
                                   #     href="data/LFQ_datasets/ubiquitin/combined_protein.tsv", 
                                   #     download="combined_protein.tsv")),
                                   # p(a("Example Manifest", target= "_blank",
                                   #     href="data/LFQ_datasets/ubiquitin/lfq_manifest.tsv", 
                                   #     download="lfq_manifest.tsv"))
                                 ),
                                 conditionalPanel(
                                   condition = "input.exp == 'TMT'",
                                   fileInput('tmt_expr',
                                             'Upload TMT report *.tsv',
                                             accept=c('text/tsv',
                                                      'text/tab-separated-values,text/plain',
                                                      '.tsv')),
                                   fileInput('tmt_annot',
                                             'Upload sample annotation',
                                             accept=c('text/tsv',
                                                      'text/tab-separated-values,text/plain',
                                                      '.tsv')),
                                 ),
                                 conditionalPanel(
                                   condition = "input.exp == 'DIA'",
                                   fileInput('dia_expr',
                                             'Upload DIA report *.tsv',
                                             accept=c('text/tsv',
                                                      'text/tab-separated-values,text/plain',
                                                      '.tsv')),
                                   fileInput('dia_manifest',
                                             'Upload Manifest',
                                             accept=c('text/tsv',
                                                      'text/tab-separated-values,text/plain',
                                                      '.tsv',
                                                      '.fp-manifest')),
                                 ),
                 tags$hr(),
                 menuItem("Advanced Options",tabName="advanced", icon = icon("cogs"), 
                          numericInput("p", 
                                       "Adjusted p-value cutoff",
                                       min = 0.0001, max = 0.1, value = 0.05),
                          numericInput("lfc",
                                       "Log2 fold change cutoff",
                                       min = 0, max = 10, value = 0.7),
                          # checkboxInput("paired",
                          #               "Paired test", FALSE),
                          radioButtons("imputation",
                                       "Imputation type",
                                       choices = c("No imputation"="none", "Perseus-type"="man", "MLE"="MLE", "knn"="knn", "min"="min", "zero"="zero"),
                                       selected = "none"),
                          radioButtons("fdr_correction",
                                       "Type of FDR correction",
                                       choices =  c("Benjamini Hochberg"="BH",
                                                    "t-statistics-based"="fdrtool"
                                       ), selected= "BH"),
                          # checkboxInput("single_peptide",
                          #               "Include single peptide identifications", FALSE),
                          numericInput("k_number",
                                       "Number of clusters in heatmap",
                                       min = 1, max = 10, value = 3)
                 ),
               tags$hr(),
               actionButton("analyze", "Start Analysis"),
               #actionButton("load_data", "Load example data")
               tags$script(HTML("
                 $(document).ready(function() {
                    $('#analyze').on('click', function(){$(this).blur()});
                  })
                "))
                 ), tabName = 'analysis'),
                  
        # convertMenuItem(menuItem('Demo', icon=icon("eye"), tabName = "demo"), tabName = "demo"),
       convertMenuItem(menuItem('User Guide', icon=icon("question"), 
		#href = "https://monashbioinformaticsplatform.github.io/LFQ-Analyst/", 
		tabName = "info"), tabName = "info")
      )
    ), # sidebar close
    
 ################################################################ 
    ## DASHBOARD BODY
 ################################################################ 
    
    dashboardBody(
      useShinyjs(), #imp to use shinyjs functions
      # tags$head(includeScript("google_analytics.js")),
     
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
                  h3("Fragpipe-Analyst: An customized version of LFQ-Analyst for FragPipe."),
                p(HTML(paste0("Fragpipe-Analyst is an easy-to-use, interactive web application developed to perform 
                  differential expression analysis with “one click” and to visualize label-free quantitative proteomic 
                  datasets preprocessed with ",
                  a(href="https://fragpipe.nesvilab.org/", target="_blank", "FragPipe"), 
                  ". It's based on the origianl ",
                  a(href="https://bioinformatics.erc.monash.edu/apps/LFQ-Analyst/", target="_blanl", "LFQ-Analyst"),
                  ". Fragpipe-Analyst provides a wealth of user-analytic features 
                  and offers numerous publication-quality result output graphics and tables to facilitate statistical 
                  and exploratory analysis of label-free quantitative datasets. "))), 
                br(),
                # HTML('<center><img src="./LFQ_analyst.svg" width="600px"></center>'),
                br(),
                h4("Sidebar tabs"),
                tags$ul(
                tags$li(tags$b("Analysis: "),"perform your own analysis"), 
                # tags$li(tags$b("Demo: "),"familiarise yourself with Fragpipe-Analyst by browsing through pre-analysed results"), 
                tags$li(tags$b("User Guide: "), "download an in-depth manual") 
                ),
                width = 12,
                solidHeader = TRUE,
                status = "primary"
                 )#box 1 closed
               
             ) #fluidrow close
            ), # home tab close
      tabItem(tabName = "analysis",
      div(id="quickstart_info",
                    fluidPage(
                      box(
                        title = "Getting Started",
                        h3(tags$b(span("Quick Start", style="text-decoration:underline"))),
                        tags$ul(
                          tags$li("Choose which experiments you performed. Currently, LFQ, TMT, and DIA are supported"),
                          tags$li("For LFQ,",
                                  tags$ul(
                                    tags$li("Upload your ", tags$b("combined_protein.tsv "), "generated by ", tags$a(href="https://fragpipe.nesvilab.org/", target="_blank", "FragPipe")),
                                    tags$li("Upload your ", tags$b("FragPipe Manifest (.fp-manifest)"))
                                  )
                          ),
                          tags$li("For TMT,",
                                  tags$ul(
                                    tags$li("Upload your TMT-I gene report", tags$b("[abundance/ratio]_gene_[normalization].tsv"), "generated by ", tags$a(href="https://tmt-integrator.nesvilab.org/", target="_blank", "TMT-Integrator")),
                                    tags$li("Upload the combined annotation tsv file", tags$b("(.tsv)"), "automatically generated by FragPipe."),
                                    tags$li("Make sure your TMT-I gene report is without duplicated column names. If there are duplicated column names, you can change the column name to this format: [original column name]_1 and [original column name]_2. Each replicate should also have a different name/row in the annotation tsv file.")
                                  )
                          ),
                          tags$li("For DIA,",
                                  tags$ul(
                                    tags$li("Upload PG matrix generated by ", tags$a(href="https://tmt-integrator.nesvilab.org/", target="_blank", "TMT-Integrator")),
                                    tags$li("Upload your ", tags$b("FragPipe Manifest (.fp-manifest)"))
                                  )
                          ),
                          tags$li(tags$b("Optional: "),"Adjust the p-value cut-off, the log2 fold change cut-off, 
                                the imputation type, FDR correction method and/or number of clusters in heatmap
                                in the", tags$b("Advanced Options")),
                          tags$li("Press ", tags$b("'Start Analysis' ")),
                          tags$li(tags$b("Hint: "), " Use the ", tags$b("User Guide ")," tab for a detailed explanation of inputs, 
                                advanced options and outputs"),
                        ),
                        # conditionalPanel(
                        #   condition = "input.exp == 'LFQ'",
                        #   p("test conditionalPanel: LFQ")
                        # ),
                        # conditionalPanel(
                        #   condition = "input.exp == 'TMT'",
                        #   p("test conditionalPanel: TMT")
                        # ),
                        # conditionalPanel(
                        #   condition = "input.exp == 'DIA'",
                        #   p("test conditionalPanel: DIA")
                        # ),
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
                            
                            infoBoxOutput("significantBox",width = 4),
                            box(
                              column(5,uiOutput("downloadreport")), # offset for dist between buttons
                              #tags$br(),
                              #column(5,uiOutput('downloadPlots')),
                              width = 4)
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
                                width= 6,
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
                                                               value = FALSE),
                                                 checkboxInput("p_adj",
                                                               "Adjusted p values",
                                                               value = FALSE),
                                                 width = 4),
                                             tags$p("Select protein from Results Table to highlight on the plot OR 
                                                    drag the mouse on plot to show expression of proteins in Table")
                                             #Add text line
                                             # tags$p("OR"),
                                             #  tags$p("Drag the mouse on plot to show expression of proteins in Table") 
                                             ),
                                           
                                           fluidRow(
                                             shinycssloaders::withSpinner(plotOutput("volcano", height = 600,
                                                        # hover = "protein_hover"),
                                                        #),
                                                        # click = "protein_click"),
                                                        brush = "protein_brush",
                                                        click = "protein_click"), color = "#3c8dbc"),
                                             downloadButton('downloadVolcano', 'Save Highlighted Plot'),
                                             actionButton("resetPlot", "Clear Selection")
                                             #)),
                                           )),
                                  tabPanel(title= "Heatmap",
                                           fluidRow(
                                             shinycssloaders::withSpinner(plotOutput("heatmap", height = 600), color = "#3c8dbc")
                                           ),
                                           fluidRow(
                                             box(numericInput("cluster_number",
                                                              "Cluster to download",
                                                              min=1, max=6, value = 3), width = 6),
                                             box(downloadButton('downloadCluster',"Save Cluster"),
                                                 downloadButton('download_hm_svg', "Save svg"),
                                                  width = 5)
                                           ),
                                      # align save button
                                      tags$style(type='text/css', "#downloadCluster {margin-top: 25px;}"),
                                      tags$style(type='text/css', "#download_hm_svg {margin-top: 25px;}")
                                  ),
                                  tabPanel(title = "Protein Plot",
                                           fluidRow(
                                             box(radioButtons("type",
                                                              "Plot type",
                                                              choices = c("Box Plot"= "boxplot",
                                                                          "Violin Plot"="violin", 
                                                                          "Interaction Plot"= "interaction",
                                                                          "Intensity Plot"="dot"
                                                              ),
                                                              selected = "boxplot", 
                                                              inline = TRUE),
                                                 width = 12
                                             ),
                                             tags$p("Select one or more rows from Results Table to plot individual 
                                                    protein intesities across conditions and replicates")
                                             ),
                                           fluidRow(
                                             shinycssloaders::withSpinner(plotOutput("protein_plot"), color = "#3c8dbc"),
                                             downloadButton('downloadProtein', 'Download Plot')
                                             )
                                            )
# Abundance plot is under consideration
#                                   navbarMenu("Abundance Plot", 
#     								                 tabPanel(title = "Abundance rank",
#     								                          fluidRow(
#     								                            tags$p("Select protein from LFQ Results Table to highlight on the plot OR 
#                                                       drag the mouse on plot to show expression of proteins in Table")
#     								                          ),
#     								                          fluidRow(
#     								                            plotOutput("abundance_rank_dm",
#     								                                       height = 600,
#     								                                       brush = "protein_brush_rank_dm",
#     								                                       click = "protein_click_rank_dm"),
#     								                            downloadButton('downloadAbundance_rank_dm', 'Save Highlighted Plot'),
#     								                            actionButton("resetPlot_rank_dm", "Clear Selection")
#     								                          )
#     								                 ),
#     								                 tabPanel("Abundance comparison",
#     								                          fluidRow(
#     								                            column(uiOutput("abundance_cntrst_dm"), width = 12),
#     								                            tags$p("Select protein from LFQ Results Table to highlight on the plot OR 
#                                                       drag the mouse on plot to show expression of proteins in Table")
#     								                          ),
#     								                          fluidRow(
#     								                            plotOutput("abundance_comp_dm", 
#     								                                       height = 600,
#     								                                       brush = "protein_brush_comp_dm",
#     								                                       click = "protein_click_comp_dm"),
#     								                            downloadButton('downloadAbundance_comp_dm', 'Save Highlighted Plot'),
#     								                            actionButton("resetPlot_comp_dm", "Clear Selection")
#     								                          )
#     								                 )
# 								                ) # navbarMenu close
								  # verbatimTextOutput("protein_info"))
                              ) # tabBox end
                              ) # box or column end
         ), # result fluidRow close
        
        ## QC Box
        fluidRow(
          id="qc_tab",
          column(
            width=6,
            tabBox(title = "QC Plots", width = 12, id="qc_tabBox", height=700,
                   tabPanel(title = "PCA Plot",
                            shinycssloaders::withSpinner(plotlyOutput("pca_plot", height=600), color = "#3c8dbc")
                                           # downloadButton('download_pca_svg', "Save svg")
                            ),
                   tabPanel(title="Sample Correlation",
                            shinycssloaders::withSpinner(plotOutput("sample_corr", height = 600), color = "#3c8dbc"),
                            downloadButton('download_corr_svg', "Save svg")
                            ),
                   tabPanel(title= "Sample CVs",
                            shinycssloaders::withSpinner(plotOutput("sample_cvs", height = 600), color = "#3c8dbc"),
                            downloadButton('download_cvs_svg', "Save svg")
                            ),
                   # conditionalPanel(condition="input.exp != 'TMT'",
                   tabPanel(title = "Protein Numbers",
                            shinycssloaders::withSpinner(plotOutput("numbers", height = 600), color = "#3c8dbc"),
                            downloadButton('download_num_svg', "Save svg")),
                   # conditionalPanel(condition="input.exp != 'TMT'",
                   tabPanel(title = "Sample coverage", value="sample_coverage_tab",
                            shinycssloaders::withSpinner(plotOutput("coverage", height = 600), color = "#3c8dbc"),
                            downloadButton('download_cov_svg', "Save svg")),
                   tabPanel(title = "Normalization", value="norm_tab",
                            shinycssloaders::withSpinner(plotOutput("norm", height = 600), color = "#3c8dbc"),
                            downloadButton('download_norm_svg', "Save svg")
                            ),
                   # tabPanel(title = "Missing values - Quant",
                   #          plotOutput("detect", height = 600)
                   # ),
                   tabPanel(title = "Missing values - Heatmap",
                            shinycssloaders::withSpinner(plotOutput("missval", height = 600), color = "#3c8dbc"),
                            downloadButton('download_missval_svg', "Save svg")
                            ),
                   tabPanel(title = "Imputation", value="imputation_tab",
                            shinycssloaders::withSpinner(plotOutput("imputation", height = 600), color = "#3c8dbc"),
                            downloadButton('download_imp_svg', "Save svg")
                            )#,
                   # tabPanel(title = "p-value Histogram",
                   #          plotOutput("p_hist", height = 600)
                   # )
                   ) # Tab box close
            ),
          column(
            width=6,
            tabBox(title = "Enrichment", width = 12, height=600,
                   tabPanel(title="Gene Ontology",
                            fluidRow(
                              column(4,uiOutput("contrast")),
                              column(4, selectInput("go_database", "GO database:",
                                                    c("Molecular Function"="GO_Molecular_Function_2021",
                                                      "Cellular Component"="GO_Cellular_Component_2021",
                                                      "Biological Process"="GO_Biological_Process_2021"))),
                              column(2, radioButtons("go_direction",
                                                     "Direction",
                                                     choices = c("Up"="UP", "Down"="DOWN"),
                                                     selected = "UP")),
                              column(12, actionButton("go_analysis", "Run Enrichment")),
                              column(12,
                                     box(width = 12, uiOutput("spinner_go"),height = 500)
                                     ),
                              column(12,
                                     downloadButton('downloadGO', 'Download Table')
                                     )
                              )
                            ),
                   tabPanel(title= "Pathway enrichment",
                            fluidRow(
                              column(4,uiOutput("contrast_1")),
                              column(4, selectInput("pathway_database", "Pathway database:",
                                                    c("Hallmark"="MSigDB_Hallmark_2020",
                                                      "KEGG"="KEGG_2021_Human",
                                                      "Reactome"="Reactome_2022"))),
                              column(2, radioButtons("pathway_direction",
                                                     "Direction",
                                                     choices = c("Up"="UP", "Down"="DOWN"),
                                                     selected = "UP")),
                              column(12, actionButton("pathway_analysis", "Run Enrichment")),
                              column(12,
                                     box(width = 12, uiOutput("spinner_pa"),height = 500)
                              ),
                              column(12,
                                     downloadButton('downloadPA', 'Download Table')
                              )
                              )
                            )
                   ) # Tab box close
            ) # column end
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
                              tags$p("Pre-filtering Results table and Venn plot based on preference Filter Condition, 
                                                    and/or changing the Sliders of each condition/group below"),
                              br(),
                              tags$h4("Subset Results Table"),
                              shinyWidgets::prettyCheckboxGroup("filtered_condition_fragpipe",
                                                                "Filtered Condition",
                                                                choices = c('Proteins with more than two peptides'),
                                                                shape = "round",
                                                                selected = NULL,
                              ),
                              tags$hr(),
                              # tags$h4("Number of replicates present"),
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
              # conditionalPanel(condition = "input.exp == 'LFQ' ", value = "test-lfq",
              #                  tabPanel('test Conditional LFQ', textOutput("This is test conditionalpanel LFQ"))),
              # conditionalPanel(condition = "input.exp == 'TMT' ", value = "test-tmt",
              #                  tabPanel('test Conditional TMT', textOutput("This is test conditionalpanel TMT")))
            ) # panel_list close
          ) # div close
      #bookmarkButton()
        )), #analysis tab close
      
      tabItem(tabName = "info",
              fluidRow( 
               box(
                   title = "User Guide",
                   h3("Fragpipe-Analyst Documentation"),
                   div(p(HTML(paste0('Learn more about our FragPipe',
                                     a(href = 'https://monashbioinformaticsplatform.github.io/LFQ-Analyst/', target='_blank', 'here'))))),
                   div(p(HTML(paste0("The user manual of original LFQ-Analyst can be accessed",
			                            a(href = 'https://bioinformatics.erc.monash.edu/apps/LFQ-Analyst/LFQ-Analyst_manual.pdf', 
			                              target='_blank', tags$b("here.")))))),
                   h4("Contact Us"),
			p("For any feedback or question regarding Fragpipe-Analyst, please contact the 
			  Proteomics & Integrative Bioinformatics Lab, University of Michigan:"),
      tags$ul(
			tags$li("Professor Alexey Nesvizhskii: nesvi@med.umich.edu")
      ),
			
			h4("How to Cite Fragpipe-Analyst?"),
			div(p(HTML(paste0("Please cite the link and the original LFQ-Analyst: Shah AD, Goode RJA, Huang C, Powell DR, Schittenhelm RB. 
		LFQ-Analyst: An easy-to-use interactive web-platform to analyze and 
		visualize proteomics data preprocessed with MaxQuant. DOI:",
			                            a(href = 'https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00496', 
			                              target='_blank', tags$b("10.1021/acs.jproteome.9b00496")))))),


      h4("News and Updates"),
			
      tags$ul(
        tags$li("07-13-2022: Fragpipe-Analyst first created")
#       tags$li("27-12-2021: LFQ-Analyst being accessed by more than 5000 users worldwide"),
#       tags$li("30-09-2021: LFQ-Analyst being accessed by more than 4000 users worldwide"),
#       tags$li("03-05-2021: LFQ-Analyst being accessed by more than 3000 users worldwide"),
#       tags$li("24-02-2021: Correlation plot now use all protein expression data"),
#       tags$li("25-11-2020: LFQ-Analyst being accessed by more than 2000 users worldwide"),
# 			tags$li("07-04-2020: LFQ-Analyst being accessed by more than 1000 users worldwide"),
#       tags$li("03-01-2020: LFQ-Analyst manuscript published in volume 19 of JPR"),
#       tags$li("28-10-2019: LFQ-Analyst paper published online in Journal of Proteome Research (JPR)"),
#       tags$li("02-10-2019: Svg figures download feature added"),
#       tags$li("09-09-2019: Paired test support added"),
#       tags$li("09-09-2019: Option to include single peptide observations in the analysis"),
#       tags$li("19-02-2019: LFQ-Analyst made public")
      ),   
                width = 12,
                solidHeader = TRUE,
                status = "primary"
                ) #includeMarkdown("www/Info.md")
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
              tags$p("Proteomics & Integrative Bioinformatics Lab (University of Michigan) and the Monash Proteomics & Metabolomics Facility (Monash University)"),
              align = "left",
              style = "margin-left: 20px;")
              # style = "position:absolute;
              #         bottom:0;
              #         width:100%;
              #         height:50px;   /* Height of the footer */")
            )
          )
      ) # Dasbboardbody close
    ) #Dashboard page close
  )#Shiny U Close
}
