#Define server logic to read selected file ----
server <- function(input, output, session) {
  options(shiny.maxRequestSize=100*1024^2)## Set maximum upload size to 100MB
  
  #  Show elements on clicking Start analysis button
  observeEvent(start_analysis(),{
    if(input$analyze==0 | !start_analysis()){
      return()
    }
    shinyjs::hide("quickstart_info")
    shinyjs::show("panel_list")
  })
   
  observeEvent(start_analysis() ,{
    if (input$imputation == "none") {
      hideTab(inputId="qc_tabBox", target="imputation_tab")
    } else {
      showTab(inputId="qc_tabBox", target="imputation_tab")
    }
    if (input$exp == "LFQ"){
      exp <- exp_design_input()
      # Hide LFQ page if only have one replicate in each sample
      if (all(is.na(exp$replicate))) {
        showTab(inputId = "tab_panels", target = "quantification_panel")
        updateTabsetPanel(session, "tab_panels", selected = "quantification_panel")
      } else if (max(exp$replicate)==1){
        hideTab(inputId = "tab_panels", target = "quantification_panel")
      } else {
        showTab(inputId = "tab_panels", target = "quantification_panel")
        showTab(inputId="qc_tabBox", target="norm_tab")
        showTab(inputId="qc_tabBox", target="sample_coverage_tab")
        # make sure occ_panel visible after users updating their analysis
        showTab(inputId = "tab_panels", target = "occ_panel")
        updateTabsetPanel(session, "tab_panels", selected = "quantification_panel")
      }
      shinyjs::show("venn_filter")
    } else if (input$exp == "TMT") {
      hideTab(inputId = "tab_panels", target = "occ_panel")
      hideTab(inputId="qc_tabBox", target="norm_tab")
      hideTab(inputId="qc_tabBox", target="sample_coverage_tab")
      showTab(inputId = "tab_panels", target = "quantification_panel")
      updateTabsetPanel(session, "tab_panels", selected = "quantification_panel")
    } else { # DIA
      hideTab(inputId="qc_tabBox", target="norm_tab")
      showTab(inputId="qc_tabBox", target="sample_coverage_tab")
      showTab(inputId = "tab_panels", target = "occ_panel")
      updateTabsetPanel(session, "tab_panels", selected = "quantification_panel")
      shinyjs::hide("venn_filter")
    }
  })

   # observeEvent(start_analysis(),{ 
   #     if(input$analyze==0 | !start_analysis()){
   #       return()
   #     }
   #   shinyjs::show("results_tab")
   # })
   # 
   # observeEvent(start_analysis() ,{ 
   #   if(input$analyze==0 | !start_analysis()){
   #     return()
   #   }
   #   shinyjs::show("qc_tab")
   # })
   # 
   # observeEvent(start_analysis, {
   #   if(input$analyze==0 | !start_analysis()){
   #     return()
   #   }
   #   shinyjs::show("enrichment_tab")
   # })

   # observeEvent(input$analyze,{
   #   shinyjs::hide("howto")
   # })
   
  # observe({
  #   if(input$body=="info"){
      # output$howto<-renderUI({
      #   box(width=12,
      #       includeMarkdown("www/Info.Rmd")
      #   )
      # })
      # 
      # observeEvent(input$analyze,{
      #   output$howto<-renderUI({NULL})
      # })
    #   }
    # else{
    #   output$howto<-renderUI({NULL})
    #     }
    # })
 
   ## Shinyalert
   start_analysis <- eventReactive(input$analyze,{ 
     if(input$analyze==0 ){
       return(F)
     } else {
       if (input$exp == "LFQ"){
         inFile <- input$lfq_expr
         exp_design_file <- input$lfq_manifest
       } else if (input$exp == "TMT") {
         inFile <- input$tmt_expr
         exp_design_file <- input$tmt_annot
       } else if (input$exp == "DIA") {
         inFile <- input$dia_expr
         exp_design_file <- input$dia_manifest
       }
       if (is.null(inFile) | is.null(exp_design_file)) {
         shinyalert("Input file missing!", "Please checkout your input files", type="info",
                    closeOnClickOutside = TRUE,
                    closeOnEsc = TRUE,
                    timer = 5000)
         return(F)
       }
     }
     
     shinyalert("In Progress!", "Data analysis has started, wait until table and plots
                appear on the screen", type="info",
                closeOnClickOutside = TRUE,
                closeOnEsc = TRUE,
                timer = 10000) # timer in miliseconds (10 sec)
     return(T)
   })
   
   # observe({
   # if (input$tabs_selected=="demo"){
   #   shinyalert("Demo results loading!...", "Wait until table and plots
   #              appear on the screen", type="info",
   #              closeOnClickOutside = TRUE,
   #              closeOnEsc = TRUE,
   #              timer = 6000)
   # }
   # })
 
   
   ####======= Render Functions
   
   output$volcano_cntrst <- renderUI({
     if (!is.null(comparisons())) {
       df <- SummarizedExperiment::rowData(dep())
       cols <- grep("_significant$",colnames(df))
       selectizeInput("volcano_cntrst",
                      "Comparison",
                      choices = gsub("_significant", "", colnames(df)[cols]))
     }
   })
   
   ##comparisons
   output$contrast <- renderUI({
     if (!is.null(comparisons())) {
       df <- SummarizedExperiment::rowData(dep())
       cols <- grep("_significant$",colnames(df))
       selectizeInput("contrast",
                      "Comparison",
                      choices = gsub("_significant", "", colnames(df)[cols]))
     }
   })
   
   output$contrast_1 <- renderUI({
     if (!is.null(comparisons())) {
       df <- SummarizedExperiment::rowData(dep())
       cols <- grep("_significant$",colnames(df))
       selectizeInput("contrast_1",
                      "Comparison",
                      choices = gsub("_significant", "", colnames(df)[cols]))
     }
   })
   
   output$downloadTable <- renderUI({
     if(!is.null(dep())){
     selectizeInput("dataset",
                    "Choose a dataset to save" ,
                    c("Results","Original_matrix",
                      "Imputed_matrix",
                      "Full_dataset"))
     }
    })
   
   output$downloadButton <- renderUI({
     if(!is.null(dep())){
     downloadButton('downloadData', 'Save')
     }
   })
   
   output$downloadZip <- renderUI({
     if(!is.null(dep())){
     downloadButton('downloadZip1', 'Download result plots')
     }
   })
    output$downloadreport <- renderUI({
      if(!is.null(dep())){
     downloadButton('downloadReport', 'Download Report')
      }
    })
   
    output$downloadPlots <- renderUI({
      if(!is.null(dep())){
      downloadButton('downloadPlots1', 'Download Plots')
      }
    })
    
    
    ## Read input files on shiny server
    ## NOTE: have to use reactive framework, otherwise throws out error
    # observeEvent(input$analyze,{
    #   fragpipe_data<-reactive({
    #     inFile<-input$file1
    #     if(is.null(inFile))
    #       return(NULL)
    #     read.table(inFile$datapath,
    #                header = TRUE,
    #                fill= TRUE, # to fill any missing data
    #                sep = "\t"
    #     )
    #   })
    # })
    
    ## make reactive elements
    fragpipe_data_input<-reactive({NULL})
    exp_design_input<-reactive({NULL})
    exp_design_example<-reactive({NULL})
    fragpipe_data_example<-reactive({NULL})
    
    fragpipe_data_input<-eventReactive(input$analyze,{
      if (input$exp == "LFQ") {
        inFile <- input$lfq_expr
      } else if (input$exp == "TMT") {
        inFile <- input$tmt_expr
      } else if (input$exp == "DIA") {
        inFile <- input$dia_expr
      }
      if(is.null(inFile))
        return(NULL)
      temp_data <- read.table(inFile$datapath,
                 header = TRUE,
                 fill= TRUE, # to fill any missing data
                 sep = "\t",
                 quote = "",
                 comment.char = "",
                 blank.lines.skip = F,
                 check.names = F)
      colnames(temp_data) <- make.unique.2(colnames(temp_data), "_")
      # validate(maxquant_input_test(temp_data))
      if (input$exp == "TMT") {
        validate(tmt_input_test(temp_data))
        # convert columns into numeric
        mut.cols <- colnames(temp_data)[!colnames(temp_data) %in% c("Index", "NumberPSM", "ProteinID", "MaxPepProb", "ReferenceIntensity")]
        temp_data[mut.cols] <- sapply(temp_data[mut.cols], as.numeric)
      } else if (input$exp == "LFQ") {
        # handle - (dash) in experiment column
        colnames(temp_data) <- gsub("-", ".", colnames(temp_data))
        validate(fragpipe_input_test(temp_data))
        # remove contam
        temp_data <- temp_data[!grepl("contam", temp_data$Protein),]
      } else { # DIA
        validate(fragpipe_DIA_input_test(temp_data))
        # temp_data <- temp_data[!grepl("contam", temp_data$Protein),]
      }
      return(temp_data)
    })

    # observeEvent(input$analyze,{
    #   exp_design<-reactive({
    #     inFile<-input$file2
    #     if (is.null(inFile))
    #       return(NULL)
    #     temp_df<-read.table(inFile$datapath,
    #                         header = TRUE,
    #                         sep="\t",
    #                         stringsAsFactors = FALSE)
    #     exp_design_test(temp_df)
    #     temp_df$label<-as.character(temp_df$label)
    #     return(temp_df)
    #   })    
    # })

    exp_design_input<-eventReactive(input$analyze,{
      if (input$exp == "LFQ"){
        inFile <- input$lfq_manifest
      } else if (input$exp == "TMT") {
        inFile <- input$tmt_annot
      } else if (input$exp == "DIA") {
        inFile <- input$dia_manifest
      }
      if (is.null(inFile))
        return(NULL)
      if (input$exp == "TMT") {
        temp_df <- read.table(inFile$datapath,
                              header = T,
                              sep="\t",
                              stringsAsFactors = FALSE)
        # To support txt file
        if (ncol(temp_df) == 1){
          # submitting annotation.txt (not experiment_annotation.tsv) will crash here
          tryCatch({
            temp_df <- read.table(inFile$datapath,
                                  header = T,
                                  sep=" ",
                                  stringsAsFactors = FALSE)
          }, error=function(e){
            validate(need(F,
                     "Error: coudn't read the experiment_annotation.tsv. Note that experiment_annotation.tsv is not annotation.txt used to denote channel assignment in each plex set."))
          })
        }
        # change it to lower case
        colnames(temp_df) <- tolower(colnames(temp_df))
        # to support - (dash) in condition column
        temp_df$condition <- gsub("-", ".", temp_df$condition)
        validate(need(try(test_TMT_annotation(temp_df)),
                           paste0("The input annotation file should have following columns: ",
                                  "plex, channel, sample, condition, replicate, condition\n",
                                  "your current input annotation file is with following columns: ", paste(colnames(temp_df), collapse=", "))))
        temp_df$label <- temp_df$sample
        # if duplicate label exists
        if (anyDuplicated(temp_df$label)) {
          # add _number for repeat labels, but need to remove _1
          temp_df$label <- paste(temp_df$label, temp_df$replicate, sep="_")
          temp_df$label <- gsub("_1$", "", temp_df$label)
          samples_with_replicate <- temp_df$label[grepl("_", temp_df$label)]
          samples_with_replicate <- unique(gsub("_\\d+$", "", samples_with_replicate))
          temp_df[temp_df$label %in% samples_with_replicate, "label"] <- paste0(temp_df[temp_df$label %in% samples_with_replicate, "label"], "_1")
        }
      } else if (input$exp == "LFQ"){
        temp_df <- read.table(inFile$datapath,
                              header = T,
                              sep="\t",
                              stringsAsFactors = FALSE)
        # exp_design_test(temp_df)
        # temp_df$label<-as.character(temp_df$label)
        # temp_df$condition<-trimws(temp_df$condition, which = "left")

        # reuse previous logic of experiment column (no label, but experiment column from manifest at that time)
        temp_df$experiment <- temp_df$label

        # make sure replicate column is not empty
        if (!all(is.na(temp_df$replicate))) {
          # handle - (dash) in sample (experiment) column
          temp_df$sample <- gsub("-", ".", temp_df$sample)
          temp_df$label <- temp_df$sample
          if (input$lfq_type == "Intensity") {
            temp_df$label <- paste(temp_df$label, "Intensity", sep=" ")
          } else if (input$lfq_type == "MaxLFQ") {
            temp_df$label <- paste(temp_df$label, "MaxLFQ.Intensity", sep=" ")
          }  else if (input$lfq_type == "Spectral Count") {
            temp_df$label <- paste(temp_df$label, "Spectral.Count", sep=" ")
          }
        }
      } else if (input$exp == "DIA") {
        temp_df <- read.table(inFile$datapath,
                              header = T,
                              sep="\t",
                              stringsAsFactors = FALSE)
        # change it to lower case
        colnames(temp_df) <- tolower(colnames(temp_df))
        # to support - (dash) in condition column
        temp_df$condition <- gsub("-", ".", temp_df$condition)
        # make sure replicate column is not empty
        if (!all(is.na(temp_df$replicate))) {
          temp_df$label <- temp_df$file
        }
      }
      return(temp_df)
    })
   
    
### Load data from Rdata
  # observeEvent(input$load_data,{
      # example_data<-reactive({
      #   load("data/example_data.RData", envir = .GlobalEnv)
      # })
      # fragpipe_data<-reactive({example_data[1]})
      # exp_design<-reactive({example_data[2]})
   #  env<-reactive({
   #    LoadToEnvironment("data/example_data.RData", env = globalenv())
   #  })
   # 
   #   observeEvent(input$load_data,{
   # # message(env()[["exp_design"]])
   #     fragpipe_data_example<-reactive({
   #     env()[["maxquant_output"]]
   #   })
   #   })
   #   observeEvent(input$load_data,{
   #     exp_design_example<-reactive({
   #     env()[["exp_design"]]
   #   })
   #  })
   # }) ## leave this commented
   #  
   #  fragpipe_data<-eventReactive(input$load_data,{
   #    env()[['maxquant_output']]
   #  })
   # 
   # exp_design<-eventReactive(input$load_data,{
   #    env()[['exp_design']]
   #  })
   
   
### Reactive components
   processed_data<- eventReactive(start_analysis(),{
     ## check which dataset
     if(!is.null (fragpipe_data_input() )){
       fragpipe_data <- reactive({fragpipe_data_input()})
     }
     
     if(!is.null (exp_design_input() )){
       exp_design<-reactive({exp_design_input()})
     }
     
     filtered_data<-fragpipe_data()
     if (input$exp == "LFQ"){
       data_unique <- DEP::make_unique(filtered_data, "Gene","Protein ID")
       
       if (input$lfq_type == "Intensity") {
         lfq_columns <- setdiff(grep("Intensity", colnames(data_unique)),
                              grep("MaxLFQ", colnames(data_unique)))
         lfq_columns <- setdiff(lfq_columns, grep("Total Intensity", colnames(data_unique)))
         lfq_columns <- setdiff(lfq_columns, grep("Unique Intensity", colnames(data_unique)))
       } else if (input$lfq_type == "MaxLFQ") {
         lfq_columns<-grep("MaxLFQ", colnames(data_unique))
         if (length(lfq_columns) == 0) {
           stop(safeError("No MaxLFQ column available. Please make sure your files have MaxLFQ intensity columns."))
         }
       } else if (input$lfq_type == "Spectral Count") {
         lfq_columns<-grep("Spectral", colnames(data_unique))
         lfq_columns <- setdiff(lfq_columns, grep("Total Spectral Count", colnames(data_unique)))
         lfq_columns <- setdiff(lfq_columns, grep("Unique Spectral Count", colnames(data_unique)))
       }
       
       ## Check for matching columns in expression report and experiment manifest file
       test_match_lfq_column_manifest(data_unique, lfq_columns, exp_design())
       data_se <- make_se_customized(data_unique, lfq_columns, exp_design(), log2transform=T)
       return(data_se)
     } else if (input$exp == "DIA") {
       data_unique <- DEP::make_unique(filtered_data, "Genes", "Protein.Group")
       cols <- colnames(data_unique)
       selected_cols <- which(!(cols %in% c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description", "ID", "name")))
       # TODO: use DIA function
       test_match_tmt_column_design(data_unique, selected_cols, exp_design())
       data_se <- make_se_customized(data_unique, selected_cols, exp_design(), log2transform=T)
       dimnames(data_se) <- list(dimnames(data_se)[[1]], colData(data_se)$sample_name)
       colData(data_se)$label <- colData(data_se)$sample_name
       return(data_se)
     } else { # TMT
       temp_exp_design <- exp_design()
       # sample without specified condition will be removed
       temp_exp_design <- temp_exp_design[!is.na(temp_exp_design$condition), ]
       temp_exp_design <- temp_exp_design[!temp_exp_design$condition == "",]
       # temp_exp_design[is.na(temp_exp_design), "replicate"] <- 1
       # need to handle duplicate columns first
       # filtered_data <- avearrays(filtered_data)
       # filtered_data <- as.data.frame(filtered_data)
       # print(apply(filtered_data, 2, is.numeric))
       data_unique <- make_unique(filtered_data, "Index", "ProteinID")
       # handle unmatched columns
       overlapped_samples <- intersect(colnames(data_unique), temp_exp_design$label)
       data_unique <- data_unique[,colnames(data_unique) %in% c("Index", "NumberPSM", "ProteinID", "MaxPepProb", "ReferenceIntensity", "name", "ID", overlapped_samples)]
       temp_exp_design <- temp_exp_design[temp_exp_design$label %in% overlapped_samples,]
       cols <- colnames(data_unique)
       selected_cols <- which(!(cols %in% c("Index", "NumberPSM", "ProteinID", "MaxPepProb", "ReferenceIntensity", "name", "ID")))
       data_unique[selected_cols] <- apply(data_unique[selected_cols], 2, as.numeric)
       test_match_tmt_column_design(data_unique, selected_cols, temp_exp_design)
       # TMT-I report is already log2 transformed
       data_se <- make_se_customized(data_unique, selected_cols, temp_exp_design)
       return(data_se)
     }
   })
   
   filtered_data <- eventReactive(input$analyze,{
     # if (input$exp == "LFQ"){ # Check number of replicates
     #   if (input$replicate_filter){
     #     if(!is.null (exp_design_input() )){
     #       exp_design<-reactive({exp_design_input()})
     #     }
     #     if(max(exp_design()$replicate)<3){
     #       threshold<-0
     #     } else if(max(exp_design()$replicate)==3){
     #       threshold<-1
     #     } else if(max(exp_design()$replicate)<6 ){
     #       threshold<-2
     #     } else if (max(exp_design()$replicate)>=6){
     #       threshold<-trunc(max(exp_design()$replicate)/2)
     #     }
     #     filtered_se <- filter_missval_customized(processed_data(), thr = threshold)
     #     return(filtered_se)
     #   }
     # }
     filtered_se <- processed_data()
     # filter by global missingness
     if (input$min_global_appearance != 0){
       filtered_se <- global_filter(processed_data(), 100 - input$min_global_appearance)
     }
     
     # filter by checking percentage of missingness in each condition
     if (input$min_appearance_each_condition != 0){
      filtered_se <- filter_by_condition(filtered_se, input$min_appearance_each_condition)
     }
     return(filtered_se)
   })
   
   unimputed_table<-reactive({
     temp<-assay(filtered_data())
     temp1<-2^(temp)
     colnames(temp1)<-paste(colnames(temp1),"original_intensity",sep="_")
     temp1<-cbind(ProteinID=rownames(temp1),temp1) 
     #temp1$ProteinID<-rownames(temp1)
     return(as.data.frame(temp1))
   })
   
   normalised_data<-reactive({
     if (input$exp == "LFQ") {
       return(normalize_vsn(filtered_data()))
     } else {
       return(filtered_data())
     }
   })
   
   imputed_data<-eventReactive(input$analyze,{
     if (input$imputation == "none"){
       imputed <- filtered_data()
       rowData(imputed)$imputed <- apply(is.na(assay(filtered_data())), 1, any)
       rowData(imputed)$num_NAs <- rowSums(is.na(assay(filtered_data())))
     } else {
       # need a customized function here since DIA data has several slashs in the column
       # TMT report might has same issue for earlier version of FragPipe (<= 18.0)
      imputed <- impute_customized(filtered_data(),input$imputation)
     } 
     return(imputed)
   })
   
   imputed_table<-reactive({
     temp<-assay(imputed_data())
     #tibble::rownames_to_column(temp,var = "ProteinID")
     temp1<-2^(temp)
     colnames(temp1)<-paste(colnames(temp1),"imputed_intensity",sep="_")
     temp1<-cbind(ProteinID=rownames(temp1),temp1) #temp1$ProteinID<-rownames(temp1)
     return(as.data.frame(temp1))
   })
   
   diff_all<-reactive({
     test_diff_customized(imputed_data(), type = "all")
   })

   dep<-eventReactive(input$analyze,{
     # TODO: test_limma for paired
     # diff_all <- test_diff_customized(imputed_data(), type = "manual",
     #                      test = c("SampleTypeTumor"), design_formula = formula(~0+SampleType))
     if(input$fdr_correction=="BH"){
       diff_all <- test_limma_customized(imputed_data(), type='all', paired = F)
     } else { # t-statistics-based
       diff_all <- test_diff_customized(imputed_data(), type = "all")
     }
     add_rejections(diff_all,alpha = input$p, lfc= input$lfc)
   })
   
   comparisons<-reactive ({
    if (input$exp == "TMT"  | input$exp == "DIA") {
       temp<-capture.output(test_diff_customized(imputed_data(), type = "all"), type = "message")
       # temp<-capture.output(test_diff_customized(imputed_data(), type = "manual", 
       #                                           test = c("SampleTypeTumor"), design_formula = formula(~0+SampleType)),
       #                      type = "message")
       gsub(".*: ","",temp)
       ## Split conditions into character vector
       unlist(strsplit(temp,","))
       ## Remove leading and trailing spaces
       trimws(temp)
     } else if (input$exp == "LFQ") {
       temp<-capture.output(test_diff(imputed_data(),type='all'),type = "message")
       gsub(".*: ","",temp)
       ## Split conditions into character vector
       unlist(strsplit(temp,","))
       ## Remove leading and trailing spaces
       trimws(temp)
     }
   })


   ## Results plot inputs
   ## PCA Plot
   pca_input<-eventReactive(input$analyze ,{ 
     if(input$analyze==0 | !start_analysis()){
       return()
     }
     ID_col <- "label"
     if (num_total()<=500){
       if(length(levels(as.factor(colData(dep())$replicate))) <= 6){
         pca_plot<- plot_pca_plotly(dep(), n=num_total(), ID_col=ID_col, exp=input$exp)
       } else{
           pca_plot<- plot_pca_plotly(dep(), n=num_total(), indicate = "condition", ID_col=ID_col, exp=input$exp)
       }
     } else {
         if(length(levels(as.factor(colData(dep())$replicate))) <= 6){
           pca_plot<- plot_pca_plotly(dep(), ID_col=ID_col, exp=input$exp)
         }
         else{
           pca_plot<-plot_pca_plotly(dep(), indicate = "condition", ID_col=ID_col, exp=input$exp)
         }
     }
     return(pca_plot)
   })
   
   pca_static_input <- eventReactive(input$analyze ,{ 
     if(input$analyze==0 | !start_analysis()){
       return()
     }
     ID_col <- "label"
     if (num_total()<=500){
       if(length(levels(as.factor(colData(dep())$replicate))) <= 6){
         pca_plot<- plot_pca_customized(dep(), n=num_total(), ID_col=ID_col) + labs(title = "PCA Plot")
       } else{
         pca_plot<- plot_pca_customized(dep(), n=num_total(), indicate = "condition", ID_col=ID_col)  + labs(title = "PCA Plot")
       }
     } else {
       if(length(levels(as.factor(colData(dep())$replicate))) <= 6){
         pca_plot<-plot_pca_customized(dep(), ID_col=ID_col) + labs(title = "PCA Plot")
       }
       else{
         pca_plot <-plot_pca_customized(dep(), indicate = "condition", ID_col=ID_col)  + labs(title = "PCA Plot")
       }
     }
     return(pca_plot)
   })
   
   ### Heatmap for differentially expressed proteins
   heatmap_cluster<-eventReactive(input$analyze, { 
     if(input$analyze==0 | !start_analysis()){
       return()
     }
     heatmap_list <- get_cluster_heatmap(dep(),
                         type="centered", kmeans = F,
                         k=6, col_limit = 6,
                         indicate = "condition", exp=input$exp
                         )
     return(heatmap_list)
   })
   
   heatmap_input <- reactive({
     heatmap_list <- heatmap_cluster()
     heatmap_list[[1]]
   })
   
   ### Volcano Plot
    volcano_input <- reactive({
      if(!is.null(input$volcano_cntrst)) {
                    plot_volcano_new(dep(),
                    input$volcano_cntrst,
                    input$fontsize,
                    input$check_names,
                    input$p_adj)
    
      }
    })
    
    volcano_df<- reactive({
      if(!is.null(input$volcano_cntrst)) {
        get_volcano_df(dep(),
                         input$volcano_cntrst)
        
      }
    })

    
    volcano_input_selected<-reactive({
      if(!is.null(input$volcano_cntrst)){
        
        if (!is.null(input$contents_rows_selected)){
          proteins_selected<-data_result()[c(input$contents_rows_selected),]## get all rows selected
        } else if(!is.null(input$protein_brush)){
          proteins_selected <- data_result()[data_result()[["Gene Name"]] %in% protein_name_brush(), ] 
        }
        
       ## convert contrast to x and padj to y
       diff_proteins <- grep(paste("^",input$volcano_cntrst, "_log2", sep = ""),
                    colnames(proteins_selected))
        if(input$p_adj=="FALSE"){
       padj_proteins <- grep(paste("^",input$volcano_cntrst, "_p.val", sep = ""),
                                     colnames(proteins_selected))
        }
       else{
         padj_proteins <- grep(paste("^",input$volcano_cntrst, "_p.adj", sep = ""),
                               colnames(proteins_selected))
       }

       df_protein <- data.frame(x = proteins_selected[, diff_proteins],
                        y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
                        name = proteins_selected$`Gene Name`)
       #print(df_protein)
       p<-plot_volcano_new(dep(),
                    input$volcano_cntrst,
                    input$fontsize,
                    input$check_names,
                    input$p_adj)
       
       p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
         ggrepel::geom_text_repel(data = df_protein,
                                  aes(x, y, label = name),
                                  size = 4,
                                  box.padding = unit(0.1, 'lines'),
                                  point.padding = unit(0.1, 'lines'),
                                  segment.size = 0.5)## use the dataframe to plot points

       }
    })
    
    protein_input<-reactive({
      protein_selected  <- data_result()[input$contents_rows_selected,1]
      protein_selected <-as.character(protein_selected)
      if(length(levels(as.factor(colData(dep())$replicate))) <= 8){
        plot_protein(dep(), protein_selected, as.character(input$type), id="label")
      } else {
        protein_plot<-plot_protein(dep(), protein_selected, as.character(input$type), id="label")
        protein_plot + scale_color_brewer(palette = "Paired")
      }
    })
     
   ## QC Inputs
   norm_input <- reactive({
     if (input$exp == "TMT") {
       plot_normalization_customized(filtered_data(),
                                     normalised_data())
     } else if (input$exp == "DIA") {
       plot_normalization_DIA_customized(filtered_data(),
                                     normalised_data())
     } else if (input$exp == "LFQ") {
       plot_normalization_DIA_customized(filtered_data(),
                                     normalised_data())
     }
   })
   
   missval_input <- reactive({
     plot_missval_customized(filtered_data(), input$exp)
   })
   
   detect_input <- reactive({
     plot_detect(filtered_data())
   })
   
   imputation_input <- reactive({
     if (input$exp == "TMT") {
       plot_imputation_customized(filtered_data(), imputed_data())
     } else {
       plot_imputation_DIA_customized(filtered_data(), imputed_data())
     }
   })
   
   # p_hist_input <- reactive({
   #   plot_p_hist(dep())
   # })
   
   numbers_input <- reactive({
     if (input$exp == "TMT") {
       plot_numbers_by_plex_set(filtered_data())
     } else {
       plot_numbers_customized(filtered_data(), exp=input$exp)
     }
   })
   
   coverage_input <- reactive({
     plot_coverage(filtered_data())
   })
   
   correlation_input<-reactive({
     plot_cor_customized(dep(), significant=FALSE, indicate="condition", exp=input$exp)
   })
   
   cvs_input<-reactive({
     plot_cvs(dep(), id="label", check.names=F)
   })
   
   num_total<-reactive({
     dep() %>%
       nrow()
   }) 
   
   ## Enrichment inputs
    
   go_input<-eventReactive(input$go_analysis,{
     withProgress(message = 'Gene ontology enrichment is in progress',
                  detail = 'Please wait for a while', value = 0, {
                    for (i in 1:15) {
                      incProgress(1/15)
                      Sys.sleep(0.25)
                    }
                  })
     
    if(!is.null(input$contrast)){
      go_results <- test_ora_mod(dep(), databases = as.character(input$go_database), contrasts = TRUE,
                                 direction = input$go_direction, log2_threshold = input$lfc, alpha = input$p)
      null_enrichment_test(go_results)
      plot_go <- plot_enrichment(go_results, number = 10, alpha = 0.05, contrasts = input$contrast,
                                 databases = as.character(input$go_database), nrow = 2, term_size = 8)
      go_list<-list("go_result"=go_results, "plot_go"=plot_go, "go_database"=as.character(input$go_database))
      return(go_list)
    }
   })

   pathway_input<-eventReactive(input$pathway_analysis,{
     progress_indicator("Pathway Analysis is running....")
     pathway_results <- test_ora_mod(dep(), databases=as.character(input$pathway_database), contrasts = TRUE,
                                     direction = input$pathway_direction, log2_threshold = input$lfc, alpha = input$p)
     null_enrichment_test(pathway_results)
     plot_pathway <-plot_enrichment(pathway_results, number = 10, alpha = 0.05, contrasts =input$contrast_1,
                                    databases=as.character(input$pathway_database), nrow = 2, term_size = 8)
     pathway_list<-list("pa_result"=pathway_results,
                        "pathway_database"=as.character(input$pathway_database),
                        "plot_pathway"=plot_pathway)
     return(pathway_list)
   })

   
#### Interactive UI
   output$significantBox <- renderUI({
     num_total <- assay(processed_data()) %>%
       nrow()
     num_signif <- dep() %>%
       .[replace_na(SummarizedExperiment::rowData(.)$significant, F), ] %>%
       nrow()
     frac <- num_signif / num_total
     
       info_box <- 		infoBox("Significant proteins",
                             paste0(num_signif,
                                    " out of ",
                                    num_total),
                             paste0(signif(frac * 100, digits = 3),
                                    "% of proteins differentially expressed across all conditions"),
                             icon = icon("stats", lib = "glyphicon"),
                             color = "olive", width=12)
     
     return(info_box)
   })

  ##### Get results dataframe from Summarizedexperiment object
   data_result<-eventReactive(start_analysis(),{
      get_results_proteins(dep(), input$exp)
      #get_results(dep())
    })

  #### Data table
   output$contents <- DT::renderDataTable({
     df<- data_result()
     return(df)
     },
     options = list(scrollX = TRUE,
                    autoWidth=TRUE,
                    columnDefs= list(list(width = '400px', targets = c(-1)))))
  
  ## Deselect all rows button
  proxy <- dataTableProxy("contents")
  
  observeEvent(input$clear,{
    proxy %>% selectRows(NULL)
  })
  
  observeEvent(input$original,{
    output$contents <- DT::renderDataTable({
      df<- data_result()
      return(df)
    },
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
  })

  protein_name_brush<- reactive({
    protein_tmp<-brushedPoints(volcano_df(), input$protein_brush,
                               xvar = "diff", yvar = "p_values")
    return(protein_tmp$name)
  })
  
  ## Select rows dynamically
  brush <- NULL
  makeReactiveBinding("brush")
  observeEvent(input$protein_brush,{
    data <- data_result()
    proxy %>% selectRows(which(data[["Gene Name"]] %in% protein_name_brush()))
  })
 
 observeEvent(input$resetPlot,{
   session$resetBrush("protein_brush")
   brush <<- NULL
   
   proxy %>% selectRows(NULL)
 })

 
  ## Render Result Plots
  output$pca_plot<-renderPlotly({
    pca_input()
  })
  
  output$heatmap<-renderPlot({
    withProgress(message = 'Heatmap rendering is in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
   heatmap_input()
  })
 
  output$volcano <- renderPlot({
    withProgress(message = 'Volcano Plot calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    if(is.null(input$contents_rows_selected)){
      volcano_input()
    } else if(!is.null(input$volcano_cntrst)){
      volcano_input_selected()
    }# else close
  })
  
  output$protein_plot<-renderPlotly({
    if(!is.null(input$contents_rows_selected)){
      protein_input()
    }
  })
 
 
  ### QC Outputs
  output$sample_corr <-renderPlot({
    correlation_input()
  })
  
  output$sample_cvs <- renderPlot({
    cvs_input()
  })
  
  output$norm <- renderPlot({
    norm_input()
  })
  
  output$missval <- renderPlot({
    missval_input()
  })
  
  output$detect <- renderPlot({
    detect_input()
  })
  
  output$imputation <- renderPlot({
    imputation_input()
  })
  
  output$p_hist <- renderPlot({
    p_hist_input()
  })
  
  output$numbers <- renderPlot({
    numbers_input()
  })
  
  output$coverage <- renderPlot({
    coverage_input()
  })
  
  ## Enrichment Outputs
  output$spinner_go <- renderUI({
    req(input$go_analysis)
    shinycssloaders::withSpinner(plotOutput("go_enrichment"), color = "#3c8dbc")
  })
  
  output$go_enrichment<-renderPlot({
    Sys.sleep(2)
    go_input()$plot_go
  })
  
  output$spinner_pa <- renderUI({
    req(input$pathway_analysis)
    shinycssloaders::withSpinner(plotOutput("pathway_enrichment"), color = "#3c8dbc")
  })
  
  output$pathway_enrichment<-renderPlot({
    Sys.sleep(2)
    return(pathway_input()$plot_pathway)
  })
  
  ##### Download Functions
  # example data
  output$lfq_example <- downloadHandler(
    filename="combined_protein.tsv",  # desired file name on client 
    content=function(con) {
      file.copy("./data/LFQ_datasets/ubiquitin/combined_protein.tsv", con)
    }
  )
  
  output$lfq_annotation <- downloadHandler(
    filename="experiment_annotation.tsv",  # desired file name on client 
    content=function(con) {
      file.copy("./data/LFQ_datasets/ubiquitin/experiment_annotation.tsv", con)
    }
  )
  
  datasetInput <- reactive({
    switch(input$dataset,
           "Results" = get_results_proteins(dep(), input$exp),
           "Original_matrix"= unimputed_table(),
           # "significant_proteins" = get_results(dep()) %>%
           #   filter(significant) %>%
           #   select(-significant),
           "Imputed_matrix" = imputed_table(),
           "Full_dataset" = get_df_wide(dep()))
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$dataset, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(datasetInput(),
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  ### === Cluster Download ==== ####
  
  individual_cluster <- reactive({
      cluster_number <- input$cluster_number
      cluster_all <- heatmap_input()[[2]]
      df <- data_result()[cluster_all[[cluster_number]],]
      return(df)
  })
  
  # output$text1 <- renderPrint({
  #   paste(individual_cluster())
  # })
  
  output$downloadCluster <- downloadHandler(
    filename = function() { paste("Cluster_info_",input$cluster_number, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(individual_cluster(),
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  output$downloadVolcano <- downloadHandler(
    filename = function() {
      paste0("Volcano_", input$volcano_cntrst, ".pdf")
    },
    content = function(file) {
      pdf(file)
      if(is.null(input$protein_brush)){
        print(volcano_input())
        dev.off()
      }
      else{
      observeEvent(input$protein_brush,{
        print(p)
      })
      print(volcano_input_selected())
      dev.off()
      }
    }
  )
  
  
  ## Protein plot download
  output$downloadProtein <- downloadHandler(
    filename = function() {
      paste0(input$type,".pdf")
    },
    content = function(file) {
      pdf(file)
      print(protein_input())
      dev.off()
    }
  )
  
  ###### ==== DOWNLOAD GO TABLE ==== ####
  output$downloadGO <- downloadHandler(
    filename = function() { paste("GO_enrichment_",input$go_database, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(go_input()$go_result,
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  ###### ==== DOWNLOAD PATHWAY TABLE ==== ####
  output$downloadPA <- downloadHandler(
    filename = function() { paste("Pathway_enrichment_",input$pathway_database, ".csv", sep = "") }, 
    ## use = instead of <-
    content = function(file) {
      write.table(pathway_input()$pa_result,
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
output$download_hm_svg<-downloadHandler(
  filename = function() { "heatmap.svg" }, 
  ## use = instead of <-
  content = function(file) {
    heatmap_plot <- heatmap_input()
    svg(file)
    print(heatmap_plot)
    dev.off()
  }
)
  
#####===== Download Report =====#####
  output$downloadReport <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "FragPipe-Analyst_report.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "LFQ_report.Rmd")
      file.copy(paste0("./reports/", input$exp, "_report.Rmd"), tempReport, overwrite = TRUE)
      
      sig_proteins<-dep() %>%
        .[replace_na(SummarizedExperiment::rowData(.)$significant, FALSE), ] %>%
        nrow()
      
      tested_contrasts<- gsub("_p.adj", "", 
                              colnames(SummarizedExperiment::rowData(dep()))[grep("p.adj", 
                              colnames(SummarizedExperiment::rowData(dep())))])
      pg_width<- ncol(imputed_data()) / 2.5
      # Set up parameters to pass to Rmd document
      params <- list(data = processed_data,
                     alpha = input$p,
                     lfc = input$lfc,
                     num_signif= sig_proteins,
                     pg_width = pg_width,
                     tested_contrasts= tested_contrasts,
                     numbers_input= numbers_input,
                     detect_input = detect_input,
                     imputation_input = imputation_input,
                     missval_input = missval_input,
                     # p_hist_input = p_hist_input,
                     pca_input = pca_static_input,
                     coverage_input= coverage_input,
                     correlation_input =correlation_input,
                     heatmap_input = heatmap_input,
                     cvs_input = cvs_input,
                     volcano_input = volcano_input,
                     dep = dep
                     )
      
      # Knit the document, passing in the `params` list
       rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
###### ==== DOWNLOAD QC plots svg ==== ####

output$download_pca_svg<-downloadHandler(
  filename = function() { "PCA_plot.svg" }, 
  content = function(file) {
    svg(file)
    print(pca_input())
    dev.off()
  }
)

output$download_corr_svg<-downloadHandler(
  filename = function() { "Correlation_plot.svg" }, 
  content = function(file) {
    svg(file)
    print(correlation_input())
    dev.off()
  }
)

output$download_cvs_svg<-downloadHandler(
  filename = function() { "Sample_CV.svg" }, 
  content = function(file) {
    svg(file)
    print(cvs_input())
    dev.off()
  }
)

output$download_num_svg<-downloadHandler(
  filename = function() { "Proteins_plot.svg" }, 
  content = function(file) {
    svg(file)
    print(numbers_input())
    dev.off()
  }
)

output$download_cov_svg<-downloadHandler(
  filename = function() { "Coverage_plot.svg" }, 
  content = function(file) {
    svg(file)
    print(coverage_input())
    dev.off()
  }
)

output$download_norm_svg<-downloadHandler(
  filename = function() { "Normalization_plot.svg" }, 
  content = function(file) {
    svg(file)
    print(norm_input())
    dev.off()
  }
)

output$download_missval_svg<-downloadHandler(
  filename = function() { "Missing_value_heatmap.svg" }, 
  content = function(file) {
    svg(file)
    print(missval_input())
    dev.off()
  }
)

output$download_imp_svg<-downloadHandler(
  filename = function() { "Imputation_plot.svg" }, 
  content = function(file) {
    svg(file)
    print(imputation_input())
    dev.off()
  }
)

  #### Occurrence page logic ####
  data_attendance<-eventReactive(start_analysis(),{
    conditions <- condition_list()
    
    
    df <- as.data.frame(assay(processed_data()), check.names=F)
    sample_cols <- colnames(df)
   
    if (input$exp == "LFQ"){
      # MaxQuant Protein.names is Description in FragPipe output
      # MaxQuant Gene.names is Gene in FragPipe output
      df$Gene <- rowData(processed_data())$Gene
      df$Description <- rowData(processed_data())$Description
      df$Combined.Total.Peptides <- as.data.frame(rowData(processed_data())["Combined Total Peptides"])
      
      if ("" %in% df$Gene){
        df$Gene[df["Gene"]==""] <- "NoGeneNameAvailable"}
      if ("" %in% df$Description){
        df$Description[df["Description"]==""] <- "NoProteinDescriptionAvailable"}
    } else { # DIA
      # "Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description" "name"
      df$Gene <- rowData(processed_data())$Genes
      df$Description <- rowData(processed_data())$First.Protein.Description
    }
    

    # filter if all intensity are NA
    df <- df[rowSums(!is.na(df[,sample_cols])) != 0,]

    # print(conditions)
    # print(colnames(df))
    if (input$exp == "LFQ"){
      for (i in 1:length(conditions)) {
        condition <- conditions[i]
        temp <- as.data.frame(colData(processed_data()))
        temp <- temp[temp$condition==condition,]
        selected_cols <- rownames(temp)
        df[paste0("#Occurences", sep = "_", condition)] <- rowSums(!is.na(df[,selected_cols, drop=F]))
        df <- dplyr::relocate(df, paste0("#Occurences",sep = "_", condition))
        cols <- grep(paste0(condition, "$"),colnames(df))
        if (!is.null(input[[paste0("",condition)]])){
          df <- df %>%
            dplyr::filter(df[[cols]] >=input[[paste0("",condition)]][1] & df[[cols]] <=input[[paste0("",condition)]][2])
        }
      }
      df <- dplyr::relocate(df, "Gene", "Description", "Combined.Total.Peptides")
    } else { # DIA
      for (i in 1:length(conditions)) {
        condition <- conditions[i]
        temp <- as.data.frame(colData(processed_data()))
        temp <- temp[temp$condition==condition,]
        selected_cols <- temp$label
        df[paste0("#Occurences", sep = "_", condition)] <- rowSums(!is.na(df[,selected_cols, drop=F]))
        df <- dplyr::relocate(df, paste0("#Occurences",sep = "_", condition))
        if (!is.null(input[[paste0("",condition)]])){
          df <- df %>%
            dplyr::filter(df[[paste("#Occurences", condition, sep="_")]] >=input[[paste0("",condition)]][1] &
                            df[[paste("#Occurences", condition, sep="_")]] <=input[[paste0("",condition)]][2])
        }
      }
      df <- dplyr::relocate(df, "Gene", "Description")
    }
    return(df)
  })
  
  data_attendance_filtered <- reactive({
    data <- data_attendance()
    if (!is.null(input$filtered_condition_fragpipe)) {
      # if(("Reverse" %in% colnames(filtered_data)) & ('Reverse sequences' %in% input$filtered_condition_maxquant)){
      #   filtered_data<-dplyr::filter(filtered_data,Reverse!="+")
      # }
      # else{filtered_data <-filtered_data}
      # if(("Potential.contaminant" %in% colnames(filtered_data)) & ('Potential contaminants' %in% input$filtered_condition_maxquant)){
      #   filtered_data<-dplyr::filter(filtered_data,Potential.contaminant!="+")
      # }
      # else{filtered_data <-filtered_data}
  
      # if(("Only.identified.by.site" %in% colnames(filtered_data)) & ('"Only identified by site" Protein' %in% input$filtered_condition_maxquant)){
      #   filtered_data<-dplyr::filter(filtered_data,Only.identified.by.site!="+")
      # }
      # else{filtered_data <-filtered_data}
      # if(("Razor...unique.peptides" %in% colnames(filtered_data)) & ('Proteins with peptiedes < 2' %in% input$filtered_condition_maxquant)){
      #   filtered_data<-dplyr::filter(filtered_data,as.numeric(Razor...unique.peptides)>=2)
      # }
      # else{filtered_data <-filtered_data}
      if (input$exp == "LFQ") {
        if(('Proteins with more than two peptides' %in% input$filtered_condition_fragpipe)){
          data <-dplyr::filter(data,Combined.Total.Peptides>=2)
        }
      }
    }
    return(data)
  })
  
  #### Data table
  output$contents_occ <- DT::renderDataTable({
    df<- data_attendance_filtered()
    return(df)},
    options = list(scrollX = TRUE,
                   scroller = TRUE,
                   autoWidth=TRUE,
                   # columnDefs= list(list(width = "10%", targets = c(1)),
                   #                  list(width = "400px", targets = grep("Protein.names", names(df)))
                   #                  )
                   columnDefs= list(list(width = '400px', targets = c(-1)))
    )
  )
  
  make_sliderInput <- function(n= 1){
    exp_design_input <- exp_design_input()
    conditions <- exp_design_input$condition %>% unique()
    
    if (input$exp == "LFQ"){
      sliderInput(paste0("",conditions[n]),
                  label=paste0("",conditions[n]),
                  min = min(0),
                  max = max(exp_design_input$replicate),
                  value = c(0, max(exp_design_input$replicate)),
                  step = 1)
    } else {
      temp <- exp_design_input[exp_design_input$condition == conditions[n],]
      sliderInput(paste0("",conditions[n]),
                  label=paste0("",conditions[n]),
                  min = min(0),
                  max = max(nrow(temp)),
                  value = c(0, max(nrow(temp))),
                  step = 1)
    }
  }
  
  slider_bars <- reactive({
    exp_design <- exp_design_input()
    lapply(X = 1:length(unique(exp_design$condition)), FUN = make_sliderInput)
  })
  
  output$sidebar <- renderUI({
    tagList(slider_bars())
  })
  
  output$download_attendance <- downloadHandler("Occurrences_results_table.csv",
                                                content = function(file){
                                                  write.table(data_attendance_filtered(),  
                                                              file,
                                                              col.names = TRUE,
                                                              row.names = FALSE,
                                                              sep =",")
                                                },
                                                contentType = "text/csv")
  
  ## Venn plot
  condition_list <- reactive({
    if(!is.null(exp_design_input())){
      conditions <- exp_design_input()$condition %>% unique()
      return(conditions)
    }
  })
  
  observeEvent(input$analyze, {
    if (length(condition_list()) == 1){
      shinyjs::hide(id = "con_2")
      shinyjs::hide(id = "con_3")
    } 
    if (length(condition_list()) == 2){
      shinyjs::hide(id = "con_3")
    } 
  })
  
  output$condition_1 <- renderUI({
    if (!is.null(condition_list())){
      selectizeInput("condition_1",
                     "Condition 1",
                     choices = condition_list(),
                     selected = condition_list()[1])
    }
  })
  
  output$condition_2 <- renderUI({
    if (!is.null(condition_list()) & length(condition_list()) > 1){
      selectizeInput("condition_2",
                     "Condition 2",
                     # choices = condition_list(),
                     choices = condition_list()[condition_list() != input$condition_1],
                     selected = condition_list()[2])
    }
  })
  
  output$condition_3 <- renderUI({
    if (!is.null(condition_list())  & length(condition_list()) > 2){
      selectizeInput("condition_3",
                     "Condition 3",
                     # choices = condition_list(),
                     choices = c("NONE", condition_list()[condition_list() != input$condition_1 & condition_list() != input$condition_2]),
                     selected = condition_list()[3])
    }
  })
  
  venn_plot_input <- reactive({
    df<- data_attendance_filtered()
    
    if(length(condition_list()) < 2){
      stop(safeError("Venn plot should contain at least two sets"))
    } else if(length(condition_list()) < 3){
      set1 <- df$Gene[df[grep(paste0("#Occurences",sep = "_",input$condition_1),colnames(df))] != 0]
      set2 <- df$Gene[df[grep(paste0("#Occurences",sep = "_",input$condition_2),colnames(df))] != 0]
      x <- list(set1,set2)
      names(x) <- c("Condition 1", "Condition 2")
    } else {
      set1 <- df$Gene[df[grep(paste0("#Occurences",sep = "_",input$condition_1),colnames(df))] != 0]
      set2 <- df$Gene[df[grep(paste0("#Occurences",sep = "_",input$condition_2),colnames(df))] != 0]
      set3 <- df$Gene[df[grep(paste0("#Occurences",sep = "_",input$condition_3),colnames(df))] != 0]
      x <- list(set1,set2,set3)
      names(x) <- c("Condition 1", "Condition 2", "Condition 3")
      if (!is.null(input$condition_3)){
        if (input$condition_3 == "NONE"){
          x <- list(set1,set2)
          names(x) <- c("Condition 1", "Condition 2")
        }
      }
    }
    ggVennDiagram::ggVennDiagram(x,label_alpha = 0) +
      scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
  })
  
  output$venn_plot <- renderPlot({
    if (!is.null(input$condition_1) & !is.null(input$condition_2)){
        venn_plot_input()
    }
  })
  
  output$download_venn_svg<-downloadHandler(
    filename = function() { "venn_plot.svg" }, 
    content = function(file) {
      svg(file)
      print(venn_plot_input())
      dev.off()
    }
  )
 
 #### Demo logic ========== #############
 
 ####======= Render Functions
 
 # output$volcano_cntrst_dm <- renderUI({
 #   if (!is.null(comparisons_dm())) {
 #     df <- SummarizedExperiment::rowData(dep_dm())
 #     cols <- grep("_significant$",colnames(df))
 #     selectizeInput("volcano_cntrst_dm",
 #                    "Comparison",
 #                    choices = gsub("_significant", "", colnames(df)[cols]))
 #   }
 # })
 
 ## comparisons for demo
 # output$contrast_dm <- renderUI({
 #   if (!is.null(comparisons_dm())) {
 #     df <- SummarizedExperiment::rowData(dep_dm())
 #     cols <- grep("_significant$",colnames(df))
 #     selectizeInput("contrast_dm",
 #                    "Comparison",
 #                    choices = gsub("_significant", "", colnames(df)[cols]))
 #   }
 # })
 
 # output$contrast_dm_1 <- renderUI({
 #   if (!is.null(comparisons_dm())) {
 #     df <- SummarizedExperiment::rowData(dep_dm())
 #     cols <- grep("_significant$",colnames(df))
 #     selectizeInput("contrast_dm_1",
 #                    "Comparison",
 #                    choices = gsub("_significant", "", colnames(df)[cols]))
 #   }
 # })
 
 # output$downloadTable_dm <- renderUI({
 #   if(!is.null(dep_dm())){
 #     selectizeInput("dataset_dm",
 #                    "Download data table" ,
 #                    c("Results",
 #                      "Full dataset"))
 #   }
 # })
 
 # output$downloadButton_dm <- renderUI({
 #   if(!is.null(dep_dm())){
 #     downloadButton('downloadData_dm', 'Save')
 #   }
 # })
 
#  output$downloadreport_dm <- renderUI({
#    if(!is.null(dep_dm())){
#      downloadButton('downloadReport_dm', 'Download Report')
#    }
#  })
#  
#  output$downloadPlots_dm <- renderUI({
#    if(!is.null(dep_dm())){
#      downloadButton('downloadPlots1_dm', 'Download Plots')
#    }
#  })
#  
#  env_dm<-reactive({
#    LoadToEnvironment("data/demo_data.RData", env = globalenv())
#  })
#  
#  
#  
#  dep_dm<-reactive({
#    env_dm()[["dep"]]
#  })
#  
#  
# comparisons_dm<-reactive({
#   comparisons<-gsub("_p.adj", "", colnames(SummarizedExperiment::rowData(dep_dm()))[grep("p.adj", 
#                           colnames(SummarizedExperiment::rowData(dep_dm())))])
# })
 ## Results plot inputs for demo
 
 ## PCA Plot for demo

# pca_label_dm<-reactive({
#  pca_lable<-levels(as.factor(colData(dep_dm())$replicate))
# print(pca_label)
# })
# 
#  pca_input_dm<-reactive({
#    if (num_total_dm()<=500){
#      if(length(levels(as.factor(colData(dep_dm())$replicate))) <= 6){
#        pca_plot<- DEP::plot_pca(dep_dm(), n=num_total_dm(), point_size = 4)
#        pca_plot<-pca_plot + labs(title = "PCA plot")
#        return(pca_plot)
#      }
#      else{
#        pca_plot<-DEP::plot_pca(dep_dm(), n=num_total_dm(), point_size = 4, indicate = "condition")
#        pca_plot<-pca_plot + labs(title = "PCA plot")
#        return(pca_plot)
#      }
#    }
#    else{
#      if(length(levels(as.factor(colData(dep_dm())$replicate))) <= 6){
#        pca_plot<-DEP::plot_pca(dep_dm(), point_size = 4)
#        pca_plot<-pca_plot + labs(title = "PCA plot")
#        return(pca_plot)
#      }else{
#        pca_label<-SummarizedExperiment::colData(dep_dm())$replicate
#        pca_plot<-DEP::plot_pca(dep_dm(), point_size = 4, indicate = "condition")
#        #pca_plot<-pca_plot + geom_point()
#        pca_plot<-pca_plot + ggrepel::geom_text_repel(aes(label=factor(rowname)),
#                                            size = 4,
#                                            box.padding = unit(0.1, 'lines'),
#                                            point.padding = unit(0.1, 'lines'),
#                                            segment.size = 0.5)
#        pca_plot<-pca_plot + labs(title = "PCA plot")
# 	  return(pca_plot)
#      }
#    }
#    
#  })
 
 ### Heatmap Differentially expressed proteins for demo
 # heatmap_input_dm<-reactive({ 
 #   get_cluster_heatmap(dep_dm(),
 #                       type="centered",kmeans = TRUE,
 #                       k=6, col_limit = 6,
 #                       indicate = "condition"
 #   )
 # })
 
 ### Volcano Plot for demo
 # volcano_input_dm <- reactive({
 #   if(!is.null(input$volcano_cntrst_dm)) {
 #     plot_volcano_new(dep_dm(),
 #                      input$volcano_cntrst_dm,
 #                      input$fontsize_dm,
 #                      input$check_names_dm,
 #                      input$p_adj_dm)
 #   }
 # })
 # 
 # volcano_df_dm<- reactive({
 #   if(!is.null(input$volcano_cntrst_dm)) {
 #     get_volcano_df(dep_dm(),
 #                    input$volcano_cntrst_dm)
 #     
 #   }
 # })
 # 
 # 
 # volcano_input_selected_dm<-reactive({
 #   if(!is.null(input$volcano_cntrst_dm)){
 #     if (!is.null(input$contents_dm_rows_selected)){
 #       proteins_selected<-data_result_dm()[c(input$contents_dm_rows_selected),]## get all rows selected
 #     }
 #     else if(!is.null(input$protein_brush_dm)){
 #       proteins_selected<-data_result_dm()[data_result_dm()[["Gene Name"]] %in% protein_name_brush_dm(), ] 
 #     }
 #     ## convert contrast to x and padj to y
 #     diff_proteins <- grep(paste(input$volcano_cntrst_dm, "_log2", sep = ""),
 #                           colnames(proteins_selected))
 #     if(input$p_adj_dm=="FALSE"){
 #       padj_proteins <- grep(paste(input$volcano_cntrst_dm, "_p.val", sep = ""),
 #                             colnames(proteins_selected))
 #     }
 #     else{
 #       padj_proteins <- grep(paste(input$volcano_cntrst_dm, "_p.adj", sep = ""),
 #                             colnames(proteins_selected))
 #     }
 #     
 #     df_protein <- data.frame(x = proteins_selected[, diff_proteins],
 #                              y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
 #                              name = proteins_selected$`Gene Name`)
 #     #print(df_protein)
 #     p<-plot_volcano_new(dep_dm(),
 #                         input$volcano_cntrst_dm,
 #                         input$fontsize_dm,
 #                         input$check_names_dm,
 #                         input$p_adj_dm)
 #     
 #     p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
 #       ggrepel::geom_text_repel(data = df_protein,
 #                                aes(x, y, label = name),
 #                                size = 4,
 #                                box.padding = unit(0.1, 'lines'),
 #                                point.padding = unit(0.1, 'lines'),
 #                                segment.size = 0.5)## use the dataframe to plot points
 #     
 #   }
 # })
 # 
 # protein_input_dm<-reactive({ 
 #   
 #   protein_selected  <- data_result_dm()[input$contents_dm_rows_selected,1]
 #  
 #   #protein<-row_selected$name
 #   if(length(levels(as.factor(colData(dep_dm())$replicate))) <= 8){
 #     plot_protein(dep_dm(), protein_selected, input$type_dm)
 #   }
 #   else{
 #     protein_plot<-plot_protein(dep_dm(), protein_selected, input$type_dm)
 #     protein_plot + scale_color_brewer(palette = "Paired")
 #    }
 # })
 
 ## Get processed data for demo
 
 # processed_data_dm<-reactive({
 #   env_dm()[["data_filter"]]
 # })
 # normalised_data_dm<-reactive({
 #   DEP::normalize_vsn(processed_data_dm())
 # })
 # 
 # 
 # imputed_data_dm<-reactive({
 #   DEP::impute(processed_data_dm(),input$imputation)
 # })
 # 
 # 
 # diff_all_dm<-reactive({
 #   if (input$exp == "LFQ") {
 #     test_diff(imputed_data_dm(),type = 'all')
 #   } else if (input$exp == "TMT") {
 #     # test_diff_customized(imputed_data_dm(), type = "manual", 
 #     #                      test = c("SampleTypeTumor"), design_formula = formula(~0+SampleType))
 #     test_diff_customized(imputed_data_dm(), type = "all")
 #   }
 # })
	
 ## QC Inputs for demo
 # norm_input_dm <- reactive({
 #   plot_normalization(processed_data_dm(),
 #                      normalised_data_dm())
 # })
 # 
 # missval_input_dm <- reactive({
 #   plot_missval(processed_data_dm())
 # })
 # 
 # detect_input_dm <- reactive({
 #   plot_detect(processed_data_dm())
 # })
 # 
 # imputation_input_dm <- reactive({
 #   plot_imputation(normalised_data_dm(),
 #                   diff_all_dm())
 # })
 # 
 # p_hist_input_dm <- reactive({
 #   plot_p_hist(dep_dm())
 # })
 # 
 # numbers_input_dm <- reactive({
 #   if (input$exp == "LFQ"){
 #     plot_numbers(normalised_data_dm())
 #   } else if (input$exp == "TMT") {
 #     plot_numbers(normalised_data_dm())
 #   }
 # })
 # 
 # coverage_input_dm <- reactive({
 #   plot_coverage(normalised_data_dm())
 # })
 # 
 # correlation_input_dm<-reactive({
 #   plot_cor(dep_dm(), significant=TRUE, indicate="condition")
 # })
 # 
 # cvs_input_dm<-reactive({
 #   if (input$exp == "LFQ") {
 #     id <- "ID"
 #   } else if (input$exp == "TMT") {
 #     id <- "label"
 #   }
 #   plot_cvs(dep_dm(), id)
 # })
 # 
 # num_total_dm<-reactive({
 #   dep_dm() %>%
 #     nrow()
 # }) 
 
 ## Enrichment inputs for demo
 # go_input_dm<-eventReactive(input$go_analysis_dm,{
 #   withProgress(message = 'Gene ontology enrichment is in progress',
 #                detail = 'Please wait for a while', value = 0, {
 #                  for (i in 1:15) {
 #                    incProgress(1/15)
 #                    Sys.sleep(0.25)
 #                  }
 #                })
 #   
 #   if(!is.null(input$contrast_dm)){
 #     enrichment_output_test(dep_dm(), as.character(input$go_database_dm))
 #     go_results<- test_ora_mod(dep_dm(), databases = as.character(input$go_database_dm), contrasts = TRUE)
 #     plot_go<- plot_enrichment(go_results, number = 10, alpha = 0.05, contrasts =input$contrast_dm,
 #                               databases = as.character(input$go_database_dm), nrow = 2, term_size = 8) + 
 #       aes(stringr::str_wrap(Term, 60)) +
 #       xlab(NULL)
 #     go_list<-list("go_result"=go_results, "plot_go"=plot_go)
 #     return(go_list)
 #   }
 # })
 # 
 # pathway_input_dm<-eventReactive(input$pathway_analysis_dm,{
 #   withProgress(message = 'Pathway enrichment is in progress',
 #                detail = 'Please wait for a while', value = 0, {
 #                  for (i in 1:15) {
 #                    incProgress(1/15)
 #                    Sys.sleep(0.25)
 #                  }
 #                })
 #   enrichment_output_test(dep_dm(), as.character(input$pathway_database_dm))
 #   pathway_results<- test_ora_mod(dep_dm(), databases=as.character(input$pathway_database_dm), contrasts = TRUE)
 #   plot_pathway<-plot_enrichment(pathway_results, number = 10, alpha = 0.05, contrasts =input$contrast_dm_1,
 #                                 databases=as.character(input$pathway_database_dm), nrow = 3, term_size = 8) + 
 #     aes(stringr::str_wrap(Term, 30)) +
 #     xlab(NULL)
 #   pathway_list<-list("pa_result"=pathway_results, "plot_pa"=plot_pathway)
 #   return(pathway_list)
 #   
 # })
 
 
 #### Interactive UI for demo
 # output$significantBox_dm <- renderInfoBox({
 #   num_total <- dep_dm() %>%
 #     nrow()
 #   num_signif <- dep_dm() %>%
 #     .[SummarizedExperiment::rowData(.)$significant, ] %>%
 #     nrow()
 #   frac <- num_signif / num_total
 #     
 #     info_box <- 		infoBox("Significant proteins",
 #                           paste0(num_signif,
 #                                  " out of ",
 #                                  num_total),
 #                           paste0(signif(frac * 100, digits = 3),
 #                                  "% of proteins differentially expressed across all conditions"),
 #                           icon = icon("stats", lib = "glyphicon"),
 #                           color = "olive",
 #                           # fill = TRUE,
 #                           width = 4)
 #     
 #     return(info_box)
 #   })
   
 
 ##### Get results dataframe from Summarizedexperiment object for demo
 # data_result_dm<-reactive({
 #   get_results_proteins(dep_dm())
 #   #get_results(dep())
 # })
 
 #### Data table
#  output$contents_dm <- DT::renderDataTable({
#    df<- data_result_dm()
#    return(df)
#  },
#  options = list(scrollX = TRUE,
# 	autoWidth=TRUE,
#                 columnDefs= list(list(width = '400px', targets = c(-1))))
#  )
 
 ## Deselect all rows button for demo
 # proxy <- dataTableProxy("contents_dm")
 
#  observeEvent(input$clear_dm,{
#    proxy %>% selectRows(NULL)
#  })
#  
#  observeEvent(input$original_dm,{
#    output$contents_dm <- DT::renderDataTable({
#      df<- data_result_dm()
#      return(df)
#    },
#    options = list(scrollX = TRUE,
# 	autoWidth=TRUE,
#                 columnDefs= list(list(width = '400px', targets = c(-1))))
#    )
#  })
#  
#  protein_name_brush_dm<- reactive({
#    #protein_tmp<-nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
#    protein_tmp<-brushedPoints(volcano_df_dm(), input$protein_brush_dm, 
#                               xvar = "diff", yvar = "p_values")
#    protein_selected<-protein_tmp$name
#  }) 
#  
#  
#  ## Select rows dynamically
#  observeEvent(input$protein_brush_dm,{
#    output$contents_dm <- DT::renderDataTable({
#      df<- data_result_dm()[data_result_dm()[["Gene Name"]] %in% protein_name_brush_dm(), ]
#      return(df)
#    },
#    options = list(scrollX= TRUE,
# autoWidth=TRUE,
#                 columnDefs= list(list(width = '400px', targets = c(-1))))
#    )
#    
#    proteins_selected<-data_result_dm()[data_result_dm()[["Gene Name"]] %in% protein_name_brush_dm(), ] #
#    # get all rows selected
#    ## convert contrast to x and padj to y
#    diff_proteins <- grep(paste(input$volcano_cntrst_dm, "_log2", sep = ""),
#                          colnames(proteins_selected))
#    if(input$p_adj=="FALSE"){
#      padj_proteins <- grep(paste(input$volcano_cntrst_dm, "_p.val", sep = ""),
#                            colnames(proteins_selected))
#    }
#    else{
#      padj_proteins <- grep(paste(input$volcano_cntrst_dm, "_p.adj", sep = ""),
#                            colnames(proteins_selected))
#    }
#    df_protein <- data.frame(x = proteins_selected[, diff_proteins],
#                             y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
#                             name = proteins_selected$`Gene Name`)
#    #print(df_protein)
#    
#    p<-plot_volcano_new(dep_dm(),
#                        input$volcano_cntrst_dm,
#                        input$fontsize_dm,
#                        input$check_names_dm,
#                        input$p_adj_dm)
#    
#    p<- p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
#      ggrepel::geom_text_repel(data = df_protein,
#                               aes(x, y, label = name),
#                               size = 4,
#                               box.padding = unit(0.1, 'lines'),
#                               point.padding = unit(0.1, 'lines'),
#                               segment.size = 0.5)
#    
#    output$volcano_dm <- renderPlot({
#      withProgress(message = 'Volcano Plot calculations are in progress',
#                   detail = 'Please wait for a while', value = 0, {
#                     for (i in 1:15) {
#                       incProgress(1/15)
#                       Sys.sleep(0.25)
#                     }
#                   })
#      p
#    })
#    return(p)
#  })
#  
#  observeEvent(input$resetPlot_dm,{
#    session$resetBrush("protein_brush_dm")
#    brush <<- NULL
#    
#    output$contents_dm <- DT::renderDataTable({
#      df<- data_result_dm()
#      return(df)
#    },
#    options = list(scrollX = TRUE,
#                   autoWidth=TRUE,
#                   columnDefs= list(list(width = '400px', targets = c(-1))))
#    )
#  })
#  
#  ## Render Result Plots
#  output$pca_plot_dm<-renderPlot({
#    pca_input_dm()
#  })
#  output$heatmap_dm<-renderPlot({
#    withProgress(message = 'Heatmap rendering is in progress',
#                 detail = 'Please wait for a while', value = 0, {
#                   for (i in 1:15) {
#                     incProgress(1/15)
#                     Sys.sleep(0.25)
#                   }
#                 })
#    heatmap_input_dm()
#  })
#  
#  output$volcano_dm <- renderPlot({
#    withProgress(message = 'Volcano Plot calculations are in progress',
#                 detail = 'Please wait for a while', value = 0, {
#                   for (i in 1:15) {
#                     incProgress(1/15)
#                     Sys.sleep(0.25)
#                   }
#                 })
#    if(is.null(input$contents_dm_rows_selected)){
#      volcano_input_dm()
#    }
#    else if(!is.null(input$volcano_cntrst_dm)){
#      volcano_input_selected_dm()
#    } # else close
#  })
#  
#  output$protein_plot_dm<-renderPlot({
#    if(!is.null(input$contents_dm_rows_selected)){
#      protein_input_dm()
#    }
#  })
 
 
 ### QC Outputs for demo
 # output$sample_corr_dm <-renderPlot({
 #   correlation_input_dm()
 # })
 # 
 # output$sample_cvs_dm <- renderPlot({
 #   cvs_input_dm()
 # })
 # 
 # output$norm_dm <- renderPlot({
 #   norm_input_dm()
 # })
 # 
 # output$missval_dm <- renderPlot({
 #   missval_input_dm()
 # })
 # 
 # output$detect_dm <- renderPlot({
 #   detect_input_dm()
 # })
 # 
 # output$imputation_dm <- renderPlot({
 #   imputation_input_dm()
 # })
 # 
 # output$p_hist <- renderPlot({
 #   p_hist_input_dm()
 # })
 # 
 # output$numbers_dm <- renderPlot({
 #   numbers_input_dm()
 # })
 # 
 # output$coverage_dm <- renderPlot({
 #   coverage_input_dm()
 # })
 # 
 # ## Enrichment Outputs
 # output$go_enrichment_dm<-renderPlot({
 #   go_input_dm()$plot_go
 # })
 # 
 # output$pathway_enrichment_dm<-renderPlot({
 #   pathway_input_dm()$plot_pa
 # })
 
 ##### Download Functions for demo
 # datasetInput_dm <- reactive({
 #   switch(input$dataset_dm,
 #          "Results" = get_results_proteins(dep_dm()),
 #          "Full dataset" = get_df_wide(dep_dm()))
 # })
 # 
 # output$downloadData_dm <- downloadHandler(
 #   filename = function() { paste(input$dataset_dm, ".csv", sep = "") }, ## use = instead of <-
 #   content = function(file) {
 #     write.table(datasetInput_dm(),
 #                 file,
 #                 col.names = TRUE,
 #                 row.names = FALSE,
 #                 sep =",") }
 # )
 
 ### === Cluster Download for demo ==== ####
 
 # individual_cluster_dm <- reactive({
 #   cluster_number <- input$cluster_number_dm
 #   cluster_all <- heatmap_input_dm()
 #   data_result_dm()[cluster_all[[cluster_number]],]
 # })
 # 
 # 
 # 
 # output$downloadCluster_dm <- downloadHandler(
 #   filename = function() { paste("Cluster_info_",input$cluster_number_dm, ".csv", sep = "") }, ## use = instead of <-
 #   content = function(file) {
 #     write.table(individual_cluster_dm(),
 #                 file,
 #                 col.names = TRUE,
 #                 row.names = FALSE,
 #                 sep =",") }
 # )
 # 
 # output$downloadVolcano_dm <- downloadHandler(
 #   filename = function() {
 #     paste0("Volcano_", input$volcano_cntrst_dm, ".pdf")
 #   },
 #   content = function(file) {
 #     pdf(file)
 #     print(volcano_input_selected_dm())
 #     dev.off()
 #   }
 # )
 
 
 ## Protein plot download for demo
 # output$downloadProtein_dm <- downloadHandler(
 #   filename = function() {
 #     paste0(input$type_dm,".pdf")
 #   },
 #   content = function(file) {
 #     pdf(file)
 #     print(protein_input_dm())
 #     dev.off()
 #   }
 # )
 
 ###### ==== DOWNLOAD GO TABLE for demo ==== ####
 # output$downloadGO_dm <- downloadHandler(
 #   filename = function() { paste("GO_enrichment_",input$go_database_dm, ".csv", sep = "") }, ## use = instead of <-
 #   content = function(file) {
 #     write.table(go_input_dm()$go_result,
 #                 file,
 #                 col.names = TRUE,
 #                 row.names = FALSE,
 #                 sep =",") }
 # )
 
 ###### ==== DOWNLOAD PATHWAY TABLE for demo ==== ####
 # output$downloadPA_dm <- downloadHandler(
 #   filename = function() { paste("Pathway_enrichment_",input$pathway_database_dm, ".csv", sep = "") }, 
 #   ## use = instead of <-
 #   content = function(file) {
 #     write.table(pathway_input_dm()$pa_result,
 #                 file,
 #                 col.names = TRUE,
 #                 row.names = FALSE,
 #                 sep =",") }
 # )
 
 #####===== Download Report for demo =====#####
 # output$downloadReport_dm <- downloadHandler(
 #   
 #   filename = "LFQ-Analyst_report.pdf",
 #   content = function(file) {
 #     file.copy("www/LFQ-Analyst_report.pdf",file)
 #     # tempReport <- file.path(tempdir(), "LFQ_report.Rmd")
 #     # file.copy("./templates/LFQ_report.Rmd", tempReport, overwrite = TRUE)
 #     # 
 #     # sig_proteins_dm<-dep_dm() %>%
 #     #   .[SummarizedExperiment::rowData(.)$significant, ] %>%
 #     #   nrow()
 #     # 
 #     # tested_contrasts_dm<- gsub("_p.adj", "", 
 #     #                            colnames(SummarizedExperiment::rowData(dep()))[grep("p.adj", 
 #     #    colnames(SummarizedExperiment::rowData(dep_dm())))])
 #     # pg_width_dm<- ncol(processed_data_dm()) / 2.5
 #     # # Set up parameters to pass to Rmd document
 #     # params <- list(data = processed_data_dm,
 #     #                alpha = 0.05,
 #     #                lfc = 1,
 #     #                num_signif= sig_proteins_dm,
 #     #                tested_contrasts= tested_contrasts_dm,
 #     #                pg_width = pg_width_dm,
 #     #                numbers_input= numbers_input_dm,
 #     #                pca_input = pca_input_dm,
 #     #                coverage_input= coverage_input_dm,
 #     #                correlation_input =correlation_input_dm,
 #     #                heatmap_input = heatmap_input_dm,
 #     #                dep = dep_dm
 #     # )
 #     # 
 #     # # Knit the document, passing in the `params` list
 #     # rmarkdown::render(tempReport, output_file = file,
 #     #                   params = params,
 #     #                   envir = new.env(parent = globalenv())
 #     # )
 #   }
 # )
}
