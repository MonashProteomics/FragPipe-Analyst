#Define server logic to read selected file ----
server <- function(input, output, session) {
  options(shiny.maxRequestSize=100*1024^2)## Set maximum upload size to 100MB
  
#  Show elements on clicking Start analysis button
   observeEvent(input$analyze ,{ 
     if(input$analyze==0){
       return()
     }
    shinyjs::hide("quickstart_info")
    shinyjs::show("panel_list")
    })
   
   # Hide LFQ page if only have one replicate in each sample
   observeEvent(input$analyze ,{ 
     exp <- exp_design_input()
     if (max(exp$replicate)==1){
       hideTab(inputId = "tab_panels", target = "lfq_panel")
     } else {
       showTab(inputId = "tab_panels", target = "lfq_panel")
       updateTabsetPanel(session, "tab_panels", selected = "lfq_panel")
     }
     
   })
   
   ## Shinyalert
   # observeEvent(input$analyze ,{ 
   #   if(input$analyze==0 ){
   #     return()
   #   }
   #   
   #   shinyalert("In Progress!", "Data analysis has started, wait until table and plots
   #              appear on the screen", type="info",
   #              closeOnClickOutside = TRUE,
   #              closeOnEsc = TRUE,
   #              timer = 10000) # timer in miliseconds (10 sec)
   # })
   
   observeEvent(input$tab_panels,{ 
     if(input$tab_panels==0 ){
       return()
     }
     
     shinyalert("In Progress!", "Data analysis has started, wait until table and plots
                appear on the screen", type="info",
                closeOnClickOutside = TRUE,
                closeOnEsc = TRUE,
                timer = 10000) # timer in miliseconds (10 sec)
   })
   
   observe({
   if (input$tabs_selected=="demo"){
     shinyalert("Demo results loading!...", "Wait until table and plots
                appear on the screen", type="info",
                closeOnClickOutside = TRUE,
                closeOnEsc = TRUE,
                timer = 6000)
   }
   })
 
   
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
       selectizeInput("contrast",
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
    
    ## make reactive elements
    maxquant_data_input<-reactive({NULL})
    exp_design_input<-reactive({NULL})
    exp_design_example<-reactive({NULL})
    maxquant_data_example<-reactive({NULL})
    
    maxquant_data_input<-eventReactive(input$analyze,{
      inFile<-input$file1
      if(is.null(inFile))
        return(NULL)
      temp_data<-read.delim(inFile$datapath,
                 header = TRUE,
                 fill= TRUE, # to fill any missing data
                 sep = "\t"
      )
      # change inconsistency column names
      colnames(temp_data)[names(temp_data) == "Description"] <- "Protein.names"
      colnames(temp_data)[names(temp_data) == "Gene" | names(temp_data) == "Gene.Names" ] <- "Gene.names"
      colnames(temp_data)[names(temp_data) == "Protein.ID"] <- "Protein.IDs"
      colnames(temp_data)[names(temp_data) == "Combined.Total.Peptides"] <- "Unique.Stripped.Peptides" # fragpipe different outputs
      # validate(maxquant_input_test(temp_data))
      return(temp_data)
    })
    
    # Creates handsontable template table
    exp_data1 <- reactive({
      req(input$file1)
      df <- read.delim(input$file1$datapath,
                       header = TRUE,
                       fill= TRUE, # to fill any missing data
                       sep = "\t")
      
      tempTable =  get_exp_design(df)
      rhandsontable(tempTable) %>% 
        hot_col("label", readOnly = T) 
    })
    
    # Outputs the template table
    output$exp_protein<- renderRHandsontable({exp_data1()})
    
    # Changes the handsontable back into a dataframe 
    exp_data2<-reactive({NULL})
    exp_data2 <- eventReactive(input$save_exp, {
      hot_to_r(input$exp_protein)
    })
    
    observeEvent(input$save_exp, {
      output$save_message <- renderText({
        if (sum(is.na(exp_data2())) != 0) {
          stop(safeError("Warning: Cells can not be empty"))
        }
        else {
          "Experiment design table saved"
        }
      })
    })
    
    # download edited exp_design table
    output$download_exp <- downloadHandler("LFQ-Analyst_experimental_design.txt",
                                           content = function(file){
                                             write.table(exp_data2(), file, 
                                                         sep = "\t", 
                                                         row.names = F,
                                                         quote = F)
                                           },
                                           contentType = "text/csv")
    
    # rewrite exp_design table
    observeEvent(input$showTable ,{
      if(input$showTable==0){
        return()
      }
      shinyjs::show("quickstart_info")
      shinyjs::hide("panel_list")
    })

    observeEvent(input$original_exp,{
      output$exp_protein<- renderRHandsontable({exp_data1()})
      output$save_message <- renderText({""})
      })
    
    # experiment design file
    exp_design_input<-eventReactive(input$analyze,{
      inFile<-input$file2
      if (is.null(inFile) || input$save_exp>0) {
        temp_df <- exp_data2()
      }
      else {
        temp_df<-read.delim(inFile$datapath,
                            header = TRUE,
                            sep="\t",
                            stringsAsFactors = FALSE)
        exp_design_test(temp_df)
        temp_df$label<-as.character(temp_df$label)
        temp_df$condition<-trimws(temp_df$condition, which = "left")
      }
      return(temp_df)
    })
   
### Reactive components
   processed_data<- reactive({
     ## check which dataset
     if(!is.null (maxquant_data_input() )){
       maxquant_data <- reactive({maxquant_data_input()})
     }
     
     if(!is.null (exp_design_input() )){
       exp_design<-reactive({exp_design_input()})
     }
     
     message(exp_design())
     
     #check maxquant columns    
     if(length(grep("^LFQ.", colnames(maxquant_data()))) !=0){
       if(any(grepl('+',maxquant_data()$Reverse))){
         filtered_data<-dplyr::filter(maxquant_data(),Reverse!="+")
       }
       else{filtered_data<-maxquant_data()}
       if(any(grepl('+',filtered_data$Potential.contaminant))){
         filtered_data<-dplyr::filter(filtered_data,Potential.contaminant!="+")
       }
       if(any(grepl('+',filtered_data$Only.identified.by.site))){
         filtered_data<-dplyr::filter(filtered_data,Only.identified.by.site!="+") 
       }
       if(input$single_peptide==TRUE){
         filtered_data <-filtered_data
       }
       else{filtered_data<-dplyr::filter(filtered_data,as.numeric(Razor...unique.peptides)>=2)}
       
       filtered_data<-ids_test(filtered_data)
       
       data_unique<- DEP::make_unique(filtered_data,"Gene.names","Protein.IDs",delim=";")
       lfq_columns<-grep("^LFQ.", colnames(data_unique))
       intensity_names <- colnames( data_unique[,lfq_columns]) %>% gsub("LFQ.intensity.", "", .)
     }
     else{
       if(input$single_peptide==TRUE){
         filtered_data <-maxquant_data()
       }
       # filtered out Stripped Peptides less than 2
       if( any(colnames(maxquant_data()) %in% "Unique.Stripped.Peptides")){
         filtered_data <-dplyr::filter(maxquant_data(),as.numeric(Unique.Stripped.Peptides)>=2)
       }
       else{
         filtered_data<-maxquant_data()
       }
       
       filtered_data<-ids_test(filtered_data)
       
       data_unique<- DEP::make_unique(filtered_data,"Gene.names","Protein.IDs",delim=";")
       
       # Compatible for Fragpipe datasets without "MaxLFQ" columns
       if(length(grep(".MaxLFQ.Intensity", colnames(data_unique))) !=0){
         lfq_columns<-grep(".MaxLFQ.Intensity", colnames(data_unique))
         intensity_names <- colnames( data_unique[,lfq_columns]) %>% gsub(".MaxLFQ.Intensity", "", .)
       } 
       else {
         # lfq_columns<-grep(".Total.Intensity", colnames(data_unique))
         # intensity_names <- colnames( data_unique[,lfq_columns]) %>% gsub(".Total.Intensity", "", .)
         all_intentisy_cols <- grep(".Intensity", colnames(data_unique))
         remove_cols <- grep("Unique.Intensity|Total.Intensity",colnames(data_unique))
         lfq_columns <- setdiff(all_intentisy_cols, remove_cols)
         intensity_names <- colnames( data_unique[,lfq_columns]) %>%
           gsub(".Intensity", "", .)  %>% 
           gsub(".MaxLFQ", "", .) %>% 
           gsub(".Razor", "",.) 
       }
       
       # rename intensity column names
       names(data_unique)[lfq_columns] <- c(intensity_names)
     }
     
     ### remove columns of samples not in experiment design table
     # remove_columns <- intensity_names[make.names(delete_prefix(intensity_names)) %in% make.names(delete_prefix(exp_design()$label)) ==FALSE]
     if (any(make.names(intensity_names) %in% make.names(exp_design()$label))){
       remove_columns <- intensity_names[make.names(intensity_names )%in% make.names(exp_design()$label) ==FALSE]
     } else {
       remove_columns <- intensity_names[make.names(delete_prefix(intensity_names) )%in% make.names(delete_prefix(exp_design()$label)) ==FALSE]
     }
     
     if (identical(remove_columns, character(0)) == FALSE){
       lfq_columns <- lfq_columns[-c(which(intensity_names %in% remove_columns))]
     } else {
       lfq_columns <- lfq_columns
     }
     
     ## Check for matching columns in maxquant and experiment design file
     test_match_lfq_column_design(data_unique,lfq_columns, exp_design())
     data_se<-DEP:::make_se(data_unique,lfq_columns,exp_design())
  
     # Check number of replicates
     if(max(exp_design()$replicate)<3){
       threshold<-0
     } else  if(max(exp_design()$replicate)==3){
       threshold<-1
     } else if(max(exp_design()$replicate)<6 ){
       threshold<-2
     } else if (max(exp_design()$replicate)>=6){
       threshold<-trunc(max(exp_design()$replicate)/2)
     }
     
     filter_missval(data_se,thr = threshold)
   })
   
   unimputed_table<-reactive({
     temp<-assay(processed_data())
     temp1<-2^(temp)
     colnames(temp1)<-paste(colnames(temp1),"original_intensity",sep="_")
     temp1<-cbind(ProteinID=rownames(temp1),temp1) 
     #temp1$ProteinID<-rownames(temp1)
     return(as.data.frame(temp1))
   })
   
   normalised_data<-reactive({
     normalize_vsn(processed_data())
   })
   
   imputed_data<-reactive({
     DEP::impute(processed_data(),input$imputation)
   })
   
   imputed_table<-reactive({
     temp<-assay(imputed_data())
     temp1<-2^(temp)
     colnames(temp1)<-paste(colnames(temp1),"imputed_intensity",sep="_")
     temp1<-cbind(ProteinID=rownames(temp1),temp1) #temp1$ProteinID<-rownames(temp1)
     return(as.data.frame(temp1))
   })
   
   diff_all<-reactive({
     test_diff(imputed_data(),type = 'all')
   })

   dep<-reactive({
     if(input$fdr_correction=="BH"){
       diff_all<-test_limma(imputed_data(),type='all', paired = input$paired)
       add_rejections(diff_all,alpha = input$p, lfc= input$lfc)
     }
     else{
       diff_all<-test_diff(imputed_data(),type='all')
       add_rejections(diff_all,alpha = input$p, lfc= input$lfc)
     }
     
   })
   
   comparisons<-reactive ({
     temp<-capture.output(DEP::test_diff(imputed_data(),type='all'),type = "message")
     gsub(".*: ","",temp)
     ## Split conditions into character vector
     unlist(strsplit(temp,","))
     ## Remove leading and trailing spaces
     trimws(temp)
    })

   ## Results plot inputs
   
   ## PCA Plot
   pca_input<-eventReactive(input$analyze ,{ 
     if(input$analyze==0 ){
       return()
     }
     if (num_total()<=500){
       if(length(levels(as.factor(colData(dep())$replicate))) <= 6){
         pca_plot<-DEP::plot_pca(dep(), n=num_total(), point_size = 4)
         pca_plot<-pca_plot + labs(title = "PCA Plot")
         return(pca_plot)
       }
       else{
         pca_plot<-DEP::plot_pca(dep(), n=num_total(), point_size = 4, indicate = "condition") 
         pca_plot<-pca_plot + labs(title = "PCA Plot")
         return(pca_plot)
       }
     }
     else{
       if(length(levels(as.factor(colData(dep())$replicate))) <= 6){
         pca_plot<-DEP::plot_pca(dep(), point_size = 4)
         pca_plot<-pca_plot + labs(title = "PCA Plot")
         return(pca_plot)
       }
       else{
	     #pca_label<-SummarizedExperiment::colData(dep())$replicate
       pca_plot<-DEP::plot_pca(dep(), point_size = 4, indicate = "condition")
       #pca_plot<-pca_plot + geom_point()
       pca_plot<-pca_plot + ggrepel::geom_text_repel(aes(label=factor(rowname)),
                                           size = 4,
                                           box.padding = unit(0.1, 'lines'),
                                           point.padding = unit(0.1, 'lines'),
                                           segment.size = 0.5)
       pca_plot<-pca_plot + labs(title = "PCA Plot")
       return(pca_plot)
       }
     }
     
   })
   
   ### Heatmap Differentially expressed proteins
   heatmap_cluster <- eventReactive(input$analyze ,{ 
     if(input$analyze==0 ){
       return()
     }
     heatmap_list <- get_cluster_heatmap(dep(),
                                         type="centered",kmeans = TRUE,
                                         k=input$k_number, col_limit = 6,
                                         indicate = "condition"
     )
     return(heatmap_list)
   })
   
   heatmap_input <- reactive({
     heatmap_list <- heatmap_cluster()
     heatmap_list[[1]]
   })
   
   # heatmap_input<-eventReactive(input$analyze ,{ 
   #   if(input$analyze==0 ){
   #     return()
   #   }
   #   # get_cluster_heatmap(dep(),
   #   #                     type="centered",kmeans = TRUE,
   #   #                     k=input$k_number, col_limit = 6,
   #   #                     indicate = "condition"
   #   #                     )
   #   
   #   heatmap_list <- get_cluster_heatmap(dep(),
   #                                       type="centered",kmeans = TRUE,
   #                                       k=input$k_number, col_limit = 6,
   #                                       indicate = "condition"
   #   )
   #   heatmap_list[[1]]
   # })
   
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
        }
        else if(!is.null(input$protein_brush)){
        proteins_selected<-data_result()[data_result()[["Gene Name"]] %in% protein_name_brush(), ] 
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
                        name = proteins_selected$`Gene Name`,
                        percent_imputation = paste0("(",round((proteins_selected$num_NAs/length(colData(dep())$label))*100,1),"%)"))
       df_protein$percent_imputation[df_protein$percent_imputation=="(0%)"] <- ""
       #print(df_protein)
       p<-plot_volcano_new(dep(),
                    input$volcano_cntrst,
                    input$fontsize,
                    input$check_names,
                    input$p_adj)
       
       p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
         ggrepel::geom_text_repel(data = df_protein,
                                  aes(x, y, label = paste(name,percent_imputation)),
                                  size = 4,
                                  box.padding = unit(0.1, 'lines'),
                                  point.padding = unit(0.1, 'lines'),
                                  segment.size = 0.5)## use the dataframe to plot points

       }
    })
    
    protein_input<-reactive({ 
      
      protein_selected  <- data_result()[input$contents_rows_selected,1]
      
      if(length(levels(as.factor(colData(dep())$replicate))) <= 8){
        plot_protein(dep(), protein_selected, input$type)
      }
      else{
        protein_plot<-plot_protein(dep(), protein_selected, input$type)
        protein_plot + scale_color_brewer(palette = "Paired")
      }
      
    })
    
     
   ## QC Inputs
   norm_input <- reactive({
     plot_normalization(processed_data(),
                        normalised_data())
   })
   
   missval_input <- reactive({
     plot_missval(processed_data())
   })
   
   detect_input <- reactive({
     plot_detect(processed_data())
   })
   
   imputation_input <- reactive({
     plot_imputation(normalised_data(),
                     diff_all())
   })
   
   p_hist_input <- reactive({
     plot_p_hist(dep())
   })
   
   numbers_input <- reactive({
     plot_numbers(normalised_data())
   })
   
   coverage_input <- reactive({
     plot_coverage(normalised_data())
   })
   
   correlation_input<-reactive({
     plot_cor(dep(),significant = FALSE)
   })
   
   cvs_input<-reactive({
     plot_cvs(dep())
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
    enrichment_output_test(dep(), input$go_database)
    go_results<- test_gsea_mod(dep(), databases = input$go_database, contrasts = TRUE)
    null_enrichment_test(go_results)
    plot_go<- plot_enrichment(go_results, number = 5, alpha = 0.05, contrasts =input$contrast,
    databases = input$go_database, nrow = 2, term_size = 8) + aes(stringr::str_wrap(Term, 60)) +
     xlab(NULL)
    go_list<-list("go_result"=go_results, "plot_go"=plot_go)
    return(go_list)
    }
   })

   pathway_input<-eventReactive(input$pathway_analysis,{
     progress_indicator("Pathway Analysis is running....")
     enrichment_output_test(dep(), input$pathway_database)
     pathway_results<- test_gsea_mod(dep(), databases=input$pathway_database, contrasts = TRUE)
     null_enrichment_test(pathway_results)
     plot_pathway<-plot_enrichment(pathway_results, number = 5, alpha = 0.05, contrasts =input$contrast_1,
               databases=input$pathway_database, nrow = 3, term_size = 8) + aes(stringr::str_wrap(Term, 30)) +
       xlab(NULL)
     pathway_list<-list("pa_result"=pathway_results, "plot_pa"=plot_pathway)
     return(pathway_list)
     })

   
#### Interactive UI
   output$significantBox <- renderInfoBox({
     num_total <- dep() %>%
       nrow()
     num_signif <- dep() %>%
       .[SummarizedExperiment::rowData(.)$significant, ] %>%
       nrow()
     frac <- num_signif / num_total
     
       info_box <- 		infoBox("Significant proteins",
                             paste0(num_signif,
                                    " out of ",
                                    num_total),
                             paste0(signif(frac * 100, digits = 3),
                                    "% of proteins differentially expressed across all conditions"),
                             icon = icon("stats", lib = "glyphicon"),
                             color = "olive",
                            # fill = TRUE,
                             width = 4)
     
     return(info_box)
   })

  ##### Get results dataframe from Summarizedexperiment object
    data_result<-reactive({
      get_results_proteins(dep())
      #get_results(dep())
    })
    
    
 #### Data table
  output$contents <- DT::renderDataTable({
    df<- data_result() %>% dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
    return(df)
  },
  options = list(scrollX = TRUE,
autoWidth=TRUE,
                columnDefs= list(list(width = '400px', targets = c(-1))))
  ) 
  
  ## Deselect all rows button
  proxy <- dataTableProxy("contents")
  
  observeEvent(input$clear,{
    proxy %>% selectRows(NULL)
  })
  
  observeEvent(input$original,{
    output$contents <- DT::renderDataTable({
      df<- data_result() %>% dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
      return(df)
    },
    options = list(scrollX = TRUE,
     autoWidth=TRUE,
                columnDefs= list(list(width = '400px', targets = c(-1))))
    )
  })

   protein_name_brush<- reactive({
     #protein_tmp<-nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
     protein_tmp<-brushedPoints(volcano_df(), input$protein_brush, 
                                xvar = "diff", yvar = "p_values")
     protein_selected<-protein_tmp$name
     }) 
   protein_name_click<- reactive({
     protein_tmp<-nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
    # protein_tmp<-brushedPoints(volcano_df(), input$protein_brush, 
                                #xvar = "diff", yvar = "p_values")
     protein_selected<-protein_tmp$name
   }) 
  
  ## Select rows dynamically
   
   brush <- NULL
   makeReactiveBinding("brush")
   
  observeEvent(input$protein_brush,{
    output$contents <- DT::renderDataTable({
      df<- data_result()[data_result()[["Gene Name"]] %in% protein_name_brush(), ] %>% 
        dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
      return(df)
    },
    options = list(scrollX= TRUE)
    )
    
    # proteins_selected<-data_result()[data_result()[["Gene Name"]] %in% protein_name_brush(), ] ## get all rows selected
    # ## convert contrast to x and padj to y
    # diff_proteins <- grep(paste("^",input$volcano_cntrst, "_log2", sep = ""),
    #                       colnames(proteins_selected))
    # if(input$p_adj=="FALSE"){
    #   padj_proteins <- grep(paste("^",input$volcano_cntrst, "_p.val", sep = ""),
    #                         colnames(proteins_selected))
    # }
    # else{
    #   padj_proteins <- grep(paste("^",input$volcano_cntrst, "_p.adj", sep = ""),
    #                         colnames(proteins_selected))
    # }
    # df_protein <- data.frame(x = proteins_selected[, diff_proteins],
    #                          y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
    #                          name = proteins_selected$`Gene Name`)
    # #print(df_protein)
    # 
    # p<-plot_volcano_new(dep(),
    #                     input$volcano_cntrst,
    #                     input$fontsize,
    #                     input$check_names,
    #                     input$p_adj)
    # 
    #   p<- p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
    #   ggrepel::geom_text_repel(data = df_protein,
    #                            aes(x, y, label = name),
    #                            size = 4,
    #                            box.padding = unit(0.1, 'lines'),
    #                            point.padding = unit(0.1, 'lines'),
    #                            segment.size = 0.5)
    # 
    # output$volcano <- renderPlot({
    #   withProgress(message = 'Volcano Plot calculations are in progress',
    #                detail = 'Please wait for a while', value = 0, {
    #                  for (i in 1:15) {
    #                    incProgress(1/15)
    #                    Sys.sleep(0.25)
    #                  }
    #                })
    #  p
    # })
    # return(p)
  })
 
 observeEvent(input$resetPlot,{
   session$resetBrush("protein_brush")
   brush <<- NULL
   
   output$contents <- DT::renderDataTable({
     df<- data_result() %>% dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
     return(df)
   },
   options = list(scrollX = TRUE,
                  autoWidth=TRUE,
                  columnDefs= list(list(width = '400px', targets = c(-1))))
   )
 })

 observeEvent(input$protein_click,{
   output$contents <- DT::renderDataTable({
     df<- data_result()[data_result()[["Gene Name"]] %in% protein_name_click(), ] %>% 
       dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
     return(df)
   },
   options = list(scrollX= TRUE,
   autoWidth=TRUE,
                columnDefs= list(list(width = '400px', targets = c(-1))))
   )
 })
 
  ## Render Result Plots
  output$pca_plot<-renderPlot({
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
    if(is.null(input$contents_rows_selected) & is.null(input$protein_brush)){
      volcano_input()
      }
    else if(!is.null(input$volcano_cntrst)){
      volcano_input_selected()
      }# else close
  })
  
  output$protein_plot<-renderPlot({
    if(!is.null(input$contents_rows_selected)){
    protein_input()
    }
  })
  
  ### Abundance plot panel ####
  #comparisons
  output$abundance_cntrst <- renderUI({
    if (!is.null(comparisons())) {
      df <- SummarizedExperiment::rowData(dep())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("abundance_cntrst",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  ## abundance rank plot
  #abundance rank plot brush
  protein_name_brush_rank <- reactive({
    protein_tmp<-brushedPoints(data_result(), input$protein_brush_rank,
                               xvar = "rank", yvar = "mean_abundance")
    protein_selected<-protein_tmp$`Gene Name`
  })
  protein_name_click_rank <- reactive({
    protein_tmp<-nearPoints(data_result(), input$protein_click_rank, maxpoints = 1)
    protein_selected<-protein_tmp$`Gene Name`
  })
  
  abundance_rank_input <- reactive({
    df <- SummarizedExperiment::rowData(dep())
    cols <- grep("_significant$",colnames(df))
    contrast <- gsub("_significant", "", colnames(df)[cols])[1]
    p_list <- plot_abundance(data_result(), contrast)
    p_list[[1]]
  })
  
  abundance_rank_input_selected<-reactive({
    if (!is.null(input$contents_rows_selected)){
      proteins_selected<-data_result()[c(input$contents_rows_selected),]## get all rows selected
    }
    else if(!is.null(input$protein_brush_rank)){
      proteins_selected<-data_result()[data_result()[["Gene Name"]] %in% protein_name_brush_rank(), ] 
    }
    
    df_protein <- data.frame(x = proteins_selected$rank,
                             y = proteins_selected$mean_abundance,
                             name = proteins_selected$`Gene Name`)
    
    df <- SummarizedExperiment::rowData(dep())
    cols <- grep("_significant$",colnames(df))
    contrast <- gsub("_significant", "", colnames(df)[cols])[1]
    p_list <- plot_abundance(data_result(),
                             contrast)
    
    p_list[[1]] + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
      ggrepel::geom_label_repel(data = df_protein,
                                aes(x, y, label = name),
                                nudge_y = 0.5,
                                size = 4,
                                box.padding = unit(0.1, 'lines'),
                                point.padding = unit(0.1, 'lines'),
                                segment.size = 0.5)## use the dataframe to plot points
    
  })
  
  output$abundance_rank <- renderPlot({
    withProgress(message = 'Abundance Plot calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    if(is.null(input$contents_rows_selected) & is.null(input$protein_brush_rank)){
      abundance_rank_input()
    }
    else {
      abundance_rank_input_selected()
    }
  })
  
  # abundance comparison plot
  #abundance rank plot brush
  protein_name_brush_comp <- reactive({
    if(!is.null(input$abundance_cntrst)){
      contrast1 <- input$abundance_cntrst %>% gsub("_vs.*", "",.)
      contrast2 <- input$abundance_cntrst %>% gsub("^.*vs_", "",.)
        
      protein_tmp<-brushedPoints(data_result(), input$protein_brush_comp,
                                 xvar = paste("mean", contrast1, sep = "_"), yvar = paste("mean", contrast2, sep = "_"))
      protein_selected<-protein_tmp$`Gene Name`
    }
  })
  protein_name_click_comp <- reactive({
    if(!is.null(input$abundance_cntrst)){
      contrast1 <- input$abundance_cntrst %>% gsub("_vs.*", "",.)
      contrast2 <- input$abundance_cntrst %>% gsub("^.*vs_", "",.)
      
      
      protein_tmp<-nearPoints(data_result(), input$protein_click_comp, 
                              xvar = paste("mean", contrast1, sep = "_"), yvar = paste("mean", contrast2, sep = "_"),
                              maxpoints = 1)
      protein_selected<-protein_tmp$`Gene Name`
    }
  })

  abundance_comp_input <- reactive({
    if(!is.null(input$abundance_cntrst)) {
      p_list <- plot_abundance(data_result(),
                               input$abundance_cntrst)
      p_list[[2]]
    }
  })

  abundance_comp_input_selected<-reactive({
    if(!is.null(input$abundance_cntrst)){

      if (!is.null(input$contents_rows_selected)){
        proteins_selected<-data_result()[c(input$contents_rows_selected),]## get all rows selected
      }
      else if(!is.null(input$protein_brush_comp)){
        proteins_selected<-data_result()[data_result()[["Gene Name"]] %in% protein_name_brush_comp(), ]
      }
      
      contrast1 <- input$abundance_cntrst %>% gsub("_vs.*", "",.)
      contrast2 <- input$abundance_cntrst %>% gsub("^.*vs_", "",.)
      
      df_protein <- data.frame(x = proteins_selected[,grep(paste("^mean", contrast1, sep = "_"), colnames(proteins_selected))],
                                   y = proteins_selected[,grep(paste("^mean", contrast2, sep = "_"), colnames(proteins_selected))],
                                   name = proteins_selected$`Gene Name`)

      p_list <- plot_abundance(data_result(),
                               input$abundance_cntrst)

      p_list[[2]] + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
        ggrepel::geom_label_repel(data = df_protein,
                                 aes(x, y, label = name),
                                 nudge_y = 0.5,
                                 size = 4,
                                 box.padding = unit(0.1, 'lines'),
                                 point.padding = unit(0.1, 'lines'),
                                 segment.size = 0.5)## use the dataframe to plot points
    }
  })
  
  output$abundance_comp <- renderPlot({
    withProgress(message = 'Abundance Plot calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    if(is.null(input$contents_rows_selected) & is.null(input$protein_brush_comp)){
      abundance_comp_input()
    }
    else if(!is.null(input$abundance_cntrst)){
      abundance_comp_input_selected()
    }
  })
  
  output$downloadAbundance_rank <- downloadHandler(
    filename = function() {
      "Abundance_rank.svg"
    },
    content = function(file) {
      if(is.null(input$contents_rows_selected)){
        p <- abundance_rank_input()
      }
      else{
        p <- abundance_rank_input_selected()
      }
      svg(file, width = 8, height = 8)
      print(p)
      dev.off()
    }
  )
  
  output$downloadAbundance_comp <- downloadHandler(
    filename = function() {
      paste0("Abundance_", input$abundance_cntrst, ".svg")
    },
    content = function(file) {
      if(is.null(input$contents_rows_selected)){
        p <- abundance_comp_input()
      }
      else{
        p <- abundance_comp_input_selected()
      }
      svg(file, width = 8, height = 8)
      print(p)
      dev.off()
    }
  )
 
  ## Select rows dynamically in abundance rank plot
  observeEvent(input$protein_brush_rank,{
    output$contents <- DT::renderDataTable({
      df<- data_result()[data_result()[["Gene Name"]] %in% protein_name_brush_rank(), ]  %>% 
        dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
      return(df)
    },
    options = list(scrollX= TRUE)
    )
  })
  
  observeEvent(input$protein_click_rank,{
    output$contents <- DT::renderDataTable({
      df<- data_result()[data_result()[["Gene Name"]] %in% protein_name_click_rank(), ]  %>% 
        dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
      return(df)
    },
    options = list(scrollX= TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
  })
  
  ## Select rows dynamically in abundance comparison plot
  
  brush_comp <- NULL
  makeReactiveBinding("brush_comp")
  
  observeEvent(input$protein_brush_comp,{
    output$contents <- DT::renderDataTable({
      df<- data_result()[data_result()[["Gene Name"]] %in% protein_name_brush_comp(), ]  %>% 
        dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
      return(df)
    },
    options = list(scrollX= TRUE)
    )
  })
  
  observeEvent(input$protein_click_comp,{
    output$contents <- DT::renderDataTable({
      df<- data_result()[data_result()[["Gene Name"]] %in% protein_name_click_comp(), ]  %>% 
        dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
      return(df)
    },
    options = list(scrollX= TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
  })
  
  # reset abundance plots
  observeEvent(input$resetPlot_rank,{
    session$resetBrush("protein_brush_rank")
    brush <<- NULL
    
    output$contents <- DT::renderDataTable({
      df<- data_result() %>% dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
      return(df)
    },
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
  })
  
  observeEvent(input$resetPlot_comp,{
    session$resetBrush("protein_brush_comp")
    brush_comp <<- NULL
    
    output$contents <- DT::renderDataTable({
      df<- data_result() %>% dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
      return(df)
    },
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
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
  output$go_enrichment<-renderPlot({
    go_input()$plot_go
  })
  
  output$pathway_enrichment<-renderPlot({
   pathway_input()$plot_pa
  })
  
  ##### Download Functions
  datasetInput <- reactive({
    switch(input$dataset,
           "Results" = get_results_proteins(dep()),
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
  
  # individual_cluster <- reactive({
  #     cluster_number <- input$cluster_number
  #     heatmap_list <- get_cluster_heatmap(dep(),
  #                                         type="centered",kmeans = TRUE,
  #                                         k=input$k_number, col_limit = 6,
  #                                         indicate = "condition"
  #     )
  #     cluster_all <-heatmap_list[[2]]
  #     data_result()[cluster_all[[cluster_number]],]
  #   })
  individual_cluster <- reactive({
    cluster_number <- input$cluster_number
    cluster_all <-heatmap_cluster()[[2]]
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
    content = function(file) {
      # heatmap_plot<-DEP::plot_heatmap(dep(),"centered", k=6, indicate = "condition")
      svg(file)
      print(heatmap_input())
      dev.off()
    }
  )
  
#####===== Download Report =====#####
  output$downloadReport <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "LFQ-Analyst_report.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "LFQ_report.Rmd")
      file.copy("LFQ_report.Rmd", tempReport, overwrite = TRUE)
      
      sig_proteins<-dep() %>%
        .[SummarizedExperiment::rowData(.)$significant, ] %>%
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
                     p_hist_input = p_hist_input,
                     pca_input = pca_input,
                     coverage_input= coverage_input,
                     correlation_input =correlation_input,
                     heatmap_input = heatmap_input,
                     cvs_input= cvs_input,
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

 
 
 #### Demo logic ========== #############
 
 ####======= Render Functions
 
 output$volcano_cntrst_dm <- renderUI({
   if (!is.null(comparisons_dm())) {
     df <- SummarizedExperiment::rowData(dep_dm())
     cols <- grep("_significant$",colnames(df))
     selectizeInput("volcano_cntrst_dm",
                    "Comparison",
                    choices = gsub("_significant", "", colnames(df)[cols]))
   }
 })
 
 ##comparisons
 output$contrast_dm <- renderUI({
   if (!is.null(comparisons_dm())) {
     df <- SummarizedExperiment::rowData(dep_dm())
     cols <- grep("_significant$",colnames(df))
     selectizeInput("contrast_dm",
                    "Comparison",
                    choices = gsub("_significant", "", colnames(df)[cols]))
   }
 })
 
 output$contrast_dm_1 <- renderUI({
   if (!is.null(comparisons_dm())) {
     df <- SummarizedExperiment::rowData(dep_dm())
     cols <- grep("_significant$",colnames(df))
     selectizeInput("contrast_dm_1",
                    "Comparison",
                    choices = gsub("_significant", "", colnames(df)[cols]))
   }
 })
 
 output$downloadTable_dm <- renderUI({
   if(!is.null(dep_dm())){
     selectizeInput("dataset_dm",
                    "Download data table" ,
                    c("Results",
                      "Full dataset"))
   }
 })
 
 output$downloadButton_dm <- renderUI({
   if(!is.null(dep_dm())){
     downloadButton('downloadData_dm', 'Save')
   }
 })
 
 
 output$downloadreport_dm <- renderUI({
   if(!is.null(dep_dm())){
     downloadButton('downloadReport_dm', 'Download Report')
   }
 })
 
 output$downloadPlots_dm <- renderUI({
   if(!is.null(dep_dm())){
     downloadButton('downloadPlots1_dm', 'Download Plots')
   }
 })
 
 env_dm<-reactive({
   LoadToEnvironment("data/demo_data.RData", env = globalenv())
 })
 
 
 
 dep_dm<-reactive({
   env_dm()[["dep"]]
 })
 
 
comparisons_dm<-reactive({
  comparisons<-gsub("_p.adj", "", colnames(SummarizedExperiment::rowData(dep_dm()))[grep("p.adj", 
                          colnames(SummarizedExperiment::rowData(dep_dm())))])
})
 ## Results plot inputs
 
 ## PCA Plot

pca_label_dm<-reactive({
 pca_lable<-levels(as.factor(colData(dep_dm())$replicate))
print(pca_label)
})

 pca_input_dm<-reactive({
   if (num_total_dm()<=500){
     if(length(levels(as.factor(colData(dep_dm())$replicate))) <= 6){
       pca_plot<- DEP::plot_pca(dep_dm(), n=num_total_dm(), point_size = 4)
       pca_plot<-pca_plot + labs(title = "PCA plot")
       return(pca_plot)
     }
     else{
       pca_plot<-DEP::plot_pca(dep_dm(), n=num_total_dm(), point_size = 4, indicate = "condition")
       pca_plot<-pca_plot + labs(title = "PCA plot")
       return(pca_plot)
     }
   }
   else{
     if(length(levels(as.factor(colData(dep_dm())$replicate))) <= 6){
       pca_plot<-DEP::plot_pca(dep_dm(), point_size = 4)
       pca_plot<-pca_plot + labs(title = "PCA plot")
       return(pca_plot)
     }else{
       pca_label<-SummarizedExperiment::colData(dep_dm())$replicate
       pca_plot<-DEP::plot_pca(dep_dm(), point_size = 4, indicate = "condition")
       #pca_plot<-pca_plot + geom_point()
       pca_plot<-pca_plot + ggrepel::geom_text_repel(aes(label=factor(rowname)),
                                           size = 4,
                                           box.padding = unit(0.1, 'lines'),
                                           point.padding = unit(0.1, 'lines'),
                                           segment.size = 0.5)
       pca_plot<-pca_plot + labs(title = "PCA plot")
	  return(pca_plot)
     }
   }
   
 })
 
 ### Heatmap Differentially expressed proteins
 heatmap_input_dm<-reactive({ 
   get_cluster_heatmap(dep_dm(),
                       type="centered",kmeans = TRUE,
                       k=6, col_limit = 6,
                       indicate = "condition"
   )
 })
 
 ### Volcano Plot
 volcano_input_dm <- reactive({
   if(!is.null(input$volcano_cntrst_dm)) {
     plot_volcano_new(dep_dm(),
                      input$volcano_cntrst_dm,
                      input$fontsize_dm,
                      input$check_names_dm,
                      input$p_adj_dm)
     
   }
 })
 
 volcano_df_dm<- reactive({
   if(!is.null(input$volcano_cntrst_dm)) {
     get_volcano_df(dep_dm(),
                    input$volcano_cntrst_dm)
     
   }
 })
 
 
 volcano_input_selected_dm<-reactive({
   if(!is.null(input$volcano_cntrst_dm)){
     if (!is.null(input$contents_dm_rows_selected)){
       proteins_selected<-data_result_dm()[c(input$contents_dm_rows_selected),]## get all rows selected
     }
     else if(!is.null(input$protein_brush_dm)){
       proteins_selected<-data_result_dm()[data_result_dm()[["Gene Name"]] %in% protein_name_brush_dm(), ] 
     }
     ## convert contrast to x and padj to y
     diff_proteins <- grep(paste(input$volcano_cntrst_dm, "_log2", sep = ""),
                           colnames(proteins_selected))
     if(input$p_adj_dm=="FALSE"){
       padj_proteins <- grep(paste(input$volcano_cntrst_dm, "_p.val", sep = ""),
                             colnames(proteins_selected))
     }
     else{
       padj_proteins <- grep(paste(input$volcano_cntrst_dm, "_p.adj", sep = ""),
                             colnames(proteins_selected))
     }
     
     # df_protein <- data.frame(x = proteins_selected[, diff_proteins],
     #                          y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
     #                          name = proteins_selected$`Gene Name`)
     df_protein <- data.frame(x = proteins_selected[, diff_proteins],
                              y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
                              name = proteins_selected$`Gene Name`,
                              percent_imputation = paste0("(",round((proteins_selected$num_NAs/length(colData(dep_dm())$label))*100,1),"%)"))
     df_protein$percent_imputation[df_protein$percent_imputation=="(0%)"] <- ""
     #print(df_protein)
     p<-plot_volcano_new(dep_dm(),
                         input$volcano_cntrst_dm,
                         input$fontsize_dm,
                         input$check_names_dm,
                         input$p_adj_dm)
     
     p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
       ggrepel::geom_text_repel(data = df_protein,
                                aes(x, y, label = paste(name,percent_imputation)),
                                size = 4,
                                box.padding = unit(0.1, 'lines'),
                                point.padding = unit(0.1, 'lines'),
                                segment.size = 0.5)## use the dataframe to plot points
     
   }
 })
 
 protein_input_dm<-reactive({ 
   
   protein_selected  <- data_result_dm()[input$contents_dm_rows_selected,1]
   
   #protein<-row_selected$name
   if(length(levels(as.factor(colData(dep_dm())$replicate))) <= 8){
   plot_protein(dep_dm(), protein_selected, input$type_dm)
   }
   else{
     protein_plot<-plot_protein(dep_dm(), protein_selected, input$type_dm)
     protein_plot + scale_color_brewer(palette = "Paired")
     }
   
 })
 
 ## Get processed data
 
 processed_data_dm<-reactive({
   env_dm()[["data_filter"]]
 })
 normalised_data_dm<-reactive({
   DEP::normalize_vsn(processed_data_dm())
 })
	
	
imputed_data_dm<-reactive({
   DEP::impute(processed_data_dm(),input$imputation)
 })
 
 
 diff_all_dm<-reactive({
   test_diff(imputed_data_dm(),type = 'all')
 })
	
 ## QC Inputs
 norm_input_dm <- reactive({
   plot_normalization(processed_data_dm(),
                      normalised_data_dm())
 })
 
 missval_input_dm <- reactive({
   plot_missval(processed_data_dm())
 })
 
 detect_input_dm <- reactive({
   plot_detect(processed_data_dm())
 })
 
 imputation_input_dm <- reactive({
   plot_imputation(normalised_data_dm(),
                   diff_all_dm())
 })
 
 p_hist_input_dm <- reactive({
   plot_p_hist(dep_dm())
 })
 
 numbers_input_dm <- reactive({
   plot_numbers(normalised_data_dm())
 })
 
 coverage_input_dm <- reactive({
   plot_coverage(normalised_data_dm())
 })
 
 correlation_input_dm<-reactive({
   plot_cor(dep_dm())
 })
 
 cvs_input_dm<-reactive({
   plot_cvs(dep_dm())
 })
 
 num_total_dm<-reactive({
   dep_dm() %>%
     nrow()
 }) 
 
 ## Enrichment inputs
 
 go_input_dm<-eventReactive(input$go_analysis_dm,{
   withProgress(message = 'Gene ontology enrichment is in progress',
                detail = 'Please wait for a while', value = 0, {
                  for (i in 1:15) {
                    incProgress(1/15)
                    Sys.sleep(0.25)
                  }
                })
   
   if(!is.null(input$contrast_dm)){
     enrichment_output_test(dep_dm(), input$go_database_dm)
     go_results<- test_gsea_mod(dep_dm(), databases = input$go_database_dm, contrasts = TRUE)
     plot_go<- plot_enrichment(go_results, number = 5, alpha = 0.05, contrasts =input$contrast_dm,
                               databases = input$go_database_dm, nrow = 2, term_size = 8) + 
       aes(stringr::str_wrap(Term, 60)) +
       xlab(NULL)
     go_list<-list("go_result"=go_results, "plot_go"=plot_go)
     return(go_list)
   }
 })
 
 pathway_input_dm<-eventReactive(input$pathway_analysis_dm,{
   withProgress(message = 'Pathway enrichment is in progress',
                detail = 'Please wait for a while', value = 0, {
                  for (i in 1:15) {
                    incProgress(1/15)
                    Sys.sleep(0.25)
                  }
                })
   enrichment_output_test(dep_dm(), input$pathway_database_dm)
   pathway_results<- test_gsea_mod(dep_dm(), databases=input$pathway_database_dm, contrasts = TRUE)
   plot_pathway<-plot_enrichment(pathway_results, number = 5, alpha = 0.05, contrasts =input$contrast_dm_1,
                                 databases=input$pathway_database_dm, nrow = 3, term_size = 8) + 
     aes(stringr::str_wrap(Term, 30)) +
     xlab(NULL)
   pathway_list<-list("pa_result"=pathway_results, "plot_pa"=plot_pathway)
   return(pathway_list)
   
 })
 
 
 #### Interactive UI
 output$significantBox_dm <- renderInfoBox({
   num_total <- dep_dm() %>%
     nrow()
   num_signif <- dep_dm() %>%
     .[SummarizedExperiment::rowData(.)$significant, ] %>%
     nrow()
   frac <- num_signif / num_total
     
     info_box <- 		infoBox("Significant proteins",
                           paste0(num_signif,
                                  " out of ",
                                  num_total),
                           paste0(signif(frac * 100, digits = 3),
                                  "% of proteins differentially expressed across all conditions"),
                           icon = icon("stats", lib = "glyphicon"),
                           color = "olive",
                           # fill = TRUE,
                           width = 4)
     
     return(info_box)
   })
   
 
 ##### Get results dataframe from Summarizedexperiment object
 data_result_dm<-reactive({
   get_results_proteins(dep_dm())
   #get_results(dep())
 })
 
 
 #### Data table
 output$contents_dm <- DT::renderDataTable({
   df<- data_result_dm() %>% dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
   return(df)
 },
 options = list(scrollX = TRUE,
	autoWidth=TRUE,
                columnDefs= list(list(width = '400px', targets = c(-1))))
 )
 
 ## Deselect all rows button
 proxy <- dataTableProxy("contents_dm")
 
 observeEvent(input$clear_dm,{
   proxy %>% selectRows(NULL)
 })
 
 observeEvent(input$original_dm,{
   output$contents_dm <- DT::renderDataTable({
     df<- data_result_dm() %>% dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
     return(df)
   },
   options = list(scrollX = TRUE,
	autoWidth=TRUE,
                columnDefs= list(list(width = '400px', targets = c(-1))))
   )
 })
 
 protein_name_brush_dm<- reactive({
   #protein_tmp<-nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
   protein_tmp<-brushedPoints(volcano_df_dm(), input$protein_brush_dm, 
                              xvar = "diff", yvar = "p_values")
   protein_selected<-protein_tmp$name
 }) 
 protein_name_click_dm<- reactive({
   protein_tmp<-nearPoints(volcano_df_dm(), input$protein_click_dm, maxpoints = 1)
   # protein_tmp<-brushedPoints(volcano_df(), input$protein_brush, 
   #xvar = "diff", yvar = "p_values")
   protein_selected<-protein_tmp$name
 }) 
 
 
 ## Select rows dynamically
 observeEvent(input$protein_brush_dm,{
   output$contents_dm <- DT::renderDataTable({
     df<- data_result_dm()[data_result_dm()[["Gene Name"]] %in% protein_name_brush_dm(), ] %>% 
       dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
     return(df)
   },
   options = list(scrollX= TRUE,
                  autoWidth=TRUE,
                  columnDefs= list(list(width = '400px', targets = c(-1))))
   )
   
   # proteins_selected<-data_result_dm()[data_result_dm()[["Gene Name"]] %in% protein_name_brush_dm(), ] #
   # # get all rows selected
   # ## convert contrast to x and padj to y
   # diff_proteins <- grep(paste(input$volcano_cntrst_dm, "_log2", sep = ""),
   #                       colnames(proteins_selected))
   # if(input$p_adj=="FALSE"){
   #   padj_proteins <- grep(paste(input$volcano_cntrst_dm, "_p.val", sep = ""),
   #                         colnames(proteins_selected))
   # }
   # else{
   #   padj_proteins <- grep(paste(input$volcano_cntrst_dm, "_p.adj", sep = ""),
   #                         colnames(proteins_selected))
   # }
   # df_protein <- data.frame(x = proteins_selected[, diff_proteins],
   #                          y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
   #                          name = proteins_selected$`Gene Name`)
   # #print(df_protein)
   # 
   # p<-plot_volcano_new(dep_dm(),
   #                     input$volcano_cntrst_dm,
   #                     input$fontsize_dm,
   #                     input$check_names_dm,
   #                     input$p_adj_dm)
   # 
   # p<- p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
   #   ggrepel::geom_text_repel(data = df_protein,
   #                            aes(x, y, label = name),
   #                            size = 4,
   #                            box.padding = unit(0.1, 'lines'),
   #                            point.padding = unit(0.1, 'lines'),
   #                            segment.size = 0.5)
   # 
   # output$volcano_dm <- renderPlot({
   #   withProgress(message = 'Volcano Plot calculations are in progress',
   #                detail = 'Please wait for a while', value = 0, {
   #                  for (i in 1:15) {
   #                    incProgress(1/15)
   #                    Sys.sleep(0.25)
   #                  }
   #                })
   #   p
   # })
   # return(p)
 })
 
 observeEvent(input$resetPlot_dm,{
   session$resetBrush("protein_brush_dm")
   brush <<- NULL
   
   output$contents_dm <- DT::renderDataTable({
     df<- data_result_dm() %>% dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
     return(df)
   },
   options = list(scrollX = TRUE,
                  autoWidth=TRUE,
                  columnDefs= list(list(width = '400px', targets = c(-1))))
   )
 })
 
 observeEvent(input$protein_click_dm,{
   output$contents_dm <- DT::renderDataTable({
     df<- data_result_dm()[data_result_dm()[["Gene Name"]] %in% protein_name_click_dm(), ] %>% 
       dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
     return(df)
   },
   options = list(scrollX= TRUE)
   )
 })
 ## Render Result Plots
 output$pca_plot_dm<-renderPlot({
   pca_input_dm()
 })
 output$heatmap_dm<-renderPlot({
   withProgress(message = 'Heatmap rendering is in progress',
                detail = 'Please wait for a while', value = 0, {
                  for (i in 1:15) {
                    incProgress(1/15)
                    Sys.sleep(0.25)
                  }
                })
   heatmap_input_dm()
 })
 
 output$volcano_dm <- renderPlot({
   withProgress(message = 'Volcano Plot calculations are in progress',
                detail = 'Please wait for a while', value = 0, {
                  for (i in 1:15) {
                    incProgress(1/15)
                    Sys.sleep(0.25)
                  }
                })
   if(is.null(input$contents_dm_rows_selected) & is.null(input$protein_brush_dm)){
     volcano_input_dm()
   }
   else if(!is.null(input$volcano_cntrst_dm)){
     volcano_input_selected_dm()
   } # else close
 })
 
 output$protein_plot_dm<-renderPlot({
   if(!is.null(input$contents_dm_rows_selected)){
     protein_input_dm()
   }
 })
 
 ### Demo Abundance plot panel ####
 #comparisons
 output$abundance_cntrst_dm <- renderUI({
   if (!is.null(comparisons_dm())) {
     df <- SummarizedExperiment::rowData(dep_dm())
     cols <- grep("_significant$",colnames(df))
     selectizeInput("abundance_cntrst_dm",
                    "Comparison",
                    choices = gsub("_significant", "", colnames(df)[cols]))
   }
 })
 
 ## abundance rank plot
 #abundance rank plot brush
 protein_name_brush_rank_dm <- reactive({
   protein_tmp<-brushedPoints(data_result_dm(), input$protein_brush_rank_dm,
                              xvar = "rank", yvar = "mean_abundance")
   protein_selected<-protein_tmp$`Gene Name`
 })
 protein_name_click_rank_dm <- reactive({
   protein_tmp<-nearPoints(data_result_dm(), input$protein_click_rank_dm, maxpoints = 1)
   protein_selected<-protein_tmp$`Gene Name`
 })
 
 abundance_rank_input_dm <- reactive({
   df <- SummarizedExperiment::rowData(dep_dm())
   cols <- grep("_significant$",colnames(df))
   contrast <- gsub("_significant", "", colnames(df)[cols])[1]
   p_list <- plot_abundance(data_result_dm(), contrast)
   p_list[[1]]
 })
 
 abundance_rank_input_selected_dm<-reactive({
   if (!is.null(input$contents_dm_rows_selected)){
     proteins_selected<-data_result_dm()[c(input$contents_dm_rows_selected),]## get all rows selected
   }
   else if(!is.null(input$protein_brush_rank_dm)){
     proteins_selected<-data_result_dm()[data_result_dm()[["Gene Name"]] %in% protein_name_brush_rank_dm(), ] 
   }
   
   df_protein <- data.frame(x = proteins_selected$rank,
                            y = proteins_selected$mean_abundance,
                            name = proteins_selected$`Gene Name`)
   
   df <- SummarizedExperiment::rowData(dep_dm())
   cols <- grep("_significant$",colnames(df))
   contrast <- gsub("_significant", "", colnames(df)[cols])[1]
   p_list <- plot_abundance(data_result_dm(),
                            contrast)
   
   p_list[[1]] + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
     ggrepel::geom_label_repel(data = df_protein,
                               aes(x, y, label = name),
                               nudge_y = 0.5,
                               size = 4,
                               box.padding = unit(0.1, 'lines'),
                               point.padding = unit(0.1, 'lines'),
                               segment.size = 0.5)## use the dataframe to plot points
   
 })
 
 output$abundance_rank_dm <- renderPlot({
   withProgress(message = 'Abundance Plot calculations are in progress',
                detail = 'Please wait for a while', value = 0, {
                  for (i in 1:15) {
                    incProgress(1/15)
                    Sys.sleep(0.25)
                  }
                })
   if(is.null(input$contents_dm_rows_selected) & is.null(input$protein_brush_rank_dm)){
     abundance_rank_input_dm()
   }
   else {
     abundance_rank_input_selected_dm()
   }
 })
 
 # abundance comparison plot
 #abundance rank plot brush
 protein_name_brush_comp_dm <- reactive({
   if(!is.null(input$abundance_cntrst_dm)){
     contrast1 <- input$abundance_cntrst_dm %>% gsub("_vs.*", "",.)
     contrast2 <- input$abundance_cntrst_dm %>% gsub("^.*vs_", "",.)
     
     protein_tmp<-brushedPoints(data_result_dm(), input$protein_brush_comp_dm,
                                xvar = paste("mean", contrast1, sep = "_"), yvar = paste("mean", contrast2, sep = "_"))
     protein_selected<-protein_tmp$`Gene Name`
   }
 })
 protein_name_click_comp_dm <- reactive({
   if(!is.null(input$abundance_cntrst_dm)){
     contrast1 <- input$abundance_cntrst_dm %>% gsub("_vs.*", "",.)
     contrast2 <- input$abundance_cntrst_dm %>% gsub("^.*vs_", "",.)
     
     
     protein_tmp<-nearPoints(data_result_dm(), input$protein_click_comp_dm, 
                             xvar = paste("mean", contrast1, sep = "_"), yvar = paste("mean", contrast2, sep = "_"),
                             maxpoints = 1)
     protein_selected<-protein_tmp$`Gene Name`
   }
 })
 
 abundance_comp_input_dm <- reactive({
   if(!is.null(input$abundance_cntrst_dm)) {
     p_list <- plot_abundance(data_result_dm(),
                              input$abundance_cntrst_dm)
     p_list[[2]]
   }
 })
 
 abundance_comp_input_selected_dm<-reactive({
   if(!is.null(input$abundance_cntrst_dm)){
     
     if (!is.null(input$contents_dm_rows_selected)){
       proteins_selected<-data_result_dm()[c(input$contents_dm_rows_selected),]## get all rows selected
     }
     else if(!is.null(input$protein_brush_comp_dm)){
       proteins_selected<-data_result_dm()[data_result_dm()[["Gene Name"]] %in% protein_name_brush_comp_dm(), ]
     }
     
     contrast1 <- input$abundance_cntrst_dm %>% gsub("_vs.*", "",.)
     contrast2 <- input$abundance_cntrst_dm %>% gsub("^.*vs_", "",.)
     
     df_protein <- data.frame(x = proteins_selected[,grep(paste("^mean", contrast1, sep = "_"), colnames(proteins_selected))],
                              y = proteins_selected[,grep(paste("^mean", contrast2, sep = "_"), colnames(proteins_selected))],
                              name = proteins_selected$`Gene Name`)
     
     p_list <- plot_abundance(data_result_dm(),
                              input$abundance_cntrst_dm)
     
     p_list[[2]] + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
       ggrepel::geom_label_repel(data = df_protein,
                                 aes(x, y, label = name),
                                 nudge_y = 0.5,
                                 size = 4,
                                 box.padding = unit(0.1, 'lines'),
                                 point.padding = unit(0.1, 'lines'),
                                 segment.size = 0.5)## use the dataframe to plot points
   }
 })
 
 output$abundance_comp_dm <- renderPlot({
   withProgress(message = 'Abundance Plot calculations are in progress',
                detail = 'Please wait for a while', value = 0, {
                  for (i in 1:15) {
                    incProgress(1/15)
                    Sys.sleep(0.25)
                  }
                })
   if(is.null(input$contents_dm_rows_selected) & is.null(input$protein_brush_comp_dm)){
     abundance_comp_input_dm()
   }
   else if(!is.null(input$abundance_cntrst_dm)){
     abundance_comp_input_selected_dm()
   }
 })
 
 output$downloadAbundance_rank_dm <- downloadHandler(
   filename = function() {
     "Abundance_rank.svg"
   },
   content = function(file) {
     if(is.null(input$contents_dm_rows_selected) & is.null(input$protein_brush_rank_dm) ){
       p <- abundance_rank_input_dm()
     }
     else{
       p <- abundance_rank_input_selected_dm()
     }
     svg(file, width = 8, height = 8)
     print(p)
     dev.off()
   }
 )
 
 output$downloadAbundance_comp_dm <- downloadHandler(
   filename = function() {
     paste0("Abundance_", input$abundance_cntrst_dm, ".svg")
   },
   content = function(file) {
     if(is.null(input$contents_dm_rows_selected) & is.null(input$protein_brush_comp_dm)){
       p <- abundance_comp_input_dm()
     }
     else{
       p <- abundance_comp_input_selected_dm()
     }
     svg(file, width = 8, height = 8)
     print(p)
     dev.off()
   }
 )
 
 ## Select rows dynamically in abundance rank plot
 
 # brush_rank_dm <- NULL
 # makeReactiveBinding("brush_rank_dm")
 
 observeEvent(input$protein_brush_rank_dm,{
   output$contents_dm <- DT::renderDataTable({
     df<- data_result_dm()[data_result_dm()[["Gene Name"]] %in% protein_name_brush_rank_dm(), ]  %>% 
       dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
     return(df)
   },
   options = list(scrollX= TRUE)
   )
 })
 
 observeEvent(input$protein_click_rank_dm,{
   output$contents_dm <- DT::renderDataTable({
     df<- data_result_dm()[data_result_dm()[["Gene Name"]] %in% protein_name_click_rank_dm(), ]  %>% 
       dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
     return(df)
   },
   options = list(scrollX= TRUE,
                  autoWidth=TRUE,
                  columnDefs= list(list(width = '400px', targets = c(-1))))
   )
 })
 
 ## Select rows dynamically in abundance comparison plot
 
 brush_comp_dm <- NULL
 makeReactiveBinding("brush_comp_dm")
 
 observeEvent(input$protein_brush_comp_dm,{
   output$contents_dm <- DT::renderDataTable({
     df<- data_result_dm()[data_result_dm()[["Gene Name"]] %in% protein_name_brush_comp_dm(), ]  %>% 
       dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
     return(df)
   },
   options = list(scrollX= TRUE)
   )
 })
 
 observeEvent(input$protein_click_comp_dm,{
   output$contents_dm <- DT::renderDataTable({
     df<- data_result_dm()[data_result_dm()[["Gene Name"]] %in% protein_name_click_comp_dm(), ]  %>% 
       dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
     return(df)
   },
   options = list(scrollX= TRUE,
                  autoWidth=TRUE,
                  columnDefs= list(list(width = '400px', targets = c(-1))))
   )
 })
 
 # reset abundance plots
 observeEvent(input$resetPlot_rank_dm,{
   session$resetBrush("protein_brush_rank_dm")
   brush <<- NULL
   
   output$contents_dm <- DT::renderDataTable({
     df<- data_result_dm() %>% dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
     return(df)
   },
   options = list(scrollX = TRUE,
                  autoWidth=TRUE,
                  columnDefs= list(list(width = '400px', targets = c(-1))))
   )
 })
 
 observeEvent(input$resetPlot_comp_dm,{
   session$resetBrush("protein_brush_comp_dm")
   brush_comp_dm <<- NULL
   
   output$contents_dm <- DT::renderDataTable({
     df<- data_result_dm() %>% dplyr::select(-dplyr::starts_with("mean"),-"rank") # drop mean abundance columns
     return(df)
   },
   options = list(scrollX = TRUE,
                  autoWidth=TRUE,
                  columnDefs= list(list(width = '400px', targets = c(-1))))
   )
 })
 
 
 
 
 
 
 ### QC Outputs
 output$sample_corr_dm <-renderPlot({
   correlation_input_dm()
 })
 
 output$sample_cvs_dm <- renderPlot({
   cvs_input_dm()
 })
 
 output$norm_dm <- renderPlot({
   norm_input_dm()
 })
 
 output$missval_dm <- renderPlot({
   missval_input_dm()
 })
 
 output$detect_dm <- renderPlot({
   detect_input_dm()
 })
 
 output$imputation_dm <- renderPlot({
   imputation_input_dm()
 })
 
 output$p_hist <- renderPlot({
   p_hist_input_dm()
 })
 
 output$numbers_dm <- renderPlot({
   numbers_input_dm()
 })
 
 output$coverage_dm <- renderPlot({
   coverage_input_dm()
 })
 
 ## Enrichment Outputs
 output$go_enrichment_dm<-renderPlot({
   go_input_dm()$plot_go
 })
 
 output$pathway_enrichment_dm<-renderPlot({
   pathway_input_dm()$plot_pa
 })
 
 ##### Download Functions
 datasetInput_dm <- reactive({
   switch(input$dataset_dm,
          "Results" = get_results_proteins(dep_dm()),
          "Full dataset" = get_df_wide(dep_dm()))
 })
 
 output$downloadData_dm <- downloadHandler(
   filename = function() { paste(input$dataset_dm, ".csv", sep = "") }, ## use = instead of <-
   content = function(file) {
     write.table(datasetInput_dm(),
                 file,
                 col.names = TRUE,
                 row.names = FALSE,
                 sep =",") }
 )
 
 ### === Cluster Download ==== ####
 
 individual_cluster_dm <- reactive({
   cluster_number <- input$cluster_number_dm
   cluster_all <- heatmap_input_dm()
   data_result_dm()[cluster_all[[cluster_number]],]
 })
 
 
 
 output$downloadCluster_dm <- downloadHandler(
   filename = function() { paste("Cluster_info_",input$cluster_number_dm, ".csv", sep = "") }, ## use = instead of <-
   content = function(file) {
     write.table(individual_cluster_dm(),
                 file,
                 col.names = TRUE,
                 row.names = FALSE,
                 sep =",") }
 )
 
 output$downloadVolcano_dm <- downloadHandler(
   filename = function() {
     paste0("Volcano_", input$volcano_cntrst_dm, ".pdf")
   },
   content = function(file) {
     pdf(file)
     print(volcano_input_selected_dm())
     dev.off()
   }
 )
 
 
 ## Protein plot download
 output$downloadProtein_dm <- downloadHandler(
   filename = function() {
     paste0(input$type_dm,".pdf")
   },
   content = function(file) {
     pdf(file)
     print(protein_input_dm())
     dev.off()
   }
 )
 
 ###### ==== DOWNLOAD GO TABLE ==== ####
 output$downloadGO_dm <- downloadHandler(
   filename = function() { paste("GO_enrichment_",input$go_database_dm, ".csv", sep = "") }, ## use = instead of <-
   content = function(file) {
     write.table(go_input_dm()$go_result,
                 file,
                 col.names = TRUE,
                 row.names = FALSE,
                 sep =",") }
 )
 
 ###### ==== DOWNLOAD PATHWAY TABLE ==== ####
 output$downloadPA_dm <- downloadHandler(
   filename = function() { paste("Pathway_enrichment_",input$pathway_database_dm, ".csv", sep = "") }, 
   ## use = instead of <-
   content = function(file) {
     write.table(pathway_input_dm()$pa_result,
                 file,
                 col.names = TRUE,
                 row.names = FALSE,
                 sep =",") }
 )
 
 
 
 #####===== Download Report =====#####
 output$downloadReport_dm <- downloadHandler(
   
   filename = "LFQ-Analyst_report.pdf",
   content = function(file) {
     file.copy("www/LFQ-Analyst_report.pdf",file)
     # tempReport <- file.path(tempdir(), "LFQ_report.Rmd")
     # file.copy("./templates/LFQ_report.Rmd", tempReport, overwrite = TRUE)
     # 
     # sig_proteins_dm<-dep_dm() %>%
     #   .[SummarizedExperiment::rowData(.)$significant, ] %>%
     #   nrow()
     # 
     # tested_contrasts_dm<- gsub("_p.adj", "", 
     #                            colnames(SummarizedExperiment::rowData(dep()))[grep("p.adj", 
     #    colnames(SummarizedExperiment::rowData(dep_dm())))])
     # pg_width_dm<- ncol(processed_data_dm()) / 2.5
     # # Set up parameters to pass to Rmd document
     # params <- list(data = processed_data_dm,
     #                alpha = 0.05,
     #                lfc = 1,
     #                num_signif= sig_proteins_dm,
     #                tested_contrasts= tested_contrasts_dm,
     #                pg_width = pg_width_dm,
     #                numbers_input= numbers_input_dm,
     #                pca_input = pca_input_dm,
     #                coverage_input= coverage_input_dm,
     #                correlation_input =correlation_input_dm,
     #                heatmap_input = heatmap_input_dm,
     #                dep = dep_dm
     # )
     # 
     # # Knit the document, passing in the `params` list
     # rmarkdown::render(tempReport, output_file = file,
     #                   params = params,
     #                   envir = new.env(parent = globalenv())
     # )
   }
 )
 
 
 ### ===== Download Plots ===== #####
 # output$downloadPlots1_dm <- downloadHandler(
 #   filename = "Plots.pdf",
 #   content = function(file) {
 #     tempReport <- file.path(tempdir(), "Plots.Rmd")
 #     file.copy("templates/Plots.Rmd", tempReport, overwrite = TRUE)
 #     
 #     tested_contrasts_dm<- gsub("_p.adj", "", 
 #                                colnames(SummarizedExperiment::rowData(dep_dm()))[grep("p.adj", 
 #                                  colnames(SummarizedExperiment::rowData(dep_dm())))])
 #     pg_width_dm<- ncol(processed_data_dm()) / 2.5
 #     
 #     # Set up parameters to pass to Rmd document
 #     params <- list(tested_contrasts= tested_contrasts_dm,
 #                    pg_width = pg_width_dm,
 #                    numbers_input= numbers_input_dm,
 #                    detect_input = detect_input_dm,
 #                    imputation_input = imputation_input_dm,
 #                    missval_input = missval_input_dm,
 #                    p_hist_input = p_hist_input_dm,
 #                    pca_input = pca_input_dm,
 #                    coverage_input= coverage_input_dm,
 #                    correlation_input =correlation_input_dm,
 #                    heatmap_input = heatmap_input_dm,
 #                    cvs_input= cvs_input_dm,
 #                    dep = dep_dm
 #     )
 #     
 #     # Knit the document, passing in the `params` list
 #     rmarkdown::render(tempReport, output_file = file,
 #                       params = params,
 #                       envir = new.env(parent = globalenv())
 #     )
 #   }
 # ) ## Download plot close
 
 #### Occurrence page logic ####
 data_attendance<-reactive({
   protein_data <- maxquant_data_input()
   
   #check maxquant columns    
   if(length(grep("^LFQ.", colnames(protein_data))) !=0){
     lfq_columns<-grep("LFQ.", colnames(protein_data))
   }
   else {
     # Compatible for Fragpipe datasets without "MaxLFQ" columns
     if(length(grep(".MaxLFQ.Intensity", colnames(protein_data))) !=0){
       lfq_columns<-grep(".MaxLFQ.Intensity", colnames(protein_data))
     } 
     else {
       # lfq_columns<-grep(".Total.Intensity", colnames(protein_data))
       all_intentisy_cols <- grep(".Intensity", colnames(protein_data))
       remove_cols <- grep("Unique.Intensity|Total.Intensity",colnames(protein_data))
       lfq_columns <- setdiff(all_intentisy_cols, remove_cols)
     }
     
   }
   intensity_names <- colnames( protein_data[,lfq_columns]) 
   
   needed_cols <- c("Protein.IDs", "Gene.names",intensity_names,"Protein.names", 
                    "Reverse", "Potential.contaminant", "Only.identified.by.site", "Razor...unique.peptides", "Unique.Stripped.Peptides")
   df <- protein_data[,colnames(protein_data) %in% needed_cols]
   df <- dplyr::relocate(df, "Protein.names", .after = last_col())

   intensity_names_1 <- intensity_names %>% 
     gsub("LFQ.intensity[.]", "", .)  %>% 
     gsub(".Intensity", "", .)  %>% 
     gsub(".MaxLFQ", "", .) %>% 
     gsub(".Razor", "",.)
   
   colnames(df)[colnames(df) %in% intensity_names] <- intensity_names_1
   
   # get conditions
   exp_design <- exp_design_input()
   conditions <- exp_design$condition %>% unique()
   
   # replace intensity column names
   replace_protein <- paste("LFQ_intensity", exp_design$condition, exp_design$replicate,sep = "_") %>% unique()
   
   if (any(make.names(intensity_names_1) %in% make.names(exp_design$label))){
     colnames(df)[colnames(df) %in% intensity_names_1] <- replace_protein[match(make.names(intensity_names_1), 
                                                                                make.names(exp_design$label), nomatch = 0)]
   } else {
     colnames(df)[colnames(df) %in% intensity_names_1] <- replace_protein[match(make.names(delete_prefix(intensity_names_1)), 
                                                                                make.names(delete_prefix(exp_design$label)), nomatch = 0)]
   }
   
   # remove intensity columns not in experimental design file.
   df <- df[!is.na(names(df))]
   # filter if all intensity are 0
   df <- df[rowSums(df[,grep("^LFQ_", colnames(df))]) != 0,]
   
   for (i in 1:length(conditions)) {
     condition <- conditions[i]
     pattern <- paste(condition,"[[:digit:]]",sep = "_")
     df[paste0("#Occurences",sep = "_",condition)] <- rowSums(df %>% select(grep(pattern, colnames(df))) != 0)
     
     # change column order
     df <- dplyr::relocate(df, paste0("#Occurences",sep = "_",condition), .before = paste("LFQ_intensity", conditions[1],"1",sep = "_"), .after = NULL)
     # print(colnames(df))
     cols <- grep(paste0(condition, "$"),colnames(df))
     
     if (!is.null(input[[paste0("",condition)]])){
       df <- df %>%
         dplyr::filter(df[[cols]] >=input[[paste0("",condition)]][1] & df[[cols]] <=input[[paste0("",condition)]][2])
     }
   }
   
   if ("" %in% df$Gene.names){
     df$Gene.names[df["Gene.names"]==""] <- "NoGeneNameAvailable"}
   if ("" %in% df$Protein.names){
     df$Protein.names[df["Protein.names"]==""] <- "NoProteinNameAvailable"}
   return(df)
 })
 
 data_attendance_filtered <- reactive({
   filtered_data <- data_attendance()
   if (is.null(input$filtered_condition_maxquant) & is.null(input$filtered_condition_fragpipe)) {
     filtered_data <- filtered_data
   }
   else {
     if(("Reverse" %in% colnames(filtered_data)) & ('Reverse sequences' %in% input$filtered_condition_maxquant)){
       filtered_data<-dplyr::filter(filtered_data,Reverse!="+")
     }
     else{filtered_data <-filtered_data}
     if(("Potential.contaminant" %in% colnames(filtered_data)) & ('Potential contaminants' %in% input$filtered_condition_maxquant)){
       filtered_data<-dplyr::filter(filtered_data,Potential.contaminant!="+")
     }
     else{filtered_data <-filtered_data}
     if(("Only.identified.by.site" %in% colnames(filtered_data)) & ('"Only identified by site" Protein' %in% input$filtered_condition_maxquant)){
       filtered_data<-dplyr::filter(filtered_data,Only.identified.by.site!="+")
     }
     else{filtered_data <-filtered_data}
     if(("Razor...unique.peptides" %in% colnames(filtered_data)) & ('Proteins with peptiedes < 2' %in% input$filtered_condition_maxquant)){
       filtered_data<-dplyr::filter(filtered_data,as.numeric(Razor...unique.peptides)>=2)
     }
     else{filtered_data <-filtered_data}
     if(("Unique.Stripped.Peptides" %in% colnames(filtered_data)) & ('Proteins with peptiedes < 2' %in% input$filtered_condition_fragpipe)){
       filtered_data <-dplyr::filter(filtered_data,as.numeric(Unique.Stripped.Peptides)>=2)
     }
     else{filtered_data <-filtered_data}
   }
   
   colnames(filtered_data)[names(filtered_data) == "Razor...unique.peptides"] <- "#Peptides"
   colnames(filtered_data)[names(filtered_data) == "Unique.Stripped.Peptides"] <- "#Peptides"
   
   drop_cols <- c("Reverse", "Potential.contaminant", "Only.identified.by.site")
   filtered_data<- filtered_data[, !(colnames(filtered_data) %in% drop_cols)]
   return(filtered_data)
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
   
   sliderInput(paste0("",conditions[n]),
               label=paste0("",conditions[n]),
               min = min(0), 
               max = max(exp_design_input$replicate),
               value = c(0, max(exp_design_input$replicate)),
               step = 1)}
 
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
     exp_design <- exp_design_input()
     conditions <- exp_design$condition %>% unique() %>% sort()
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
     set1 <- df$Protein.IDs[df[grep(paste0("#Occurences",sep = "_",input$condition_1),colnames(df))] != 0]
     set2 <- df$Protein.IDs[df[grep(paste0("#Occurences",sep = "_",input$condition_2),colnames(df))] != 0]
     x <- list(set1,set2)
     names(x) <- c("Condition 1", "Condition 2")
   } else {
     set1 <- df$Protein.IDs[df[grep(paste0("#Occurences",sep = "_",input$condition_1),colnames(df))] != 0]
     set2 <- df$Protein.IDs[df[grep(paste0("#Occurences",sep = "_",input$condition_2),colnames(df))] != 0]
     set3 <- df$Protein.IDs[df[grep(paste0("#Occurences",sep = "_",input$condition_3),colnames(df))] != 0]
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

 #### Occurrence demo page logic ####
 env_raw <- reactive({
   LoadToEnvironment("data/raw_demo_data.RData", env = globalenv())
 })
  
 exp_design_dm <- reactive({
   env_raw()[["demo_exp_design"]]
 }) 
 raw_demo_data <- reactive({
   env_raw()[["demo_input"]]
   })
   
 data_attendance_dm<-reactive({
   protein_data <- raw_demo_data()

   #check maxquant columns
   if(length(grep("^LFQ.", colnames(protein_data))) !=0){
     lfq_columns<-grep("LFQ.", colnames(protein_data))
   }
   else {
     # Compatible for Fragpipe datasets without "MaxLFQ" columns
     if(length(grep(".MaxLFQ.Intensity", colnames(protein_data))) !=0){
       lfq_columns<-grep(".MaxLFQ.Intensity", colnames(protein_data))
     }
     else {
       # lfq_columns<-grep(".Total.Intensity", colnames(protein_data))
       all_intentisy_cols <- grep(".Intensity", colnames(data_unique))
       remove_cols <- grep("Unique.Intensity|Total.Intensity",colnames(data_unique))
       lfq_columns <- setdiff(all_intentisy_cols, remove_cols)
     }

   }
   intensity_names <- colnames( protein_data[,lfq_columns])

   needed_cols <- c("Protein.IDs", "Gene.names",intensity_names,"Protein.names",
                    "Reverse", "Potential.contaminant", "Only.identified.by.site", "Razor...unique.peptides", "Unique.Stripped.Peptides")
   # df <- protein_data %>% select(Protein.IDs, Gene.names, all_of(lfq_columns), Protein.names)
   df <- protein_data[,colnames(protein_data) %in% needed_cols]
   # colnames(df)[names(df) == "Peptide.counts..razor.unique."] <- "Peptide_counts(razor+unique)"
   df <- dplyr::relocate(df, "Protein.names", .after = last_col())
   # rename intensity names
   intensity_names_1 <- intensity_names %>% 
     gsub("LFQ.intensity[.]", "", .)  %>% 
     gsub(".Intensity", "", .)  %>% 
     gsub(".MaxLFQ", "", .) %>% 
     gsub(".Razor", "",.)
   colnames(df)[colnames(df) %in% intensity_names] <- intensity_names_1
   
   # get conditions
   exp_design <- exp_design_dm()
   conditions <- exp_design$condition %>% unique()

   # replace intensity column names
   replace_protein <- paste("LFQ_intensity", exp_design$condition, exp_design$replicate,sep = "_") %>% unique()
   colnames(df)[colnames(df) %in% intensity_names_1] <- replace_protein[match(colnames(df), exp_design$label, nomatch = 0)]

   # remove intensity columns not in experimental design file.
   df <- df[!is.na(names(df))]
   # filter if all intensity are 0
   df <- df[rowSums(df[,grep("^LFQ_", colnames(df))]) != 0,]

   for (i in 1:length(conditions)) {
     condition <- conditions[i]
     pattern <- paste(condition,"[[:digit:]]",sep = "_")
     df[paste0("#Occurences",sep = "_",condition)] <- rowSums(df %>% select(grep(pattern, colnames(df))) != 0)

     # change column order
     df <- dplyr::relocate(df, paste0("#Occurences",sep = "_",condition), .before = paste("LFQ_intensity", conditions[1],"1",sep = "_"), .after = NULL)
     # print(colnames(df))
     cols <- grep(paste0(condition, "$"),colnames(df))

     if (!is.null(input[[paste0("",condition)]])){
       df <- df %>%
         dplyr::filter(df[[cols]] >=input[[paste0("",condition)]][1] & df[[cols]] <=input[[paste0("",condition)]][2])
     }
   }
   
   if ("" %in% df$Gene.names){
     df$Gene.names[df["Gene.names"]==""] <- "NoGeneNameAvailable"}
   if ("" %in% df$Protein.names){
     df$Protein.names[df["Protein.names"]==""] <- "NoProteinNameAvailable"}
   return(df)
 })

 data_attendance_filtered_dm <- reactive({
   filtered_data <- data_attendance_dm()
   if (is.null(input$filtered_condition_maxquant_dm) & is.null(input$filtered_condition_fragpipe_dm)) {
     filtered_data <- filtered_data
   }
   else {
     if(("Reverse" %in% colnames(filtered_data)) & ('Reverse sequences' %in% input$filtered_condition_maxquant_dm)){
       filtered_data<-dplyr::filter(filtered_data,Reverse!="+")
     }
     else{filtered_data <-filtered_data}
     if(("Potential.contaminant" %in% colnames(filtered_data)) & ('Potential contaminants' %in% input$filtered_condition_maxquant_dm)){
       filtered_data<-dplyr::filter(filtered_data,Potential.contaminant!="+")
     }
     else{filtered_data <-filtered_data}
     if(("Only.identified.by.site" %in% colnames(filtered_data)) & ('"Only identified by site" Protein' %in% input$filtered_condition_maxquant_dm)){
       filtered_data<-dplyr::filter(filtered_data,Only.identified.by.site!="+")
     }
     else{filtered_data <-filtered_data}
     if(("Razor...unique.peptides" %in% colnames(filtered_data)) & ('Proteins with peptiedes < 2' %in% input$filtered_condition_maxquant_dm)){
       filtered_data<-dplyr::filter(filtered_data,as.numeric(Razor...unique.peptides)>=2)
     }
     else{filtered_data <-filtered_data}
     if(("Unique.Stripped.Peptides" %in% colnames(filtered_data)) & ('Proteins with peptiedes < 2' %in% input$filtered_condition_fragpipe_dm)){
       filtered_data <-dplyr::filter(filtered_data,as.numeric(Unique.Stripped.Peptides)>=2)
     }
     else{filtered_data <-filtered_data}
   }

   colnames(filtered_data)[names(filtered_data) == "Razor...unique.peptides"] <- "#Peptides"
   colnames(filtered_data)[names(filtered_data) == "Unique.Stripped.Peptides"] <- "#Peptides"

   drop_cols <- c("Reverse", "Potential.contaminant", "Only.identified.by.site")
   filtered_data<- filtered_data[, !(colnames(filtered_data) %in% drop_cols)]
   return(filtered_data)
 })

 #### Data table
 output$contents_occ_dm <- DT::renderDataTable({
   df<- data_attendance_filtered_dm()
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

 make_sliderInput_dm <- function(n= 1){
   exp_design_input <- exp_design_dm()
   conditions <- exp_design_input$condition %>% unique()

   sliderInput(paste0("",conditions[n]),
               label=paste0("",conditions[n]),
               min = min(0),
               max = max(exp_design_input$replicate),
               value = c(0, max(exp_design_input$replicate)),
               step = 1)}

 slider_bars_dm <- reactive({
   exp_design <- exp_design_dm()
   lapply(X = 1:length(unique(exp_design$condition)), FUN = make_sliderInput_dm)
 })

 output$sidebar_dm <- renderUI({
   tagList(slider_bars_dm())
 })

 output$download_attendance_dm <- downloadHandler("Occurrences_results_table.csv",
                                               content = function(file){
                                                 write.table(data_attendance_filtered_dm(),
                                                             file,
                                                             col.names = TRUE,
                                                             row.names = FALSE,
                                                             sep =",")
                                               },
                                               contentType = "text/csv")
 
 ## Demo Venn plot
 output$condition_1_dm <- renderUI({
   selectizeInput("condition_1_dm",
                  "Condition 1",
                  choices = c("Benign", "Malignant"),
                  selected = "Benign")
 })
 
 output$condition_2_dm <- renderUI({
   selectizeInput("condition_2_dm",
                  "Condition 2",
                  choices = c("Benign", "Malignant")[c("Benign", "Malignant")!= input$condition_1_dm],
                  selected = "Malignant")
 })
 
 venn_plot_input_dm <- reactive({
   df<- data_attendance_filtered_dm()
   set1 <- df$Protein.IDs[df[grep(paste0("#Occurences",sep = "_",input$condition_1_dm),colnames(df))] != 0]
   set2 <- df$Protein.IDs[df[grep(paste0("#Occurences",sep = "_",input$condition_2_dm),colnames(df))] != 0]
   x <- list(set1,set2)
   names(x) <- c("Condition 1", "Condition 2")
   ggVennDiagram::ggVennDiagram(x,label_alpha = 0) +
     scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
 })
 
 output$venn_plot_dm <- renderPlot({
   venn_plot_input_dm()
 })
 
 output$download_venn_svg_dm<-downloadHandler(
   filename = function() { "venn_plot.svg" }, 
   content = function(file) {
     svg(file)
     print(venn_plot_input_dm())
     dev.off()
   }
 )
 
  
}
