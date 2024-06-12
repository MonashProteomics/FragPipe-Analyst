enrichr_mod <- function(genes, databases = NULL) {
  # check gene type
  if (length(genes) != 0){
    if (all(startsWith(genes, "ENSG"))) {
      genes_map <- ensembldb::select(EnsDb.Hsapiens.v86,
                                 keys= genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
      genes <- genes_map$SYMBOL
    }
    httr::set_config(httr::config(ssl_verifypeer = 0L))
    cat("Uploading data to Enrichr... ")
    if (is.vector(genes) & ! all(genes == "") & length(genes) != 0) {
      temp <- POST(url="http://maayanlab.cloud/Enrichr/enrich",
                   body=list(list=paste(genes, collapse="\n")))
    } else if (is.data.frame(genes)) {
      temp <- POST(url="http://maayanlab.cloud/Enrichr/enrich",
                   body=list(list=paste(paste(genes[,1], genes[,2], sep=","),
                                        collapse="\n")))
    } else {
      warning("genes must be a non-empty vector of gene names or a dataframe with genes and score.")
    }
    GET(url="http://maayanlab.cloud/Enrichr/share")
    cat("Done.\n")
    dbs <- as.list(databases)
    dfSAF <- options()$stringsAsFactors
    options(stringsAsFactors = FALSE)
    result <- lapply(dbs, function(x) {
      cat("  Querying ", x, "... ", sep="")
      r <- GET(url="http://maayanlab.cloud/Enrichr/export",
               query=list(file="API", backgroundType=x))
      r <- gsub("&#39;", "'", intToUtf8(r$content))
      tc <- textConnection(r)
      r <- read.table(tc, sep = "\t", header = TRUE, quote = "", comment.char="")
      close(tc)
      cat("Done.\n")
      return(r)
    })
    options(stringsAsFactors = dfSAF)
    cat("Parsing results... ")
    names(result) <- dbs
    cat("Done.\n")
  } 
  else { # no genes provided
    result <- data.frame(Term=character(),
                         Overlap=character(),
                         P.value=double(),
                         Adjusted.P.value=double(),
                         Old.P.value=double(),
                         Old.Adjusted.P.value=double(),
                         Odds.Ratio=double(),
                         Combined.Score=double(),
                         Genes=character())
  }
  return(result)
}


###### ========= Test_gsea new


test_ora_mod <- function(dep,
                         databases,
                         contrasts = TRUE, direction="UP", log2_threshold=0.7, alpha=0.05) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(databases),
                          is.logical(contrasts),
                          length(contrasts) == 1)
  
  
  row_data <- rowData(dep, use.names = FALSE)
  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }
  
  
  
  # Run background list
  message("Background")
  if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "protein") {
    background <- unique(row_data$Gene)
  } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "gene") {
    background <- unique(row_data$ID)
  } else if (metadata(dep)$level == "protein") {
      background <- unique(gsub("[.].*", "", row_data$name))
  } else if (metadata(dep)$exp == "LFQ" & metadata(dep)$level == "peptide") {
    background <- unique(row_data$Gene)
  } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "peptide") {
    background <- unique(row_data$Gene)
  } else if (metadata(dep)$exp == "DIA"  & metadata(dep)$level == "peptide") {
    background <- unique(row_data$Genes)
  }
  
  background_enriched <- enrichr_mod(background, databases)
  df_background <- NULL
  for(database in databases) {
    temp <- background_enriched[database][[1]] %>%
      mutate(var = database)
    df_background <- rbind(df_background, temp)
  }
  df_background$contrast <- "background"
  df_background$n <- length(background)
  
  OUT <- df_background %>%
    mutate(bg_IN = as.numeric(gsub("/.*", "", Overlap)),
           bg_OUT = n - bg_IN) %>%
    select(Term, bg_IN, bg_OUT)
  
  if(contrasts) {
    # Get gene symbols
    
    df <- row_data %>%
      as.data.frame() %>%
      # select(name, ends_with("_significant")) %>%
      mutate(name = gsub("[.].*", "", name))

    constrast_columns <- df %>% select(ends_with("_significant")) %>% colnames()
    constrasts <- gsub("_significant", "", constrast_columns)

    # Run enrichR for every contrast
    df_enrich <- NULL
    for(contrast in constrast_columns) {
      message(gsub("_significant", "", contrast))
      # contrast column might have NA
      df[is.na(df[[contrast]]),contrast] <- F
      significant <- df
      if (direction == "UP"){
        significant <- significant[(significant[gsub("_significant", "_diff", contrast)] > log2_threshold),]
      } else if (direction == "DOWN") {
        significant <- significant[significant[gsub("_significant", "_diff", contrast)] < -log2_threshold,]
      }
      significant <- significant[significant[gsub("_significant", "_p.adj", contrast)] < alpha,]
      
      if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "protein") {
        genes <- unique(significant$Gene)
      } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "gene") {
        genes <- unique(significant$ID)
      } else if (metadata(dep)$level == "protein") {
        genes <- significant$name
      } else if (metadata(dep)$exp == "LFQ" & metadata(dep)$level == "peptide") {
        genes <- unique(significant$Gene)
      } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "peptide") {
        genes <- unique(significant$Gene)
      } else if (metadata(dep)$exp == "DIA" & metadata(dep)$level == "peptide") {
        genes <- unique(significant$Genes)
      }
      
      message(paste0(length(genes), " genes are submitted"))
      if (length(genes) != 0){
        enriched <- enrichr_mod(genes, databases)
        # print(colnames(enriched[["KEGG_2021_Human"]]))
        # Tidy output
        contrast_enrich <- NULL
        for(database in databases) {
          temp <- enriched[database][[1]] %>%
            mutate(var = database)
          contrast_enrich <- rbind(contrast_enrich, temp)
        }
        if (nrow(contrast_enrich) != 0) { # has enrichment
          contrast_enrich$contrast <- contrast
          contrast_enrich$n <- length(genes)
          # Background correction
          cat("Background correction... ")
          contrast_enrich <- contrast_enrich %>%
            mutate(IN = as.numeric(gsub("/.*", "", Overlap)),
                   OUT = n - IN) %>%
            select(-n) %>%
            left_join(OUT, by = "Term") %>%
            mutate(log_odds = log2((IN * bg_OUT) / (OUT * bg_IN)))
          cat("Done.")
        }
        df_enrich <- rbind(df_enrich, contrast_enrich) %>%
          mutate(contrast = gsub("_significant", "", contrast))
      } else {
        cat("No significant genes for enrichment analysis")
      }
    }
  } else {
    # Get gene symbols
    significant <- row_data %>%
      as.data.frame() %>%
      select(name, significant) %>%
      filter(significant) %>%
      mutate(name = gsub("[.].*", "", name))
    
    # Run enrichR
    if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "protein") {
      genes <- unique(significant$Gene)
    } else if (metadata(dep)$exp == "TMT" & metadata(dep)$level == "gene") {
      genes <- unique(significant$ID)
    } else if (metadata(dep)$level == "protein") {
      genes <- significant$name
    } else if (metadata(dep)$level == "peptide") {
      genes <- unique(significant$Gene)
    }

    enriched <- enrichr_mod(genes, databases)
    
    # Tidy output
    df_enrich <- NULL
    for(database in databases) {
      temp <- enriched[database][[1]] %>%
        mutate(var = database)
      df_enrich <- rbind(df_enrich, temp)
    }
    df_enrich$contrast <- "significant"
    df_enrich$n <- length(genes)
    
    # Background correction
    cat("Background correction... ")
    df_enrich <- df_enrich %>%
      mutate(IN = as.numeric(gsub("/.*", "", Overlap)),
             OUT = n - IN) %>%
      select(-n) %>%
      left_join(OUT, by = "Term") %>%
      mutate(log_odds = log2((IN * bg_OUT) / (OUT * bg_IN)))
    cat("Done.")
  }
  
  if (nrow(df_enrich) != 0) {
    df_enrich$p_hyper = phyper(q=(df_enrich$IN-1), m = df_enrich$bg_IN, n = df_enrich$bg_OUT, k = (df_enrich$IN+df_enrich$OUT),
                               lower.tail = F )
    df_enrich$p.adjust_hyper = p.adjust(df_enrich$p_hyper, method = "BH")
  }
  return(df_enrich)
}



#######################################################
## Plot Enrichment Results
#######################################################

plot_enrichment <- function(gsea_results, number = 10, alpha = 0.05,
                            contrasts = NULL, databases = NULL,  adjust=F, use_whole_proteome=F,
                            nrow = 1, term_size = 8) {
  assertthat::assert_that(is.data.frame(gsea_results),
                          is.numeric(number),
                          length(number) == 1,
                          is.numeric(alpha),
                          length(alpha) == 1,
                          is.numeric(term_size),
                          length(term_size) == 1,
                          is.numeric(nrow),
                          length(nrow) == 1)
  
  # Check gsea_results object
  if(any(!c("Term", "var",
            "contrast","Adjusted.P.value")
         %in% colnames(gsea_results))) {
    stop("'", deparse(substitute(gsea_results)),
         "' does not contain the required columns",
         "\nMake sure that HGNC gene symbols are present",
         "\n in your 'Gene Names' column of Results table",
         call. = FALSE)
  }
  
  no_enrichment_text <- paste("\n   No enrichment found.\n",
                              "       You can still download enrichment result table. \n")
  
  if(!is.null(contrasts)) {
    assertthat::assert_that(is.character(contrasts))
    
    
    valid_contrasts <- unique(gsea_results$contrast)
    
    if(!all(contrasts %in% valid_contrasts)) {
      return(ggplot() +
               annotate("text", x = 4, y = 25, size=8, label = no_enrichment_text) + 
               theme_void()
      )
    }
    if(!any(contrasts %in% valid_contrasts)) {
      contrasts <- contrasts[contrasts %in% valid_contrasts]
      message("Not all contrasts found",
              "\n Following contrasts are found: '",
              paste0(contrasts, collapse = "', '"), "'")
    }
    
    gsea_results <- filter(gsea_results, contrast %in% contrasts)
  }
  if(!is.null(databases)) {
    assertthat::assert_that(is.character(databases))
    
    valid_databases <- unique(gsea_results$var)
    
    if(all(!databases %in% valid_databases)) {
      valid_cntrsts_msg <- paste0("Valid databases are: '",
                                  paste0(valid_databases, collapse = "', '"),
                                  "'")
      stop("Not a valid database, please run `plot_gsea()`",
           "with valid databases as argument\n",
           valid_cntrsts_msg,
           call. = FALSE)
    }
    if(any(!databases %in% valid_databases)) {
      databases <- databases[databases %in% valid_databases]
      message("Not all databases found",
              "\nPlotting the following databases: '",
              paste0(databases, collapse = "', '"), "'")
    }
    
    gsea_results <- filter(gsea_results, var %in% databases)
  }
  # Get top enriched gene sets
  if (!use_whole_proteome) {
    if (adjust) {
      terms <- gsea_results %>%
        dplyr::group_by(contrast, var) %>%
        dplyr::filter(p.adjust_hyper <= alpha) %>%
        dplyr::arrange(p.adjust_hyper) %>%
        dplyr::slice(seq_len(number)) %>%
        .$Term
      subset <- gsea_results %>%
        dplyr::filter(Term %in% terms) %>%
        dplyr::arrange(var, p.adjust_hyper)
    } else {
      terms <- gsea_results %>%
        dplyr::group_by(contrast, var) %>%
        dplyr::filter(p_hyper <= alpha) %>%
        dplyr::arrange(p_hyper) %>%
        dplyr::slice(seq_len(number)) %>%
        .$Term
      subset <- gsea_results %>%
        dplyr::filter(Term %in% terms) %>%
        dplyr::arrange(var, p_hyper)
    }
  } else {
    if (adjust) {
      terms <- gsea_results %>%
        dplyr::group_by(contrast, var) %>%
        dplyr::filter(Adjusted.P.value <= alpha) %>%
        dplyr::arrange(Adjusted.P.value) %>%
        dplyr::slice(seq_len(number)) %>%
        .$Term
      subset <- gsea_results %>%
        dplyr::filter(Term %in% terms) %>%
        dplyr::arrange(var, Adjusted.P.value)
    } else {
      terms <- gsea_results %>%
        dplyr::group_by(contrast, var) %>%
        dplyr::filter(P.value <= alpha) %>%
        dplyr::arrange(P.value) %>%
        dplyr::slice(seq_len(number)) %>%
        .$Term
      subset <- gsea_results %>%
        dplyr::filter(Term %in% terms) %>%
        dplyr::arrange(var, P.value)
    }
  }
  subset$Term <- readr::parse_factor(subset$Term, levels = unique(subset$Term))
  subset$var <- readr::parse_factor(subset$var, levels = unique(subset$var))
  
  if (nrow(subset) == 0) {
    return(ggplot() +
             annotate("text", x = 4, y = 25, size=8, label = no_enrichment_text) + 
             theme_void()
    )
  } else {
    # Plot top enriched gene sets
    subset$Overlap_ratio <- sapply(subset$Overlap, function(x) eval(parse(text=x)))
    # return(ggplot(subset, aes(y = reorder(Term, Overlap_ratio), x=Overlap_ratio, size=IN, color=Adjusted.P.value)) +
    #   geom_point() +
    #   facet_wrap(~contrast, nrow = nrow) +
    #   scale_color_continuous(low="red", high="blue", name = "Adjusted.P.value",
    #                            guide=guide_colorbar(reverse=TRUE)) +
    #   labs(y = "Term") +
    #   theme_bw() +
    #   theme(legend.position = "top", legend.text = element_text(size = 9))
    # )
    if (!use_whole_proteome) {
      if (adjust){
        return(ggplot(subset, aes(y = reorder(Term, log_odds), x=log_odds, size=IN, color=p.adjust_hyper)) +
                 geom_point() +
                 facet_wrap(~contrast, nrow = nrow) +
                 scale_color_continuous(low="red", high="blue", name = "p.adjust",
                                        guide=guide_colorbar(reverse=T,
                                                             label.theme = element_text(angle = 90),
                                                             label.vjust = 0.5)) +
                 labs(y = "Term", x = "log2 Odds ratio", size = "size") +
                 theme_bw() +
                 theme(legend.position = "top", legend.text = element_text(size = 9)))
      } else {
        return(ggplot(subset, aes(y = reorder(Term, log_odds), x=log_odds, size=IN, color=p_hyper)) +
                 geom_point() +
                 facet_wrap(~contrast, nrow = nrow) +
                 scale_color_continuous(low="red", high="blue", name = "p",
                                        guide=guide_colorbar(reverse=T,
                                                             label.theme = element_text(angle = 90),
                                                             label.vjust = 0.5)) +
                 labs(y = "Term", x = "log2 Odds ratio", size = "size") +
                 theme_bw() +
                 theme(legend.position = "top", legend.text = element_text(size = 9)))
      }
    } else {
      if (adjust){
        return(ggplot(subset, aes(y = reorder(Term, Odds.Ratio), x=log2(Odds.Ratio), size=IN, color=Adjusted.P.value)) +
                 geom_point() +
                 facet_wrap(~contrast, nrow = nrow) +
                 scale_color_continuous(low="red", high="blue", name = "Adjusted.P.value",
                                        guide=guide_colorbar(reverse=T,
                                                             label.theme = element_text(angle = 90),
                                                             label.vjust = 0.5)) +
                 labs(y = "Term", x = "log2 Odds ratio", size = "size") +
                 theme_bw() +
                 theme(legend.position = "top", legend.text = element_text(size = 9)))
      } else {
        return(ggplot(subset, aes(y = reorder(Term, Odds.Ratio), x=log2(Odds.Ratio), size=IN, color=P.value)) +
                 geom_point() +
                 facet_wrap(~contrast, nrow = nrow) +
                 scale_color_continuous(low="red", high="blue", name = "P.value",
                                        guide=guide_colorbar(reverse=T,
                                                             label.theme = element_text(angle = 90),
                                                             label.vjust = 0.5)) +
                 labs(y = "Term", x = "log2 Odds ratio", size = "size") +
                 theme_bw() +
                 theme(legend.position = "top", legend.text = element_text(size = 9)))
      }
    }
  }
}

#gene_names_true<-read_table("R/gene_names.txt",col_names = F)
