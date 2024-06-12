
# make make.unique generates new names strat from 1 rather than nothing
make.unique.2 = function(x, sep='.'){
  ave(x, x, FUN=function(a){if(length(a) > 1){paste(a, 1:length(a), sep=sep)} else {a}})
}

matrixplot_modify<-function(data, mapping, pts=list(), smt=list(), ...){
    ggplot(data = data, mapping = mapping, ...) + 
        do.call(geom_point, pts) +
        do.call(geom_smooth, smt) 
}


coef_variation<-function(x){
  coef=sd(x)/mean(x)
}

#### Plot CVs

plot_cvs <- function(se, id="ID", scale=T, check.names=T) {
  
  ## backtransform data
  untransformed_intensity<- 2^(assay(se))
  exp_design<-colData(se)

  ### merge untransformed to exp design and calculate cvs
  if (id == "ID") {
    cvs_group<- untransformed_intensity %>% data.frame() %>%
      tibble::rownames_to_column() %>%
      tidyr::gather("ID", "Intensity", -rowname) %>%
      dplyr::left_join(.,data.frame(exp_design), by="ID") %>%
      dplyr::group_by(rowname,condition) %>%
      dplyr::summarise(cvs=coef_variation(Intensity)) %>%
      dplyr::group_by(condition)%>%
      dplyr::mutate(condition_median=median(cvs))
  } else {
    cvs_group<- untransformed_intensity %>% data.frame(check.names=check.names) %>%
      tibble::rownames_to_column() %>%
      tidyr::gather("ID", "Intensity", -rowname) %>%
      dplyr::left_join(.,data.frame(exp_design), by=c("ID"=id)) %>%
      dplyr::group_by(rowname,condition) %>%
      dplyr::summarise(cvs=coef_variation(Intensity)) %>%
      dplyr::group_by(condition)%>%
      dplyr::mutate(condition_median=median(cvs, na.rm = T))
  }
  if (scale) {
    p1 <- ggplot(cvs_group, aes(cvs, color=condition, fill=condition)) +
      geom_histogram(alpha=.5, bins= 20, show.legend = FALSE) +
      facet_wrap(~condition) +
      geom_vline(aes(xintercept=condition_median, group=condition),
                 color='grey40',
                 linetype="dashed") +
      scale_x_continuous(labels = scales::percent, limits=c(0, 1)) +
      labs(title= 'Sample Coefficient of Variation', x="Coefficient of Variation", y="Count") +
      theme_DEP2() +
      theme(plot.title = element_text(hjust = 0.5,face = "bold"))
  } else {
    p1 <- ggplot(cvs_group, aes(cvs, color=condition, fill=condition)) +
      geom_histogram(alpha=.5, bins= 20, show.legend = FALSE) +
      facet_wrap(~condition) +
      geom_vline(aes(xintercept=condition_median, group=condition),
                 color='grey40',
                 linetype="dashed") +
      scale_x_continuous(labels = scales::percent) +
      labs(title= 'Sample Coefficient of Variation', x="Coefficient of Variation", y="Count") +
      theme_DEP2() +
      theme(plot.title = element_text(hjust = 0.5,face = "bold"))
  }
  

  p <- p1 + geom_text(aes(x=0.9,
                          y=max(ggplot_build(p1)$data[[1]]$ymax*1.1), 
                     label=paste0("Median =",round(condition_median,2)*100,"%",by="")),
                 show.legend = FALSE, size=4)
  return(p)
}


#### Get individual clusters from heatmap
get_cluster_heatmap <- function(dep, type = c("contrast", "centered"),
                                kmeans = FALSE, k = 6,
                                col_limit = 6, indicate = NULL,
                                alpha = 0.01, lfc = 1,
                                clustering_distance = c("euclidean", "maximum", "manhattan", "canberra",
                                                        "binary", "minkowski", "pearson", "spearman", "kendall", "gower"),
                                row_font_size = 6, col_font_size = 10, plot = TRUE, exp="LFQ", show_row_names=F, ...) {
  
  # Show error if inputs are not the required classes
  if(is.integer(k)) k <- as.numeric(k)
  if(is.integer(col_limit)) col_limit <- as.numeric(col_limit)
  if(is.integer(row_font_size)) row_font_size <- as.numeric(row_font_size)
  if(is.integer(col_font_size)) col_font_size <- as.numeric(col_font_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(type),
                          is.logical(kmeans),
                          is.numeric(k),
                          length(k) == 1,
                          is.numeric(col_limit),
                          length(col_limit) == 1,
                          is.numeric(row_font_size),
                          length(row_font_size) == 1,
                          is.numeric(col_font_size),
                          length(col_font_size) == 1,
                          is.logical(plot),
                          length(plot) == 1)
  
  # Show error if inputs do not contain required columns
  type <- match.arg(type)
  clustering_distance <- match.arg(clustering_distance)
  
  # Extract row and col data
  row_data <- rowData(dep)
  col_data <- colData(dep) %>%
    as.data.frame()
  
  # Show error if inputs do not contain required columns
  if(any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '",
                deparse(substitute(dep)), "'"),
         call. = FALSE)
  }
  if(length(grep("_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }
  if(!"significant" %in% colnames(row_data)) {
    stop(paste0("'significant' column is not present in '",
                deparse(substitute(dep)),
                "'.\nRun add_rejections() to obtain the required column."),
         call. = FALSE)
  }
  
  # Heatmap annotation
  if(!is.null(indicate) & type == "contrast") {
    warning("Heatmap annotation only applicable for type = 'centered'",
            call. = FALSE)
  }
  if(!is.null(indicate) & type == "centered") {
    ha1 <- get_annotation(dep, indicate)
  } else {
    ha1 <- NULL
  }
  
  conditions <- gsub("_diff", "",colnames(row_data)[grepl("_diff", colnames(row_data))])
  cols_p <- paste0(conditions, "_p.adj")
  cols_lfc <- paste0(conditions, "_diff")
  p <- as.matrix(row_data[,cols_p]) <= alpha
  lfc <- abs(as.matrix(row_data[,cols_lfc])) >= lfc
  p[is.na(p)] <- FALSE
  lfc[is.na(p)] <- FALSE
  colnames(p) <- conditions
  colnames(lfc) <- conditions
  combined_rejections <- p
  colnames(combined_rejections) <- conditions
  combined_rejections <- p & lfc
  filtered <- dep[apply(combined_rejections, 1, any), ]
  
  if (nrow(filtered) == 0){
    return(ggplot() +
      annotate("text", x = 4, y = 25, size=8, label = "No differential expressed genes available for the heatmap") +
      theme_void())
  }

  # Check for missing values
  if(any(is.na(assay(filtered)))) {
    warning("Missing values in '", deparse(substitute(dep)), "'. ",
            "Using clustering_distance = 'gower'",
            call. = FALSE)
    clustering_distance <- "gower"
    obs_NA <- TRUE
  } else {
    obs_NA <- FALSE
  }
  
  # Get centered intensity values ('centered')
  if(type == "centered") {
    rowData(filtered)$mean <- rowMeans(assay(filtered), na.rm = TRUE)
    df <- assay(filtered) - rowData(filtered)$mean
  }
  # Get contrast fold changes ('contrast')
  if(type == "contrast") {
    df <- rowData(filtered) %>%
      data.frame(.) %>%
      column_to_rownames(var = "name") %>%
      select(dplyr::ends_with("_diff"))
    colnames(df) <-
      gsub("_diff", "", colnames(df)) %>%
      gsub("_vs_", " vs ", .)
  }
  
  # Facultative kmeans clustering
  if(kmeans & obs_NA) {
    warning("Cannot perform kmeans clustering with missing values",
            call. = FALSE)
    kmeans <- FALSE
  }
  if(kmeans & !obs_NA) {
    set.seed(1)
    df_kmeans <- kmeans(df, k)
    if(type == "centered") {
      # Order the k-means clusters according to the maximum fold change
      # in all samples averaged over the proteins in the cluster
      order <- data.frame(df) %>%
        cbind(., cluster = df_kmeans$cluster) %>%
        dplyr::mutate(row = apply(.[, seq_len(ncol(.) - 1)], 1, function(x) max(x))) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarize(index = sum(row)/n()) %>%
        dplyr::arrange(desc(index)) %>%
        dplyr::pull(cluster) %>%
        match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
    if(type == "contrast") {
      # Order the k-means clusters according to their average fold change
      order <- cbind(df, cluster = df_kmeans$cluster) %>%
        dplyr::gather(condition, diff, -cluster) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarize(row = mean(diff)) %>%
        dplyr::arrange(desc(row)) %>%
        dplyr::pull(cluster) %>%
        match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
  }
  
  if(ncol(df) == 1) {
    col_clust = FALSE
  } else {
    col_clust = TRUE
  }
  if(nrow(df) == 1) {
    row_clust = FALSE
  } else {
    row_clust = TRUE
  }
  if(clustering_distance == "gower") {
    clustering_distance <- function(x) {
      dist <- cluster::daisy(x, metric = "gower")
      dist[is.na(dist)] <- max(dist, na.rm = TRUE)
      return(dist)
    }
  }
  
  # Legend info
  legend <- ifelse(type == "contrast",
                   "log2 Fold change",
                   "log2 Centered intensity")

  # use sample name for the heatmap
  temp <- colData(filtered)
  rownames(temp) <- temp$label
  colnames(df) <- temp[colnames(df), "sample_name"]

  # Heatmap
  ht1 = Heatmap(df,
                col = circlize::colorRamp2(
                  seq(-col_limit, col_limit, (col_limit/5)),
                  rev(RColorBrewer::brewer.pal(11, "RdBu"))),
                split = if(kmeans) {df_kmeans$cluster} else {NULL},
                show_row_names = show_row_names,
                cluster_rows = col_clust,
                cluster_columns = row_clust,
                row_names_side = "left",
                column_names_side = "top",
                clustering_distance_rows = clustering_distance,
                clustering_distance_columns = clustering_distance,
                heatmap_legend_param = list(color_bar = "continuous",
                                            legend_direction = "horizontal",
                                            legend_width = unit(5, "cm"),
                                            title_position = "lefttop"),
                name = legend,
                row_names_gp = gpar(fontsize = row_font_size),
                column_names_gp = gpar(fontsize = col_font_size),
                top_annotation = ha1,
                ...)
  # return (row_order(ht1))
  # Return data.frame
  p <- draw(ht1, heatmap_legend_side = "top")
  row_clusters<- row_order(ht1)
  #mat<-as.matrix(df)
  
  # for (i in 1:length(row_clusters)){
  #   if (i==1){
  #     clu <-t(t(row.names(ht1[row_clusters[[i]],])))
  #     out <-cbind (clu, paste("cluster", i, sep=""))
  #     colnames(out)<- c("ProteinID", "Cluster")
  #   }
  #   else{
  #     clu <- t(t(row.names(ht1[row_clusters[[i]],])))
  #     clu <- cbind(clu, paste("cluster", i, sep = ""))
  #     out <- cbind(out, clu)
  #   }
  # }
  heatmap_list <- list(p, row_clusters)
  return(heatmap_list)
}

# Internal function to get ComplexHeatmap::HeatmapAnnotation object
get_annotation <- function(dep, indicate) {
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(indicate))
  
  # Check indicate columns
  col_data <- colData(dep) %>%
    as.data.frame()
  columns <- colnames(col_data)
  if(all(!indicate %in% columns)) {
    stop("'",
         paste0(indicate, collapse = "' and/or '"),
         "' column(s) is/are not present in ",
         deparse(substitute(dep)),
         ".\nValid columns are: '",
         paste(columns, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  if(any(!indicate %in% columns)) {
    indicate <- indicate[indicate %in% columns]
    warning("Only used the following indicate column(s): '",
            paste0(indicate, collapse = "', '"),
            "'")
  }
  
  # Get annotation
  anno <- dplyr::select(col_data, indicate)
  
  # Annotation color
  names <- colnames(anno)
  anno_col <- vector(mode="list", length=length(names))
  names(anno_col) <- names
  for(i in names) {
    var = anno[[i]] %>% unique() %>% sort()
    if(length(var) == 1) {
      cols <- c("black")
    } else if(length(var) == 2) {
      cols <- c("orangered", "cornflowerblue")
    } else if(length(var) < 7 & length(var) > 2) {
      cols <- RColorBrewer::brewer.pal(length(var), "Pastel1")
    } else if(length(var) <= 12) {
      cols <- RColorBrewer::brewer.pal(length(var), "Set3")
    } else {
      cols <- colorRampPalette(brewer.pal(12, "Set3"))(length(var))
    }
    names(cols) <- var
    anno_col[[i]] <-  cols
    names(cols) <- var
    anno_col[[i]] <-  cols
  }
  
  # HeatmapAnnotation object
  ComplexHeatmap::HeatmapAnnotation(df = anno,
                    col = anno_col,
                    show_annotation_name = TRUE)
}


#### ===== limma BH FDR ===== #####

test_limma <- function(se, type = c("control", "all", "manual"),
                       control = NULL, test = NULL,
                       design_formula = formula(~ 0 + condition),
                       paired = FALSE) {
  #require("dplyr", "tidyr", "purrr")
  
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type),
                          class(design_formula) == "formula")
  if (paired == FALSE){
    design_formula <- design_formula
  }else{
    design_formula<-formula(~ 0 + condition + replicate)
  }
  
  
  # Show error if inputs do not contain required columns
  type <- match.arg(type)
  
  col_data <- colData(se)
  raw <- assay(se)
  
  if(any(!c("name", "ID") %in% colnames(rowData(se)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if(any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)), "'")
  }
  
  if(!is.null(control)) {
    # Show error if control input is not valid
    assertthat::assert_that(is.character(control),
                            length(control) == 1)
    if(!control %in% unique(col_data$condition)) {
      stop("run test_diff() with a valid control.\nValid controls are: '",
           paste0(unique(col_data$condition), collapse = "', '"), "'",
           call. = FALSE)
    }
  }
  
  # variables in formula
  variables <- terms.formula(design_formula) %>%
    attr(., "variables") %>%
    as.character() %>%
    .[-1]
  
  # Throw error if variables are not col_data columns
  if(any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }
  if(variables[1] != "condition") {
    stop("first factor of 'design_formula' should be 'condition'")
  }
  
  # Obtain variable factors
  for(var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }
  
  # Make an appropriate design matrix
  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))
  
  # Generate contrasts to be tested
  # Either make all possible combinations ("all"),
  # only the contrasts versus the control sample ("control") or
  # use manual contrasts
  conditions <- as.character(unique(condition))
  if(type == "all") {
    # All possible combinations
    cntrst <- apply(utils::combn(conditions, 2), 2, paste, collapse = " - ")
    
    if(!is.null(control)) {
      # Make sure that contrast containing
      # the control sample have the control as denominator
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if(length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>%
          gsub(paste(control, "- ", sep = " "), "", .) %>%
          paste(" - ", control, sep = "")
      }
    }
    
  }
  if(type == "control") {
    # Throw error if no control argument is present
    if(is.null(control))
      stop("run test_diff(type = 'control') with a 'control' argument")
    
    # Make contrasts
    cntrst <- paste(conditions[!conditions %in% control],
                    control,
                    sep = " - ")
  }
  if(type == "manual") {
    # Throw error if no test argument is present
    if(is.null(test)) {
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(test))
    
    if(any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("run test_diff() with valid contrasts in 'test'",
           ".\nValid contrasts should contain combinations of: '",
           paste0(conditions, collapse = "', '"),
           "', for example '", paste0(conditions[1], "_vs_", conditions[2]),
           "'.", call. = FALSE)
    }
    
    cntrst <- gsub("_vs_", " - ", test)
    
  }
  # Print tested contrasts
  message("Tested contrasts: ",
          paste(gsub(" - ", "_vs_", cntrst), collapse = ", "))
  
  # Test for differential expression by empirical Bayes moderation
  # of a linear model on the predefined contrasts
  fit <- lmFit(raw, design = design)
  made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
  contrast_fit <- contrasts.fit(fit, made_contrasts)
  
  if(any(is.na(raw))) {
    for(i in cntrst) {
      covariates <- strsplit(i, " - ") %>% unlist
      single_contrast <- makeContrasts(contrasts = i, levels = design[, covariates])
      single_contrast_fit <- contrasts.fit(fit[, covariates], single_contrast)
      contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 1]
      contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 1]
    }
  }
  
  eB_fit <- eBayes(contrast_fit)
  
  # function to retrieve the results of
  # the differential expression test using 'fdrtool'
  retrieve_fun <- function(comp, fit = eB_fit){
    res <- topTable(fit, sort.by = "t", adjust.method="BH", coef = comp,
                    number = Inf, confint = TRUE)
    # res <- res[!is.na(res$t),]
    #fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    # res$qval <- res$adj.P.Value
    #res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- tibble::rownames_to_column(res)
    return(res)
  }
  
  #limma_res<- topTable(eB_fit, sort.by = 'B', adjust.method="BH", coef = cntrst, number = Inf, confint = T )
  # limma_res$comparison <- rep(cntrst, dim(limma_res)[1])
  #limma_res <- rownames_to_column(limma_res)
  # Retrieve the differential expression test results
  limma_res <- purrr::map_df(cntrst, retrieve_fun)
  
  # Select the logFC, CI and qval variables
  table <- limma_res %>%
    dplyr::select(rowname, logFC, CI.L, CI.R, P.Value, adj.P.Val, comparison) %>%
    dplyr::mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
    tidyr::gather(variable, value, -c(rowname,comparison)) %>%
    dplyr::mutate(variable = dplyr::recode(variable, logFC = "diff", P.Value = "p.val", adj.P.Val = "p.adj")) %>%
    tidyr::unite(temp, comparison, variable) %>%
    tidyr::spread(temp, value)
  rowData(se) <- merge(rowData(se), table,
                       by.x = "name", by.y = "rowname", all.x = TRUE)
  return(se)
  #return(table)
}

get_results_proteins <- function(dep) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"))
  
  row_data <- rowData(dep)
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
  
  # Obtain average protein-centered enrichment values per condition
  row_data$mean <- rowMeans(assay(dep), na.rm = TRUE)
  centered <- assay(dep) - row_data$mean
  centered <- data.frame(centered) %>%
    tibble::rownames_to_column() %>%
    tidyr::gather(ID, val, -rowname) %>%
    dplyr::left_join(., data.frame(colData(dep)), by = c("ID"="label"))

  centered <- dplyr::group_by(centered, rowname, condition) %>%
    dplyr::summarize(val = mean(val, na.rm = TRUE)) %>%
    dplyr::mutate(val = signif(val, digits = 3)) %>%
    tidyr::spread(condition, val)
  colnames(centered)[2:ncol(centered)] <-
    paste(colnames(centered)[2:ncol(centered)], "_centered", sep = "")
  
  # Obtain average enrichments of conditions versus the control condition
  ratio <- as.data.frame(row_data) %>%
    #tibble::column_to_rownames("name") %>%
    dplyr::select(dplyr::ends_with("diff")) %>%
    signif(., digits = 3) %>%
    tibble::rownames_to_column()
  colnames(ratio)[2:ncol(ratio)] <-
    gsub("_diff", "_log2 fold change", colnames(ratio)[2:ncol(ratio)])
 # df <- left_join(ratio, centered, by = "rowname")
  
  # Select the adjusted p-values and significance columns
  pval <- as.data.frame(row_data) %>%
    #tibble::column_to_rownames("name") %>%
    dplyr::select(dplyr::ends_with("p.val"),
                  dplyr::ends_with("p.adj"),
                  dplyr::ends_with("significant")) %>%
   tibble::rownames_to_column()
  pval[, grep("p.adj", colnames(pval))] <-
    pval[, grep("p.adj", colnames(pval))] %>%
    signif(digits = 3)
  pval[, grep("p.val", colnames(pval))] <-
    pval[, grep("p.val", colnames(pval))] %>%
    signif(digits = 3)
  
  # Join into a results table
  if (metadata(dep)$exp == "LFQ") {
    if(metadata(dep)$level == "protein") {
      ids <- as.data.frame(row_data) %>% dplyr::select(ID, name)
      table <- dplyr::left_join(ids,ratio, by=c("ID"="rowname"))
      table <- dplyr::left_join(table, pval, by = c("ID" = "rowname"))
      table <-as.data.frame(row_data) %>%
        dplyr::select(ID, imputed, num_NAs, Description) %>%
        dplyr::left_join(table, ., by = "ID")
      table <- table %>% dplyr::arrange(desc(significant))
      colnames(table)[1] <- c("Protein ID")
      colnames(table)[2] <- c("Gene Name")
    } else if(metadata(dep)$level == "peptide") {
      ids <- as.data.frame(row_data) %>% dplyr::select(ID, name, Gene)
      table <- dplyr::left_join(ids, ratio, by=c("ID"="rowname"))
      table <- dplyr::left_join(table, pval, by = c("ID" = "rowname"))
      table <- as.data.frame(row_data) %>%
        dplyr::select(ID, imputed, num_NAs, Description) %>%
        dplyr::left_join(table, ., by = "ID")
      table <- table %>% dplyr::arrange(desc(significant))
      colnames(table)[1] <- c("Index")
      colnames(table)[2] <- c("Protein ID")
      colnames(table)[3] <- c("Gene Name")
    }
  } else if (metadata(dep)$exp == "TMT") {
    if(metadata(dep)$level == "gene") {
      ids <- as.data.frame(row_data) %>% dplyr::select(ID, name)
      table <- dplyr::left_join(ids,ratio, by=c("ID"="rowname"))
      table <- dplyr::left_join(table, pval, by = c("ID" = "rowname"))
      table <- as.data.frame(row_data) %>%
        dplyr::select(ID, imputed, num_NAs) %>%
        dplyr::left_join(table, ., by = "ID")
      table <- table %>% dplyr::arrange(desc(significant))
      colnames(table)[1] <- c("Gene Name")
      colnames(table)[2] <- c("Protein ID")
    } else if (metadata(dep)$level == "protein") {
      ids <- as.data.frame(row_data) %>% dplyr::select(name, Gene)
      table <- dplyr::left_join(ids, ratio, by=c("name"="rowname"))
      table <- dplyr::left_join(table, pval, by = c("name" = "rowname"))
      table <- as.data.frame(row_data) %>%
        dplyr::select(name, imputed, num_NAs) %>%
        dplyr::left_join(table, ., by = "name")
      table <- table %>% dplyr::arrange(desc(significant))
      colnames(table)[1] <- c("Protein ID")
      colnames(table)[2] <- c("Gene Name")
    } else if (metadata(dep)$level == "peptide") {
      ids <- as.data.frame(row_data) %>% dplyr::select(ID, ProteinID, Gene)
      table <- dplyr::left_join(ids, ratio, by=c("ID"="rowname"))
      table <- dplyr::left_join(table, pval, by = c("ID" = "rowname"))
      table <- as.data.frame(row_data) %>%
        dplyr::select(ID, imputed, num_NAs) %>%
        dplyr::left_join(table, ., by = "ID")
      table <- table %>% dplyr::arrange(desc(significant))
      colnames(table)[1] <- c("Index")
      colnames(table)[2] <- c("Protein ID")
      colnames(table)[3] <- c("Gene Name")
    }
  } else if (metadata(dep)$exp == "DIA") {
    if (metadata(dep)$level == "protein") {
      ids <- as.data.frame(row_data) %>% dplyr::select(ID, name)
      table <- dplyr::left_join(ids,ratio, by=c("ID"="rowname"))
      table <- dplyr::left_join(table, pval, by = c("ID" = "rowname"))
      table <- as.data.frame(row_data) %>%
        dplyr::select(ID, imputed, num_NAs) %>%
        dplyr::left_join(table, ., by = "ID")
      table <- table %>% dplyr::arrange(desc(significant))
      colnames(table)[1] <- c("Protein ID")
      colnames(table)[2] <- c("Gene Name")
    } else if (metadata(dep)$level == "peptide") {
      ids <- as.data.frame(row_data) %>% dplyr::select(ID, name, Genes)
      table <- dplyr::left_join(ids,ratio, by=c("ID"="rowname"))
      table <- dplyr::left_join(table, pval, by = c("ID" = "rowname"))
      table <- as.data.frame(row_data) %>%
        dplyr::select(ID, imputed, num_NAs) %>%
        dplyr::left_join(table, ., by = "ID")
      table <- table %>% dplyr::arrange(desc(significant))
      colnames(table)[1] <- c("Index")
      colnames(table)[2] <- c("Protein ID")
      colnames(table)[3] <- c("Gene Name")
    }
  }
  return(table)
}


#### ==== get prefix function 

get_prefix <- function(words) {
  # Show error if input is not the required class
  assertthat::assert_that(is.character(words))
  
  # Show error if 'words' contains 1 or less elements
  if(length(words) <= 1) {
    stop("'words' should contain more than one element")
  }
  # Show error if 'words' contains NA
  if(any(is.na(words))) {
    stop("'words' contains NAs")
  }
  
  # Truncate words to smallest name
  minlen <- min(nchar(words))
  truncated <- substr(words, 1, minlen)
  
  # Show error if one of the elements is shorter than one character
  if(minlen < 1) {
    stop("At least one of the elements is too short")
  }
  
  # Get identifical characters
  mat <- data.frame(strsplit(truncated, ""), stringsAsFactors = FALSE)
  identical <- apply(mat, 1, function(x) length(unique(x)) == 1)
  
  # Obtain the longest common prefix
  prefix <- as.logical(cumprod(identical))
  paste(mat[prefix, 1], collapse = "")
}

#### ===== delete prefix function

delete_prefix <- function(words) {
  # Get prefix
  prefix <- get_prefix(words)
  # Delete prefix from words
  gsub(paste0("^", prefix), "", words)
}
