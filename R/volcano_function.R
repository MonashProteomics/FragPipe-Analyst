## New function for volcano plot
plot_volcano_new <- function(dep, contrast, label_size = 3, name_col = NULL,
                         add_names = TRUE, adjusted = T, lfc = 1, alpha = 0.05, plot = TRUE, show_gene = F) {
  # Show error if inputs are not the required classes
  if(is.integer(label_size)) label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(contrast),
                          length(contrast) == 1,
                          is.numeric(label_size),
                          length(label_size) == 1,
                          is.logical(add_names),
                          length(add_names) == 1,
                          is.logical(adjusted),
                          length(adjusted) == 1,
                          is.logical(plot),
                          length(plot) == 1)
  
  row_data <- rowData(dep, use.names = FALSE)
  # Show error if inputs do not contain required columns
  if (is.null(name_col)) {
    name_col <- "ID"
  }
  if (any(!c("name", "ID", name_col) %in% colnames(row_data))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun make_unique() to obtain required columns."),
         call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }
  if(length(grep("_significant", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_significant' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun add_rejections() to obtain the required columns."),
         call. = FALSE)
  }
  
  # Show error if an unvalid contrast is given
  if (length(grep(paste("^",contrast,"_diff", sep = ""),
                  colnames(row_data))) == 0) {
    valid_cntrsts <- row_data %>%
      data.frame() %>%
      select(ends_with("_diff")) %>%
      colnames(.) %>%
      gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                paste0(valid_cntrsts, collapse = "', '"),
                                "'")
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n",
         valid_cntrsts_msg,
         call. = FALSE)
  }
  
  # Generate a data.frame containing all info for the volcano plot
  diff <- grep(paste("^",contrast,"_diff", sep = ""),
               colnames(row_data))
  if(adjusted) {
    p_values <- grep(paste("^",contrast, "_p.adj", sep = ""),
                     colnames(row_data))
  } else {
    p_values <- grep(paste("^",contrast, "_p.val", sep = ""),
                     colnames(row_data))
  }
  signif <- abs(row_data[,diff]) >= lfc & row_data[, p_values] <= alpha
  if (!show_gene) {
    if (metadata(dep)$exp == "LFQ") {
      if (metadata(dep)$level != "peptide") {
        df_tmp <- data.frame(diff = row_data[, diff],
                             p_values = -log10(row_data[, p_values]),
                             signif = signif,
                             name = row_data$name,
                             ID = row_data$ID,
                             label = row_data[,name_col])
      } else {
        df_tmp <- data.frame(diff = row_data[, diff],
                             p_values = -log10(row_data[, p_values]),
                             signif = signif,
                             name = row_data$name,
                             ID = row_data$ID,
                             label = row_data[,name_col])
      }
    } else if (metadata(dep)$exp == "TMT") {
      if (metadata(dep)$level == "protein") {
        df_tmp <- data.frame(diff = row_data[, diff],
                             p_values = -log10(row_data[, p_values]),
                             signif = signif,
                             name = row_data$ID)
      } else if (metadata(dep)$level == "gene") {
        df_tmp <- data.frame(diff = row_data[, diff],
                             p_values = -log10(row_data[, p_values]),
                             signif = signif,
                             name = row_data$ID)
      } else { # peptide
        df_tmp <- data.frame(diff = row_data[, diff],
                             p_values = -log10(row_data[, p_values]),
                             signif = signif,
                             name = row_data$Index)
      }
    } else if (metadata(dep)$exp == "DIA") {
      if (metadata(dep)$level != "peptide") {
        df_tmp <- data.frame(diff = row_data[, diff],
                             p_values = -log10(row_data[, p_values]),
                             signif = signif,
                             name = row_data$ID)
      } else {
        df_tmp <- data.frame(diff = row_data[, diff],
                             p_values = -log10(row_data[, p_values]),
                             signif = signif,
                             name = row_data$Index)
      }
    }
  } else {
    if (metadata(dep)$exp == "LFQ") {
      if (metadata(dep)$level != "peptide") {
        df_tmp <- data.frame(diff = row_data[, diff],
                             p_values = -log10(row_data[, p_values]),
                             signif = signif,
                             name = row_data$Gene)
      } else {
        df_tmp <- data.frame(diff = row_data[, diff],
                             p_values = -log10(row_data[, p_values]),
                             signif = signif,
                             name = paste0(row_data$Gene, "_", row_data$Peptide.Sequence))
      }
    } else if (metadata(dep)$exp == "TMT") {
      if (metadata(dep)$level == "protein") {
        df_tmp <- data.frame(diff = row_data[, diff],
                             p_values = -log10(row_data[, p_values]),
                             signif = signif,
                             name = row_data$Gene)
      } else if (metadata(dep)$level == "gene") {
        df_tmp <- data.frame(diff = row_data[, diff],
                             p_values = -log10(row_data[, p_values]),
                             signif = signif,
                             name = row_data$ID)
      } else if (metadata(dep)$level == "peptide") {
        df_tmp <- data.frame(diff = row_data[, diff],
                             p_values = -log10(row_data[, p_values]),
                             signif = signif,
                             name = paste0(row_data$Gene, "_", row_data$Peptide))
      }
    } else if (metadata(dep)$exp == "DIA") {
      if (metadata(dep)$level != "peptide") {
        df_tmp <- data.frame(diff = row_data[, diff],
                           p_values = -log10(row_data[, p_values]),
                           signif = signif,
                           name = row_data$Genes)
      } else {
        df_tmp <- data.frame(diff = row_data[, diff],
                             p_values = -log10(row_data[, p_values]),
                             signif = signif,
                             name = paste0(row_data$Genes, "_", row_data$Stripped.Sequence))
      }
    }
  }
  df <- df_tmp %>% data.frame() %>% filter(!is.na(signif)) %>%
    arrange(signif)

  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)

  # Plot volcano with or without labels
  p <- ggplot(df, aes(diff, p_values)) +
    geom_vline(xintercept = 0) +
    geom_point(aes(col = signif)) +
    geom_text(data = data.frame(), aes(x = c(Inf, -Inf),
                                       y = c(-Inf, -Inf),
                                       hjust = c(1, 0),
                                       vjust = c(-1, -1),
                                       label = c(name1, name2),
                                       size = 5,
                                       fontface = "bold")) +
    labs(title = contrast,
         x = expression(log[2]~"Fold change")) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey"))
  if (add_names) {
    p <- p + ggrepel::geom_text_repel(data = filter(df, signif),
                                      aes(label = name),
                                      size = label_size,
                                      box.padding = unit(0.1, 'lines'),
                                      point.padding = unit(0.1, 'lines'),
                                      segment.size = 0.5)
  }
  if(adjusted) {
    p <- p + labs(y = expression(-log[10]~"Adjusted p-value"))
  } else {
    p <- p + labs(y = expression(-log[10]~"P-value"))
  }
  if(plot) {
   # return(list(p, df))
   # return(df)
    return(p)
  } else {
    df <- df %>%
      select(name, diff, p_values, signif)
    colnames(df)[c(1,2,3)] <- c("protein", "log2_fold_change", "p_value_-log10")
    if(adjusted) {
      colnames(df)[3] <- "adjusted_p_value_-log10"
    }
    return(df)
  }
}


#####====== get_volcano_df =======#######
get_volcano_df <- function(dep, contrast, adjusted = FALSE) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(contrast),
                          length(contrast) == 1)
  
  row_data <- rowData(dep, use.names = FALSE)

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(row_data))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun make_unique() to obtain required columns."),
         call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }
  if(length(grep("_significant", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_significant' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun add_rejections() to obtain the required columns."),
         call. = FALSE)
  }
  
  # Show error if an unvalid contrast is given
  if (length(grep(paste(contrast, "_diff", sep = ""),
                  colnames(row_data))) == 0) {
    valid_cntrsts <- row_data %>%
      data.frame() %>%
      select(ends_with("_diff")) %>%
      colnames(.) %>%
      gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                paste0(valid_cntrsts, collapse = "', '"),
                                "'")
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n",
         valid_cntrsts_msg,
         call. = FALSE)
  }
  
  # Generate a data.frame containing all info for the volcano plot
  diff <- grep(paste(contrast, "_diff", sep = ""),
               colnames(row_data))
  if(adjusted) {
    p_values <- grep(paste(contrast, "_p.adj", sep = ""),
                     colnames(row_data))
  } else {
    p_values <- grep(paste(contrast, "_p.val", sep = ""),
                     colnames(row_data))
  }
  signif <- grep(paste(contrast, "_significant", sep = ""),
                 colnames(row_data))
  df_tmp <- data.frame(diff = row_data[, diff],
                       p_values = -log10(row_data[, p_values]),
                       signif = row_data[, signif],
                       ID = row_data$ID,
                       name = row_data$name)
  df<- df_tmp %>% data.frame() %>% filter(!is.na(signif)) %>%
    arrange(signif)
  
 return(df)
}

### Function to plot intensities of individual proteins
plot_protein<-function(dep, protein, type, id="ID"){
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(protein),
                          is.character(type))

  subset<-dep[protein]
  
  df_reps <- data.frame(assay(subset), check.names = F) %>%
    rownames_to_column() %>%
    gather(ID, val, -rowname) %>%
    left_join(., data.frame(colData(subset)), by = c("ID"=id))
  
  df_reps$rowname <- parse_factor(as.character(df_reps$rowname), levels = protein)
  
  df_CI<- df_reps %>%
    group_by(condition, rowname) %>%
    summarize(mean = mean(val, na.rm = TRUE),
              sd = sd(val, na.rm = TRUE),
              n = n()) %>%
    mutate(error = qnorm(0.975) * sd / sqrt(n),
           CI.L = mean - error,
           CI.R = mean + error) %>%
    as.data.frame()
  df_CI$rowname <- parse_factor(as.character(df_CI$rowname), levels = protein)

  df_reps$condition <- as.factor(df_reps$condition)
  df_reps <- df_reps[!is.na(df_reps$val),]
  df_reps$replicate <- as.character(df_reps$replicate)
  if(type=="violin"){
    if (max(df_reps$replicate) == 1){
      p <- plot_ly(df_reps,
                   x = ~rowname,
                   y = ~val,
                   color = ~condition,
                   text = ~sample_name,
                   hoverinfo = "text",
                   type = "violin") %>%
        plotly::layout(violinmode = "group",
                       xaxis = list(title = ''),
                       yaxis = list(title = 'Abundance'), showlegend=T)
    } else {
      p <- plot_ly(df_reps,
                   x = ~rowname,
                   y = ~val,
                   color = ~condition,
                   text = ~sample_name,
                   hoverinfo = "text",
                   type = "violin") %>%
        plotly::layout(violinmode = "group",
                       xaxis = list(title = ''),
                       yaxis = list(title = 'Abundance'), showlegend=T)
    }
    return(p)
  } else if(type=="boxplot"){
    if (max(df_reps$replicate) == 1){
      p <- plot_ly(df_reps,
                   x = ~rowname,
                   y = ~val,
                   color = ~condition,
                   text = ~sample_name,
                   hoverinfo = "text",
                   type = "box",
                   boxpoints = "all", jitter = 0.3, pointpos = 0) %>%
        plotly::layout(boxmode = "group",
                       xaxis = list(title = ''),
                       yaxis = list(title = 'Abundance'), showlegend=T)
    } else {
      p <- plot_ly(df_reps,
                   x = ~rowname,
                   y = ~val,
                   color = ~condition,
                   text = ~sample_name,
                   hoverinfo = "text",
                   type = "box",
                   boxpoints = "all", jitter = 0.3, pointpos = 0) %>%
        plotly::layout(boxmode = "group",
                       xaxis = list(title = ''),
                       yaxis = list(title = 'Abundance'), showlegend=T)
    }
    return(p)
  }
}

plot_volcano_mod <- function(dep, contrast, label_size = 3, name_col = NULL,
                             add_names = TRUE, adjusted = FALSE, plot = TRUE) {
  # Show error if inputs are not the required classes
  if(is.integer(label_size)) label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(contrast),
                          length(contrast) == 1,
                          is.numeric(label_size),
                          length(label_size) == 1,
                          is.logical(add_names),
                          length(add_names) == 1,
                          is.logical(adjusted),
                          length(adjusted) == 1,
                          is.logical(plot),
                          length(plot) == 1)
  
  row_data <- rowData(dep, use.names = FALSE)

  # Show error if inputs do not contain required columns
  if (is.null(name_col)) {
    name_col <- "ID"
  }
  if (any(!c("name", "ID", name_col) %in% colnames(row_data))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun make_unique() to obtain required columns."),
         call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }
  if(length(grep("_significant", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_significant' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun add_rejections() to obtain the required columns."),
         call. = FALSE)
  }
  
  # Show error if an unvalid contrast is given
  if (length(grep(paste("^",contrast, "_diff", sep = ""),
                  colnames(row_data))) == 0) {
    valid_cntrsts <- row_data %>%
      data.frame() %>%
      select(ends_with("_diff")) %>%
      colnames(.) %>%
      gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                paste0(valid_cntrsts, collapse = "', '"),
                                "'")
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n",
         valid_cntrsts_msg,
         call. = FALSE)
  }
  
  # Generate a data.frame containing all info for the volcano plot
  diff <- grep(paste("^",contrast, "_diff", sep = ""),
               colnames(row_data))
  if(adjusted) {
    p_values <- grep(paste("^",contrast, "_p.adj", sep = ""),
                     colnames(row_data))
  } else {
    p_values <- grep(paste("^", contrast, "_p.val", sep = ""),
                     colnames(row_data))
  }
  signif <- grep(paste("^",contrast, "_significant", sep = ""),
                 colnames(row_data))
  df <- data.frame(x = row_data[, diff],
                   y = -log10(row_data[, p_values]),
                   p.val = row_data[, p_values],
                   significant = row_data[, signif],
                   name = row_data$name,
                   ID = row_data$ID,
                   label = row_data[,name_col]) %>%
    filter(!is.na(significant)) %>%
    arrange(significant)
  
  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)

  # Plot volcano with or without labels
  p <- ggplot(df, aes(x, y)) +
    geom_vline(xintercept = 0) +
    geom_point(aes(col = significant)) +
    geom_text(data = data.frame(), aes(x = c(Inf, -Inf),
                                       y = c(-Inf, -Inf),
                                       hjust = c(1, 0),
                                       vjust = c(-1, -1),
                                       label = c(name1, name2),
                                       size = 5,
                                       fontface = "bold")) +
    labs(title = contrast,
         x = expression(log[2]~"Fold change")) +
    theme_DEP1() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey"))
  if (add_names) {
    p <- p + ggrepel::geom_text_repel(data = filter(df, significant),
                                      aes(label = label),
                                      size = label_size,
                                      box.padding = unit(0.1, 'lines'),
                                      point.padding = unit(0.1, 'lines'),
                                      segment.size = 0.5)
  }
  if(adjusted) {
    p <- p + labs(y = expression(-log[10]~"Adjusted p-value"))
  } else {
    p <- p + labs(y = expression(-log[10]~"P-value"))
  }
  if(plot) {
    return(p)
  } else {
    df <- df %>%
      select(name, x, y, significant) %>%
      arrange(desc(x))
    colnames(df)[c(1,2,3)] <- c("protein", "log2_fold_change", "p_value_-log10")
    if(adjusted) {
      colnames(df)[3] <- "adjusted_p_value_-log10"
    }
    return(df)
  }
}

