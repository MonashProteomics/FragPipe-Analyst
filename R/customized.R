# customized functions in addition to LFQ-Analyst

# original: make_se from
# https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/functions.R
#' Data.frame to SummarizedExperiment object
#' conversion using an experimental design
#'
#' \code{make_se_customized} creates a SummarizedExperiment object
#' based on two data.frames: the protein table and experimental design.
#'
#' @param proteins_unique Data.frame,
#' Protein table with unique names annotated in the 'name' column
#' (output from \code{\link{make_unique}()}).
#' @param columns Integer vector,
#' Column numbers indicating the columns containing the assay data.
#' @param expdesign Data.frame,
#' Experimental design with 'label', 'condition'
#' and 'replicate' information.
#' See \code{\link{UbiLength_ExpDesign}} for an example experimental design.
#' @return A SummarizedExperiment object
#' with log2-transformed values.
#' @examples
#' # Load example
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#' @export
make_se_customized <- function(proteins_unique, columns, expdesign, log2transform=F, exp="LFQ", lfq_type=NULL, level=NULL) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins_unique),
                          is.integer(columns),
                          is.data.frame(expdesign))
  
  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)),
         "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
         "are not present in the experimental design",
         call. = FALSE)
  }
  if(any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric",
         "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }
  
  # If input is a tibble, convert to data.frame
  if(tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  if(tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)
  
  # Select the assay data
  rownames(proteins_unique) <- proteins_unique$ID
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  if (log2transform) {
    raw <- log2(raw)
  }
  # Generate the colData from the experimental design
  # and match these with the assay data
  # expdesign <- mutate(expdesign, condition = make.names(condition)) %>%
  #   unite(ID, condition, replicate, remove = FALSE)
  # rownames(expdesign) <- expdesign$ID
  rownames(expdesign) <- expdesign$label
  # print(expdesign)
  # print(colnames(raw))
  matched <- match(make.names(expdesign$label),
                   make.names(colnames(raw)))

  if(any(is.na(matched))) {
    stop("None of the labels in the experimental design match ",
         "with column names in 'proteins_unique'",
         "\nRun make_se() with the correct labels in the experimental design",
         "and/or correct columns specification")
  }
  
  colnames(raw)[matched] <- expdesign$label
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]

  # Select the rowData
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$ID

  # Generate the SummarizedExperiment object
  se <- SummarizedExperiment(assays = as.matrix(raw),
                             colData = expdesign,
                             rowData = row_data,
                             metadata = list("exp"=exp,
                                             "lfq_type"=lfq_type,
                                             "level"=level,
                                             "log2transform"=log2transform))
  
  return(se)
}

# original: https://github.com/arnesmits/DEP/blob/master/R/functions.R
test_diff_customized <- function(se, type = c("control", "all", "manual"),
                                 control = NULL, test = NULL,
                                 design_formula = formula(~ 0 + condition)) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type),
                          class(design_formula) == "formula")
  
  # Show error if inputs do not contain required columns
  type <- match.arg(type)
  
  col_data <- colData(se)
  raw <- assay(se)
  
  if(any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
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
  # if(variables[1] != "condition") {
  #   stop("first factor of 'design_formula' should be 'condition'")
  # }
  
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
  conditions <- as.character(unique(col_data$condition))
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
    
    # if(any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
    #   stop("run test_diff() with valid contrasts in 'test'",
    #        ".\nValid contrasts should contain combinations of: '",
    #        paste0(conditions, collapse = "', '"),
    #        "', for example '", paste0(conditions[1], "_vs_", conditions[2]),
    #        "'.", call. = FALSE)
    # }
    
    cntrst <- gsub("_vs_", " - ", test)
    
  }
  # Print tested contrasts
  message("Tested contrasts: ",
          paste(gsub(" - ", "_vs_", cntrst), collapse = ", "))
  
  # Test for differential expression by empirical Bayes moderation
  # of a linear model on the predefined contrasts
  # print(design)
  # print(cntrst)
  fit <- lmFit(raw, design = design)
  made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
  # print(made_contrasts)
  contrast_fit <- contrasts.fit(fit, made_contrasts)
  
  if (type != "manual") {
    if(any(is.na(raw))) {
      for(i in cntrst) {
        covariates <- strsplit(i, " - ") %>% unlist
        single_contrast <- makeContrasts(contrasts = i, levels = design[, covariates])
        single_contrast_fit <- contrasts.fit(fit[, covariates], single_contrast)
        contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 1]
        contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 1]
      }
    }
  }
  eB_fit <- eBayes(contrast_fit)
  
  # function to retrieve the results of
  # the differential expression test using 'fdrtool'
  retrieve_fun <- function(comp, fit = eB_fit){
    res <- topTable(fit, sort.by = "t", coef = comp,
                    number = Inf, confint = TRUE)
    res <- res[!is.na(res$t),]
    fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- rownames_to_column(res)
    return(res)
  }
  
  # Retrieve the differential expression test results
  limma_res <- map_df(cntrst, retrieve_fun)

  # Select the logFC, CI and qval variables
  table <- limma_res %>%
    select(rowname, logFC, CI.L, CI.R, P.Value, qval, comparison) %>%
    mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
    gather(variable, value, -c(rowname,comparison)) %>%
    mutate(variable = recode(variable, logFC = "diff", P.Value = "p.val", qval = "p.adj")) %>%
    unite(temp, comparison, variable) %>%
    spread(temp, value)
  rowData(se) <- as.data.frame(left_join(as.data.frame(rowData(se)), table,
                                         by=c("ID"="rowname")))
  return(se)
}

# customized from DEP's add_rejections: https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/functions.R#L1032
#' Mark significant proteins
#'
#' \code{add_rejections_customized} marks significant proteins based on defined cutoffs.
#'
#' @param diff SummarizedExperiment,
#' Proteomics dataset on which differential enrichment analysis
#' has been performed (output from \code{\link{test_diff}()}).
#' @param alpha Numeric(1),
#' Sets the threshold for the adjusted P value.
#' @param lfc Numeric(1),
#' Sets the threshold for the log2 fold change.
#' @return A SummarizedExperiment object
#' annotated with logical columns indicating significant proteins.
#' @examples
#' # Load example
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' # Filter, normalize and impute missing values
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinProb", q = 0.01)
#'
#' # Test for differentially expressed proteins
#' diff <- test_diff(imputed, "control", "Ctrl")
#' dep <- add_rejections(diff, alpha = 0.05, lfc = 1)
#' @export
add_rejections_customized <- function(diff, alpha = 0.05, lfc = 1) {
  # Show error if inputs are not the required classes
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(inherits(diff, "SummarizedExperiment"),
                          is.numeric(alpha),
                          length(alpha) == 1,
                          is.numeric(lfc),
                          length(lfc) == 1)

  row_data <- rowData(diff, use.names = FALSE) %>%
    as.data.frame()
  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(diff)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(diff)),
         "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }

  # get all columns with adjusted p-values and log2 fold changes
  cols_p <- grep("_p.adj", colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))

  # Mark differential expressed proteins by
  # applying alpha and log2FC parameters per protein
  if(length(cols_p) == 1) {
    rowData(diff)$significant <-
      row_data[, cols_p] <= alpha & abs(row_data[, cols_diff]) >= lfc
    rowData(diff)$contrast_significant <-
      rowData(diff, use.names = FALSE)$significant
    colnames(rowData(diff))[ncol(rowData(diff, use.names = FALSE))] <-
      gsub("p.adj", "significant", colnames(row_data)[cols_p])
  }
  if(length(cols_p) > 1) {
    p_reject <- row_data[, cols_p] <= alpha
    p_reject[is.na(p_reject)] <- FALSE
    diff_reject <- abs(row_data[, cols_diff]) >= lfc
    diff_reject[is.na(diff_reject)] <- FALSE
    sign_df <- p_reject & diff_reject
    sign_df <- cbind(sign_df,
                     significant = apply(sign_df, 1, function(x) any(x)))
    colnames(sign_df) <- gsub("_p.adj", "_significant", colnames(sign_df))
    sign_df <- cbind(ID = row_data$ID, as.data.frame(sign_df))
    rowData(diff) <- as.data.frame(left_join(as.data.frame(rowData(diff)), sign_df,
                                           by=c("ID"="ID")))
  }
  return(diff)
}


# similar to test_match_lfq_column_design
test_match_tmt_column_design <- function(unique_data, lfq_columns, exp_design){
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(unique_data),
                          is.integer(lfq_columns),
                          is.data.frame(exp_design))
  
  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(unique_data))) {
    stop(safeError("'Gene name' and/or 'Protein ID' columns are not present in
          protein groups input file"
    ))
  }
  
  if(any(!c("label", "condition", "replicate") %in% colnames(exp_design))) {
    stop(safeError("'label', 'condition' and/or 'replicate' columns
         are not present in the experimental design"))
  }
  
  if(any(!apply(unique_data[, lfq_columns], 2, is.numeric))) {
    stop(safeError("specified 'columns' should be numeric
         Run make_se_parse() with the appropriate columns as argument"))
  }
  
  raw <- unique_data[, lfq_columns]
  # expdesign <- mutate(exp_design, condition = make.names(condition)) %>%
  #   unite(ID, label, remove = FALSE)
  # rownames(expdesign) <- expdesign$ID
  expdesign <- exp_design
  # print(expdesign)
  rownames(expdesign) <- expdesign$label
  # print(make.names(expdesign$label))
  # print(make.names(colnames(raw)))
  matched <- match(make.names(expdesign$label),
                   make.names(colnames(raw)))
  
  # TODO: give warning message to indicate which columns are not matched
  # if(any(is.na(matched))) {
  if(all(is.na(matched))) {
    stop(safeError("The labels/'run names' in the experimental design DID NOT match
         with column names in TMT-I report.
         Please run FragPipe-Analyst with correct labels in the experimental design"))
  }
}

# similar to test_match_lfq_column_design
test_match_DIA_column_design <- function(unique_data, lfq_columns, exp_design){
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(unique_data),
                          is.integer(lfq_columns),
                          is.data.frame(exp_design))

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(unique_data))) {
    stop(safeError("'Gene name' and/or 'Protein ID' columns are not present in
          protein groups input file"
    ))
  }

  if(any(!c("label", "condition", "replicate") %in% colnames(exp_design))) {
    stop(safeError("'label', 'condition' and/or 'replicate' columns
         are not present in the experimental design"))
  }
  
  if(any(!apply(unique_data[, lfq_columns], 2, is.numeric))) {
    stop(safeError("specified 'columns' should be numeric
         Run make_se_parse() with the appropriate columns as argument"))
  }
  
  raw <- unique_data[, lfq_columns]
  # expdesign <- mutate(exp_design, condition = make.names(condition)) %>%
  #   unite(ID, label, remove = FALSE)
  # rownames(expdesign) <- expdesign$ID
  expdesign <- exp_design
  # print(expdesign)
  rownames(expdesign) <- expdesign$label
  # print(make.names(expdesign$label))
  # print(make.names(colnames(raw)))
  matched <- match(make.names(expdesign$label),
                   make.names(colnames(raw)))
  
  # TODO: give warning message to indicate which columns are not matched
  # if(any(is.na(matched))) {
  if(all(is.na(matched))) {
    stop(safeError("The labels/'run names' in the experimental design DID NOT match
         with column names in TMT-I report.
         Please run FragPipe-Analyst with correct labels in the experimental design"))
  }
}

# original: plot_coverage from
# https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/plot_functions_frequencies.R#L129
#' Plot protein coverage
#'
#' \code{plot_coverage_customized} generates a barplot
#' of the protein coverage in all samples.
#'
#' @param se SummarizedExperiment,
#' Data object for which to plot observation frequency.
#' @param plot Logical(1),
#' If \code{TRUE} (default) the barplot is produced.
#' Otherwise (if \code{FALSE}), the data which the
#' barplot is based on are returned.
#' @return Barplot of protein coverage in samples
#' (generated by \code{\link[ggplot2]{ggplot}})
#' @export
plot_coverage_customized <- function(se, plot = TRUE) {
  # Show error if input is not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.logical(plot),
                          length(plot) == 1)

  # Make a binary long data.frame (1 = valid value, 0 = missing value)
  df <- assay(se) %>%
    data.frame() %>%
    rownames_to_column() %>%
    gather(ID, bin, -rowname) %>%
    mutate(bin = ifelse(is.na(bin), 0, 1))
  # Identify the number of experiments a protein was observed
  stat <- df %>%
    group_by(rowname) %>%
    summarize(sum = sum(bin))
  # Get the frequency of the number of experiments proteins were observerd
  # and plot the cumulative sum of these numbers
  table <- table(stat$sum) %>%
    data.frame()
  p <- ggplot(table, aes(x = "all", y = Freq, fill = Var1)) +
    geom_col(col = "white") +
    scale_fill_grey(start = 0.8, end = 0.2) +
    labs(title = "Feature coverage",
         x = "",
         y = "Number of features",
         fill = "Samples") +
    theme_DEP1()
  if(plot) {
    return(p)
  } else {
    df <- as.data.frame(table)
    colnames(df) <- c("samples", "features")
    return(df)
  }
}


# original: get_results from 
# https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/functions.R
#' Generate a results table
#'
#' \code{get_results_customized} generates a results table from a proteomics dataset
#' on which differential enrichment analysis was performed.
#'
#' @param dep SummarizedExperiment,
#' Data object for which differentially enriched proteins are annotated
#' (output from \code{\link{test_diff}()} and \code{\link{add_rejections}()}).
#' @return A data.frame object
#' containing all results variables from the performed analysis.
#' @examples
#' # Load example
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' # Filter, normalize and impute missing values
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinProb", q = 0.01)
#'
#' # Test for differentially expressed proteins
#' diff <- test_diff(imputed, "control", "Ctrl")
#' dep <- add_rejections(diff, alpha = 0.05, lfc = 1)
#'
#' # Get results
#' results <- get_results(dep)
#' colnames(results)
#'
#' significant_proteins <- results[results$significant,]
#' nrow(significant_proteins)
#' head(significant_proteins)
#' @export
get_results_customized <- function(dep) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"))
  
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
  
  # Obtain average protein-centered enrichment values per condition
  row_data$mean <- rowMeans(assay(dep), na.rm = TRUE)
  centered <- assay(dep) - row_data$mean
  centered <- data.frame(centered) %>%
    rownames_to_column() %>%
    gather(ID, val, -rowname) %>%
    left_join(., data.frame(colData(dep)), by = c("ID"="label"))
  centered <- group_by(centered, rowname, condition) %>%
    summarize(val = mean(val, na.rm = TRUE)) %>%
    mutate(val = signif(val, digits = 3)) %>%
    spread(condition, val)
  colnames(centered)[2:ncol(centered)] <-
    paste(colnames(centered)[2:ncol(centered)], "_centered", sep = "")
  
  # Obtain average enrichments of conditions versus the control condition
  ratio <- as.data.frame(row_data) %>%
    column_to_rownames("name") %>%
    select(ends_with("diff")) %>%
    signif(., digits = 3) %>%
    rownames_to_column()
  colnames(ratio)[2:ncol(ratio)] <-
    gsub("_diff", "_ratio", colnames(ratio)[2:ncol(ratio)])
  df <- left_join(ratio, centered, by = "rowname")
  
  # Select the adjusted p-values and significance columns
  pval <- as.data.frame(row_data) %>%
    column_to_rownames("name") %>%
    select(ends_with("p.val"),
           ends_with("p.adj"),
           ends_with("significant")) %>%
    rownames_to_column()
  pval[, grep("p.adj", colnames(pval))] <-
    pval[, grep("p.adj", colnames(pval))] %>%
    signif(digits = 3)
  
  # Join into a results table
  ids <- as.data.frame(row_data) %>% select(name, ID)
  table <- left_join(ids, pval, by = c("name" = "rowname"))
  table <- left_join(table, df, by = c("name" = "rowname")) %>%
    arrange(desc(significant))
  return(table)
}

plot_pca_plotly <- function(dep, x = 1, y = 2, indicate = c("condition", "replicate"),
                    label = FALSE, n = 500, point_size = 8, label_size = 3, plot = TRUE, ID_col="ID", exp="LFQ", scale=F) {
  if(is.integer(x)) x <- as.numeric(x)
  if(is.integer(y)) y <- as.numeric(y)
  if(is.integer(n)) n <- as.numeric(n)
  if(is.integer(point_size)) point_size <- as.numeric(point_size)
  if(is.integer(label_size)) label_size <- as.numeric(label_size)
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.numeric(x),
                          length(x) == 1,
                          is.numeric(y),
                          length(y) == 1,
                          is.numeric(n),
                          length(n) == 1,
                          is.character(indicate),
                          is.logical(label),
                          is.numeric(point_size),
                          length(point_size) == 1,
                          is.numeric(label_size),
                          length(label_size) == 1,
                          is.logical(plot),
                          length(plot) == 1)
  
  # Check for valid x and y values
  if(x > ncol(dep) | y > ncol(dep)) {
    stop(paste0("'x' and/or 'y' arguments are not valid\n",
                "Run plot_pca() with 'x' and 'y' <= ",
                ncol(dep), "."),
         call. = FALSE)
  }
  

  # Check for valid 'indicate'
  columns <- colnames(colData(dep))
  if(!is.null(indicate)) {
    if(length(indicate) > 3) {
      stop("Too many features in 'indicate'
        Run plot_pca() with a maximum of 3 indicate features")
    }
    if(any(!indicate %in% columns)) {
      stop(paste0("'",
                  paste0(indicate, collapse = "' and/or '"),
                  "' column(s) is/are not present in ",
                  deparse(substitute(dep)),
                  ".\nValid columns are: '",
                  paste(columns, collapse = "', '"),
                  "'."),
           call. = FALSE)
    }
  }

  data <- assay(dep)
  data <- data[complete.cases(data), ]
  # Get the variance per protein and take the top n variable proteins
  var <- apply(data, 1, sd)
  # Check for valid 'n' value
  if(n > nrow(data)) {
    message(paste("'n' argument is larger than number of features availble(",
                  nrow(data), ").", nrow(data), "features will be used for PCA calculation."))
    df <- data
    n <- nrow(data)
  } else {
    df <- data[order(var, decreasing = TRUE)[seq_len(n)],]
  }
  
  # Calculate PCA
  pca <- prcomp(t(df), scale = scale)
  pca_df <- pca$x %>%
    data.frame() %>%
    rownames_to_column() %>%
    left_join(., data.frame(colData(dep)), by = c("rowname" = ID_col))
  
  # Calculate the percentage of variance explained
  percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  
  # Make factors of indicate features
  for(feat in indicate) {
    pca_df[[feat]] <- as.factor(pca_df[[feat]])
  }
  
  
  if(length(indicate) == 1) {
    if (exp == "TMT") {
    pca_df$plex <- as.factor(pca_df$plex)
    p <- plot_ly() %>%
      add_trace(data=pca_df, type = 'scatter', marker = list(size = point_size),
                mode = 'markers',
                x = ~PC1,
                y = ~PC2,
                color = as.formula(paste0('~', indicate[1])),
                text = pca_df$sample_name,
                xaxis="x",
                yaxis="y",
                legendgroup=indicate[1],
                legendgrouptitle_text=indicate[1]) %>%
      add_trace(data=pca_df, type = "scatter", marker = list(size = point_size),
                mode = 'markers',
                x = ~PC1,
                y = ~PC2,
                text = pca_df$sample_name,
                color = as.formula(paste0('~', "plex")),
                legendgroup="plex",
                legendgrouptitle_text="plex",
                xaxis="x2",
                yaxis="y2", visible = FALSE, inherit = FALSE) %>%
      plotly::layout(title = paste0('PCA plot (', n, " features used)"),
                     xaxis = list(title = paste0("PC", x, ": ", percent[x], "%")),
                     xaxis2 = list(title = paste0("PC", x, ": ", percent[x], "%"), overlaying="x", visible=F),
                     yaxis = list(title = paste0("PC", y, ": ", percent[y], "%")),
                     yaxis2 = list(title = paste0("PC", y, ": ", percent[y], "%"), overlaying="y", visible=F),
                     updatemenus = list(
                       list(
                         y = 0.8,
                         buttons = list(
                           list(method = "update",
                                args = list(list(visible=unlist(Map(rep, x = c(T, F), each = c(length(unique(pca_df$condition)),
                                                                                               length(unique(pca_df$plex)))))),
                                            list(xaxis = list(title = paste0("PC", x, ": ", percent[x], "%"),
                                                              visible = TRUE),
                                                 xaxis2 = list(overlaying = "x", visible = FALSE),
                                                 yaxis = list(title = paste0("PC", y, ": ", percent[y], "%"),
                                                              visible = TRUE),
                                                 yaxis2 = list(overlaying = "y", visible = FALSE))),
                                label = "by condition"),
                           list(method = "update",
                                args = list(list(visible=unlist(Map(rep, x = c(F, T), each = c(length(unique(pca_df$condition)),
                                                                                               length(unique(pca_df$plex)))))),
                                            list(xaxis = list(visible = F),
                                                 xaxis2 = list(title = paste0("PC", x, ": ", percent[x], "%"),
                                                               overlaying = "x", visible = T),
                                                 yaxis = list(visible = F),
                                                 yaxis2 = list(title = paste0("PC", y, ": ", percent[y], "%"),
                                                               overlaying = "y", visible = T))),
                                label = "by plex")
                         )
                       )
                     )
                  )
    } else {
      p <- plot_ly(data=pca_df, type = 'scatter', mode = 'markers', marker = list(size = point_size)) %>%
        plotly::layout(title = paste0('PCA plot (', n, " features used)"),
                       xaxis = list(title = paste0("PC", x, ": ", percent[x], "%")),
                       yaxis = list(title = paste0("PC", y, ": ", percent[y], "%"))) %>%
        add_trace(type = "scatter",
                  x = ~PC1,
                  y = ~PC2,
                  text = pca_df$sample_name,
                  color = as.formula(paste0('~', indicate[1])),
                  mode = 'markers',
                  legendgroup=indicate[1],
                  legendgrouptitle_text=indicate[1])
    }
  } else if(length(indicate) == 2) {
    if (exp == "TMT" | exp == "TMT-peptide"){
      pca_df$plex <- as.factor(pca_df$plex)
      p <- plot_ly() %>%
        #Overlay color for gears
        add_trace(data=pca_df, type = "scatter",
                  x = ~PC1,
                  y = ~PC2,
                  symbol = as.formula(paste0('~', indicate[2])),
                  marker = list(color = "grey", size = point_size + 3),
                  text = pca_df$sample_name,
                  hoverinfo = 'text',
                  mode = 'markers',
                  legendgroup=indicate[2],
                  legendgrouptitle_text=indicate[2]) %>%
        add_trace(data=pca_df, type = "scatter",
                  x = ~PC1,
                  y = ~PC2,
                  marker = list(size = point_size),
                  color = as.formula(paste0('~', indicate[1])),
                  mode = 'markers',
                  text = pca_df$sample_name,
                  hoverinfo = 'text',
                  legendgroup=indicate[1],
                  legendgrouptitle_text=indicate[1]) %>%
        add_trace(data=pca_df, type = "scatter",
                  x = ~PC1,
                  y = ~PC2,
                  marker = list(size = point_size),
                  color = as.formula(paste0('~', "plex")),
                  text = pca_df$sample_name,
                  hoverinfo = 'text',
                  mode = 'markers',
                  legendgroup="plex",
                  legendgrouptitle_text="plex",
                  xaxis="x2", yaxis="y2", visible=F) %>%
        plotly::layout(title = paste0('PCA plot (', n, " features used)"),
                       xaxis = list(title = paste0("PC", x, ": ", percent[x], "%")),
                       xaxis2 = list(title = paste0("PC", x, ": ", percent[x], "%"), overlaying="x", visible=F),
                       yaxis = list(title = paste0("PC", y, ": ", percent[y], "%")),
                       yaxis2 = list(title = paste0("PC", y, ": ", percent[y], "%"), overlaying="y", visible=F),
                       legend=list(itemclick = FALSE,
                                   itemdoubleclick = FALSE,
                                   groupclick = FALSE),
                       updatemenus = list(
                         list(
                           y = 0.8,
                           buttons = list(
                             list(method = "update",
                                  args = list(list(visible=unlist(Map(rep, x = c(T, T, F), each = c(length(unique(pca_df$condition)),
                                                                                                    length(unique(pca_df$replicate)),
                                                                                                    length(unique(pca_df$plex)))))),
                                              list(xaxis = list(title = paste0("PC", x, ": ", percent[x], "%"),
                                                                visible = TRUE),
                                                   xaxis2 = list(overlaying = "x", visible = FALSE),
                                                   yaxis = list(title = paste0("PC", y, ": ", percent[y], "%"),
                                                                visible = TRUE),
                                                   yaxis2 = list(overlaying = "y", visible = FALSE))),
                                  label = "by condition"),
                             list(method = "update",
                                  args = list(list(visible=unlist(Map(rep, x = c(F, F, T), each = c(length(unique(pca_df$condition)),
                                                                                                    length(unique(pca_df$replicate)),
                                                                                                    length(unique(pca_df$plex)))))),
                                              list(xaxis = list(visible = F),
                                                   xaxis2 = list(title = paste0("PC", x, ": ", percent[x], "%"),
                                                                 overlaying = "x", visible = T),
                                                   yaxis = list(visible = F),
                                                   yaxis2 = list(title = paste0("PC", y, ": ", percent[y], "%"),
                                                                 overlaying = "y", visible = T))),
                                  label = "by plex")
                           )
                         )
                       )
                       )
    } else {
      p <- plot_ly(data=pca_df, type = 'scatter',
                   mode = 'markers', marker = list(size = point_size), text=~rowname) %>%
        #Overlay color for gears
        add_trace(type = "scatter",
                  x = ~PC1,
                  y = ~PC2,
                  symbol = as.formula(paste0('~', indicate[2])),
                  marker = list(color = "grey", size = point_size + 3),
                  mode = 'markers',
                  text = pca_df$sample_name,
                  hoverinfo = 'text',
                  legendgroup=indicate[2],
                  legendgrouptitle_text=indicate[2]) %>%
        add_trace(type = "scatter",
                  x = ~PC1,
                  y = ~PC2,
                  color = as.formula(paste0('~', indicate[1])),
                  mode = 'markers',
                  text = pca_df$sample_name,
                  hoverinfo = 'text',
                  legendgroup=indicate[1],
                  legendgrouptitle_text=indicate[1]) %>%
        plotly::layout(title = paste0('PCA plot (', n, " features used)"),
                       xaxis = list(title = paste0("PC", x, ": ", percent[x], "%")),
                       yaxis = list(title = paste0("PC", y, ": ", percent[y], "%")),
                       legend=list(itemclick = FALSE,
                                   itemdoubleclick = FALSE,
                                   groupclick = FALSE))
    }
  }
  
  if(plot) {
    return(p)
  } else {
    df <- pca_df %>%
      select(rowname, paste0("PC", c(x, y)), match(indicate, colnames(pca_df)))
    colnames(df)[1] <- "sample"
    return(df)
  }
}

# original plot_pca from
# https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/plot_functions_explore.R
#' Plot PCA
#'
#' \code{plot_pca} generates a PCA plot using the top variable proteins.
#'
#' @param dep SummarizedExperiment,
#' Data object for which differentially enriched proteins are annotated
#' (output from \code{\link{test_diff}()} and \code{\link{add_rejections}()}).
#' @param x Integer(1),
#' Sets the principle component to plot on the x-axis.
#' @param y Integer(1),
#' Sets the principle component to plot on the y-axis.
#' @param indicate Character,
#' Sets the color, shape and facet_wrap of the plot
#' based on columns from the experimental design (colData).
#' @param label Logical,
#' Whether or not to add sample labels.
#' @param n Integer(1),
#' Sets the number of top variable proteins to consider.
#' @param point_size Integer(1),
#' Sets the size of the points.
#' @param label_size Integer(1),
#' Sets the size of the labels.
#' @param plot Logical(1),
#' If \code{TRUE} (default) the PCA plot is produced.
#' Otherwise (if \code{FALSE}), the data which the
#' PCA plot is based on are returned.
#' @return A scatter plot (generated by \code{\link[ggplot2]{ggplot}}).
#' @examples
#' # Load example
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' # Filter, normalize and impute missing values
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinProb", q = 0.01)
#'
#' # Test for differentially expressed proteins
#' diff <- test_diff(imputed, "control", "Ctrl")
#' dep <- add_rejections(diff, alpha = 0.05, lfc = 1)
#'
#' # Plot PCA
#' plot_pca(dep)
#' plot_pca(dep, indicate = "condition")
#' @export
plot_pca_customized <- function(dep, x = 1, y = 2, indicate = c("condition", "replicate"),
                     label = FALSE, n = 500, point_size = 4, label_size = 3, plot = TRUE, ID_col="ID", scale=F) {
  if(is.integer(x)) x <- as.numeric(x)
  if(is.integer(y)) y <- as.numeric(y)
  if(is.integer(n)) n <- as.numeric(n)
  if(is.integer(point_size)) point_size <- as.numeric(point_size)
  if(is.integer(label_size)) label_size <- as.numeric(label_size)
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.numeric(x),
                          length(x) == 1,
                          is.numeric(y),
                          length(y) == 1,
                          is.numeric(n),
                          length(n) == 1,
                          is.character(indicate),
                          is.logical(label),
                          is.numeric(point_size),
                          length(point_size) == 1,
                          is.numeric(label_size),
                          length(label_size) == 1,
                          is.logical(plot),
                          length(plot) == 1)
  
  # Check for valid x and y values
  if(x > ncol(dep) | y > ncol(dep)) {
    stop(paste0("'x' and/or 'y' arguments are not valid\n",
                "Run plot_pca() with 'x' and 'y' <= ",
                ncol(dep), "."),
         call. = FALSE)
  }
  
  # Check for valid 'indicate'
  columns <- colnames(colData(dep))
  if(!is.null(indicate)) {
    if(length(indicate) > 3) {
      stop("Too many features in 'indicate'
        Run plot_pca() with a maximum of 3 indicate features")
    }
    if(any(!indicate %in% columns)) {
      stop(paste0("'",
                  paste0(indicate, collapse = "' and/or '"),
                  "' column(s) is/are not present in ",
                  deparse(substitute(dep)),
                  ".\nValid columns are: '",
                  paste(columns, collapse = "', '"),
                  "'."),
           call. = FALSE)
    }
  }
  
  data <- assay(dep)
  data <- data[complete.cases(data), ]
  # Get the variance per protein and take the top n variable proteins
  var <- apply(data, 1, sd)
  # Check for valid 'n' value
  if(n > nrow(data)) {
    message(paste("'n' argument is larger than number of features availble(",
                  nrow(data), ").", nrow(data), "features will be used for PCA calculation."))
    df <- data
    n <- nrow(data)
  } else {
    df <- data[order(var, decreasing = TRUE)[seq_len(n)],]
  }
  
  # Calculate PCA
  pca <- prcomp(t(df), scale = scale)
  pca_df <- pca$x %>%
    data.frame() %>%
    rownames_to_column() %>%
    left_join(., data.frame(colData(dep)), by = c("rowname" = ID_col))
  
  # Calculate the percentage of variance explained
  percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  
  # Make factors of indicate features
  for(feat in indicate) {
    pca_df[[feat]] <- as.factor(pca_df[[feat]])
  }
  
  # Plot the PCA plot
  p <- ggplot(pca_df, aes(get(paste0("PC", x)), get(paste0("PC", y)))) +
    labs(title = paste0("PCA plot - top ", n, " variable proteins"),
         x = paste0("PC", x, ": ", percent[x], "%"),
         y = paste0("PC", y, ": ", percent[y], "%")) +
    coord_fixed() +
    theme_DEP1()
  
  if(length(indicate) == 0) {
    p <- p + geom_point(size = point_size)
  }
  if(length(indicate) == 1) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]]),
                        size = point_size) +
      labs(col = indicate[1])
  }
  if(length(indicate) == 2) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]],
                            shape = pca_df[[indicate[2]]]),
                        size = point_size) +
      labs(col = indicate[1],
           shape = indicate[2])
  }
  if(length(indicate) == 3) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]],
                            shape = pca_df[[indicate[2]]]),
                        size = point_size) +
      facet_wrap(~pca_df[[indicate[3]]])
    labs(col = indicate[1],
         shape = indicate[2])
  }
  if(label) {
    p <- p + geom_text(aes(label = rowname), size = label_size)
  }
  if(plot) {
    return(p)
  } else {
    df <- pca_df %>%
      select(rowname, paste0("PC", c(x, y)), match(indicate, colnames(pca_df)))
    colnames(df)[1] <- "sample"
    return(df)
  }
}

# Customized from DEP's plot_numbers
# https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/plot_functions_frequencies.R
#' Plot protein numbers
#'
#' \code{plot_numbers_customized} generates a barplot
#' of the number of identified proteins per sample.
#'
#' @param se SummarizedExperiment,
#' Data object for which to plot protein numbers
#' (output from \code{\link{make_se}()} or \code{\link{make_se_parse}()}).
#' @param plot Logical(1),
#' If \code{TRUE} (default) the barplot is produced.
#' Otherwise (if \code{FALSE}), the data which the
#' barplot is based on are returned.
#' @return Barplot of the number of identified proteins per sample
#' (generated by \code{\link[ggplot2]{ggplot}})

plot_numbers_customized <- function(se, plot = TRUE, exp = "LFQ") {
  # Show error if input is not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.logical(plot),
                          length(plot) == 1)
  
  # Make a binary long data.frame (1 = valid value, 0 = missing value)
  # print(rowData(se))
  df <- assay(se) %>%
    data.frame(check.names = F) %>%
    rownames_to_column() %>%
    gather(ID, bin, -rowname) %>%
    mutate(bin = ifelse(is.na(bin), 0, 1))
  # print(df)


  # Summarize the number of proteins identified
  # per sample and generate a barplot
  stat <- df %>%
    group_by(ID) %>%
    summarize(n = n(), sum = sum(bin)) %>%
    left_join(., data.frame(colData(se)), by = c("ID"="label"))

  p <- ggplot(stat, aes(x = sample_name, y = sum, fill = condition)) +
    geom_col() +
    geom_hline(yintercept = unique(stat$n)) +
    labs(title = "Features per sample", x = "", y = "Number of features") +
    theme_DEP2()

  if(plot) {
    return(p)
  } else {
    df <- as.data.frame(stat)
    colnames(df)[seq_len(3)] <- c("sample", "total_features", "features_in_sample")
    return(df)
  }
}

plot_numbers_by_plex_set <- function(se, ...) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  unique_plexes <- unique(colData(se)$plex)
  prot_v <- c()
  for(i in 1:length(unique_plexes)){
    n_prot <- assay(se[, se$plex == unique_plexes[i]]) %>%
      data.frame() %>%
      filter(if_all(everything(), ~!is.na(.))) %>%
      nrow()
    prot_v <- c(prot_v, n_prot)
  }
  df_prot <- data.frame(plex=factor(unique_plexes), num_protein=prot_v)
  return(ggplot(df_prot, aes(x = plex, y = num_protein)) +
           geom_bar(stat="identity") +
           labs(title = "Number of proteins across plex sets", x = "Plex",
                y = "Number of proteins") +
           theme_DEP2())
}

# customized from 
# https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/plot_functions_QC.R
#' Visualize normalization
#'
#' \code{plot_normalization_customized} generates boxplots
#' of all conditions for input objects, e.g. before and after normalization.
#'
#' @param se SummarizedExperiment,
#' Data object, e.g. before normalization (output from \code{\link{make_se}()}
#' or \code{\link{make_se_parse}()}).
#' @param ... Additional SummarizedExperiment object(s),
#' E.g. data object after normalization
#' (output from \code{\link{normalize_vsn}}).
#' @return Boxplots of all conditions
#' for input objects, e.g. before and after normalization
#'   (generated by \code{\link[ggplot2]{ggplot}}).
#' Adding components and other plot adjustments can be easily done
#' using the ggplot2 syntax (i.e. using '+')
#' @examples
#' # Load example
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' # Filter and normalize
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#'
#' # Plot normalization
#' plot_normalization(se, filt, norm)
#' @export
plot_normalization_customized <- function(se, ...) {
  # Get arguments from call
  call <- match.call()
  arglist <- lapply(call[-1], function(x) x)
  var.names <- vapply(arglist, deparse, character(1))
  arglist <- lapply(arglist, eval.parent, n = 2)
  names(arglist) <- var.names
  
  # Show error if inputs are not the required classes
  lapply(arglist, function(x) {
    assertthat::assert_that(inherits(x,
                                     "SummarizedExperiment"),
                            msg = "input objects need to be of class 'SummarizedExperiment'")
    if (any(!c("label", "condition", "replicate") %in% colnames(colData(x)))) {
      # ID is not required
      stop("'label', 'condition' and/or 'replicate' ",
           "columns are not present in (one of) the input object(s)",
           "\nRun make_se() or make_se_parse() to obtain the required columns",
           call. = FALSE)
    }
  })
  
  # Function to get a long data.frame of the assay data
  # annotated with sample info
  gather_join <- function(se) {
    assay(se) %>%
      data.frame() %>%
      gather(ID, val) %>%
      left_join(., data.frame(colData(se)), by = c("ID"="label"))
  }
  
  df <- map_df(arglist, gather_join, .id = "var") %>%
    mutate(var = factor(var, levels = names(arglist)))
  
  # Boxplots for conditions with facet_wrap
  # for the original and normalized values
  ggplot(df, aes(x = ID, y = val, fill = condition)) +
    geom_boxplot(notch = TRUE, na.rm = TRUE) +
    coord_flip() +
    facet_wrap(~var, ncol = 1) +
    labs(x = "", y = expression(log[2]~"Intensity")) +
    theme_DEP1()
}

plot_normalization_DIA_customized <- function(se, ...) {
  # Get arguments from call
  call <- match.call()
  arglist <- lapply(call[-1], function(x) x)
  var.names <- vapply(arglist, deparse, character(1))
  arglist <- lapply(arglist, eval.parent, n = 2)
  names(arglist) <- var.names
  
  # Show error if inputs are not the required classes
  lapply(arglist, function(x) {
    assertthat::assert_that(inherits(x,
                                     "SummarizedExperiment"),
                            msg = "input objects need to be of class 'SummarizedExperiment'")
    if (any(!c("label", "condition", "replicate") %in% colnames(colData(x)))) {
      # ID is not required
      stop("'label', 'condition' and/or 'replicate' ",
           "columns are not present in (one of) the input object(s)",
           "\nRun make_se() or make_se_parse() to obtain the required columns",
           call. = FALSE)
    }
  })
  
  # Function to get a long data.frame of the assay data
  # annotated with sample info
  gather_join <- function(se) {
    assay(se) %>%
      data.frame(check.names = F) %>%
      gather(ID, val) %>%
      left_join(., data.frame(colData(se)), by = c("ID"="label"))
  }
  
  df <- map_df(arglist, gather_join, .id = "var") %>%
    mutate(var = factor(var, levels = names(arglist)))
  
  # Boxplots for conditions with facet_wrap
  # for the original and normalized values
  df$label <- paste(df$experiment, df$replicate, sep="_")
  p <- ggplot(df, aes(x = label, y = val, fill = condition)) +
    geom_boxplot(notch = TRUE, na.rm = TRUE) +
    coord_flip() +
    facet_wrap(~var, ncol = 1) +
    labs(x = "", y = expression(log[2]~"Intensity")) +
    theme_DEP1()
  return(p)
}

# Customized from DEP's plot_imputation
# https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/plot_functions_QC.R
#' Visualize imputation
#'
#' \code{plot_density} generates density plots
#' of all conditions for input objects, e.g. before and after imputation.
#'
#' @param ses a named list of SummarizedExperiment,
#' (output from \code{\link{impute}()}).
#' @return Density plots of all conditions for input objects, e.g. before and
#' after imputation (generated by \code{\link[ggplot2]{ggplot}}).
#' @examples
#' # Plot density
#' plot_density(list("filtered data"=filt, "normalized data"=norm, "imputed data"=imputed))
#' @export
plot_density <- function(ses) {
  # Function to get a long data.frame of the assay data
  # annotated with sample info
  gather_join <- function(se) {
    assay(se) %>%
      data.frame(check.names = F) %>%
      gather(ID, val) %>%
      left_join(., data.frame(colData(se)), by = c("ID"="label"))
  }
  
  df <- map_df(ses, gather_join, .id = "var") %>%
    mutate(var = factor(var, levels = names(ses)))

  # Density plots for different conditions with facet_wrap
  # for original and imputed samles
  ggplot(df, aes(val, col = condition)) +
    geom_density(na.rm = TRUE) +
    facet_wrap(~var, ncol = 1) +
    labs(x = expression(log[2]~"Intensity"), y = "Density") +
    theme_DEP1()
}

plot_density_spectral_count <- function(ses) {
  # Function to get a long data.frame of the assay data
  # annotated with sample info
  gather_join <- function(se) {
    assay(se) %>%
      data.frame(check.names = F) %>%
      gather(ID, val) %>%
      left_join(., data.frame(colData(se)), by = c("ID"="label"))
  }
  
  df <- map_df(ses, gather_join, .id = "var") %>%
    mutate(var = factor(var, levels = names(ses)))
  
  # Density plots for different conditions with facet_wrap
  # for original and imputed samles
  ggplot(df, aes(val, col = condition)) +
    geom_density(na.rm = TRUE) +
    facet_wrap(~var, ncol = 1) +
    labs(x = "Spectral count", y = "Density") +
    theme_DEP1()
}

test_limma_customized <- function(se, type = c("control", "all", "manual"),
                       control = NULL, test = NULL,
                       design_formula = formula(~ 0 + condition),
                       paired = FALSE) {
  
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
  
  if(any(!c("name") %in% colnames(rowData(se)))) {
    stop("'name' column is not present in '",
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
  rowData(se) <- as.data.frame(left_join(as.data.frame(rowData(se)), table,
                                         by=c("ID"="rowname")))
  return(se)
}

# https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/functions.R
#' Imputation by random draws from a manually defined distribution
#'
#' \code{manual_impute_customized} imputes missing values in a proteomics dataset
#' by random draws from a manually defined distribution.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}). It is adviced to first remove
#' proteins with too many missing values using \code{\link{filter_missval}()}
#' and normalize the data using \code{\link{normalize_vsn}()}.
#' @param shift Numeric(1),
#' Sets the left-shift of the distribution (in standard deviations) from
#' the median of the original distribution.
#' @param scale Numeric(1),
#' Sets the width of the distribution relative to the
#' standard deviation of the original distribution.
#' @return An imputed SummarizedExperiment object.
#' @examples
#' # Load example
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' # Filter and normalize
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#'
#' # Impute missing values manually
#' imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
#' @export
manual_impute_customized <- function(se, scale = 0.3, shift = 1.8) {
  if(is.integer(scale)) scale <- is.numeric(scale)
  if(is.integer(shift)) shift <- is.numeric(shift)
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.numeric(scale),
                          length(scale) == 1,
                          is.numeric(shift),
                          length(shift) == 1)
  
  se_assay <- assay(se)
  
  # Show error if there are no missing values
  if(!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)), "'",
         call. = FALSE)
  }
  
  # Get descriptive parameters of the current sample distributions
  stat <- se_assay %>%
    data.frame(check.names = F) %>%
    rownames_to_column() %>%
    gather(samples, value, -rowname) %>%
    filter(!is.na(value))  %>%
    group_by(samples) %>%
    summarise(mean = mean(value),
              median = median(value),
              sd = sd(value),
              n = n(),
              infin = nrow(se_assay) - n)
  # Impute missing values by random draws from a distribution
  # which is left-shifted by parameter 'shift' * sd and scaled by parameter 'scale' * sd.
  for (a in seq_len(nrow(stat))) {
    
    set.seed(123)
    assay(se)[is.na(assay(se)[, stat$samples[a]]), stat$samples[a]] <-
      rnorm(stat$infin[a],
            mean = stat$median[a] - shift * stat$sd[a],
            sd = stat$sd[a] * scale)
  }
  return(se)
}


# https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/functions.R
#' Impute missing values
#'
#' \code{impute_customized} imputes missing values in a proteomics dataset.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}). It is adviced to first remove
#' proteins with too many missing values using \code{\link{filter_missval}()}
#' and normalize the data using \code{\link{normalize_vsn}()}.
#' @param fun "bpca", "knn", "QRILC", "MLE", "MinDet",
#' "MinProb", "man", "min", "zero", "mixed" or "nbavg",
#' Function used for data imputation based on \code{\link{manual_impute}}
#' and \code{\link[MSnbase:impute-methods]{impute}}.
#' @param ... Additional arguments for imputation functions as depicted in
#' \code{\link{manual_impute}} and \code{\link[MSnbase:impute-methods]{impute}}.
#' @return An imputed SummarizedExperiment object.
#' @examples
#' # Load example
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' # Filter and normalize
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#'
#' # Impute missing values using different functions
#' imputed_MinProb <- impute(norm, fun = "MinProb", q = 0.05)
#' imputed_QRILC <- impute(norm, fun = "QRILC")
#'
#' imputed_knn <- impute(norm, fun = "knn", k = 10, rowmax = 0.9)
#' imputed_MLE <- impute(norm, fun = "MLE")
#'
#' imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
#' @export
impute_customized <- function(se, fun = c("bpca", "knn", "QRILC", "MLE",
                               "MinDet", "MinProb", "man", "min", "zero",
                               "mixed", "nbavg"), ...) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(fun))

  # Show error if inputs do not contain required columns
  fun <- match.arg(fun)
  
  if(any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  
  # Annotate whether or not there are missing values and how many
  rowData(se)$imputed <- apply(is.na(assay(se)), 1, any)
  rowData(se)$num_NAs <- rowSums(is.na(assay(se)))
  
  # Don't impute rows with all missing values
  se <- se[!rowData(se)$num_NAs == dim(se)[2],]
  
  # Show error if there are no missing values
  if(!any(is.na(assay(se)))) {
    warning("No missing values in '", deparse(substitute(se)), "'. ",
            "Returning the unchanged object.",
            call. = FALSE)
    return(se)
  }
  
  # if the "man" function is selected, use the manual imputation method
  if(fun == "man") {
    se <- manual_impute_customized(se, ...)
  }
  # else use the MSnSet::impute function
  else {
    MSnSet_data <- as(se, "MSnSet")
    MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun, ...)
    assay(se) <- MSnbase::exprs(MSnSet_imputed)
  }
  return(se)
}

# modified from LFQ-Analyst's plot_volcano_new
plot_volcano_customized <- function(dep, contrast, label_size = 3,
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
  # if(any(!c("name", "ID") %in% colnames(row_data))) {
  #   stop(paste0("'name' and/or 'ID' columns are not present in '",
  #               deparse(substitute(dep)),
  #               "'.\nRun make_unique() to obtain required columns."),
  #        call. = FALSE)
  # }
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
  signif <- grep(paste("^",contrast, "_significant", sep = ""),
                 colnames(row_data))
  df_tmp <- data.frame(diff = row_data[, diff],
                       p_values = -log10(row_data[, p_values]),
                       signif = row_data[, signif],
                       name = row_data$name)
  df<- df_tmp %>% data.frame() %>% filter(!is.na(signif)) %>%
    arrange(signif)
  
  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  #return(df)
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
    theme_DEP1() +
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
      select(name, diff, p_value, signif) %>%
      arrange(desc(x))
    colnames(df)[c(1,2)] <- c("protein", "log2_fold_change")
    if(adjusted) {
      colnames(df)[3] <- "adjusted_p_value_-log10"
    } else {
      colnames(df)[3] <- "p_value_-log10"
    }
    return(df)
  }
}


# modified from DEP's plot_cor
# https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/plot_functions_explore.R
plot_cor_customized <- function(dep, significant = TRUE, lower = -1, upper = 1,
                                pal = "PRGn", pal_rev = FALSE, indicate = NULL,
                                font_size = 12, plot = TRUE, exp="LFQ", ...) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.logical(significant),
                          length(significant) == 1,
                          is.numeric(lower),
                          length(lower) == 1,
                          is.numeric(upper),
                          length(upper) == 1,
                          is.character(pal),
                          length(pal) == 1,
                          is.logical(pal_rev),
                          length(pal_rev) == 1,
                          is.numeric(font_size),
                          length(font_size) == 1,
                          is.logical(plot),
                          length(plot) == 1)
  
  # Check for valid lower and upper values
  if(!(lower >= -1  & upper >= -1 & lower <= 1 & upper <= 1)) {
    stop("'lower' and/or 'upper' arguments are not valid
         Run plot_pca() with 'lower' and 'upper' between -1 and 1",
         call. = FALSE)
  }
  
  # Check for valid pal
  pals <- RColorBrewer::brewer.pal.info %>%
    rownames_to_column() %>%
    filter(category != "qual")
  if(!pal %in% pals$rowname) {
    stop("'", pal,"' is not a valid color panel",
         " (qualitative panels also not allowed)\n",
         "Run plot_pca() with one of the following 'pal' options: ",
         paste(pals$rowname, collapse = "', '"), "'",
         call. = FALSE)
  }
  
  # if(any(is.na(assay(dep)))) {
  #   stop("Missing values in '", deparse(substitute(dep)), "'. Use plot_dist() instead")
  # }
  
  # Heatmap annotation
  if(!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))
    
    col_data <- colData(dep) %>%
      as.data.frame()
    columns <- colnames(col_data)
    if(any(!indicate %in% columns)) {
      stop("'",
           paste0(indicate, collapse = "' and/or '"),
           "' column(s) is/are not present in ",
           deparse(substitute(dep)),
           ".\nValid columns are: '",
           paste(columns, collapse = "', '"),
           "'.",
           call. = FALSE)
    }
    
    # Get annotation
    anno <- colData(dep) %>%
      data.frame() %>%
      select(indicate)
    
    # Annotation color
    names <- colnames(anno)
    anno_col <- vector(mode="list", length=length(names))
    names(anno_col) <- names
    for(i in names) {
      var <- anno[[i]] %>%
        unique() %>%
        sort()
      if (length(var) == 1) {
        cols <- c("black")
      } else if (length(var) == 2) {
        cols <- c("orangered", "cornflowerblue")
      } else if (length(var) < 7 & length(var) > 2) {
        cols <- brewer.pal(length(var), "Pastel1")
      } else if (length(var) <= 12) {
        cols <- brewer.pal(length(var), "Set3")
      } else {
        cols <- colorRampPalette(brewer.pal(12, "Set3"))(length(var))
      }
      names(cols) <- var
      anno_col[[i]] <-  cols
    }

    # HeatmapAnnotation object
    ha1 = HeatmapAnnotation(df = anno,
                            col = anno_col,
                            show_annotation_name = TRUE)
  } else {
    ha1 <- NULL
  }
  
  # Filter for significant proteins
  if(significant) {
    
    # Check for significant column
    if(!"significant" %in% colnames(rowData(dep, use.names = FALSE))) {
      stop("'significant' column is not present in '",
           deparse(substitute(dep)),
           "'\nRun add_rejections() to obtain the required column",
           call. = FALSE)
    }
    
    dep <- dep[replace_na(rowData(dep, use.names = FALSE)$significant, FALSE), ]
  }
  
  # Calculate correlation matrix
  data <- assay(dep)

  # use sample name for the heatmap
  temp <- colData(dep)
  rownames(temp) <- temp$label
  colnames(data) <- temp[colnames(data), "sample_name"]

  cor_mat <- cor(data, use="complete.obs")
  lower <- min(cor_mat)
  upper <- max(cor_mat)

  # Plot heatmap
  ht1 = Heatmap(cor_mat,
                col = circlize::colorRamp2(
                  seq(lower, upper, ((upper-lower)/7)),
                  if(pal_rev) {
                    rev(RColorBrewer::brewer.pal(8, pal))
                  } else {
                    RColorBrewer::brewer.pal(8, pal)
                  }),
                heatmap_legend_param = list(
                  color_bar = "continuous",
                  legend_direction = "horizontal",
                  legend_width = unit(5, "cm"),
                  title_position = "topcenter"),
                name = "Pearson correlation",
                column_names_gp = gpar(fontsize = font_size),
                row_names_gp = gpar(fontsize = font_size),
                top_annotation = ha1,
                ...)
  if(plot) {
    draw(ht1, heatmap_legend_side = "top")
  } else {
    df <- as.data.frame(cor_mat)
    return(df)
  }
}

# similar to plot_cor but generate Jaccard similarity matrix for protein identification results
plot_Jaccard <- function(dep, lower = -1, upper = 1,
                                pal = "PRGn", pal_rev = FALSE, indicate = NULL,
                                font_size = 12, plot = TRUE, exp="LFQ", ...) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.numeric(lower),
                          length(lower) == 1,
                          is.numeric(upper),
                          length(upper) == 1,
                          is.character(pal),
                          length(pal) == 1,
                          is.logical(pal_rev),
                          length(pal_rev) == 1,
                          is.numeric(font_size),
                          length(font_size) == 1,
                          is.logical(plot),
                          length(plot) == 1)

  # Check for valid lower and upper values
  if(!(lower >= -1  & upper >= -1 & lower <= 1 & upper <= 1)) {
    stop("'lower' and/or 'upper' arguments are not valid
         Run plot_pca() with 'lower' and 'upper' between -1 and 1",
         call. = FALSE)
  }

  # Check for valid pal
  pals <- RColorBrewer::brewer.pal.info %>%
    rownames_to_column() %>%
    filter(category != "qual")
  if(!pal %in% pals$rowname) {
    stop("'", pal,"' is not a valid color panel",
         " (qualitative panels also not allowed)\n",
         "Run plot_pca() with one of the following 'pal' options: ",
         paste(pals$rowname, collapse = "', '"), "'",
         call. = FALSE)
  }

  # if(any(is.na(assay(dep)))) {
  #   stop("Missing values in '", deparse(substitute(dep)), "'. Use plot_dist() instead")
  # }

  # Heatmap annotation
  if(!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))

    col_data <- colData(dep) %>%
      as.data.frame()
    columns <- colnames(col_data)
    if(any(!indicate %in% columns)) {
      stop("'",
           paste0(indicate, collapse = "' and/or '"),
           "' column(s) is/are not present in ",
           deparse(substitute(dep)),
           ".\nValid columns are: '",
           paste(columns, collapse = "', '"),
           "'.",
           call. = FALSE)
    }

    # Get annotation
    anno <- colData(dep) %>%
      data.frame() %>%
      select(indicate)

    # Annotation color
    names <- colnames(anno)
    anno_col <- vector(mode="list", length=length(names))
    names(anno_col) <- names
    for(i in names) {
      var = anno[[i]] %>% unique() %>% sort()
      if(length(var) == 1)
        cols <- c("black")
      if(length(var) == 2)
        cols <- c("orangered", "cornflowerblue")
      if(length(var) < 7 & length(var) > 2)
        cols <- RColorBrewer::brewer.pal(length(var), "Pastel1")
      if(length(var) >= 7)
        cols <- RColorBrewer::brewer.pal(length(var), "Set3")
      names(cols) <- var
      anno_col[[i]] <-  cols
    }
    
    # HeatmapAnnotation object
    ha1 = HeatmapAnnotation(df = anno,
                            col = anno_col,
                            show_annotation_name = TRUE)
  } else {
    ha1 <- NULL
  }
  
  # Calculate Jaccard index matrix
  data <- assay(dep)
  
  # use sample name for the heatmap
  temp <- colData(dep)
  rownames(temp) <- temp$label
  colnames(data) <- temp[colnames(data), "sample_name"]
  
  cor_mat <- 1 - as.matrix(vegdist(t(data), method="jaccard", na.rm=T))
  
  lower <- min(cor_mat)
  upper <- max(cor_mat)
  
  # Plot heatmap
  ht1 = Heatmap(cor_mat,
                col = circlize::colorRamp2(c(lower, (upper+lower)/2, upper), c("blue", "lightyellow", "red")),
                # circlize::colorRamp2(
                #   seq(lower, upper, ((upper-lower)/7)),
                #   if(pal_rev) {
                #     rev(RColorBrewer::brewer.pal(8, pal))
                #   } else {
                #     RColorBrewer::brewer.pal(8, pal)
                #   })
                heatmap_legend_param = list(
                  color_bar = "continuous",
                  legend_direction = "horizontal",
                  legend_width = unit(5, "cm"),
                  title_position = "topcenter"),
                name = "Jaccard similarity",
                column_names_gp = gpar(fontsize = font_size),
                row_names_gp = gpar(fontsize = font_size),
                top_annotation = ha1,
                ...)
  if(plot) {
    draw(ht1, heatmap_legend_side = "top")
  } else {
    df <- as.data.frame(cor_mat)
    return(df)
  }
}

# customized from https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/functions.R
#' Filter on missing values
#'
#' \code{filter_missval_customized} filters a proteomics dataset based on missing values.
#' The dataset is filtered for proteins that have a maximum of 'thr' missing values in at least one condition.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}).
#' @param thr Integer(1),
#' Sets the threshold for the allowed number of missing values
#' in at least one condition.
#' @return A filtered SummarizedExperiment object.
#' @examples
#' # Load example
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' # Filter
#' stringent_filter <- filter_missval_customized(se, thr = 0)
#' less_stringent_filter <- filter_missval_customized(se, thr = 1)
#' @export
filter_missval_customized <- function(se, thr = 0) {
  # Show error if inputs are not the required classes
  if(is.integer(thr)) thr <- as.numeric(thr)
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.numeric(thr),
                          length(thr) == 1)
  
  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(colData(se)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  max_repl <- max(colData(se)$replicate)
  if(thr < 0 | thr > max_repl) {
    stop("invalid filter threshold applied",
         "\nRun filter_missval() with a threshold ranging from 0 to ",
         max_repl)
  }
  
  # Make assay values binary (1 = valid value)
  bin_data <- assay(se)
  idx <- is.na(assay(se))
  bin_data[!idx] <- 1
  bin_data[idx] <- 0
  
  # print(colData(se))
  # Filter se on the maximum allowed number of
  # missing values per condition (defined by thr)
  keep <- bin_data %>%
    data.frame(check.names = F) %>%
    rownames_to_column() %>%
    gather(label, value, -rowname) %>%
    left_join(., data.frame(colData(se)), by = "label", check.names=F) %>%
    group_by(rowname, condition) %>%
    summarize(miss_val = n() - sum(value)) %>%
    filter(miss_val <= thr) %>%
    spread(condition, miss_val)
  se_fltrd <- se[keep$rowname, ]
  return(se_fltrd)
}

# customized from https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/plot_functions_QC.R

#' Plot a heatmap of proteins with missing values
#'
#' \code{plot_missval_customized} generates a heatmap of proteins
#' with missing values to discover whether values are missing by random or not.
#'
#' @param se SummarizedExperiment,
#' Data object with missing values.
#' @return A heatmap indicating whether values are missing (0) or not (1)
#' (generated by \code{\link[ComplexHeatmap]{Heatmap}}).

plot_missval_customized <- function(se) {
  # Show error if input is not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  
  se_assay <- assay(se)
  # Show error if there are no missing values
  if(!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)), "'",
         call. = FALSE)
  }
  
  # Make assay data binary (1 = valid value, 0 = missing value)
  df <- se_assay %>% data.frame(.,check.names = F)

  missval <- df[apply(df, 1, function(x) any(is.na(x))), ]
  missval <- ifelse(is.na(missval), 0, 1)

  # use sample name for the heatmap
  temp <- colData(se)
  rownames(temp) <- temp$label
  colnames(missval) <- temp[colnames(missval), "sample_name"]

  # Plot binary heatmap
  if (dim(missval)[1] >= 65536) {
    dist <- factoextra::get_dist(missval, "euclidean")
    mat.hc <- fastcluster::hclust(dist, method = "complete")
    mat.dend <- as.dendrogram(mat.hc)
    ht2 <- Heatmap(missval,
                   col = c("white", "black"),
                   cluster_rows = mat.dend,
                   column_names_side = "top",
                   show_row_names = F,
                   show_column_names = T,
                   show_row_dend = F,
                   name = paste0("Missing values pattern (", dim(missval)[1], " proteins )"),
                   column_names_gp = gpar(fontsize = 16),
                   heatmap_legend_param = list(at = c(0, 1),
                                               labels = c("Missing value", "Valid value")))
  } else {
    ht2 <- Heatmap(missval,
                   col = c("white", "black"),
                   column_names_side = "top",
                   show_row_names = FALSE,
                   show_column_names = TRUE,
                   name = paste0("Missing values pattern (", dim(missval)[1], " proteins )"),
                   column_names_gp = gpar(fontsize = 16),
                   heatmap_legend_param = list(at = c(0, 1),
                                               labels = c("Missing value", "Valid value")))
  }

  draw(ht2, heatmap_legend_side = "top")
}

# make_unique from DEP
# https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/functions.R
#' \code{make_unique} generates unique identifiers
#' for a proteomics dataset based on "name" and "id" columns.
#' @export
make_unique <- function(proteins, names, ids, delim = ";") {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins),
                          is.character(names),
                          length(names) == 1,
                          is.character(ids),
                          length(ids) == 1,
                          is.character(delim),
                          length(delim) == 1)
  
  col_names <- colnames(proteins)
  # Show error if inputs do not contain required columns
  if(!names %in% col_names) {
    stop("'", names, "' is not a column in '",
         deparse(substitute(proteins)), "'",
         call. = FALSE)
  }
  if(!ids %in% col_names) {
    stop("'", ids, "' is not a column in '",
         deparse(substitute(proteins)), "'",
         call. = FALSE)
  }
  
  # If input is a tibble, convert to data.frame
  if(tibble::is_tibble(proteins))
    proteins <- as.data.frame(proteins)
  
  # Select the name and id columns, and check for NAs
  double_NAs <- apply(proteins[,c(names, ids)], 1, function(x) all(is.na(x)))
  if(any(double_NAs)) {
    stop("NAs in both the 'names' and 'ids' columns")
  }
  
  # Take the first identifier per row and make unique names.
  # If there is no name, the ID will be taken.
  proteins_unique <- proteins %>%
    mutate(name = gsub(paste0(delim, ".*"), "", get(names)),
           ID = gsub(paste0(delim, ".*"), "", get(ids)),
           name = ifelse(name == "" | is.na(name), ID, name))
  return(proteins_unique)
}


# normalize_vsn from DEP
# https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/functions.R
#' Normalization using vsn
#' @export
normalize_vsn <- function(se) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  
  # Variance stabilization transformation on assay data
  se_vsn <- se
  vsn.fit <- vsn::vsnMatrix(2 ^ assay(se_vsn))
  assay(se_vsn) <- vsn::predict(vsn.fit, 2 ^ assay(se_vsn))
  return(se_vsn)
}

# test_diff from DEP
# https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/functions.R
#' performs a differential enrichment test based on
#' protein-wise linear models and empirical Bayes
#' statistics using \pkg{limma}. False Discovery Rates are estimated
#' using \pkg{fdrtool}.
#' @export
test_diff <- function(se, type = c("control", "all", "manual"),
                      control = NULL, test = NULL,
                      design_formula = formula(~ 0 + condition)) {
  
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type),
                          class(design_formula) == "formula")
  
  # Show error if inputs do not contain required columns
  type <- match.arg(type)
  
  col_data <- colData(se)
  raw <- assay(se)
  
  if(any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
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
    res <- topTable(fit, sort.by = "t", coef = comp,
                    number = Inf, confint = TRUE)
    res <- res[!is.na(res$t),]
    fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- rownames_to_column(res)
    return(res)
  }
  
  # Retrieve the differential expression test results
  limma_res <- map_df(cntrst, retrieve_fun)
  
  # Select the logFC, CI and qval variables
  table <- limma_res %>%
    select(rowname, logFC, CI.L, CI.R, P.Value, qval, comparison) %>%
    mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
    gather(variable, value, -c(rowname,comparison)) %>%
    mutate(variable = recode(variable, logFC = "diff", P.Value = "p.val", qval = "p.adj")) %>%
    unite(temp, comparison, variable) %>%
    spread(temp, value)
  rowData(se) <- merge(rowData(se, use.names = FALSE), table,
                       by.x = "name", by.y = "rowname", all.x = TRUE, sort=FALSE)
  return(se)
}

# theme_DEP1 from DEP https://github.com/arnesmits/DEP/blob/b425d8d0db67b15df4b8bcf87729ef0bf5800256/R/functions.R
#' @export
theme_DEP1 <- function() {
  # Use theme_bw() as default
  basesize <- 12
  theme <- ggplot2::theme_bw(base_size = basesize)
  
  # Change plot title appearance
  theme$plot.title$face <- "bold"
  theme$plot.title$size <- basesize + 2
  theme$plot.title$hjust <- 0.5
  
  # Change axis title appearance
  theme$axis.title.x$size <- basesize + 2
  
  theme$axis.title.y$size <- basesize + 2
  
  # Change axis text appearance
  theme$axis.text$size <- basesize
  theme$axis.text$colour <- "black"
  
  # Change legend title appearance
  theme$legend.title$size <- basesize + 2
  
  # Change legend text appearance
  theme$legend.text$size <- basesize
  
  # Change strip text (facet headers) appearance
  theme$strip.text$face <- "bold"
  theme$strip.text$size <- basesize + 2
  theme$strip.text$colour <- "black"
  
  return(theme)
}

theme_DEP2 <- function() {
  # Get vertical x-axis labels
  theme <- theme_DEP1()
  theme$axis.text.x$angle <- 90
  theme$axis.text.x$hjust <- 1
  theme$axis.text.x$vjust <- 0.5
  return(theme)
}

#' Generate a wide data.frame from a SummarizedExperiment (from DEP)
#'
#' \code{get_df_wide} generate a wide data.frame from a SummarizedExperiment.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}).
#' @return A data.frame object
#' containing all data in a wide format,
#' where each row represents a protein.
#' @examples
#' @export
get_df_wide <- function(se) {
  # Show error if inputs are not the required classes
  assert_that(inherits(se, "SummarizedExperiment"))
  
  # Show error if inputs do not contain required columns
  if (!"name" %in% colnames(rowData(se, use.names = FALSE))) {
    stop("'name' column is not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  
  # Extract row data
  row_data <- rowData(se, use.names = FALSE) %>%
    data.frame()
  # Extract assay data
  assay_data <- assay(se) %>%
    data.frame() %>%
    rownames_to_column()
  colnames(assay_data)[1] <- "name"
  
  # Merge row and assay data into a wide data.frame
  wide <- full_join(assay_data, row_data, by = "name")
  
  return(wide)
}