# customized functions to make LFQ-FP functional

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
make_se_customized <- function(proteins_unique, columns, expdesign, log2transform=F) {
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
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  if (log2transform) {
    raw <- log2(raw)
  }
  # Generate the colData from the experimental design
  # and match these with the assay data
  expdesign <- mutate(expdesign, condition = make.names(condition)) %>%
    unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  
  matched <- match(make.names(delete_prefix(expdesign$label)),
                   make.names(delete_prefix(colnames(raw))))
  if(any(is.na(matched))) {
    stop("None of the labels in the experimental design match ",
         "with column names in 'proteins_unique'",
         "\nRun make_se() with the correct labels in the experimental design",
         "and/or correct columns specification")
  }
  
  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
  
  # Select the rowData
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name
  
  # Generate the SummarizedExperiment object
  se <- SummarizedExperiment(assays = as.matrix(raw),
                             colData = expdesign,
                             rowData = row_data)
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
  print(made_contrasts)
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
  rowData(se) <- merge(rowData(se, use.names = FALSE), table,
                       by.x = "name", by.y = "rowname", all.x = TRUE, sort=FALSE)
  return(se)
}