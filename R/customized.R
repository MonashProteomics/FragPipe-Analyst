# customized function from FragPipe-Analyst (https://github.com/hsiaoyi0504/FragPipe-Analyst)

plot_pca_plotly <- function(dep, x = 1, y = 2, indicate = c("condition", "replicate"),
                            label = FALSE, n = 500, point_size = 8, label_size = 3, plot = TRUE, ID_col="ID") {
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
  
  # Check for valid 'n' value
  if(n > nrow(dep)) {
    stop(paste0("'n' argument is not valid.\n",
                "Run plot_pca() with 'n' <= ",
                nrow(dep),
                "."),
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
  
  # Get the variance per protein and take the top n variable proteins
  var <- apply(assay(dep), 1, sd)
  df <- assay(dep)[order(var, decreasing = TRUE)[seq_len(n)],]
  
  # Calculate PCA
  pca <- prcomp(t(df), scale = FALSE)
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
    p <- plot_ly(data=pca_df, type = 'scatter', mode = 'markers', marker = list(size = point_size)) %>%
      plotly::layout(title = 'PCA plot', xaxis = list(title = paste0("PC", x, ": ", percent[x], "%")), yaxis = list(title = paste0("PC", y, ": ", percent[y], "%"))) %>%
      add_trace(type = "scatter",
                x = ~PC1, 
                y = ~PC2,
                color = as.formula(paste0('~', indicate[1])),
                mode = 'markers',
                legendgroup=indicate[1],
                legendgrouptitle_text=indicate[1])
  }
  if(length(indicate) == 2) {
    p <- plot_ly(data=pca_df, type = 'scatter',
                 mode = 'markers', marker = list(size = point_size), text=~rowname) %>%
      #Overlay color for gears
      add_trace(type = "scatter",
                x = ~PC1, 
                y = ~PC2,
                symbol = as.formula(paste0('~', indicate[2])),
                marker = list(color = "grey", size = point_size + 3),
                mode = 'markers',
                legendgroup=indicate[2],
                legendgrouptitle_text=indicate[2]) %>%
      add_trace(type = "scatter",
                x = ~PC1, 
                y = ~PC2,
                color = as.formula(paste0('~', indicate[1])),
                mode = 'markers',
                legendgroup=indicate[1],
                legendgrouptitle_text=indicate[1]) %>%
      plotly::layout(title = 'PCA plot',
                     xaxis = list(title = paste0("PC", x, ": ", percent[x], "%")),
                     yaxis = list(title = paste0("PC", y, ": ", percent[y], "%")),
                     legend=list(itemclick = FALSE,
                                 itemdoubleclick = FALSE,
                                 groupclick = FALSE))
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