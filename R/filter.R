global_filter <- function(se, percentage=50){
  percentage <- percentage / 100
  ridx <- rowSums(is.na(assay(se))) / ncol(assay(se)) <= percentage
  se <- se[ridx,]
  return(se)
}

filter_by_condition <- function(se, min_percentage=50) {
  min_percentage <- min_percentage / 100
  conditions <- unique(colData(se)$condition)
  row_ids <- rep(0, nrow(assay(se)))
  for (c in conditions){
    se_c <- se[,colData(se)$condition == c]
    ridx <- rowSums(!is.na(assay(se_c))) / ncol(assay(se_c)) >= min_percentage
    row_ids <- row_ids + ridx
  }
  se <- se[row_ids > 0,]
  return(se)
}