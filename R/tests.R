### Test if column names are proper in experiment design file

exp_design_test<-function(exp_design){
  col_names<-colnames(exp_design)
  ## 
  if(!"label" %in% col_names){
    stop(safeError("The column 'label'(case sensitive) is not found in the Experimental Design File"))
  }
  
  else if (!"condition" %in% col_names){
    stop(safeError("The column 'condition' (case sensitive) is not found in the Experimental Design File"))
  }
  
  else if (!"replicate" %in% col_names){
    stop(safeError("The column 'replicate' (case sensitive) is not found in the Experimental Design File"))
  }
  
}

# mimic the original maxquant part (exp_design_test)
# temp2 <- read.table("../test-LFQ-Analyst/lfq_manifest.fp-manifest", sep="\t", stringsAsFactors = F)
fragpipe_manifest_test<-function(fragpipe_manifest){
  col_names<-colnames(fragpipe_manifest)
  # if(!"label" %in% col_names){
  #   stop(safeError("The column 'label'(case sensitive) is not found in the Experimental Design File"))
  # }
  # 
  # else if (!"condition" %in% col_names){
  #   stop(safeError("The column 'condition' (case sensitive) is not found in the Experimental Design File"))
  # }
  # 
  # else if (!"replicate" %in% col_names){
  #   stop(safeError("The column 'replicate' (case sensitive) is not found in the Experimental Design File"))
  # }
}

# mimic the original maxquant part (maxquant_input_test)
# temp <- read.csv("../test-LFQ-Analyst/combined_protein.tsv", sep="\t", stringsAsFactors = F)
fragpipe_input_test<-function(fragpipe_input){
  col_names<-colnames(fragpipe_input)
  if(!"Gene" %in% col_names){
    print(col_names)
    stop(safeError("The column 'Gene' is not found in the FragPipe combined_protein File"))
  }
  
  else if (any(grepl("Intensity", col_names))==FALSE){
    stop(safeError("Columns ending with 'Intensity' are not found in the FragPipe combined_protein File"))
  }
  
  # More like razor proteins
  # else if (!"Protein.IDs" %in% col_names){
  #   stop(safeError("The column 'Protein IDs' is not found in the MaxQuant proteinGroups File"))
  # }
  
  # else if (!"Reverse" %in% col_names){
  #   stop(safeError("The column 'Reverse' is not found in the MaxQuant proteinGroups File"))
  # }
  
  # else if (!"Potential.contaminant" %in% col_names){
  #   stop(safeError("The column 'Potential contaminant' is not found in the MaxQuant proteinGroups File"))
  # }
  
  # else if (!"Only.identified.by.site" %in% col_names){
  #   stop(safeError("The column 'Only identified by site' is not found in the MaxQuant proteinGroups File"))
  # }
  
  # else if (!"Razor...unique.peptides" %in% col_names){
  #   stop(safeError("The column 'Razor + unique peptides' is not found in the MaxQuant proteinGroups File"))
  # }
  
  else if (!"Description" %in% col_names){
    stop(safeError("The column 'Description' is not found in the FragPipe combined_protein.tsv"))
  }
  
}

fragpipe_DIA_input_test<-function(fragpipe_input){
  col_names<-colnames(fragpipe_input)
  if(!"Genes" %in% col_names){
    stop(safeError("The column 'Genes' is not found in the uploaded PG matrix file"))
  }
  # else if (!"Description" %in% col_names){
  #   stop(safeError("The column 'Description' is not found in the FragPipe combined_protein File"))
  # }

}


tmt_input_test<-function(tmt_input){
  col_names<-colnames(tmt_input)
  if(!"Index" %in% col_names){
    stop(safeError("The column 'Index' is not found in the TMT-I gene report."))
  }
  else if(!"ProteinID" %in% col_names & !"Gene" %in% col_names){
    stop(safeError("The column 'ProteinID' or 'Gene' is not found in the TMT-I gene report."))
  }
  else if(!"NumberPSM" %in% col_names){
    stop(safeError("The column 'NumberPSM' is not found in the TMT-I gene report."))
  }
}

### Test if experimental design names and LFQ column names match
test_match_lfq_column_manifest <-function(unique_data, lfq_columns, exp_design){
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
  expdesign <- exp_design
  expdesign$ID <- expdesign$label
  rownames(expdesign) <- expdesign$ID
  # print(make.names(expdesign$label))
  # print(make.names(colnames(raw)))
  matched <- match(make.names(expdesign$label),
                   make.names(colnames(raw)))
  
  if(any(is.na(matched))) {
    stop(safeError("The labels/'run names' in the experimental design DID NOT match
         with lfq column names in maxquants proteinGroups file
         Run LFQ-Analyst with correct labels in the experimental design"))
  }
}

null_enrichment_test <- function(gsea_result, alpha=0.05){
  if(nrow(gsea_result)==0){
    stop(safeError("The proteins/genes you submitted are not found in the database.
                   Please choose other database or lower thresholds and try again."))
  }
  gsea_df <- gsea_result %>% group_by(contrast, var) %>% dplyr::filter(Adjusted.P.value <= alpha)
  if(nrow(gsea_df)==0){
    stop(safeError("No enriched term found at FDR cutoff 0.05. 
                   Enrichment plot could not be displayed. 
                   However, the results (non-significant hits) can still be accessed 
                   through 'Download table' button."))
  }
}

ids_test<-function(filtered_data){
  if("Evidence.IDs" %in% colnames(filtered_data)){
    filtered_data$`Evidence.IDs`<-stringr::str_trunc(as.character(filtered_data$`Evidence.IDs`), 25000)
  }
  if("MS.MS.IDs" %in% colnames(filtered_data)){
    filtered_data$`MS.MS.IDs`<-stringr::str_trunc(as.character(filtered_data$`MS.MS.IDs`), 25000)
  }
  
  return(filtered_data)
  
}

test_TMT_annotation <- function(df) {
  required_columns <- c("channel", "plex", "sample", "condition", "replicate", "sample_name")
  if (any(!required_columns %in% colnames(df))) {
    return(FALSE)
  }
  return(TRUE)
}

