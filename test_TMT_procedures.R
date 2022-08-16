source("./R/tests.R")
source("./R/customized.R")
source("./R/functions.R")
library(DEP)

temp_data <- read.table("./data/TMT_datasets/luad_ratio_gene_MD_phospho_glyco.tsv",
                      header = TRUE,
                      fill= TRUE, # to fill any missing data
                      sep = "\t",
                      check.names = F)

temp_df <- read.table("./data/TMT_datasets/luad.annotation.updated_20220722.tsv",
                      header = T,
                      sep="\t",
                      stringsAsFactors = FALSE)

temp_exp_design <- temp_df[!is.na(temp_df$condition), ]
temp_exp_design <- temp_exp_design[!temp_exp_design$condition == "",]

data_unique <- DEP::make_unique(temp_data, "Index", "ProteinID")

# handle unmatched columns
overlapped_samples <- intersect(colnames(data_unique), temp_exp_design$label)
data_unique <- data_unique[,colnames(data_unique) %in% c("Index", "NumberPSM", "ProteinID", "MaxPepProb", "ReferenceIntensity", "name", "ID", overlapped_samples)]


temp_exp_design <- temp_exp_design[temp_exp_design$label %in% overlapped_samples,]
cols <- colnames(data_unique)
selected_cols <- which(!(cols %in% c("Index", "NumberPSM", "ProteinID", "MaxPepProb", "ReferenceIntensity", "name", "ID")))
data_unique[selected_cols] <- apply(data_unique[selected_cols], 2, as.numeric)
test_match_tmt_column_design(data_unique, selected_cols, temp_exp_design)
data_se <- make_se_customized(data_unique, selected_cols, temp_exp_design)

plot_numbers_by_plex_set(data_se)

# impute
imputed_data <- DEP::impute(data_se, "man")

DE_result <- test_limma_customized(imputed_data, type = "all")
dep <- add_rejections(DE_result, alpha=0.05, lfc=log2(1.5))
data_results <- get_results_proteins(dep, "TMT")
data_results[data_results["Gene Name"] == "CA9",]

DE_result2 <- test_diff_customized(imputed_data, type = "all")
dep2 <- add_rejections(DE_result2, alpha=0.05, lfc=log2(1.5))
data_results2 <- get_results_proteins(dep2, "TMT")
# data_results <- get_results_customized(dep)
data_results2[data_results2["Gene Name"] == "CA9",]

plot_pca_plotly(DE_result, n=500, indicate = "condition", ID_col="label")

tested_contrasts<- gsub("_p.adj", "", 
                        colnames(SummarizedExperiment::rowData(dep2))[grep("p.adj",
                                                                           colnames(SummarizedExperiment::rowData(dep2)))])

plot_volcano_customized(dep2, contrast = tested_contrasts[1],label_size = 2, add_names = F)
