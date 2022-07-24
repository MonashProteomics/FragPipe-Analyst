source("./R/tests.R")
source("./R/customized.R")
source("./R/functions.R")
library(DEP)

temp_data <- read.table("./data/luad_ratio_gene_MD_phospho_glyco.tsv",
                      header = TRUE,
                      fill= TRUE, # to fill any missing data
                      sep = "\t",
                      check.names = F)

temp_df <- read.table("./data/luad.annotation.updated_20220722.tsv",
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

plot_numbers_customized(data_se)

# DE analysis
# data_diff_manual <- test_diff_customized(data_se, type = "manual", 
#                               test = c("SampleTypeTumor"), design_formula = formula(~0+SampleType))
data_diff_manual <- test_diff_customized(data_se, type = "all")
dep <- add_rejections(data_diff_manual, alpha=0.05, lfc=log2(1.5))
data_results <- get_results_customized(dep)

# impute
imputed_data <- DEP::impute(data_se, "man")


