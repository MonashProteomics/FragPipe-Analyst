source("./R/tests.R")
source("./R/customized.R")
source("./R/functions.R")
library(DEP)

temp_data <- read.table("./data/DIA_datasets/diann-output.pg_matrix.tsv",
                      header = TRUE,
                      fill= TRUE, # to fill any missing data
                      sep = "\t",
                      check.names = F)

temp_df <- read.table("./data/DIA_datasets/dia_manifest.fp-manifest",
                      header = F,
                      sep="\t",
                      stringsAsFactors = FALSE)
colnames(temp_df) <- c("path", "condition", "replicate", "Data.type")
temp_df$label <- temp_df$path


data_unique <- DEP::make_unique(temp_data, "Genes", "Protein.Group")
cols <- colnames(data_unique)
selected_cols <- which(!(cols %in% c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description", "ID", "name")))
test_match_tmt_column_design(data_unique, selected_cols, temp_df)
data_se <- make_se_customized(data_unique, selected_cols, temp_df, log2transform=T)

imputed <- impute_customized(data_se, "man")

DE_result <- test_limma_customized(imputed, type = "all")
dep <- add_rejections(DE_result, alpha=0.05, lfc=log2(1.5))
data_results <- get_results_proteins(dep, "DIA")
data_results[data_results["Gene Name"] == "CA9",]

DE_result2 <- test_diff_customized(imputed, type = "all")
dep2 <- add_rejections(DE_result2, alpha=0.05, lfc=log2(1.5))
data_results2 <- get_results_proteins(dep2, "DIA")
# data_results <- get_results_customized(dep)
data_results2[data_results2["Gene Name"] == "CA9",]
