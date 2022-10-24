source("./R/tests.R")
source("./R/customized.R")
source("./R/functions.R")
library(DEP)
library(plotly)

temp_data <- read.table("C:\\Users\\yihsiao\\Documents\\workspace\\FragPipe_analyst_Datasets\\DIA_datasets\\updated\\diann-output.pg_matrix.tsv",
                      header = TRUE,
                      fill= TRUE, # to fill any missing data
                      sep = "\t",
                      check.names = F)

temp_df <- read.table("C:\\Users\\yihsiao\\Documents\\workspace\\FragPipe_analyst_Datasets\\DIA_datasets\\updated\\dia_manifest.fp-manifest",
                      header = F,
                      sep="\t",
                      stringsAsFactors = FALSE)
colnames(temp_df) <- c("path", "experiment", "replicate", "Data.type")
temp_df$condition <- gsub("_.*", "", temp_df$experiment)
temp_df$label <- temp_df$path
temp_df$label <- temp_df$path


data_unique <- DEP::make_unique(temp_data, "Genes", "Protein.Group")
cols <- colnames(data_unique)
selected_cols <- which(!(cols %in% c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description", "ID", "name")))
test_match_tmt_column_design(data_unique, selected_cols, temp_df)
data_se <- make_se_customized(data_unique, selected_cols, temp_df, log2transform=T)

# non-imputed version
imputed <- data_se
rowData(imputed)$imputed <- apply(is.na(assay(data_se)), 1, any)
rowData(imputed)$num_NAs <- rowSums(is.na(assay(data_se)))
normalised_data <- normalize_vsn(data_se)
plot_numbers_customized(normalised_data, exp="DIA")

DE_result <- test_limma_customized(imputed, type = "all")
dep <- add_rejections(DE_result, alpha=0.05, lfc=log2(1.5))

plot_cvs(dep, id="label", check.names=F)
plot_cor_customized(dep, significant=FALSE, indicate="condition", exp="DIA")
plot_normalization_DIA_customized(data_se,
                                  normalised_data)
plot_missval_customized(data_se, "DIA")

data_results <- get_results_proteins(dep, "DIA")
get_cluster_heatmap(dep,
                    type="centered" ,kmeans = TRUE,
                    k=3, col_limit = 6,
                    indicate = "condition",
                    exp="DIA"
)



# imputed version
imputed <- impute_customized(data_se, "man")
DE_result <- test_limma_customized(imputed, type = "all")
dep2 <- add_rejections(DE_result, alpha=0.05, lfc=log2(1.5))
data_results2 <- get_results_proteins(dep2, "DIA")

