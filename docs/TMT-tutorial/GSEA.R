library(clusterProfiler)
library(msigdbr)
library(enrichplot)

df <- read.csv("./Results.csv", header=TRUE)
gene_list <- df$Tumor_vs_Normal_log2.fold.change
names(gene_list) <- df$Gene.Name
gene_list <- na.omit(gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

# Get the hallmark gene set
H <- msigdbr(species = "Homo sapiens", category = "H")
H.select <- dplyr::select(H, gs_name, gene_symbol)

gsea_result <- GSEA(
  geneList = gene_list,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  TERM2GENE = H.select,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

dotplot(gsea_result, x="NES", showCategory=10)
