library(readxl)

mapping <- read_excel("./data/CPTAC-Pancan-Data_Case_sample_ID_mapping.xlsx", sheet="Tab1.CPTAC3_DNA_RNA_Proteome")
mapping <- mapping[mapping$`Tumor Type` == "LUAD" & mapping$`Sample Type` == "Tissue for Protein",]
mapping <- mapping[!is.na(mapping$`Case ID`),]
mapping <- as.data.frame(mapping)
rownames(mapping) <- mapping$`Aliquot ID`

data <- read.csv("./data/luad.annotation.raw.tsv", sep="\t", stringsAsFactors = F)
data$condition <- mapping[data$ID, "Category"]
write.table(data, "./data/luad.annotation.updated.tsv", sep="\t", quote = F, na = "", row.names = F)
