###############################
#' TODO: Data preparation
#' RNA-seq of TCGA solid tumor
#' driver SNV mutation status profile of TCGA solid tumor
#' Cytolitic activity calculation  
#' 
###############################################################################
setwd("/~")

# Working directory
data_home <- "./data/"
###################################################################

###
# 1. SNV profile --------------------
#' https://gdc.cancer.gov/about-data/publications/panimmune
#'###############################################################################
mc3_maf <- read.csv("./data/TCGA_panimmune/mc3.v0.2.8.PUBLIC.maf", sep = "\t")
dim(mc3_maf) # 3600963     114



###
# 2. driver list --------------------
#' Cancer and non-cancer drivers damaged in TCGA.
#'###############################################################################
load("./data/TCGA/combined.solidTumor.driver.RData")



###
# 3. driver SNV profile AND mutation status matrix --------------------
#' 
#'###############################################################################
length(intersect(combined.solidTumor.driver$'Gene', mc3_maf$Hugo_Symbol)) # 487 
mc3_maf_driver <- mc3_maf[which(mc3_maf$Hugo_Symbol %in% combined.solidTumor.driver$'Gene'),]
# mc3_maf_common$Variant_Classification


# Turn to 01 mutation status matrix
mutation_status <- mc3_maf_driver[, c("Hugo_Symbol", "Tumor_Sample_Barcode")]
mutation_status <- unique(mutation_status) # Multiple mutations in same gene are treated as one (count once)
mutation_status$mutation_status <- "1"


library(tidyr)
mutation_status_mat <- spread(mutation_status, key = "Tumor_Sample_Barcode", value = "mutation_status") # 487 9851
rownames(mutation_status_mat) <- mutation_status_mat$Hugo_Symbol
mutation_status_mat <- mutation_status_mat[, -1]
mutation_status_mat <- as.matrix(mutation_status_mat)
mutation_status_mat[which(is.na(mutation_status_mat))] <- "0"



###
# 4. quality-controlled RNA-seq samples --------------------
#' https://gdc.cancer.gov/about-data/publications/panimmune
#'###############################################################################
### sample_annotations
sample_annotations <- read.csv("./data/TCGA_panimmune/merged_sample_quality_annotations.tsv", sep = "\t") # 79286    12

# quality-controlled RNA-seq samples
sample_annotations <- subset(sample_annotations, platform == "IlluminaHiSeq_RNASeqV2" & Do_not_use == "False") # 10258

# solid tumor types
# length(unique(sample_annotations$'cancer.type')) # 33 
ref_cancer_types <- unique(sample_annotations$'cancer.type')
sample_annotations <- subset(sample_annotations, cancer.type %in% setdiff(ref_cancer_types, c("LAML", "DLBC"))) # 10037 12 31 cancer types


### intersection with RNA-seq data
# 20531 11070
RNA_data <- read.csv("./data/TCGA_panimmune/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", sep = "\t", check.names = FALSE, row.names = 1)
common_samples <- intersect(sample_annotations$aliquot_barcode, colnames(RNA_data)) # 9777
RNA_data <- RNA_data[, common_samples]



###
# 5. mutation status matrix of quality-controlled solid tumor RNA-seq samples --------------------
#' 
#'###############################################################################

# 8223 common samples with RNA-seq
colnames(mutation_status_mat) <- substr(colnames(mutation_status_mat), 1, 16)
colnames(RNA_data) <- substr(colnames(RNA_data), 1, 16)

common_samples <- intersect(colnames(mutation_status_mat), colnames(RNA_data)); length(common_samples) # 8223
mutation_status_mat <- mutation_status_mat[, common_samples] # 487 8223
RNA_data_common <- RNA_data[, common_samples] # 20531 8223



save(mutation_status_mat, file = paste0(data_home, "Solid_tumor_mutation_status_mat.RData")) # 487 8223
write.csv(mutation_status_mat, file = paste0(data_home, "Solid_tumor_mutation_status_mat.csv"), quote = FALSE,
          row.names = TRUE)
save(RNA_data_common, file = paste0(data_home, "Solid_tumor_RNA_data.RData"))

a <- mutation_status_mat
storage.mode(a) <- "numeric"
b <- rowSums(a)
summary(b) # Each gene was mutated in at least 11 samples
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 11.0   104.0   169.0   232.1   280.0  3242.0



###
# 6. CYT score calculation --------------------
#' geometric mean
#'###############################################################################
library(psych)

CYT_gene <- c("GZMA|3001", "PRF1|5551")
CYT_exp <- RNA_data_common[CYT_gene, ]
CYT_exp <- log2(CYT_exp + 0.01)


CYT_score <- apply(CYT_exp, 2, function(x){
  idx <- which(x >= 0)
  if(length(idx) == 0){
    return(0)
  } else {
    return(geometric.mean(x[idx]))
  }
})

CYT_score_df <- data.frame(PatientID = names(CYT_score), CYT_score = CYT_score)
save(CYT_score_df, file = paste0(data_home, "Solid_tumor_CYT_score_df.RData"))

write.csv(CYT_score_df, file = paste0(data_home, "Solid_tumor_CYT_score_df.csv"), quote = FALSE,
          row.names = FALSE)
