###############################
#' TODO: Prioritization of CYT regulators based on 
#' weights in the fully connected layer of the Transformer and regression coefficient
#' 
#'###############################################################################
setwd("/~")

###
# 1. Data preparation --------------------
#' 
#'###############################################################################
mutation_status_mat <- read.csv("./data/Solid_tumor_mutation_status_mat.csv", row.names = 1, check.names = FALSE) # 487 8223
CYT_score_df <- read.csv("./data/Solid_tumor_CYT_score_df.csv", row.names = 1) # 8223 1

### sample cancer type annotations 
sample_annotations <- read.csv("./data/merged_sample_quality_annotations.tsv", sep = "\t") # 79286    12
sample_annotations$sample_barcode <- substr(sample_annotations$aliquot_barcode, 1, 16) 
CYT_score_df$cancer_type <- sample_annotations$'cancer.type'[match(rownames(CYT_score_df), sample_annotations$sample_barcode)]

### data seperate
CYT_score_df_list <- split(CYT_score_df, CYT_score_df$cancer_type)
load(file = "./data/mutation_status_mat_list.RData")


### identified CYT regulators
load(file = "./Results/nonzero_coef_return.RData") # regression coefficient beta
load(file = "./Results/non_zero_variables.RData")
load(file = "./Results/embeddings_glmnet_fit_beta_list.RData")


### immune pathways
# 17 immune pathways from ImmPort
ImmPort_GeneList <- read.table("./data/ImmunePathways/ImmPort_GeneList.txt", sep = "\t", header = TRUE) # 2013 6
# ImmuneSigDB_GeneList <- clusterProfiler::read.gmt("./data/ImmunePathways/c7.all.v2023.1.Hs.symbols.gmt") # 2013 6
ImmuneSigDB_GeneList <- clusterProfiler::read.gmt("./data/ImmunePathways/c7.immunesigdb.v2023.1.Hs.symbols.gmt") # 2013 6
# CYS_signature
load(file = "./data/ImmunePathways/CYS_signature.RData")


CYT_Benchmark <- read.table("./data/CYT_Benchmark.txt", sep = "\t", header = TRUE) # 50 5
Experimental_TIME_driver <- read.table("./data/Experimental TIME driver.txt", sep = "\t", header = TRUE)



###
# 2. benchmark map rank --------------------
#' 
#'###############################################################################
### weights in the fully connected layer of the Transformer
fc_weight <- unlist(read.csv("./code_by_jby_for_attention/Results/fc1weight.csv", header = FALSE)) # 1  487
names(fc_weight) <- rownames(mutation_status_mat)
fc_weight <- fc_weight[order(fc_weight, decreasing = TRUE)]

scale_fcweight <- (fc_weight-min(fc_weight))/(max(fc_weight)-min(fc_weight))
scale_fcweight <- scale_fcweight[order(scale_fcweight, decreasing = TRUE)]


idx <- match(CYT_Benchmark$'Experimental.TIME.driver', names(scale_fcweight))
# names(scale_fcweight)[sort(na.omit(idx))]

non_zero_variables_beta <- data.frame(Gene = names(nonzero_coef_return$beta), 
                                      Coef = nonzero_coef_return$beta,
                                      abs_Coef = abs(nonzero_coef_return$beta))
non_zero_variables_beta <- non_zero_variables_beta[-which(non_zero_variables_beta$Gene %in% c("cp1", "cp2")),]
non_zero_variables_beta$rank <- rank(non_zero_variables_beta$abs_Coef)


idx <- match(CYT_Benchmark$'Experimental.TIME.driver', non_zero_variables_beta$Gene)
non_zero_variables_beta$rank[na.omit(idx)]
# [1] 239   1 224  39 179 238 153 most at tail 
# [1] "MEN1"   "CREBBP" "UBR5"   "SCAF4"  "FAS"    "TGFBR2" "B2M" 



###
# 3. pancancer top-bottom rank comparison --------------------
#'   PPI
#'###############################################################################

### scale_fcweight + absolute beta
integrate_weight_add <- scale_fcweight[non_zero_variables_beta$Gene] + non_zero_variables_beta$abs_Coef
integrate_weight_add <- integrate_weight_add[order(integrate_weight_add, decreasing = TRUE)]

### Analysis in STRING database
cutoff <- round(length(integrate_weight_add) * 0.05)
write.table(c(names(head(integrate_weight_add, n = cutoff)), "PRF1", "GZMA"), 
            file = "./Results/PPI/pancancer_integrate_weight_top0.05.txt",
            sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(c(names(tail(integrate_weight_add, n = cutoff)), "PRF1", "GZMA"),
            file = "./Results/PPI/pancancer_integrate_weight_bottom0.05.txt",
            sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)


### barplot
integrate_weight_df <- c(head(integrate_weight_add, n = cutoff), tail(integrate_weight_add, n = cutoff))
integrate_weight_df <- data.frame(Gene = names(integrate_weight_df), 
                                  integrate_weight = integrate_weight_df)

p1 <- ggplot(integrate_weight_df, aes(x = reorder(Gene, integrate_weight, decreasing = TRUE), y = integrate_weight, fill = integrate_weight, color = integrate_weight)) +
  geom_bar(stat = "identity", color = "white") + 
  labs(x = "", y = "Integrated weight") + 
  scale_fill_gradient(low = "white", high = '#f87669') + 
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1, size = 12))

integrate_weight_df$col_ID <- 24:1

p2 <- ggplot(integrate_weight_df[1:15,], aes(x = reorder(Gene, col_ID, decreasing = TRUE), y = col_ID, fill = col_ID, color = col_ID)) +
  geom_bar(stat = "identity", color = "white") + 
  labs(x = "", y = "Integrated weight") + 
  scale_fill_gradient(low = "white", high = 'firebrick2') + 
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1, size = 14))

p3 <- ggplot(integrate_weight_df[1:15,], aes(x = reorder(Gene, col_ID, decreasing = TRUE), y = col_ID, fill = col_ID, color = col_ID)) +
  geom_bar(stat = "identity", color = "white") + 
  labs(x = "", y = "Integrated weight") + 
  scale_fill_gradient(low = "white", high = 'steelblue') + 
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1, size = 14))

library(patchwork)
pdf("./Results/pancancer_twotail/integrate_weight_twotail_genes.pdf", height = 12, width = 9)
p1 / p2 / p3
dev.off()
