###############################
#' TODO: Identification of CYT regulators based on mutantion embeddings profile
#'      
#'###############################################################################
setwd("/~")

###
# 1. Data preparation --------------------
#' 
#'###############################################################################
embeddings <- read.csv("./Results/embedding.csv", header = FALSE) # 8223  487
mutation_status_mat <- read.csv("./data/Solid_tumor_mutation_status_mat.csv", row.names = 1, check.names = FALSE) # 487 8223

rownames(embeddings) <- colnames(mutation_status_mat)
colnames(embeddings) <- rownames(mutation_status_mat)

CYT_score_df <- read.csv("./data/Solid_tumor_CYT_score_df.csv", row.names = 1) # 8223 1
CYT_Benchmark <- read.table("./data/CYT_Benchmark.txt", sep = "\t", header = TRUE) # 50 5
Experimental_TIME_driver <- read.table("./data/Experimental TIME driver.txt", sep = "\t", header = TRUE)


###
# 2. Lasso‑regularised ordinal regression --------------------
#' 
#'###############################################################################

#'######
# (1) overall --------------------
#'######
library(tidyverse)
CYT_score_df <- CYT_score_df %>% mutate(CYT_score_sperate = cut(CYT_score_df$CYT_score, breaks = quantile(CYT_score_df$CYT_score, probs = seq(0, 1, 1/3)), 
                                       labels = c("low", "medium", "high"), include.lowest = TRUE))

library(glmnetcr)
glmnet_fit <- glmnetcr(x = embeddings, y = CYT_score_df$CYT_score_sperate)
# a0: Intercept sequence of length ‘length(lambda)’
# beta: a ‘nvars x length(lambda)’ matrix of coefficients
# lambda: The actual sequence of ‘lambda’ values used

### nonzero_coef_variables
AIC <- select.glmnetcr(glmnet_fit, which = "AIC") # A numeric value of length one representing the step number
nonzero_coef_return <- nonzero.glmnetcr(glmnet_fit, s = AIC) 
save(nonzero_coef_return, file = "./Results/nonzero_coef_return.RData")
# a0: intercept estimate
# beta: non-zero estimates for variables and ordinal thresholds
non_zero_variables <- setdiff(names(nonzero_coef_return$beta), c("cp1", "cp2"))  # 250
save(non_zero_variables, file = "./Results/non_zero_variables.RData")


#'######
# (2) single cancer type --------------------
#'######
### sample cancer type annotations 
sample_annotations <- read.csv("./data/merged_sample_quality_annotations.tsv", sep = "\t") # 79286    12
sample_annotations$sample_barcode <- substr(sample_annotations$aliquot_barcode, 1, 16) 
CYT_score_df$cancer_type <- sample_annotations$'cancer.type'[match(rownames(CYT_score_df), sample_annotations$sample_barcode)]

### data seperate
CYT_score_df_list <- split(CYT_score_df, CYT_score_df$cancer_type)
save(CYT_score_df_list, file = "./data/CYT_score_df_list.RData")

embeddings_mat_list <- lapply(names(CYT_score_df_list), function(ca){
    idx <- match(rownames(CYT_score_df_list[[ca]]), rownames(embeddings))
    tmp_embeddings_mat <- embeddings[idx, ]

    # mutation frequency
    freq_count <- colSums(tmp_embeddings_mat)
    if(length(which(freq_count == 0)) > 0) tmp_embeddings_mat <- tmp_embeddings_mat[, -which(freq_count == 0)]

        return(as.data.frame(tmp_embeddings_mat))
})
names(embeddings_mat_list) <- names(CYT_score_df_list)
sapply(embeddings_mat_list, ncol) # check number of features

### CYT score turn to ordered factor
library(tidyverse)
CYT_score_df_list <- lapply(CYT_score_df_list, function(x){
    x <- x %>% mutate(CYT_score_sperate = cut(x$CYT_score, breaks = quantile(x$CYT_score, probs = seq(0, 1, 1/3)), 
                                       labels = c("low", "medium", "high"), include.lowest = TRUE)) 
    return(x)
})

### Lasso‑regularised ordinal regression
library(glmnetcr)
embeddings_glmnet_fit_beta_list <- lapply(names(CYT_score_df_list), function(ca){
    glmnet_fit <- glmnetcr(x = embeddings_mat_list[[ca]], y = CYT_score_df_list[[ca]]$CYT_score_sperate)
    # a0: Intercept sequence of length ‘length(lambda)’
    # beta: a ‘nvars x length(lambda)’ matrix of coefficients
    # lambda: The actual sequence of ‘lambda’ values used

    ### nonzero_coef_variables
    AIC <- select.glmnetcr(glmnet_fit, which = "AIC") # A numeric value of length one representing the step number
    nonzero_coef_return <- nonzero.glmnetcr(glmnet_fit, s = AIC) 
    # a0: intercept estimate
    # beta: non-zero estimates for variables and ordinal thresholds
    
    non_zero_variables <- setdiff(names(nonzero_coef_return$beta), c("cp1", "cp2"))
    if(length(non_zero_variables) > 0){
        return(data.frame(Gene = non_zero_variables, Coef = nonzero_coef_return$beta[non_zero_variables]))
    }else{
        return(NULL)
    }
})
names(embeddings_glmnet_fit_beta_list) <- names(CYT_score_df_list)
save(embeddings_glmnet_fit_beta_list, file = "./Results/embeddings_glmnet_fit_beta_list.RData")

### Number of samples in each cancer type
sample_number <- data.frame(table(CYT_score_df$cancer_type))
sample_number <- sample_number[order(sample_number$Freq, decreasing = TRUE),]
sample_number$Var1 <- factor(sample_number$Var1, levels = sample_number$Var1)
save(sample_number, file = "./Results/sample_number.RData")

p <- ggplot(sample_number, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "lightskyblue1") + 
  labs(x = "", y = "Number of samples") + 
  geom_text(aes(label = Freq), size = 5, vjust = -0.5) + 
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1, size = 14))
    
pdf("./Results/Number of samples.pdf", height = 6, width = 12)
p
dev.off()



###
# 3. check results --------------------
#' 
#'###############################################################################
#'######
# (1) Overlap with confirmed CYT regulatory drivers --------------------
#'######

### pan-cancer
length(unique(non_zero_variables)) # 250
length(intersect(non_zero_variables, CYT_Benchmark$'Experimental.TIME.driver')) # 7
# [1] "B2M"    "CREBBP" "FAS"    "MEN1"   "SCAF4"  "TGFBR2" "UBR5" 

### single cancer type
embeddings_glmnet_fit_beta_list <- lapply(names(embeddings_glmnet_fit_beta_list), function(ca){
    beta_df <- embeddings_glmnet_fit_beta_list[[ca]]
    if(length(beta_df) > 0){
        beta_df$cancer_type <- ca
        return(beta_df)
    }else{
        return(NULL)
    }
})

imm_reg_driver_df <- do.call(rbind.data.frame, embeddings_glmnet_fit_beta_list) # 1848    3
length(unique(imm_reg_driver_df$Gene)) # 484
length(intersect(CYT_Benchmark$'Experimental.TIME.driver', imm_reg_driver_df$Gene)) # 15
length(intersect(CYT_Benchmark$'Experimental.TIME.driver', rownames(mutation_status_mat))) # 15 

sort(table(imm_reg_driver_df$cancer_type))
# KIRC MESO PCPG TGCT  LGG READ THCA  UVM PAAD THYM STAD UCEC BRCA COAD HNSC CESC 
#    1    1    1    1    2    2    2    3    4    6    7    9   11   30   78  187 
# LIHC BLCA LUAD LUSC SKCM
#  253  269  288  338  355

confirmed_imm_reg_driver_df <- imm_reg_driver_df[which(imm_reg_driver_df$Gene %in% CYT_Benchmark$'Experimental.TIME.driver'), ]
sort(table(confirmed_imm_reg_driver_df$cancer_type))
# COAD HNSC CESC BLCA LUSC LIHC LUAD SKCM 
#    1    2    6    9    9   10   10   12


# barplot
library(ggpubr)
confirmed_imm_reg_num_df <- as.data.frame(sort(table(confirmed_imm_reg_driver_df$cancer_type)))
colnames(confirmed_imm_reg_num_df)[1] <- "Cancer_type"

p <- ggplot(confirmed_imm_reg_num_df, aes(x = reorder(Cancer_type, Freq, decreasing = TRUE), y = Freq, fill = Cancer_type)) + 
  geom_bar(stat = "identity", width = 0.6) + 
  # scale_fill_manual(values = cell_colors) + 
  geom_text(aes(label = Freq), size = 4, vjust = -0.5) +  
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), legend.position = "none") + # 坐标轴倾斜
  labs(y = "Number of immune regulators", x = "Cancer type", 
       title = "Distribution of confirmed immune regulators identified by immregInformer")
pdf("./Results/Distribution of confirmed immune regulators.pdf", height = 5, width = 6)
p
dev.off()

# dotplot
# ordered
embeddings_ordered_beta_list <- lapply(names(embeddings_glmnet_fit_beta_list), function(ca){
  beta_df <- embeddings_glmnet_fit_beta_list[[ca]]
  if(length(beta_df) > 0){
    # beta_df$cancer_type <- ca
    beta_df <- beta_df[order(beta_df$Coef, decreasing = TRUE),]
    return(beta_df)
  }else{
    return(NULL)
  }
})
names(embeddings_ordered_beta_list) <- names(embeddings_glmnet_fit_beta_list)

# CYT benckmark regulators map
CYT_benckmark_regs_rank <- lapply(embeddings_ordered_beta_list, function(x){
  idx <- match(CYT_Benchmark$'Experimental.TIME.driver', x$Gene)
  # return(x[sort(na.omit(idx)),])
  return(sort(na.omit(idx)))
})
match_num <- sapply(CYT_benckmark_regs_rank, length)

gene_cancer_beta_df <- do.call(rbind.data.frame, lapply(embeddings_ordered_beta_list[which(match_num > 0)], function(x){
  idx <- match(CYT_Benchmark$'Experimental.TIME.driver', x$Gene)
  return(x[na.omit(idx),])
}))
tmp <- sort(table(gene_cancer_beta_df$Gene))
gene_cancer_beta_df$gene_order <- tmp[gene_cancer_beta_df$Gene]
tmp <- sort(table(gene_cancer_beta_df$cancer_type))
gene_cancer_beta_df$cancer_order <- tmp[gene_cancer_beta_df$cancer_type]

library(ggpubr)
p1 <- ggplot(gene_cancer_beta_df, aes(y = reorder(cancer_type, cancer_order, decreasing = TRUE), 
                                      x = reorder(Gene, gene_order, decreasing = TRUE))) + 
  # geom_point(aes(fill = Coef), color = "white", shape = 21, alpha = 0.5, stroke = 2) + 
  geom_point(aes(fill = Coef), color = "black", shape = 21, size = 6) + 
  scale_fill_gradient2(midpoint = 0,  
                       low = '#2fa1dd',
                       mid = "white",
                       high = '#f87669') + 
  labs(y = "Cancer type", x = "Driver genes") +  
  theme_classic2() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1), legend.position = 'top')

pdf(file = "./Results/Confirmed immune regulators in each cancer.pdf", height = 5, width = 6)
p1
dev.off()



#'######
# (2) Overlap with confirmed immune regulatory drivers --------------------
#'######
### pan-cancer
length(intersect(non_zero_variables, Experimental_TIME_driver$'Experimental.TIME.driver')) # 44 

pos_drivers <- setdiff(names(nonzero_coef_return$beta)[which(nonzero_coef_return$beta > 0)], c("cp1", "cp2"))
neg_drivers <- setdiff(names(nonzero_coef_return$beta)[which(nonzero_coef_return$beta < 0)], c("cp1", "cp2"))

sapply(list(pos_drivers, neg_drivers), function(x){
    # return(length(intersect(x, Experimental_TIME_driver$'Experimental.TIME.driver')))
    return(intersect(x, Experimental_TIME_driver$'Experimental.TIME.driver'))
})



#'######
# (3) Frequency of identified immune regulatory drivers in each cancer --------------------
#'######
imm_reg_drivers_freq <- as.data.frame(table(imm_reg_driver_df$Gene)) 
summary(imm_reg_drivers_freq$Freq)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.000   3.000   4.000   3.818   5.000   7.000 
sapply(2:6, function(n){
    return(length(which(imm_reg_drivers_freq$Freq > n)))
})
# 413 301 139  36   6
imm_reg_drivers_freq[which(imm_reg_drivers_freq$Freq > 6), 1]


# driver mutations in over 2 cancer types
sapply(embeddings_mat_list, ncol) # check number of features
driver_mut_count <- unlist(sapply(embeddings_mat_list, function(x) return(colnames(x))))
driver_mut_freq <- as.data.frame(table(driver_mut_count))
summary(driver_mut_freq$Freq)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.00   20.00   22.00   22.03   24.00   30.00

sapply(2:6, function(n){
    return(length(which(imm_reg_drivers_freq$Freq > n)))
})



#'######
# (4) proportion of identified immune regulators --------------------
#'######
sapply(embeddings_mat_list, ncol) # check number of features
# sapply(mutation_status_mat_list, ncol)
imm_reg_drivers_in_each_cancer <- sapply(embeddings_glmnet_fit_beta_list, nrow)
imm_reg_drivers_in_each_cancer <- sapply(imm_reg_drivers_in_each_cancer, function(x){
    if(is.null(x)) {
        return(0)
    }else{
        return(x)
    }
})

### 
signif(sort(imm_reg_drivers_in_each_cancer/sapply(embeddings_mat_list, ncol)), 3)



###
# 4. Mutation frequency comparison --------------------
#' identified immune regulatory drivers and non-regulators
#'###############################################################################
#'######
# (1) pancancer --------------------
#'######
pancancer_mutation_freq <- rowSums(mutation_status_mat)/ncol(mutation_status_mat)
pancancer_immReg_freq <- data.frame(Gene = names(pancancer_mutation_freq), mut_freq = pancancer_mutation_freq,
                                    CYT_reg = sapply(names(pancancer_mutation_freq), function(x) {
                                            return(is.element(x, non_zero_variables))
                                            }))
## boxplot
library(ggpubr)
library(patchwork)

p <- ggviolin(pancancer_immReg_freq, x = "CYT_reg", y = "mut_freq", color = "CYT_reg",
               fill = "CYT_reg", 
               add = "boxplot", add.params = list(fill = "white", color = "black")) + 
        stat_compare_means() + # Add global p-value
          labs(y = "Mutation frequency", x = "") + 
          theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.5), legend.position = "none")
pdf(file = "./Results/Mutation frequency comparison in pancancer.pdf", height = 4, width = 3)
p
dev.off()