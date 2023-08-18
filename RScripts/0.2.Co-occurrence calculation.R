###############################
#' TODO: Mutation co-occurrence analysis in pan-cancer samples
#'       
###############################################################################
setwd("/~")

# Working directory
data_home <- "./data/"
###################################################################

###
# 1. Data load --------------------
#' 
#'###############################################################################
# mutation status matrix of quality-controlled solid tumor RNA-seq samples
mutation_status_mat <- get(load(file = paste0(data_home, "Solid_tumor_mutation_status_mat.RData"))) # 487 8223
storage.mode(mutation_status_mat) <- "numeric"
  
# CYT_score_df
CYT_score_df <- get(load(file = paste0(data_home, "Solid_tumor_CYT_score_df.RData")))

# identical(colnames(mutation_status_mat), rownames(CYT_score_df))



###
# 2. mutation co-occurrence analysis --------------------
#' 
#'###############################################################################
source("./RCodes/Function/getChisqTest.R")
source("./RCodes/Function/get_mut_inter_fisher.R")

mutation_cooccurrence_fisher_res <- get_mut_inter_fisher(mutationMat = mutation_status_mat, outfile = NULL) # 118341      2


mutation_cooccurrence_fisher_res$p_adjust <- p.adjust(mutation_cooccurrence_fisher_res$pvalue, method = "BH")
### check results
length(which(mutation_cooccurrence_fisher_res$p_adjust < 0.001))/nrow(mutation_cooccurrence_fisher_res) # 93%
save(mutation_cooccurrence_fisher_res, file = paste0(data_home, "Solid_tumor_mutation_status_cooccurrenceFisher.RData"))
length(which(mutation_cooccurrence_fisher_res$p_adjust > 0.001))



###
# 3. Turn to mutation cooccurrence matrix --------------------
#' 
#'###############################################################################
### Set the oddsratio = 0 when not significant 
idx1 <- which(mutation_cooccurrence_fisher_res$p_adjust > 0.05)
idx2 <- which(mutation_cooccurrence_fisher_res$oddsratio == 0)
mutation_cooccurrence_fisher_res$oddsratio[which(mutation_cooccurrence_fisher_res$p_adjust > 0.05)] <- 0
mutation_cooccurrence_fisher_res$norm_OR <- (mutation_cooccurrence_fisher_res$oddsratio - 0)/max(mutation_cooccurrence_fisher_res$oddsratio)


gene_pairs_df <- do.call(rbind.data.frame, lapply(rownames(mutation_cooccurrence_fisher_res), function(x){
  genes <- unlist(strsplit(x, "_"))
  gene_pairs <- data.frame(gene1 = genes[1], gene2 = genes[2])
  return(gene_pairs)
}))
mutation_cooccurrence_fisher_res <- data.frame(gene_pairs_df, mutation_cooccurrence_fisher_res) 


library(tidyr)
mutation_cooccurrence_mat <- spread(mutation_cooccurrence_fisher_res[, c("gene1", "gene2", "norm_OR")], key = "gene2", value = "norm_OR")
rownames(mutation_cooccurrence_mat) <- mutation_cooccurrence_mat$gene1
mutation_cooccurrence_mat <- mutation_cooccurrence_mat[, -1]
identical(rownames(mutation_cooccurrence_mat)[-1], colnames(mutation_cooccurrence_mat)[-486])



# square matrix 
mutation_cooccurrence_mat <- rbind(cbind(A1CF = NA, mutation_cooccurrence_mat), ZNRF3 = NA)
mutation_cooccurrence_mat <- as.matrix(mutation_cooccurrence_mat)
diag(mutation_cooccurrence_mat) <- 1

all(is.na(lower.tri(mutation_cooccurrence_mat, diag = FALSE)))

upper.tri(mutation_cooccurrence_mat, diag = FALSE)
all(is.na(mutation_cooccurrence_mat[lower.tri(mutation_cooccurrence_mat, diag = FALSE)]))

a <- mutation_cooccurrence_mat
t_a <- t(a)
a[lower.tri(a, diag = FALSE)] <- t_a[lower.tri(t_a, diag = FALSE)]
identical(a[1,], a[, 1])
save(a, file = paste0(data_home, "Solid_tumor_mutation_cooccurrence_mat.RData"))
write.csv(a, file = paste0(data_home, "Solid_tumor_mutation_cooccurrence_mat.csv"), quote = FALSE,
          row.names = TRUE)




###
# 4. Turn to mutation cooccurrence p.adjust matrix --------------------
#' 
#'###############################################################################
gene_pairs_df <- do.call(rbind.data.frame, lapply(rownames(mutation_cooccurrence_fisher_res), function(x){
  genes <- unlist(strsplit(x, "_"))
  gene_pairs <- data.frame(gene1 = genes[1], gene2 = genes[2])
  return(gene_pairs)
}))
mutation_cooccurrence_fisher_res <- data.frame(gene_pairs_df, mutation_cooccurrence_fisher_res) 


library(tidyr)
mutation_cooccurrence_mat <- spread(mutation_cooccurrence_fisher_res[, c("gene1", "gene2", "p_adjust")], key = "gene2", value = "p_adjust")
rownames(mutation_cooccurrence_mat) <- mutation_cooccurrence_mat$gene1
mutation_cooccurrence_mat <- mutation_cooccurrence_mat[, -1]
identical(rownames(mutation_cooccurrence_mat)[-1], colnames(mutation_cooccurrence_mat)[-486])


# square matrix 
mutation_cooccurrence_mat <- rbind(cbind(A1CF = NA, mutation_cooccurrence_mat), ZNRF3 = NA)
mutation_cooccurrence_mat <- as.matrix(mutation_cooccurrence_mat)
diag(mutation_cooccurrence_mat) <- 1

all(is.na(lower.tri(mutation_cooccurrence_mat, diag = FALSE)))

upper.tri(mutation_cooccurrence_mat, diag = FALSE)
all(is.na(mutation_cooccurrence_mat[lower.tri(mutation_cooccurrence_mat, diag = FALSE)]))

a <- mutation_cooccurrence_mat
t_a <- t(a)
a[lower.tri(a, diag = FALSE)] <- t_a[lower.tri(t_a, diag = FALSE)]
identical(a[1,], a[, 1])
save(a, file = paste0(data_home, "Solid_tumor_mutation_cooccurrence_p_adjust.RData"))

