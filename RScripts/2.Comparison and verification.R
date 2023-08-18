###############################
#' TODO: Comparison of function potential of the
#' identified immune regulatory drivers and non-regulators
#' 
#'###############################################################################
setwd("/~")

###
# 1. Data preparation --------------------
#' 
#'###############################################################################
### identified CYT regulators
load(file = "./Results/non_zero_variables.RData")
mutation_status_mat <- read.csv("./data/Solid_tumor_mutation_status_mat.csv", row.names = 1, check.names = FALSE) # 487 8223
load(file = "./Results/embeddings_glmnet_fit_beta_list.RData")
load(file = "./data/mutation_status_mat_list.RData")

### 17 immune pathways from ImmPort
ImmPort_GeneList <- read.table("./data/ImmunePathways/ImmPort_GeneList.txt", sep = "\t", header = TRUE) # 2013 6
# ImmuneSigDB_GeneList <- clusterProfiler::read.gmt("./data/ImmunePathways/c7.all.v2023.1.Hs.symbols.gmt") # 2013 6
ImmuneSigDB_GeneList <- clusterProfiler::read.gmt("./data/ImmunePathways/c7.immunesigdb.v2023.1.Hs.symbols.gmt") # 2013 6

CYS_signature <- unique(c("APBA2", "APOL3", "CTSW", "DUSP2", "GNLY", "GZMA", "GZMH", "KLRB1", "KLRD1", "KLRF1", "KLRK1", "NKG7", "RORA", "RUNX3", "SIGIRR", "WHAMMP3", "ZBTB16", 
                ImmPort_GeneList$Symbol[which(ImmPort_GeneList$Category == "NaturalKiller_Cell_Cytotoxicity")]))
save(CYS_signature, file = "./data/ImmunePathways/CYS_signature.RData")


### PPI
load(file = "./data/PPI_STRING/gene_inter_data.RData")
gene_inter_data1 <- unique(gene_inter_data[, c("gene1", "gene2")]) # 



###
# 2. Comparison of number of involving PPI --------------------
#' 
#'###############################################################################
#'######
# (1) pan-cancer --------------------
#'######
write.table(non_zero_variables, file = "./Results/non_zero_variables.txt",
            sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(setdiff(rownames(mutation_status_mat), non_zero_variables), file = "./Results/non_reg_variables.txt",
            sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(unique(c(non_zero_variables, CYS_signature)),  
            file = "./Results/reg_CYS.txt",
            sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unique(c(setdiff(rownames(mutation_status_mat), non_zero_variables), CYS_signature)), 
            file = "./Results/non_reg_CYS.txt",
            sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)



###
# 2. Comparison of PPI attributes --------------------
#' Results of the PPI analysis from STRING
#' https://version-12-0.string-db.org/
#'###############################################################################

### PPI attributes
PPI_attributes <- data.frame(driver_class = c("CYT_regulators", "non_regulators"), edges_num = c(6054, 4395), 
            avg_degree = c(37.4, 28), avg_cc = c(0.555, 0.509))

p_list <- lapply(setdiff(colnames(PPI_attributes), "driver_class"), function(at){
    p <- ggplot(PPI_attributes, aes(x = driver_class, y = PPI_attributes[, at], fill = driver_class, color = driver_class)) +
    geom_bar(stat = "identity") + 
    labs(x = "", y = at) + 
    geom_text(label = PPI_attributes[, at], position = position_stack(1.1)) +
    theme_classic2 () 

    return(p)
})

library(patchwork)
pdf("./Results/PPI/PPI_attributes.pdf", width = 8)
p_list[[1]] + p_list[[2]] + p_list[[3]] + plot_layout(guides = 'collect')
dev.off()



### nodes degree
node_degrees_reg_CYS <- read.csv("./Results/PPI/string_node_degrees_reg_CYS.tsv", sep = "\t")
node_degrees_reg_CYS$node_class <- "CYT_reg"

node_degrees_nonreg_CYS <- read.csv("./Results/PPI/string_node_degrees_nonreg_CYS.tsv", sep = "\t")
node_degrees_nonreg_CYS$node_class <- "non_reg"

node_degrees_CYS <- rbind.data.frame(node_degrees_reg_CYS, node_degrees_nonreg_CYS)
p <- ggviolin(node_degrees_CYS, x = "node_class", y = "node_degree", color = "node_class",
               fill = "node_class", 
               add = "boxplot", add.params = list(fill = "white", color = "black")) + 
        stat_compare_means() + # Add global p-value
          labs(y = "Degree of nodes in PPI(STRING)", x = "") + 
          theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5), legend.position = "none")
pdf("./Results/PPI/node_degrees_CYS_boxplot1.pdf", width = 5)
p
dev.off()