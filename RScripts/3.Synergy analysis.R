###############################
#' TODO: Mutation co-occurence analysis
#'       
#'###############################################################################
setwd("/~")

###
# 1. Data preparation --------------------
#' 
#'###############################################################################
selfattention_cor <- read.csv("./Results/selfattention.csv", header = FALSE) # 487  487
mutation_status_mat <- read.csv("./data/Solid_tumor_mutation_status_mat.csv", row.names = 1, check.names = FALSE) # 487 8223
rownames(selfattention_cor) <- colnames(selfattention_cor) <- rownames(mutation_status_mat)
load(file = "./Results/non_zero_variables.RData")

mutation_cooccurrence_mat <- read.csv("./data/Solid_tumor_mutation_cooccurrence_mat.csv", row.names = 1, check.names = FALSE) # 487 8223
mutation_cooccurrence_pad_mat <- get(load(file = "./data/Solid_tumor_mutation_cooccurrence_p_adjust.RData"))


###
# 3. Density of attention --------------------
#' 
#'###############################################################################
attention_df <- do.call(rbind.data.frame, lapply(1:nrow(selfattention_cor), function(i){
    df <- data.frame(source = rownames(selfattention_cor)[i], 
                     target = colnames(selfattention_cor),
                     attention = as.numeric(selfattention_cor[i, ]))
    return(df)
}))
summary(attention_df$attention)

### density
attention_df$CYT_reg <- sapply(attention_df$source, function(x) return(is.element(x, non_zero_variables)))
attention_df$target_reg <- sapply(attention_df$target, function(x) return(is.element(x, non_zero_variables)))


p <- ggplot(attention_df, aes(x = scale_attention, color = CYT_reg)) +
    geom_density() +
    theme_classic2 ()
pdf("./Results/co_occurence/scale_attention_density_regdiff.pdf", height = 3)
p
dev.off()


attention_df$scale_attention <- (attention_df$attention-min(attention_df$attention))/max(attention_df$attention)
summary(attention_df$scale_attention)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.0892  0.1045  0.1062  0.1172  0.9912 

attention_df$attention_rank <- rank(attention_df$attention)
save(attention_df, file = "./Results/attention_df.RData")



###
# 4. Correlation of attention and Odd Ratio --------------------
#' dotplot
#'###############################################################################
attention_df$class <- "reg_intra"
attention_df$class[which(attention_df$CYT_reg == FALSE & attention_df$target_reg == FALSE)] <- "nonreg_intra"
attention_df$class[which(attention_df$CYT_reg == TRUE & attention_df$target_reg == FALSE)] <- "reg_nonreg"
attention_df$class[which(attention_df$CYT_reg == FALSE & attention_df$target_reg == TRUE)] <- "nonreg_reg"

sub_attention_df <- subset(attention_df, class == "reg_intra" & scale_attention > quantile(attention_df$scale_attention, 0.95))
sub_attention_df$co_occurence_OR <- sapply(1:nrow(sub_attention_df), function(i){
  return(mutation_cooccurrence_mat[sub_attention_df$source[i], sub_attention_df$target[i]])
})
sub_attention_df$co_occurence_pad <- sapply(1:nrow(sub_attention_df), function(i){
  return(mutation_cooccurrence_pad_mat[sub_attention_df$source[i], sub_attention_df$target[i]])
})

length(which(sub_attention_df$co_occurence_pad < 0.05))/nrow(sub_attention_df) # 96.3% 
sub_attention_df$co_occurence_signif <- ifelse(sub_attention_df$co_occurence_pad < 0.05, "adjusted P < 0.05", "adjusted P >= 0.05")


idx <- c(which(sub_attention_df$co_occurence_OR > 0.75),
        which(sub_attention_df$'co_occurence_pad' >= 0.05 & sub_attention_df$scale_attention > 0.3))
sub_attention_df$label = ""
sub_attention_df$label[idx] <- paste(sub_attention_df$source[idx], sub_attention_df$target[idx], sep = "_")

library(ggrepel)
p6 <- ggplot(sub_attention_df, aes(x = scale_attention, y = co_occurence_OR, fill = co_occurence_signif, color = co_occurence_signif)) +
    geom_point(shape = 21) +
    geom_abline(intercept = 0, slope = 1) +
    # geom_smooth(method = "lm") + 
    geom_text_repel(aes(label = label), size = 3) + 
    theme_classic2 () 

pdf("./Results/co_occurence/attention_OR_dotplot.pdf", width = 9)
p6
dev.off()


### top0.01 heatmap
sub_attention_df_top <- subset(sub_attention_df, scale_attention > quantile(attention_df$scale_attention, 0.99))
sub_selfattention_cor <- selfattention_cor[unique(sub_attention_df_top$source), unique(sub_attention_df_top$target)]
pdf("./Results/co_occurence/attention_top0.01_cor.pdf", height = 8)
pheatmap(sub_selfattention_cor, fontsize = 8, border = "NA")
dev.off()
