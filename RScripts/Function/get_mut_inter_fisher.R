#' 使用列联表检验来识别突变之间的共发生或互斥关系
#'
#' @param mutationMat 突变矩阵
#' @param outfile 热图的输出路径
#' @param min_count 设置基因的最小突变数目
#'
#' @return 检验结果，包括 odds ratio 和 p value

get_mut_inter_fisher <- function(mutationMat=NULL, outfile=NULL, min_count=2){
    # 提取突变数目大于min_count的基因 
    ex_idx <- which(rowSums(mutationMat)>min_count)
    if(length(ex_idx)>0){
        mutationMat_sub <- mutationMat[ex_idx, ]
        
        test_genes <- rownames(mutationMat_sub)
        test_gene_combn <- combn(test_genes, 2)
        res_list <- list()
        for(comb_set_idx in 1:dim(test_gene_combn)[2]){
            # 构建四格表
            test_table <- matrix(NA, nrow=2, ncol=2)
            rownames(test_table) <- paste0(test_gene_combn[1, comb_set_idx], c("_MUT", "_WT"))
            colnames(test_table) <- paste0(test_gene_combn[2, comb_set_idx], c("_MUT", "_WT"))
            
            test_table_tmp <- table(rowSums(t(mutationMat_sub[test_gene_combn[, comb_set_idx], ])))
            test_table[1, 1] <- test_table_tmp["2"]
            test_table[2, 2] <- test_table_tmp["0"]
            
            test_table_tmp2 <- table(mutationMat_sub[test_gene_combn[1, comb_set_idx], ]-mutationMat_sub[test_gene_combn[2, comb_set_idx], ])
            test_table[1, 2] <- test_table_tmp2["1"]
            test_table[2, 1] <- test_table_tmp2["-1"]
            # 将 NA转化为 0
            na_idx <- which(is.na(test_table))
            if(length(na_idx)>0)  test_table[na_idx] <- 0
			
            # 卡方检验，获得odd ratio 和 p value
            res_tmp <- getChisqTest(testTable=test_table)
            res_list[[comb_set_idx]] <- as.data.frame(res_tmp[c("oddsratio", "pvalue")])
        }
        res <- Reduce("rbind", res_list)
        rownames(res) <- apply(test_gene_combn, 2, paste, collapse="_")
		
        ## 是否输出热图
        if(!is.null(outfile)){
          library(pheatmap)
            # 提取基因名
            gene_names <- unique(unlist(strsplit(rownames(res), split="_")))
            # 建立热图矩阵以及匹配的P值矩阵
            heatmap_arr <- matrix(NA, nrow=length(gene_names), ncol=length(gene_names))
            pvalue_arr <- matrix(NA, nrow=length(gene_names), ncol=length(gene_names))
            rownames(heatmap_arr) <- colnames(heatmap_arr) <- gene_names
            rownames(pvalue_arr) <- colnames(pvalue_arr) <- gene_names
            
            heatmap_arr_idx <- Reduce("rbind", strsplit(rownames(res), split="_"))
            for(i in 1:dim(heatmap_arr_idx)[1]){
                heatmap_arr[heatmap_arr_idx[i, 1], heatmap_arr_idx[i, 2]] <- res[i, "oddsratio"]
                pvalue_arr[heatmap_arr_idx[i, 1], heatmap_arr_idx[i, 2]] <- res[i, "pvalue"]
                heatmap_arr[heatmap_arr_idx[i, 1], heatmap_arr_idx[i, 1]] <- 1
                pvalue_arr[heatmap_arr_idx[i, 1], heatmap_arr_idx[i, 1]] <- 1
            }
            #  log2 转化
            log2heatmap_arr <- log2(heatmap_arr)
            # 替换Inf值
            if(length(which(log2heatmap_arr==Inf))>0)
                log2heatmap_arr[which(log2heatmap_arr==Inf)] <- 1000
            if(length(which(log2heatmap_arr==(-Inf)))>0)
                log2heatmap_arr[which(log2heatmap_arr==(-Inf))] <- -1000
            # 设置颜色
            heat_colors <- colorRampPalette(c("blue", "white", "white", "red"))(30)
            color_breaks <- seq(-4, 4, length=31)
            # 输出
            pdf(outfile)
            print(pheatmap::pheatmap(mat=t(log2heatmap_arr), color=heat_colors, breaks=color_breaks, cluster_rows=FALSE, cluster_cols=FALSE))
            dev.off()
        }
        return(res)
    }else{
        return(NA)
    }
}


