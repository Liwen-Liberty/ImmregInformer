#' 列联表检验
#'
#' @param testTable 四格表
#'
#' @return 检验结果，包括 odds ratio,  p value 和 检验方法

getChisqTest <- function(testTable=NULL){
	testTable <- as.matrix(testTable)
	# 计算理论频次
	thFreq <- function(o1, o2, total){
		(max(o1, o2)/total)*min(o1, o2)
	}
	n1 <- thFreq(sum(testTable[1,]), sum(testTable[,1]), sum(testTable))
	n2 <- thFreq(sum(testTable[1,]), sum(testTable[,2]), sum(testTable))
	n3 <- thFreq(sum(testTable[2,]), sum(testTable[,1]), sum(testTable))
	n4 <- thFreq(sum(testTable[2,]), sum(testTable[,2]), sum(testTable))
	theorFreq <- c(n1, n2, n3, n4)
	testRes <- NULL
	# 根据数据的大小，对不同类型的列联表进行检验
	if(sum(testTable)>40){
		if(!any(theorFreq<5)){
			testRes <- chisq.test(testTable, correct=FALSE)
		}
		if(sum(theorFreq<5)==1){
			testRes <- chisq.test(testTable, correct=TRUE)
		}
		if(sum(theorFreq<5)>1 || any(theorFreq<1)){
			testRes <- fisher.test(testTable)
		}
	}else{
		testRes <- fisher.test(testTable)
	}
	#library(fmsb)
	# 计算 odds ratio
	oddsratioRes <- fmsb::oddsratio(testTable)
	
	resList <- list(oddsratioRes$estimate, testRes$p.value, oddsratioRes$conf.int, oddsratioRes$method, testRes$method)
	names(resList) <- c("oddsratio", "pvalue", "oddsratio_95CI", "oddsratio_method", "test_method")
	return(resList)
}



