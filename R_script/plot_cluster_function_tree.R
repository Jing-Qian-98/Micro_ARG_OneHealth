.libPaths(c("/PERSONALBIO/work/microbio/m13/miniconda3/envs/R4.3/lib/R/library/", "/PERSONALBIO/work/microbio/m13/miniconda3/envs/R4.2/lib/R/library/"))
suppressMessages(library(vegan))
suppressMessages(library(phangorn))
suppressMessages(library(ape))

argv <- commandArgs(TRUE)
in_func <- argv[1]
in_method <- argv[2]
category <- argv[3]
##距离矩阵算法： "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", or "robust.aitchison"

# 读取功能丰度表及预处理
funcAbd <- read.table(in_func,sep = "\t",header = T, row.names = 1)
# 如果最后一列是字符串，则需要去除
if(!is.numeric(funcAbd[, ncol(funcAbd)])){
funcAbd <- funcAbd[-ncol(funcAbd)]
}


# 转换相对丰度并将行列转置
funcAbd_rel <- apply(funcAbd,2,function(x){x / sum(x)})
funcAbd_rel <- data.frame(t(funcAbd_rel))

# 计算距离矩阵并进行upgma聚类
dfupgma <- upgma(vegdist(funcAbd_rel, method=in_method))
#dfupgma <- hclust(vegdist(funcAbd_rel, method=in_method), method = "average")

# draw.pictures
if(dfupgma$Nnode >= 30){
  pdf(paste0("upgma_", category, ".pdf"), height= 15)
}else{
  pdf(paste0("upgma_", category, ".pdf"))
}
plot.phylo(dfupgma)
edgelabels(round(dfupgma$edge.length,3), cex = .6, adj = c(0.5, -0.25), bg = "white", frame = "none")
dev.off()


