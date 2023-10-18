# 以等级树的形式展示样本间的相似度，通过聚类树的分枝长度衡量聚类效果的好坏。与 MDS 分析相同，聚类分析可以采用任何距离评价样本之间的相似度

.libPaths(c("/PERSONALBIO/work/microbio/m13/miniconda3/envs/R4.3/lib/R/library/", "/PERSONALBIO/work/microbio/m13/miniconda3/envs/R4.2/lib/R/library/"))
suppressMessages(library(vegan))
suppressMessages(library(phangorn))
suppressMessages(library(ape))
#suppressMessages(library(stats))

argv <- commandArgs(TRUE)
in_taxa <- argv[1]
in_method <- argv[2]
in_levels <- argv[3]

##距离矩阵算法： "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", or "robust.aitchison"

# 读取物种丰度表及预处理
taxaAbd <- read.table(in_taxa,sep = "\t",header = T, row.names = 1)
taxaAbd <- taxaAbd[ , -ncol(taxaAbd)]

# 转换相对丰度并将行列转置
taxaAbd_rel <- apply(taxaAbd,2,function(x){x / sum(x)})
taxaAbd_rel <- data.frame(t(taxaAbd_rel))

# 计算距离矩阵并进行upgma聚类
dfupgma <- upgma(vegdist(taxaAbd_rel, method=in_method))
#dfupgma <- hclust(vegdist(taxaAbd_rel, method=in_method), method = "average")

# draw.pictures
if(dfupgma$Nnode >= 30){
  pdf(paste0("upgma_", in_levels, ".pdf"), height= 15)
}else{
  pdf(paste0("upgma_", in_levels, ".pdf"))
}
plot.phylo(dfupgma)
edgelabels(round(dfupgma$edge.length,3), cex = .6, adj = c(0.5, -0.25), bg = "white", frame = "none")

dev.off()


