# ​NMDS 分析 (Nonmetric Multidimensional Scaling) 与上述 PCoA 分析类似，也是一种基于样本距离矩阵的 MDS 分析方法，通过降维处理简化数据结构，在新的低维坐标系中对样本重新排序，从而在特定距离尺度下描述样本的分布特征。与 PCoA 分析不同，NMDS 分析不依赖于特征根和特征向量的计算，而是通过对样本距离进行等级排序，使样本在低维空间中的排序尽可能符合彼此之间的距离远近关系 (而非确切的距离数值)。因此，NMDS 分析不受样本距离的数值影响，仅考虑彼此之间的大小关系，对于结构复杂的数据，排序结果可能更稳定。

library(vegan)
library(ggplot2)
library(grid)

# 传参及定义变量
argv <- commandArgs(TRUE)
AbdDir <- argv[1]
map_file <- argv[2]
dist_method <- argv[3]
levels <- argv[4]

# 读取功能丰度表及预处理
taxaAbd <- read.table(AbdDir,sep = "\t",header = T, row.names = 1)
# 如果最后一列是字符串，则需要去除
if(!is.numeric(taxaAbd[, ncol(taxaAbd)])){
taxaAbd <- taxaAbd[-ncol(taxaAbd)]
}

# 读取分组表及预处理
map <- read.table(map_file, sep = "\t", header = F, stringsAsFactors = F, quote = "", colClasses = "character")
sum_group <- length(unique(map$V4))

# 转换相对丰度并将行列转置
taxaAbd_rel <- apply(taxaAbd,2,function(x){x / sum(x)})
taxaAbd_rel <- data.frame(t(taxaAbd_rel))

# 计算NMDS值及画图数据处理
dfNmds <- metaMDS(taxaAbd_rel, distance=dist_method, k = 2)
df <- data.frame(dfNmds$points)
df$Group<-map[match(rownames(df),map$V1),4]
write.table(df, "nmds.txt", sep = "\t")


# draw.pictures
if(sum_group <20){
  para <- 1
}else{
  para <-2
}

q <- ggplot(df, aes(MDS1, MDS2, colour = Group,shape = Group,label=rownames(df)))+geom_point(size=5)+theme_bw()+
           guides(col = guide_legend(ncol = para))+theme(legend.key = element_rect(linetype='blank',fill = 'white'))+
           geom_text(show.legend=F,hjust=0.5,cex=3,vjust=-1.5,alpha=0.8)+labs(title="NMDS") + theme(plot.title = element_text(hjust = 0.5)) +
                     geom_hline(data = df, yintercept = 0, linetype = "dashed", colour = "black", size = 0.25) +
                     geom_vline(data = df, xintercept = 0, linetype = "dashed", colour = "black", size = 0.25)
if(sum_group <11){q <- q+scale_colour_manual(values = c("red","blue","darkgreen","hotpink","blueviolet","brown","cadetblue","chocolate1","slateblue","green","black","cyan"))
}
if (sum_group > 5){
  q <- q + scale_shape_manual(values =c(rep(0:25,sum_group %/% 26 + 1)[1:sum_group]))
}
q <- q + labs(color="Group",shape="Group")
q2 <- ggplot_gtable(ggplot_build(q))
q2$layout$clip[q2$layout$name == "panel"] <- "off"


## 输出PDF
pdf(paste0("nmds.",levels,".pdf"), width =10, height =8, bg="white", onefile=F)
grid.draw(q2)
dev.off()

## 输出svg
svg(paste0("nmds.",levels,".svg"), width =10, height =8, bg="white", onefile=F)
grid.draw(q2)
dev.off()

## 输出PNG
png(paste0("nmds.",levels,".png"), width =10, height =8, units="in", res=100, bg="white")
grid.draw(q2)
dev.off()

