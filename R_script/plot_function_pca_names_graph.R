suppressMessages(library(permute))
suppressMessages(library(vegan))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))

# 传参及定义变量
argv <- commandArgs(TRUE)
if(length(argv) > 0){
  if (length(argv) != 3) stop('The script needs two files')
  funcAbdDir <- argv[1]
  map_file <- argv[2]
  category <- argv[3]
}

# 读取物种丰度表及预处理
funcAbd <- read.table(funcAbdDir, sep = "\t", header = T, row.names = 1, check.names = F, quote = "")
# 如果最后一列是字符串，则需要去除
if(!is.numeric(funcAbd[, ncol(funcAbd)])){
funcAbd <- funcAbd[-ncol(funcAbd)]
}


# 转换相对丰度
funcAbd_rel <- apply(funcAbd,2,function(x){x / sum(x)})
funcAbd_rel <- as.data.frame(funcAbd_rel)

#重命名并行列转置
rownames(funcAbd_rel) <- paste("X", as.character(1:nrow(funcAbd_rel)), sep = "")
df <- t(funcAbd_rel)
df <- df[order(rownames(df)),]


# 读取map表及预处理
map <- read.table(map_file, sep = "\t", header = T, comment.char = "@", stringsAsFactors = F, colClasses = "character", quote = "")
## 保留样品名和分组名并按样品名排序
map <- map[,c(-2,-3,-(ncol(map))), order(map[,1])]
#map <- map[,c(-2,-3,-(ncol(map)))]
#map <- map[order(map[,1]),]
group<-map[,1]
map <- map[, -1, drop = F]


# draw.picture
## 可同时生成多个分组表对应的结果
for (i in 1:ncol(map)){
  dft <- df[map[,i] != "Undefined",]
  result <- rda(dft)
  pca <- summary(result)$sites
  Proportion_Explained <- summary(result)$cont$importance[2, 1:(ncol(pca))]
  xl <- paste("PC1(" ,as.character(round(Proportion_Explained[1] * 100, 2)), "%)", 
              sep = "")
  yl <- paste("PC2(" ,as.character(round(Proportion_Explained[2] * 100, 2)), "%)", 
              sep = "")
  pca <- as.data.frame(pca)
  pca$Group <- map[map[,i] != "Undefined", i]
  pca$names<-group
  q <- ggplot(pca, aes(PC1, PC2,colour = factor(Group,levels=unique(map[map[,i] != "Undefined", i])), 
                       shape = factor(Group,levels=unique(map[map[,i] != "Undefined", i])))) +
    geom_point(size = 5) +theme(legend.text=element_text(size=12))+
    geom_hline(data = pca, yintercept = 0, linetype = "dashed", colour = "black", size = 0.25) +
    geom_vline(data = pca, xintercept = 0, linetype = "dashed", colour = "black", size = 0.25) +
    theme_bw() + geom_text(aes(label=names),show.legend =F,hjust=0.5,vjust= -1.5,alpha=0.8,cex=3)+
    theme(legend.key = element_rect(linetype='blank',fill = 'white'))+
    labs(x = xl, y = yl)+ggtitle("PCA") + theme(plot.title = element_text(hjust = 0.5))
  if (length(unique(pca$Group)) > 6){
    q <- q + scale_shape_manual(values = c(rep(0:25,(length(unique(pca$Group))) %/% 26 + 1)[1:(length(unique(pca$Group)))])) +
      scale_colour_hue(l = 50, c = 150) +scale_fill_distiller(values=pca$Group)
  }
  q <- q + labs(color="Group",shape="Group")
  q2 <- ggplot_gtable(ggplot_build(q))
  q2$layout$clip[q2$layout$name == "panel"] <- "off"
  
  ## 输出PDF
  pdf(paste0("PCA_", category, "_", colnames(map)[i], ".pdf"),
      width =9,
      height =8,
      bg="white", onefile=F)
  grid.draw(q2) 
  dev.off()
  
  ## 输出svg
  svg(paste0("PCA_", category, "_", colnames(map)[i], ".svg"),
      width =9,
      height =8,
      bg="white", onefile=F)
  grid.draw(q2) 
  dev.off()
  
  ## 输出PNG
  png(paste0("PCA_", category, "_", colnames(map)[i], ".png"),
      width =9, 
      height =8,
      units="in",
      res=100,
      bg="white")
  grid.draw(q2) 
  dev.off()
  
  ## 输出PCA结果表
  pca <- rbind(pca, Proportion_Explained)
  rownames(pca)[nrow(pca)] <- "Proportion"
  pca[nrow(pca), ncol(pca)] <- NA
  sample <- rownames(pca)
  pca <- cbind(sample, pca)
  pca$Group <- NULL
  write.table(pca, paste("PCA_", category, "_", colnames(map)[i], ".txt", sep = ""), 
              sep = "\t", quote = F, row.names = F)
}
