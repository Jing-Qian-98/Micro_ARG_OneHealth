# PCoA分析是一种经典的 MDS 分析方法。它可以基于任意距离尺度评价样本之间的相似度。MDS分析同样通过对样本距离矩阵进行降维分解，从而简化数据结构，展现样本在某种特定距离尺度下的自然分布

suppressMessages(library(vegan))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(plyr))
suppressMessages(library(ape))
suppressMessages(library(RColorBrewer))

# 传参及定义变量
argv <- commandArgs(TRUE)
in_taxa <- argv[1]
in_method <- argv[2]
in_file <- argv[3]
##距离矩阵算法： "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", or "robust.aitchison"

# 读取物种丰度表及预处理
taxaAbd <- read.table(in_taxa,sep = "\t",header = T, row.names = 1)
taxaAbd <- taxaAbd[ , -ncol(taxaAbd)]

# 读取分组文件
map <- read.table(in_file, row.names = 1)
Group=map$V4
Group=factor(as.character(map$V4),levels=unique(map$V4))

# 转换相对丰度并将行列转置
taxaAbd_rel <- apply(taxaAbd,2,function(x){x / sum(x)})
taxaAbd_rel <- data.frame(t(taxaAbd_rel))

# bray_curtis距离计算
dist <- vegdist(taxaAbd_rel, method=in_method)
dist <- as.matrix(dist)

# PCoA计算
pcoa <- pcoa(dist, correction="lingoes",rn=rownames(dist))
if(pcoa$correction[2]==1){name_eig="Relative_eig";axis.1=pcoa$vectors[,1];axis.2=pcoa$vectors[, 2];rn=map[rownames(pcoa$vectors),3]};
if(pcoa$correction[2]==2){name_eig="Rel_corr_eig";axis.1=pcoa$vectors.cor[,1];axis.2=pcoa$vectors.cor[, 2];rn=map[rownames(pcoa$vectors.cor),3]};

# 根据样品数量选取颜色以及点的形状
if (length(levels(factor(Group)))> 20){
col=rep(c("magenta4","tomato2","darkolivegreen1","royalblue1","firebrick1","palegoldenrod","grey30",           
       "hotpink", "goldenrod4","forestgreen","thistle1","cyan2","maroon4","aquamarine","yellow",
       "grey60","darkred","lightpink","peru","lawngreen","indianred3","mediumseagreen","mediumorchid3",
       "orangered","chocolate2","slateblue","grey80","deepskyblue2","green","gray10"),t=100)[1:length(levels(factor(Group)))] 
}else{
col=c(c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[-4][-9][-9][-11],"#B87333", "#8B00FF", "#CCB38C", "#FF2400")[1:length(levels(factor(Group)))]    #Insufficient color
}

shapes=c(15,16,17,18,25)
if (length(levels(factor(Group)))<6) {
  shape=shapes[1:length(levels(factor(Group)))]
}else{
  shape=rep(shapes[2],t=length(levels(factor(Group))))
}

# 画图数据处理
plot.data <- data.frame(Axis.1=axis.1, Axis.2=axis.2,Group=Group)
plot.var_explained <- round(100*pcoa$values[1:2, name_eig]/sum(pcoa$values[, name_eig]), digits=1);
plot.hulls <- ddply(plot.data, "Group", function(df) df[chull(df$Axis.1, df$Axis.2), ]);
plot.centroids <- ddply(plot.data, "Group", function(df) c(mean(df$Axis.1), mean(df$Axis.2)));
plot.df <- merge(plot.data, plot.centroids, by="Group");
plot.hulls$Group=factor(as.character(plot.hulls$Group),levels=unique(map$V4))
plot.df$Group=factor(as.character(plot.df$Group),levels=unique(map$V4))
plot.data$Group=factor(as.character(plot.data$Group),levels=unique(map$V4))

A=max(plot.data$Axis.1)-min(plot.data$Axis.1)
B=max(plot.data$Axis.2)-min(plot.data$Axis.2)

# draw.picture
## 带数据边界
q1 <- ggplot(data = plot.df, aes_string(x="Axis.1", y="Axis.2")) +
  geom_segment(data=plot.df, aes(x=Axis.1, y=Axis.2, xend=V1, yend=V2), color="grey",linetype="solid", size=0.4,alpha=0.9) +
  geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.2) +
  geom_point(aes(color = Group,fill = Group,shape=Group),size=2.5)+
  theme_bw()+
  xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
  ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep="")) +
  theme(axis.text.x=element_text(color = "black",size = 10),
        axis.text.y=element_text(color = "black", size = 10),
        axis.title.x=element_text(color = "black",size = 12),
        axis.title.y=element_text(color = "black", size = 12),
        axis.ticks.length = unit(.1, "cm"),
        axis.ticks=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.justification = c(0.5,0.5),
        legend.position="right",
        legend.text=element_text(face="plain",size=12),
        legend.title = element_text(face="plain",size=12),
        legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(face = "bold")
  )+ 
  scale_x_continuous()+ 
  scale_y_continuous()+ 
  scale_color_manual(values=col)+
  scale_fill_manual(values=col)+
  scale_shape_manual(values=shape)

xmin=min(plot.data$Axis.1)
xmax=max(plot.data$Axis.1)
ymin=min(plot.data$Axis.2)
ymax=max(plot.data$Axis.2)
## 调整图片比例
if (((xmax-xmin)/(ymax-ymin) < 10 && (xmax-xmin)/(ymax-ymin) > 1) | ((xmax-xmin)/(ymax-ymin)<1 && (xmax-xmin)/(ymax-ymin) > 0.1)) {
if (xmin < ymin & xmax < ymax ){
 q1 =q1+ coord_fixed(1,xlim=c(xmin-abs(xmax-ymax)/2,xmax+abs(xmax-ymax)/2),ylim=c(ymin-abs(xmin-ymin)/2,ymax+abs(xmin-ymin)/2),expand =T)
}
if (xmin > ymin & xmax < ymax ){
  q1 =q1+coord_fixed(1,xlim=c(xmin-abs(ymax-ymin-xmax+xmin)/2,xmax+abs(ymax-ymin-xmax+xmin)/2),ylim=c(ymin,ymax),expand =T)
}
if (xmin < ymin & xmax > ymax ){
  q1 =q1+coord_fixed(1,xlim=c(xmin,xmax),ylim=c(ymin-abs(ymax-ymin-xmax+xmin)/2,ymax+abs(ymax-ymin-xmax+xmin)/2),expand =T)
}
if (xmin > ymin & xmax > ymax ){
  q1 =q1+coord_fixed(1,xlim=c(xmin-abs(xmin-ymin)/2,xmax+abs(xmin-ymin)/2),ylim=c(ymin-abs(xmax-ymax)/2,ymax+abs(xmax-ymax)/2),expand =T)
}
}
## 带置信区间
q2 <- ggplot(data = plot.df, aes_string(x="Axis.1", y="Axis.2")) +
  geom_segment(data=plot.df, aes(x=Axis.1, y=Axis.2, xend=V1, yend=V2), color="grey",linetype="solid", size=0.4,alpha=0.9) +
  stat_ellipse(level=0.95, segments=101, alpha=0.2,aes(color=Group,fill=Group),size=0.6,linetype=2, geom = "polygon") +
  geom_point(aes(color = Group,fill = Group,shape=Group),size=2.5)+
  theme_bw()+
  xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
  ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep="")) +
  theme(axis.text.x=element_text(color = "black",size = 10),
        axis.text.y=element_text(color = "black", size = 10),
        axis.title.x=element_text(color = "black",size = 12),
        axis.title.y=element_text(color = "black", size = 12),
        axis.ticks.length = unit(.1, "cm"),
        axis.ticks=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.justification = c(0.5,0.5),
        legend.position="right",
        legend.text=element_text(face="plain",size=12),
        legend.title = element_text(face="plain",size=12),
        legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(face = "bold")
  )+ 
  scale_x_continuous()+ 
  scale_y_continuous()+ 
  scale_color_manual(values=col)+
  scale_fill_manual(values=col)+
  scale_shape_manual(values=shape)
## 调整图片比例
if (((xmax-xmin)/(ymax-ymin) < 10 && (xmax-xmin)/(ymax-ymin) > 1) | ((xmax-xmin)/(ymax-ymin)<1 && (xmax-xmin)/(ymax-ymin) > 0.1)) {
if (xmin < ymin & xmax < ymax ){
 q2 =q2+ coord_fixed(1,xlim=c(xmin-abs(xmax-ymax)/2,xmax+abs(xmax-ymax)/2),ylim=c(ymin-abs(xmin-ymin)/2,ymax+abs(xmin-ymin)/2),expand =T)
}
if (xmin > ymin & xmax < ymax ){
  q2 =q2+coord_fixed(1,xlim=c(xmin-abs(ymax-ymin-xmax+xmin)/2,xmax+abs(ymax-ymin-xmax+xmin)/2),ylim=c(ymin,ymax),expand =T)
}
if (xmin < ymin & xmax > ymax ){
  q2 =q2+coord_fixed(1,xlim=c(xmin,xmax),ylim=c(ymin-abs(ymax-ymin-xmax+xmin)/2,ymax+abs(ymax-ymin-xmax+xmin)/2),expand =T)
}
if (xmin > ymin & xmax > ymax ){
  q2 =q2+coord_fixed(1,xlim=c(xmin-abs(xmin-ymin)/2,xmax+abs(xmin-ymin)/2),ylim=c(ymin-abs(xmax-ymax)/2,ymax+abs(xmax-ymax)/2),expand =T)
}
}


q3 <- ggplot(data = plot.data, aes_string(x="Axis.1", y="Axis.2")) +
  geom_point(aes(color = Group,fill = Group,shape=Group),size=3)+
  theme_bw()+
  xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
  ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep="")) +
  theme(axis.text.x=element_text(color = "black",size = 10),
        axis.text.y=element_text(color = "black", size = 10),
        axis.title.x=element_text(color = "black",size = 12),
        axis.title.y=element_text(color = "black", size = 12),
        axis.ticks.length = unit(.1, "cm"),
        axis.ticks=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.justification = c(0.5,0.5),
        legend.position="right",
        legend.text=element_text(face="plain",size=12),
        legend.title = element_text(face="plain",size=12),
        legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(face = "bold")
  )+ 
  scale_x_continuous()+ 
  scale_y_continuous()+ 
  scale_color_manual(values=col)+
  scale_fill_manual(values=col)+
  scale_shape_manual(values=shape)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)
## 调整图片比例
if (((xmax-xmin)/(ymax-ymin) < 10 && (xmax-xmin)/(ymax-ymin) > 1) | ((xmax-xmin)/(ymax-ymin)<1 && (xmax-xmin)/(ymax-ymin) > 0.1)) {
if (xmin < ymin & xmax < ymax ){
 q3 =q3+ coord_fixed(1,xlim=c(xmin-abs(xmax-ymax)/2,xmax+abs(xmax-ymax)/2),ylim=c(ymin-abs(xmin-ymin)/2,ymax+abs(xmin-ymin)/2),expand =T)
}
if (xmin > ymin & xmax < ymax ){
  q3 =q3+coord_fixed(1,xlim=c(xmin-abs(ymax-ymin-xmax+xmin)/2,xmax+abs(ymax-ymin-xmax+xmin)/2),ylim=c(ymin,ymax),expand =T)
}
if (xmin < ymin & xmax > ymax ){
  q3 =q3+coord_fixed(1,xlim=c(xmin,xmax),ylim=c(ymin-abs(ymax-ymin-xmax+xmin)/2,ymax+abs(ymax-ymin-xmax+xmin)/2),expand =T)
}
if (xmin > ymin & xmax > ymax ){
  q3 =q3+coord_fixed(1,xlim=c(xmin-abs(xmin-ymin)/2,xmax+abs(xmin-ymin)/2),ylim=c(ymin-abs(xmax-ymax)/2,ymax+abs(xmax-ymax)/2),expand =T)
}
}


# hull
## 输出PDF
pdf(paste0("PCoA.", in_method,".hull.pdf"), width =6, height =5, bg="white", onefile=F)
print(q1)
dev.off()

## 输出svg
svg(paste0("PCoA.", in_method, ".hull.svg"),width =6,height =5,bg="white", onefile=F)
print(q1)
dev.off()

## 输出PNG
png(paste0("PCoA.", in_method,".hull.png"),width =6,height =5,units="in",res=100,bg="white")
print(q1)
dev.off()

# ellipse
## 输出PDF
pdf(paste0("PCoA.", in_method,".ellipse.pdf"),width =6, height =5, bg="white", onefile=F)
print(q2)
dev.off()

## 输出svg
svg(paste0("PCoA.", in_method,".ellipse.svg"), width =6, height =5, bg="white", onefile=F)
print(q2)
dev.off()

## 输出PNG
png(paste0("PCoA.", in_method,".ellipse.png"), width =6, height =5, units="in", res=100, bg="white")
print(q2)
dev.off()

## 输出PDF
pdf(paste0("PCoA.", in_method,".common.pdf"), width =6,  height =5, bg="white", onefile=F)
print(q3)
dev.off()

## 输出svg
svg(paste0("PCoA.", in_method,".common.svg"), width =6, height =5, bg="white", onefile=F)
print(q3)
dev.off()

## 输出PNG
png(paste0("PCoA.", in_method,".common.png"), width =6, height =5, units="in", res=100, bg="white")
print(q3)
dev.off()
