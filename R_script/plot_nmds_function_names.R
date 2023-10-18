library(vegan)
library(ggplot2)
library(grid)

# 传参及定义变量
argv <- commandArgs(TRUE)
funcAbdDir <- argv[1]
map_file <- argv[2]
dist_method <- argv[3]
category <- argv[4]

# 读取功能丰度表及预处理
funcAbd <- read.table(funcAbdDir,sep = "\t",header = T, row.names = 1)
# 如果最后一列是字符串，则需要去除
if(!is.numeric(funcAbd[, ncol(funcAbd)])){
funcAbd <- funcAbd[-ncol(funcAbd)]
}

# 读取分组表及预处理
map <- read.table(map_file, sep = "\t", header = F, stringsAsFactors = F, quote = "", colClasses = "character")
sum_group <- length(unique(map$V4))

# 转换相对丰度并将行列转置
funcAbd_rel <- apply(funcAbd,2,function(x){x / sum(x)})
funcAbd_rel <- data.frame(t(funcAbd_rel))

# 计算NMDS值及画图数据处理
dfNmds <- metaMDS(funcAbd_rel, distance=dist_method, k = 2)
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
pdf(paste0("nmds.",category,".pdf"), width =10, height =8, bg="white", onefile=F)
grid.draw(q2)
dev.off()

## 输出svg
svg(paste0("nmds.",category,".svg"), width =10, height =8, bg="white", onefile=F)
grid.draw(q2)
dev.off()

## 输出PNG
png(paste0("nmds.",category,".png"), width =10, height =8, units="in", res=100, bg="white")
grid.draw(q2)
dev.off()

