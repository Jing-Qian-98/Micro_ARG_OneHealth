suppressMessages(library(pheatmap))
suppressMessages(library(vegan))
suppressMessages(library(RColorBrewer))

# 传递参数以及定义变量
argv <- commandArgs(TRUE);
if (length(argv) > 0){
  if (length(argv) >3){
    stop('The script needs 2 file and 1 num:otu_modified_L6.txt,map.txt ,the num of the samples,also can have a width defaulted 7.')
  }
  in_file <- argv[1]
  in_group <- argv[2]
  in_num <- argv[3]
}

# 读取丰度文件以及预处理
df <- read.table(in_file, sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
# 求和并排序
df$sum <- apply(df[,2:ncol(df)], 1, sum)
df <- df[order(df$sum, decreasing = T),]
df<-subset(df,sum>0)
df$sum <- NULL

# 读取分组文件以及预处理
group <- read.table(in_group, sep = "\t", header = T, comment.char = "@", stringsAsFactors = F, check.names = F, quote = "")
colnames(group)[4]="Group"

# 读取组名字符长度便于画图修改参数
maxL=1
for (i in 1:length(factor(group$Group))){
  if(maxL<length(strsplit(as.character(group$Group),"")[[i]])){
    maxL=length(strsplit(as.character(group$Group),"")[[i]])
  }
}

group$Group[which(group$Group==group$Group[1])]=paste(group$Group[1],
                                                      strrep("   ",t=1+maxL-length(strsplit(as.character(group$Group[1]),"")[[1]])),
                                                      strrep("  ",t=length(strsplit(as.character(group$Group[1]),"")[[1]])),sep="")
group$Group=factor(as.character(group$Group),levels=unique(as.character(group$Group)))


# 选取Top
data <- df
if(nrow(data) < as.numeric(in_num)){
  t <- nrow(data)
}else{
  t <- as.numeric(in_num)
}
data <- data[1:t,]
rownames(data) <- data[,1]
data <- data[,-1]

# 画图数据处理
annotation_col =data.frame(group=factor(group$Group))
rownames(annotation_col)=colnames(data)
cols=c(c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[-4][-9][-9][-11],"#B87333", "#8B00FF", "#CCB38C", "#FF2400")[1:length(levels(factor(group$Group)))]
names(cols)=as.character(levels(factor(group$Group)))
ann_colors = list(
  group=cols
)

# 调整画图参数
if (maxL > 9) {
  angle_col <- 45
}else{
  angle_col <- 90 
} #字符长度大于9则倾斜字体

cellwidth <- 15
cellheight <- 12
fontsize_col <- 12
fontsize_row <- 12

## 聚类树的高度
treeheight_row <- 10+nrow(data)
treeheight_col <- 10+ncol(data)

## 画布大小
width = (max(nchar(colnames(data)))%/%2 + ncol(data) + treeheight_row%/%2 + 3)*0.2
height = (max(nchar(rownames(data)))%/%2 + nrow(data) + treeheight_col%/%2 )*0.2

# draw.pictures
## 未分组
heatmap.both_clustered <- pheatmap(data,
                                   cluster_rows=T,
                                   cluster_cols=T,
                                   cellwidth = cellwidth,
                                   cellheight = cellheight,
                                   scale="row",
                                   border_color="grey",
                                   fontsize_row =fontsize_row,
                                   fontsize_col =fontsize_col,
                                   treeheight_row = treeheight_row,
                                   treeheight_col = treeheight_col,
                                   color = colorRampPalette(rev(c("red","white","blue")))(200),
                                   clustering_method = "average",
                                   clustering_distance_rows = "correlation",
                                   #clustering_distance_cols = dcols,
                                   angle_col=angle_col)

pdf('heatmap.both_clustered.pdf',onefile=F,width=width,height=height,bg="white")
heatmap.both_clustered
dev.off()

svg('heatmap.both_clustered.svg',onefile=F,width=width,height=height,bg="white")
heatmap.both_clustered
dev.off()
dev.off()
system("gs -dSAFER -dBATCH -dNOPAUSE -r300 -sDEVICE=pngalpha -sOutputFile=heatmap.both_clustered.png heatmap.both_clustered.pdf")


heatmap.taxa_clustered <- pheatmap(data,
                                   cluster_rows=T,
                                   cluster_cols=F,
                                   cellwidth = cellwidth,
                                   cellheight = cellheight,
                                   scale="row",
                                   border_color="grey",
                                   fontsize_row =fontsize_row,
                                   fontsize_col =fontsize_col,
                                   treeheight_row = treeheight_row,
                                   treeheight_col = treeheight_col,
                                   color = colorRampPalette(rev(c("red","white","blue")))(200),
                                   clustering_method = "average",
                                   clustering_distance_rows = "correlation",
                                   angle_col=angle_col,
) 

pdf('heatmap.taxa_clustered.pdf',onefile=F,width=width,height=height,bg="white")
heatmap.taxa_clustered
dev.off()

svg('heatmap.taxa_clustered.svg',onefile=F,width=width,height=height,bg="white")
heatmap.taxa_clustered
dev.off()
dev.off()
system("gs -dSAFER -dBATCH -dNOPAUSE -r300 -sDEVICE=pngalpha -sOutputFile=heatmap.taxa_clustered.png heatmap.taxa_clustered.pdf")


## 分组
width = (max(nchar(colnames(data)))%/%2 + ncol(data) + treeheight_row%/%2 + 9)*0.2
height = (max(nchar(rownames(data)))%/%2 + nrow(data) + treeheight_col%/%2 )*0.2

if (length(levels(factor(group$Group))) != ncol(data)) 
{
  
  if (length(levels(factor(group$Group)))<21) 
  {
   
    pdf('heatmap_g.both_clustered.pdf',onefile=F,width=width,height=height,bg="white")
    pheatmap(data,
             cluster_rows=T,
             cluster_cols=T,
             cellwidth = cellwidth,
             cellheight = cellheight,
             scale="row",
             border_color="grey",
             fontsize_row =fontsize_row,
             fontsize_col =fontsize_col,
             treeheight_row = treeheight_row,
             treeheight_col = treeheight_col,
             color = colorRampPalette(rev(c("red","white","blue")))(200),
             clustering_method = "average",
             clustering_distance_rows = "correlation",
             #clustering_distance_cols = dcols,
             angle_col=angle_col,
             annotation_col = annotation_col,
             annotation_colors = ann_colors 
    )
    dev.off()
    
    svg('heatmap_g.both_clustered.svg',onefile=F,width=width,height=height,bg="white")
    pheatmap(data,
             cluster_rows=T,
             cluster_cols=T,
             cellwidth = cellwidth,
             cellheight = cellheight,
             scale="row",
             border_color="grey",
             fontsize_row =fontsize_row,
             fontsize_col =fontsize_col,
             treeheight_row = treeheight_row,
             treeheight_col = treeheight_col,
             color = colorRampPalette(rev(c("red","white","blue")))(200),
             clustering_method = "average",
             clustering_distance_rows = "correlation",
             #clustering_distance_cols = dcols,
             angle_col=angle_col,
             annotation_col = annotation_col,
             annotation_colors = ann_colors 
    )
    dev.off()
    system("gs -dSAFER -dBATCH -dNOPAUSE -r300 -sDEVICE=pngalpha -sOutputFile=heatmap_g.both_clustered.png heatmap_g.both_clustered.pdf")
    
    ##
    pdf('heatmap_g.taxa_clustered.pdf',onefile=F,width=width,height=height,bg="white")
    pheatmap(data,
             cluster_rows=T,
             cluster_cols=F,
             cellwidth = cellwidth,
             cellheight = cellheight,
             scale="row",
             border_color="grey",
             fontsize_row =fontsize_row,
             fontsize_col =fontsize_col,
             treeheight_row = treeheight_row,
             treeheight_col = treeheight_col,
             color = colorRampPalette(rev(c("red","white","blue")))(200),
             clustering_method = "average",
             clustering_distance_rows = "correlation",
             angle_col=angle_col,
             annotation_col = annotation_col,
             annotation_colors = ann_colors
    ) 
    dev.off()
    
    svg('heatmap_g.taxa_clustered.svg',onefile=F,width=width,height=height,bg="white")
    pheatmap(data,
             cluster_rows=T,
             cluster_cols=F,
             cellwidth = cellwidth,
             cellheight = cellheight,
             scale="row",
             border_color="grey",
             fontsize_row =fontsize_row,
             fontsize_col =fontsize_col,
             treeheight_row = treeheight_row,
             treeheight_col = treeheight_col,
             color = colorRampPalette(rev(c("red","white","blue")))(200),
             clustering_method = "average",
             clustering_distance_rows = "correlation",
             angle_col=angle_col,
             annotation_col = annotation_col,
             annotation_colors = ann_colors
    )     
    dev.off()
    system("gs -dSAFER -dBATCH -dNOPAUSE -r300 -sDEVICE=pngalpha -sOutputFile=heatmap_g.taxa_clustered.png heatmap_g.taxa_clustered.pdf")
  }
  else
  {
    
    ###
    pdf('heatmap_g.both_clustered.pdf',onefile=F,width=width,height=height,bg="white")
    pheatmap(data,
             cluster_rows=T,
             cluster_cols=T,
             cellwidth = cellwidth,
             cellheight = cellheight,
             scale="row",
             border_color="grey",
             fontsize_row =fontsize_row,
             fontsize_col =fontsize_col,
             treeheight_row = treeheight_row,
             treeheight_col = treeheight_col,
             color = colorRampPalette(rev(c("red","white","blue")))(200),
             clustering_method = "average",
             clustering_distance_rows = "correlation",
             #clustering_distance_cols = dcols,
             angle_col=angle_col
    )
    dev.off()
    
    svg('heatmap_g.both_clustered.svg',onefile=F,width=width,height=height,bg="white")
    pheatmap(data,
             cluster_rows=T,
             cluster_cols=T,
             cellwidth = cellwidth,
             cellheight = cellheight,
             scale="row",
             border_color="grey",
             fontsize_row =fontsize_row,
             fontsize_col =fontsize_col,
             treeheight_row = treeheight_row,
             treeheight_col = treeheight_col,
             color = colorRampPalette(rev(c("red","white","blue")))(200),
             clustering_method = "average",
             clustering_distance_rows = "correlation",
             #clustering_distance_cols = dcols,
             angle_col=angle_col
    )
    dev.off()
    system("gs -dSAFER -dBATCH -dNOPAUSE -r300 -sDEVICE=pngalpha -sOutputFile=heatmap_g.both_clustered.png heatmap_g.both_clustered.pdf")
    
    ####
    pdf('heatmap_g.taxa_clustered.pdf',onefile=F,width=width,height=height,bg="white")
    pheatmap(data,
             cluster_rows=T,
             cluster_cols=F,
             cellwidth = cellwidth,
             cellheight = cellheight,
             scale="row",
             border_color="grey",
             fontsize_row =fontsize_row,
             fontsize_col =fontsize_col,
             treeheight_row = treeheight_row,
             treeheight_col = treeheight_col,
             color = colorRampPalette(rev(c("red","white","blue")))(200),
             clustering_method = "average",
             clustering_distance_rows = "correlation",
             angle_col=angle_col,
    ) 
    dev.off()
    
    svg('heatmap_g.taxa_clustered.svg',onefile=F,width=width,height=height,bg="white")
    pheatmap(data,
             cluster_rows=T,
             cluster_cols=F,
             cellwidth = cellwidth,
             cellheight = cellheight,
             scale="row",
             border_color="grey",
             fontsize_row =fontsize_row,
             fontsize_col =fontsize_col,
             treeheight_row = treeheight_row,
             treeheight_col = treeheight_col,
             color = colorRampPalette(rev(c("red","white","blue")))(200),
             clustering_method = "average",
             clustering_distance_rows = "correlation",
             angle_col=angle_col,
    )
    dev.off()
    
  }
  system("gs -dSAFER -dBATCH -dNOPAUSE -r300 -sDEVICE=pngalpha -sOutputFile=heatmap_g.taxa_clustered.png heatmap_g.taxa_clustered.pdf")
  
}



