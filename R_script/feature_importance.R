library(RColorBrewer)
library(ggplot2)

# # 传参及定义变量
argv <- commandArgs(TRUE)
importance_file <- argv[1]
num_taxa <- argv[2]

#读取贡献度文件
imp = read.table(importance_file, header=T, sep="\t")

# 选取贡献度top N
if(nrow(imp) > num_taxa){
  imp_top <- imp[1:num_taxa, , drop = F]}

# draw.pictures
col_source <- rep(c(c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[-4][-9][-9][-11],"#B87333", "#8B00FF", "#CCB38C", "#FF2400"),t=100)
col <- c(rev(col_source[1:num_taxa]),"grey30")

p=ggplot(data = imp_top, mapping = aes(x=feature,y=importance,fill=feature)) + 
       xlab("Feature Taxa") +
       ylab("importance") +
       scale_color_manual(values = col) +
       scale_fill_manual(values = col,guide = guide_legend(nrow =min(35,nrow(imp_top)))) + 
       geom_bar(stat="identity")+coord_flip()+theme_bw()


# 输出pdf
pdf(paste0("feature_impotance.pdf"), width = 10, height =10)
print(p)
dev.off()
