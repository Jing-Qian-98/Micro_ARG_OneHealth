# RDA分析：基于线性模型，检测环境因子、样本、菌群三者之间的关系或者两两之间的关系是否显著，与PCA分析的区别在于是否有除样品和菌群之外的第三个变量。

library(vegan)
library(ggplot2)
library(ggrepel)

# 传参及定义变量
argv <- commandArgs(TRUE)
in_taxa <- argv[1]
env_data <- argv[2]
map_file <- argv[3]
in_levels <- argv[4]

# 读取物种丰度表及预处理
taxaAbd <- read.table(in_taxa,sep = "\t",header = T, row.names = 1)
taxaAbd <- taxaAbd[ , -ncol(taxaAbd)]

# 转换相对丰度并将行列转置
taxaAbd_rel <- apply(taxaAbd,2,function(x){x / sum(x)})
phylum <- t(taxaAbd_rel)

# 读取分组文件
group <- read.table(map_file, row.names = 1, header = F, sep = "\t")
colnames(group)[3] <- "group"

##读取环境因子及表达量转化
env <- read.delim(env_data, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
env=log10(env)

###对响应变量做hellinger转化
phylum1=decostand(phylum,method = "hellinger")

# 计算RDA值
rda_tb<-rda(phylum1~.,env)
result<-summary(rda_tb)

#coef() 提取RDA结果典范系数
rda_coef <- coef(rda_tb)

# R2校正
r2 <- RsquareAdj(rda_tb)
rda_noadj <- r2$r.squared #原始 R2
rda_adj <- r2$adj.r.squared   #校正后的 R2

#所有约束轴的置换检验，基于 999 次置换检验
rda_tb_test <- anova(rda_tb, permutations = 999)

#各约束轴逐一检验，基于 999 次置换检验
rda_tb_test_axis <- anova(rda_tb, by = 'axis', permutations = 999)

# p值校正（Bonferroni 为例）
rda_tb_test_axis$`Pr(>F)` <- p.adjust(rda_tb_test_axis$`Pr(>F)`, method = 'bonferroni')

# 提取残差特征值
pca_eig <- rda_tb$CA$eig

# Kaiser-Guttman准则
pca_eig[pca_eig > mean(pca_eig)]

# 断棍模型
n <- length(pca_eig)
bsm <- data.frame(j=seq(1:n), p = 0)
bsm$p[1] <- 1/n
for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
bsm$p <- 100*bsm$p/n

# 计算方差膨胀因子
vif.cca(rda_tb)

#画图数据处理：ordiR2step() 前向选择，基于 999 次置换检验
rda_tb_forward_r <- ordiR2step(rda(phylum1~1, env, scale = FALSE), scope = formula(rda_tb), R2scope = rda_adj, direction = 'forward', permutations = 999)
rda_tb_forward_r.scaling1 <- summary(rda_tb_forward_r, scaling = 1)
rda_tb_forward_r.site <- data.frame(rda_tb_forward_r.scaling1$sites)[1:2]
rda_tb_forward_r.env <- data.frame(rda_tb_forward_r.scaling1$biplot)[1:2]
rda_tb_forward_r.site$sample <- rownames(rda_tb_forward_r.site)
rda_tb_forward_r.site <- merge(rda_tb_forward_r.site, group, by = 'row.names')
rda_tb_forward_r.env$sample <- rownames(rda_tb_forward_r.env)
colnames(rda_tb_forward_r.site)[7] <- "group"

# draw.pictures
color=c( "#3C5488B2","#00A087B2", 
         "#F39B7FB2","#91D1C2B2", 
         "#8491B4B2", "#DC0000B2", 
         "#7E6148B2","yellow", 
         "darkolivegreen1", "lightskyblue", 
         "darkgreen", "deeppink", "khaki2", 
         "firebrick", "brown1", "darkorange1", 
         "cyan1", "royalblue4", "darksalmon", 
         "darkgoldenrod1", "darkseagreen", "darkorchid")

p <- ggplot(rda_tb_forward_r.site, aes(RDA1, RDA2)) +
  geom_point(aes(color = group,shape = group)) +
  stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE, linetype = 2) +
  scale_color_manual(values = color[1:length(unique(group$group))]) +
  theme(panel.grid.major = element_line(color = 'gray', size = 0.1), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'RDA1', y = 'RDA2') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = rda_tb_forward_r.env, aes(x = 0, y = 0, xend = RDA1,yend = RDA2), arrow = arrow(length = unit(0.1, 'cm')), size = 0.3, color = 'blue') +
  geom_text(data = rda_tb_forward_r.env, aes(RDA1 * 1.1, RDA2 * 1.1, label = sample), color = 'blue', size = 3)+
  geom_label_repel(aes(label =sample, color = group), size = 3, box.padding = unit(0, 'lines'), show.legend = FALSE)

## 输出PDF
pdf(paste0("RDA.", in_levels,".pdf"), width =12, height =10, bg="white", onefile=F)
print(p)
dev.off()
