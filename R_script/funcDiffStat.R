#!/usr/bin/env Rscript

set.seed(666)
## Load R packages
suppressWarnings(suppressMessages(library(tidytable)))
suppressWarnings(suppressMessages(library(tibble)))
suppressWarnings(suppressMessages(library(broom)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(PMCMRplus)))
suppressWarnings(suppressMessages(library(ggsignif)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(patchwork)))
suppressWarnings(suppressMessages(library(car)))
suppressWarnings(suppressMessages(library(purrr)))
suppressWarnings(suppressMessages(library(data.table)))

# Load functions
# Filter based on sample frequencies
remove_rare <- function(table, cutoff_pro)
{
    row2keep <- c()
    cutoff <- ceiling(cutoff_pro * ncol(table))
    for (i in 1:nrow(table)) {
        row_nonzero <- length(which(table[i, ] > 0))
        if (row_nonzero >= cutoff) {
            row2keep <- c(row2keep, i)
        }
    }
    return(table[row2keep, , drop = FALSE])
}

# Load mapping data
load_map <- function(map)
{
    meta <- fread(map, sep = "\t", header = TRUE, stringsAsFactors = FALSE, colClasses = "character", data.table = FALSE)
    meta <- as.data.frame(meta)
    meta <- meta[, -c(2, 3)]
    meta <- meta[, c(1, 3, 2)]
    return(meta)
}

# Load abundance data
load_data <- function(data, map)
{
    df <- fread(data, sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, quote = "", data.table = FALSE)
    df <- as.data.frame(df)
    rownames(df) <- df[, 1]
    df <- df[, -1]
    df = remove_rare(df, freq)
    delUnanno <- FALSE
    if (delUnanno && length(rownames(df)) > 0) {
        df = df[!grepl("[Uu]nclassified|[Oo]ther|[Uu]nidentified|[Uu]ncultured|[Uu]ncultivated|[Ii]ncertae", gsub(".*__", "", rownames(as.matrix(df))),)]
    }
    # Handle cases where filtering by sample frequency might result in an empty dataset; prompt for adjustments
    if (length(rownames(df)) > 0) {
        # If the last column is a string, remove it
        if (!is.numeric(df[, ncol(df)])) {
            df <- df[-ncol(df)]
        }
        # If a mapping table exists, match sample names with abundance data
        if (!is.null(map)) {
            df = df[, match(map[, 1], colnames(df))]
            colnames(df) = map[match(colnames(df), map[, 1]), 2]
        }
        # Filter by top N
        df = get_topn_data(df, map, topN)
        df = t(df)
        df = cbind(rownames(df), df)
        df = as.data.frame(df)
        if (!is.null(map)) {
            df$Group = factor(as.character(map[, 3]), levels = unique(as.character(map[, 3])))
            df[, 2:(ncol(df) - 1)] = apply(df[, 2:(ncol(df) - 1), drop = FALSE], 2, function(x) { as.numeric(x) })
        } else {
            df[, 2:ncol(df)] = apply(df[, 2:ncol(df), drop = FALSE], 2, function(x) { as.numeric(x) })
        }
        # Take the logarithm of abundance data
        if (transform == "log10") {
            df[, 2:(ncol(df) - 1)] = apply(df[, 2:(ncol(df) - 1), drop = FALSE], 2, function(x) { log10(x + 1) })
        } else if (transform == "log2") {
            df[, 2:(ncol(df) - 1)] = apply(df[, 2:(ncol(df) - 1), drop = FALSE], 2, function(x) { log2(x + 1) })
        } else {
            df = df
        }
    }
    return(df)
}

# Get top N: by group or sample
get_topn_data <- function(data, map, topn)
{
    if (process == "group" && !is.null(map)) {
        group = as.character(unique(map[, 3]))
        data_group = list()
        data_group = data[, map[which(map[, 3] == group[1]), 1], drop = FALSE]
        data_g = data.frame(rowMeans(as.data.frame(data_group)))
        colnames(data_g) = group[1]
        if (length(group) > 1) {
            for (i in 2:length(group)) {
                data_group = data[, map[which(map[, 3] == group[i]), 1]]
                data_tmp = data.frame(rowMeans(as.data.frame(data_group)))
                colnames(data_tmp) = group[i]
                data_g = cbind(data_g, data_tmp)
            }
        }
        data_g = data.frame(rowSums(data_g) / length(group))
        colnames(data_g) = "Overall"
        processed_data = as.data.frame(sort(rowMeans(as.data.frame(data_g)), decreasing = TRUE))
        colnames(processed_data) = "Overall"
    } else {
        processed_data <- sort(rowMeans(as.data.frame(data)), decreasing = TRUE)
        processed_data <- as.data.frame(processed_data)
    }
    topn = ifelse(topn > nrow(processed_data), nrow(processed_data), topn)
    top_tax = rownames(as.data.frame(processed_data))[1:topn]
    newdf = data[match(top_tax, rownames(data)), , drop = FALSE]
    return(newdf)
}


## 主题样式
get_main_theme <- function(justification = c(0.5,0.5),position="right",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),sub_family = "sans",sub_face = "plain",sub_colour = "black",sub_size = 14,sub_hjust = 0.5, sub_vjust = 0.5,sub_angle = 0,fill="white",grid_color="white",family="sans",face = "plain",colour="black",size=16,hjust = 0.5,vjust = 0.5,angle=0)
{
 return(theme(legend.justification = justification,legend.position=position,
                      plot.margin = margin,
					  plot.background=element_rect(fill='transparent'),
                      strip.text = element_text(family = sub_family,                     ## 子标题字体
                      face = sub_face,                                                         ## 子标题加粗、倾斜等格式
                      colour = sub_colour,                                              ## 子标题颜色
                      size = sub_size,                                                 ## 子标题字号
                      hjust = sub_hjust,                                                            ## 其他
                      vjust = sub_vjust,                                                            ## 其他
                      angle = sub_angle),
                      panel.background=element_rect(fill=fill),                            ##背景色
                      panel.grid=element_line(color=grid_color), ## 图表样式网格线
                      plot.title=element_text(family = family,                       ## 主标题字体
                      face = face,                                                         ## 主标题加粗、倾斜等格式
                      colour = colour,                                              ## 主标题颜色
                      size = size,                                                 ## 主标题字号
                      hjust = hjust,                                                            ## 其他
                      vjust = vjust,                                                            ## 其他
                      angle = angle))) 
}

get_x_theme <- function(title_family="sans", title_face = "plain", title_colour = "black", title_size=14,  title_hjust = 0.5, title_vjust = 0.5, title_angle=0, 
	text_family="sans",text_face = "plain",text_colour = "black",text_size=12,text_hjust = 0.5,text_vjust = 0.5,text_angle=0)
{
    return(theme(axis.title.x = element_text(family = title_family,                          ## x轴标题字体
                      face = title_face,                                                         ## x轴标题加粗、倾斜等格式
                      colour = title_colour,                                                  ## x轴标题颜色
                      size = title_size,                                                     ## x轴标题字号
                      hjust = title_hjust,                                                              ## 其他
                      vjust = title_vjust,                                                   ## 其他
					  angle = title_angle),
                      axis.text.x=element_text(family = text_family,                          ## x轴标题字体
                      face = text_face,                                                         ## x轴标题加粗、倾斜等格式
                      colour = text_colour,                                                  ## x轴标题颜色
                      size = text_size,                                                     ## x轴标题字号
                      hjust = text_hjust,                                                              ## 其他
                      vjust = text_vjust,                                                   ## 其他
					  angle = text_angle)))                                                            ## 其他
}

get_y_theme <- function(title_family="sans", title_face = "plain", title_colour = "black", title_size=14,  title_hjust = 0.5, title_vjust = 0.5, title_angle=90,
                        text_family="sans", text_face = "plain", text_colour = "black", text_size=12, text_hjust = 0.5, text_vjust = 0.5, text_angle=0)
{
    return(theme(#axis.title.y = element_blank(),
	            axis.title.y=element_text(family = title_family,                          ## x轴标题字体
                      face = title_face,                                                         ## x轴标题加粗、倾斜等格式
                      colour = title_colour,                                                  ## x轴标题颜色
                      size = title_size,                                                     ## x轴标题字号
                      hjust = title_hjust,                                                              ## 其他
                      vjust = title_vjust,                                                            ## 其他
					  angle = title_angle),
                axis.text.y=element_text(family = text_family,                          ## x轴标题字体
                      face = text_face,                                                         ## x轴标题加粗、倾斜等格式
                      colour = text_colour,                                                  ## x轴标题颜色
                      size = text_size,                                                     ## x轴标题字号
                      hjust = text_hjust,                                                              ## 其他
                      vjust = text_vjust,                                                            ## 其他
					  angle = text_angle)))                                                            ## 其他
}

get_legend_theme <- function(title_family = "sans",title_face = "plain",title_colour = "black",title_size = 14,title_hjust = 0.5,title_vjust = 0.5,title_angle = 0,text_family = "sans",text_colour = "black",text_size = 14,text_hjust = 0,text_vjust = 0.5,text_angle = 0)
{
    return(theme(legend.key.width=unit(0.8,'cm'),legend.key.height=unit(0.8,'cm'),
	                  legend.title=element_text(family = title_family,                   ## y轴标题字体
                      face = title_face,                                                         ## y轴标题加粗、倾斜等格式
                      colour = title_colour,                                            ## y轴标题颜色
                      size = title_size,                                               ## y轴标题字号
                      hjust = title_hjust,                                                              ## 其他
                      vjust = title_vjust,                                                            ## 其他
                      angle = title_angle),                                                             ## 其他
                      legend.text=element_text(family = text_family,                     ## y轴标题字体
                      colour = text_colour,                                             ## y轴标题颜色
                      size = text_size,                                                ## y轴标题字号
                      hjust = text_hjust,                                                              ## 其他
                      vjust = text_vjust,                                                            ## 其他
                      angle = text_angle)))                                                             ## 其他
}

# 事后检验方法
posthoc_test = function(formula, data, posthoc_method){
	if (posthoc_method == "Steel_dwass") {
		dunn= dscfAllPairsTest(formula,data=data)
	} else if (posthoc_method == "Nemenyi") {
		dunn= kwAllPairsNemenyiTest(formula,data=data)
	} else if (posthoc_method == "Conover-iman") {
		dunn= PMCMRplus::kwAllPairsConoverTest(formula,data=data)
	} else{
		dunn= kwAllPairsDunnTest(formula,data=data)
	}
	r1=melt(round(dunn$p.value,6))
	r1=r1[which(!is.na(r1[,3])),]
}

# 获取盒型图数据
boxstats_per=function(data,group)
{
	data=as.numeric(data)
	group=as.character(group)
	Groups=unique(group)
	result=c()
	for(i in 1:length(Groups)){
		tmp=data[which(group==Groups[i])]
		r=boxplot.stats(tmp)$stats
		result=rbind(result,data.frame("group"=Groups[i],"pos"=c("bottom whisker", "first quartile", "median", "third quartile", "top whisker"),"value"=r))
	}
	return(result)
}

## 箱线图、提琴图、柱状图
alpha_plot=function(data,nindex,ngroup,main_theme,x_theme,y_theme,legend_theme,col,map,index_order, add_background)
{
	data[,2:(ncol(data)-1)]=apply(data[,2:(ncol(data)-1), drop=FALSE],2,function(x){as.numeric(as.character(x))})
	data=t(data[,2:(ncol(data)-1)])
	data=as.data.frame(data)
	data=data[index_order,]
	data2=as.data.frame(melt(cbind(rownames(data),data),id=1))
	data2[,4]=map[data2[,2],3]
	colnames(data2)=c("Var1","sample","value","Group")

	if(plot_type == "barplot"){
		data2 <- data2 %>% 
			group_by(Var1, Group) %>%
			summarise(value = mean(value), .groups="keep")
		data2$Var1=factor(as.character(data2$Var1),levels=index_order)
		data2$Group=factor(as.character(data2$Group),levels=unique(map[,3]))
		data2$Y=factor(as.character(paste0(data2$Var1,as.numeric(data2$Group))),
					levels=paste0(rep(index_order,each=2+length(levels(data2$Group))),rep(0:(length(unique(data2$Group))+1),t=length(index_order))))
		# 画柱状图			
		p <- ggplot(data=data2, aes(x=value,y=Y))
		if(add_background){
			for (i in 1:(length(index_order))){
				p <- p + annotate('rect', ymin = 2+(i-1)*(9+(length(unique(map[,3]))-2)*2-2), ymax=2+(i-1)*(9+(length(unique(map[,3]))-2)*2-2)+6+(length(unique(map[,3]))-2)*2-2, xmin=-Inf, xmax=Inf, fill='#D3D3D3')
			}
		}
		p <- p + geom_bar(aes(fill=Group,color=Group), stat="identity", width=2)
	}else{
		data2$Var1=factor(as.character(data2$Var1),levels=index_order)
		data2$Group=factor(as.character(data2$Group),levels=unique(map[,3]))
		data2$Y=factor(as.character(paste0(data2$Var1,as.numeric(data2$Group))),
					levels=paste0(rep(index_order,each=2+length(levels(data2$Group))),rep(0:(length(unique(data2$Group))+1),t=length(index_order))))
		violin <- geom_violin(size=0.5,width=1.4,aes(fill=Group,color=Group),alpha=0.5)
		boxplot <- geom_boxplot(size=0.5,width=1.4,aes(fill=Group,color=Group),alpha=0,outlier.colour = "black",outlier.shape = 19,outlier.size = 1.5)
		p <- ggplot(data=data2, aes(x=value,y=Y))
		if(add_background){
			for (i in 1:(length(index_order))){
				p <- p + annotate('rect', ymin = 2+(i-1)*(9+(length(unique(map[,3]))-2)*2-2), ymax=2+(i-1)*(9+(length(unique(map[,3]))-2)*2-2)+6+(length(unique(map[,3]))-2)*2-2, xmin=-Inf, xmax=Inf, fill='#D3D3D3')
			}
		}
		if(grepl("violin", plot_type)){
			# 画提琴图部分
			p <- p + violin
		}
		if(grepl("boxplot", plot_type)){
			# 画盒型图部分
			p <- p + boxplot
		}
	}
	add_grid = ifelse(add_background, FALSE, add_grid)
	# 设置标题
	switch(transform,
		"log10" = {xt1 <- paste0("Log10-transformed ", abundance, "+1")},
		"log2" = {xt1 <- paste0("Log2-transformed ", abundance, "+1")},
		{xt1 <- abundance}
	)
	x_title1 = ifelse(x_title1 == "", xt1, x_title1)
	y_title = ifelse(y_title == "", "", y_title)
	p <- p+xlab(x_title1)+ 
	ylab(y_title)+
	scale_y_discrete(
		limits = paste0(rep(index_order,each=(length(levels(data2$Group))+1)*2+1),seq(from=0, to=length(unique(data2$Group))+1, by=0.5)),
		breaks= paste0(index_order,(length(levels(data2$Group))+1)/2),
		labels= index_order,
		expand = expansion(mult =0,add = 4)
	)+
	scale_color_manual(values = col,guide = "legend")+
	theme_bw()+                                                                                      ## 清空主题
	# labs(title=main_title)+                                                                          ## 主标题内容
	# main_theme+                                                                                      ## 添加整体主题
	scale_fill_manual(values=col,guide = "legend")+                                                  ## 填充颜色
	guides(fill=guide_legend(title=legend_title),
		   color=guide_legend(title=legend_title,nrow =min(ifelse((nindex%in%c(4,6)),20,8),ngroup)))+               ## 图例标题内容
	y_theme+                                                                                         ## 添加y轴主题
	x_theme+                                                                                         ## 添加y轴主题
	legend_theme+
	theme(
		plot.background = element_blank(),
		panel.border=element_rect(fill='transparent',color="black"),
		panel.background=element_rect(fill='transparent',color="black"),
		title = element_blank(),
		panel.grid=element_line(color=ifelse(add_grid, "grey90", "white")),
		legend.background = element_blank(),
		axis.text.y = element_text(hjust = 1)
	)
	return(p)
}

## 差异比较p值
add_signif1=function(data,ngroup,index_order,x_theme)
{
	data[,2:(ncol(data)-1)]=apply(data[,2:(ncol(data)-1), drop=FALSE],2,function(x){as.numeric(x)})
	data = data[, -1]

	normal_test = data.frame()
	bartlett_test = data.frame()
	for(index in 1:(ncol(data)-1)){
		tmp_data = data[, c(index, ncol(data))]
		# 检验每个物种数据正态性
		tryCatch(
			{
				# 数据做正态性检验
				tmp_levene_result <- tmp_data %>% select_if(is.numeric) %>% map_dfr(~ broom::tidy(shapiro.test(melt(.)$value)),.id='var')
				normal_test[colnames(tmp_data)[1], "variable"] <- colnames(tmp_data)[1]
				normal_test[colnames(tmp_data)[1], "statistic"] <- round(tmp_levene_result$statistic, 6)
				normal_test[colnames(tmp_data)[1], "p.value"] <- round(tmp_levene_result$p.value, 6)
				normal_test[colnames(tmp_data)[1], "method"] <- "Shapiro-Wilk normality test"
			},
			error = function(e){
				# 数据都相同时不能做正态性检验，设P值为1，不通过检验
				normal_test[colnames(tmp_data)[1], "variable"] <- colnames(tmp_data)[1]
				normal_test[colnames(tmp_data)[1], "statistic"] <- "-"
				normal_test[colnames(tmp_data)[1], "p.value"] <- 1
				normal_test[colnames(tmp_data)[1], "method"] <- "Shapiro-Wilk normality test"
				normal_test <<- normal_test
			}
		)
		# 检验每个物种数据方差齐性
		tryCatch(
			{
				bartlett_result <- tmp_data %>% select_if(is.numeric) %>% map_dfr(~ broom::tidy(bartlett.test(. ~ Group,data=data)),.id='var')
				bartlett_test[colnames(tmp_data)[1], "variable"] <- colnames(tmp_data)[1]
				bartlett_test[colnames(tmp_data)[1], "statistic"] <- round(bartlett_result$statistic, 6)
				bartlett_test[colnames(tmp_data)[1], "p.value"] <- round(bartlett_result$p.value, 6)
				bartlett_test[colnames(tmp_data)[1], "method"] <- "Bartlett test of homogeneity of variances"
			},
			error = function(e){
				bartlett_test[colnames(tmp_data)[1], "variable"] <- colnames(tmp_data)[1]
				bartlett_test[colnames(tmp_data)[1], "statistic"] <- "-"
				bartlett_test[colnames(tmp_data)[1], "p.value"] <- 1
				bartlett_test[colnames(tmp_data)[1], "method"] <- "Bartlett test of homogeneity of variances"
				bartlett_test <<- bartlett_test
			}
		)
	}
	write.table(normal_test,paste0(outdir,"/Normality_P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)
	write.table(bartlett_test,paste0(outdir,"/bartlett_P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)
	# 差异检验
	if(ngroup == 2){
		if(test_method == "Student's test"){
			# 根据每个物种方差齐性逐一设置var.equal
			tmp_data = data[, c(1, ncol(data))]
			tmp_bartlett_value <- bartlett_test[colnames(tmp_data)[1], "p.value"]
			tmp_equal <- ifelse(tmp_bartlett_value <= 0.05 || is.na(tmp_bartlett_value), FALSE, TRUE)
			diff <- tmp_data %>% select_if(is.numeric) %>% map_dfr(~ broom::tidy(t.test(. ~ Group,data=tmp_data, var.equal=tmp_equal)),.id='var')
			diff <- as.data.frame(diff[,c("var","p.value")])
			if(ncol(data)-1 > 1){
				for(index in 2:(ncol(data)-1)){
					tmp_data = data[, c(index, ncol(data))]
					tmp_bartlett_value <- bartlett_test[colnames(tmp_data)[1], "p.value"]
					tmp_equal <- ifelse(tmp_bartlett_value <= 0.05 || is.na(tmp_bartlett_value), FALSE, TRUE)
					# tmp_levene_result <- tmp_data %>% select_if(is.numeric) %>% map_dfr(~ broom::tidy(bartlett.test(. ~ Group,data=tmp_data)),.id='var')
					# tmp_equal <- ifelse(tmp_levene_result$p.value > 0.05, TRUE, FALSE)
					tmp_diff <- tmp_data %>% select_if(is.numeric) %>% map_dfr(~ broom::tidy(t.test(. ~ Group,data=tmp_data, var.equal=tmp_equal)),.id='var')
					diff = rbind(diff, as.data.frame(tmp_diff[,c("var","p.value")]))
				}
			}
		}else{
			diff <- data %>% select_if(is.numeric) %>% map_dfr(~ broom::tidy(wilcox.test(. ~ Group,data=data)),.id='var')
		}
	}else{
		if(test_method == "ANOVA"){
			# diff <- data %>% select_if(is.numeric) %>% map_dfr(~ broom::tidy(oneway.test(. ~ Group,data=data)),.id='var')
			# 根据每个物种方差齐性逐一设置var.equal
			tmp_data = data[, c(1, ncol(data))]
			tmp_bartlett_value <- bartlett_test[colnames(tmp_data)[1], "p.value"]
			tmp_equal <- ifelse(tmp_bartlett_value <= 0.05 || is.na(tmp_bartlett_value), FALSE, TRUE)
			diff <- tmp_data %>% select_if(is.numeric) %>% map_dfr(~ broom::tidy(oneway.test(. ~ Group,data=tmp_data, var.equal=tmp_equal)),.id='var')
			diff <- as.data.frame(diff[,c("var","p.value")])
			if(ncol(data)-1 > 1){
				for(index in 2:(ncol(data)-1)){
					tmp_data = data[, c(index, ncol(data))]
					tmp_bartlett_value <- bartlett_test[colnames(tmp_data)[1], "p.value"]
					tmp_equal <- ifelse(tmp_bartlett_value <= 0.05 || is.na(tmp_bartlett_value), FALSE, TRUE)
					# tmp_levene_result <- tmp_data %>% select_if(is.numeric) %>% map_dfr(~ broom::tidy(bartlett.test(. ~ Group,data=tmp_data)),.id='var')
					# tmp_equal <- ifelse(tmp_levene_result$p.value > 0.05, TRUE, FALSE)
					tmp_diff <- tmp_data %>% select_if(is.numeric) %>% map_dfr(~ broom::tidy(oneway.test(. ~ Group,data=tmp_data, var.equal=tmp_equal)),.id='var')
					diff = rbind(diff, as.data.frame(tmp_diff[,c("var","p.value")]))
				}
			}
		}else{
			diff <- data %>% select_if(is.numeric) %>% map_dfr(~ broom::tidy(kruskalTest(. ~ Group,data=data)),.id='var')
		}
	}
	diff_result = as.data.frame(diff[,c("var","p.value")])
	diff_result$adj.p.value = p.adjust(diff_result$p.value,method=adjust_method)
	filter_result <- diff_result %>% filter(round(p.value, 6) <= filter_pvalue)
	rownames(diff_result)=diff_result$var
	rownames(filter_result)=filter_result$var
	write.table(diff_result,paste0(outdir,"/P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)
	diff_result$var=paste0(diff_result$var,rep(c(1:length(unique(data$Group)))[round((ngroup+1)/2)],t=length(index_order)))
	diff_result$var=paste0(index_order,(length(levels(data$Group))+1)/2)
	index_order = intersect(index_order,rownames(filter_result))
	diff_result=diff_result[index_order,]
	diff_result$p.value = round(diff_result$p.value, 6)
	diff_result$adj.p.value = round(diff_result$adj.p.value, 6)

	if(length(index_order) == 0){
		p= ggplot(diff_result,aes(var,p.value))
	}else{
		p= ggplot(diff_result,aes(var,p.value)) +
			geom_text(aes(x=0,y=var),label=round(diff_result$p.value,6), size = 5,hjust=0.5,vjust=0,size=6)+
			scale_y_discrete(limits = paste0(rep(index_order,each=(length(levels(data$Group))+1)*2+1),seq(from=0, to=length(unique(data$Group))+1, by=0.5)),
						breaks= paste0(index_order,(length(levels(data$Group))+1)/2),
						labels= index_order,expand = expansion(mult =0,add = 4)
						)+
			xlab(x_title3)+
			x_theme+
			theme(plot.background=element_rect(fill='transparent'),
				panel.border=element_blank(),
				panel.background=element_rect(fill='transparent',color="black"),
				panel.grid=element_blank(),
				axis.line=element_blank(),
				axis.ticks=element_blank(),
				axis.text=element_blank(),
				axis.text.x=element_blank(),
				axis.text.y=element_blank(),
				axis.title.y=element_blank()
			)
	}
	return (list(p, index_order))
}

## 事后检验
add_signif2=function(data,ngroup,index_order,x_theme)
{
    data[,2:(ncol(data)-1)]=apply(data[,2:(ncol(data)-1), drop=FALSE],2,function(x){as.numeric(x)})
	data=data.frame(V1=data[,1],data[,match(index_order, colnames(data))],Group=data[,ncol(data)], check.names=F)

	tmp <- data %>% select_if(is.numeric) %>% map_dfr(~ posthoc_test(. ~ Group, data=data, posthoc_method = posthoc_method),.id='var')
	tmp$var=stringr::str_replace_all(tmp$var,"\\."," ")
	write.table(tmp,paste0(outdir,"/posthoc_test-P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)
    dff=tmp[which(tmp$value!="NA" & tmp$value <= 0.05),]
	# write.table(dff,paste0(outdir,"/posthoc_test-P.value-sig.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)

	if(nrow(dff)>0){
        # limitsX = paste0(rep(index_order,each=2+ngroup),rep(0:(ngroup+1),t=length(index_order)))
		limitsX = paste0(rep(index_order,each=(length(levels(data$Group))+1)*2+1),seq(from=0, to=length(unique(data$Group))+1, by=0.5))
		dff$Var1=factor(as.character(dff$Var1),levels=unique(data$Group))
		dff$Var2=factor(as.character(dff$Var2),levels=unique(data$Group))
        p=ggplot(data=dff)
		tax=dff$var[1]
		h=1.5
		for (j in 1:nrow(dff)){
	        # xmin=as.numeric(factor(paste0(dff$var[j],as.numeric(dff$Var2[j])),levels=c(limitsX)))
			# xmax=as.numeric(factor(paste0(dff$var[j],as.numeric(dff$Var1[j])),levels=c(limitsX)))
			xmin=paste0(dff$var[j],as.numeric(dff$Var2[j]))
			xmax=paste0(dff$var[j],as.numeric(dff$Var1[j]))
			text_xmin = as.numeric(factor(paste0(dff$var[j],as.numeric(dff$Var2[j])),levels=c(limitsX)))
			text_xmax = as.numeric(factor(paste0(dff$var[j],as.numeric(dff$Var1[j])),levels=c(limitsX)))
            if(tax==dff$var[j]){h=h+0.25}else{h=1.5; tax=dff$var[j] }
			if(dff$value[j]<=0.05 & dff$value[j] > 0.01) {
				p = p+geom_signif(data=dff,tip_length=0,xmin=xmin,xmax=xmax,
					# y_position=h,aes(annotations="*"),inherit.aes=F,manual = T)+
					y_position=h, aes(annotations=""), inherit.aes=T,manual = T)+
					# geom_text(text="*", x=(text_xmin+text_xmax)/2, y=h+0.05, label="*", angle=-90)+
					annotate("text", x=(text_xmin+text_xmax)/2, y=h+0.05, label="*", angle=-90)
					# coord_flip() 
			}
			if(dff$value[j]<=0.01 & dff$value[j] > 0.001) {
				p = p+geom_signif(data=dff,tip_length=0,xmin=xmin,xmax=xmax,
					# y_position=h,aes(annotations="**"),inherit.aes=T,manual = T, map_signif_level=T)+
					y_position=h,aes(annotations=""),inherit.aes=T,manual = T)+
					# geom_text(text="**", x=(text_xmin+text_xmax)/2, y=h+0.05, label="**", angle=-90)+
					annotate("text", x=(text_xmin+text_xmax)/2, y=h+0.05, label="**", angle=-90)
					# coord_flip() 
			}
			if(dff$value[j]<=0.001 ) {
				p = p+geom_signif(data=dff,tip_length=0,xmin=xmin,xmax=xmax,
					# y_position=h,aes(annotations="***"),inherit.aes=T,manual = T, map_signif_level=T) +
					y_position=h,aes(annotations=""),inherit.aes=T,manual = T) +
					# geom_text(text="***", x=(text_xmin+text_xmax)/2, y=h+0.05, label="***", angle=-90)+
					annotate("text", x=(text_xmin+text_xmax)/2, y=h+0.05, label="***", angle=-90)
					# coord_flip() 
			}
		}
	
		x_title2 <- ifelse(x_title2 == "", paste0(posthoc_method, " test"), x_title2)
		p = p+
			coord_flip()+
			scale_x_discrete(
				limits = paste0(rep(index_order,each=(length(levels(data$Group))+1)*2+1),seq(from=0, to=length(unique(data$Group))+1, by=0.5)),
				breaks= paste0(index_order,(length(levels(data$Group))+1)/2),
				labels= index_order,
				expand = expansion(mult =0,add = 4)
			)+
			# scale_x_continuous(limits = c(-3,length(limitsX)+4),expand = expansion(mult =0,add = 0)) +
			ylab(x_title2) + 
			# scale_y_continuous(expand =expand_scale(mult=c(0.1,0.2)))+
			scale_y_continuous(expand =expand_scale(mult=0.5))+
			x_theme +
			theme(plot.background=element_rect(fill='transparent'),
				panel.border=element_rect(fill='transparent',color="black"),
				panel.background=element_rect(fill='transparent',color="black"),
				panel.grid=element_blank(),
				axis.line=element_blank(),
				axis.ticks=element_blank(),
				axis.text=element_blank(),
				axis.text.y=element_blank(),
				axis.text.x=element_blank(),
				axis.title.y=element_blank(),
			)
	}
}

# 画errorbar
add_errorbar=function(data,map,ngroup,index_order,x_theme, col, add_background)
{
    data[,2:(ncol(data)-1)]=apply(data[,2:(ncol(data)-1), drop=FALSE],2,function(x){as.numeric(x)})
    data=data.frame(V1=data[,1],data[,match(index_order, colnames(data))],Group=data[,ncol(data)], check.names=F)

	# tmp <- data %>% select_if(is.numeric) %>% map_dfr(~ Nemenyi_test(. ~ Group, data=data),.id='var')
	tmp <- data %>% select_if(is.numeric) %>% map_dfr(~ posthoc_test(. ~ Group, data=data, posthoc_method = posthoc_method),.id='var')
	tmp$var=stringr::str_replace_all(tmp$var,"\\."," ")
	write.table(tmp,paste0(outdir,"/posthoc_test-P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)

	test_result <- data %>% select_if(is.numeric) %>% map_dfr(~ broom::tidy(t.test(. ~ Group,data=data)),.id='var')
	diff_result<- test_result[, c("var", "estimate", "conf.low", "conf.high")]
	diff_result$Group <- c(ifelse(diff_result$estimate > 0, levels(data$Group)[1], levels(data$Group)[2]))
	# diff_result$varY = paste0(diff_result$var,rep(c(1:length(unique(diff_result$Group)))[round((ngroup+1)/2)],t=length(index_order)))
	diff_result$varY = paste0(diff_result$var,(ngroup+1)/2)
	
	p <- ggplot(diff_result, aes(varY, estimate, fill=Group))+
		x_theme+
		scale_x_discrete(
			limits = paste0(rep(index_order,each=5+length(levels(data$Group))),seq(from=0, to=length(unique(data$Group))+1, by=0.5)),
			breaks= paste0(index_order,rep(c(1:length(unique(data$Group)))[round((ngroup+1)/2)],t=length(index_order))),
			labels= index_order,
			expand = expansion(mult =0,add = 4)
		)+ 
		coord_flip()
	if(add_background){
		for (i in 1:(length(index_order))){
			p <- p + annotate('rect', xmin = 2+(i-1)*(9+(length(unique(map[,3]))-2)*2-2), xmax=2+(i-1)*(9+(length(unique(map[,3]))-2)*2-2)+6+(length(unique(map[,3]))-2)*2-2, ymin=-Inf, ymax=Inf, fill='#D3D3D3')
		}
	}
	p <- p + geom_errorbar(aes(ymin = conf.low, ymax=conf.high), width=0.1, size=0.1)+
		ylab(x_title2) + 
		geom_point(shape=21, size=3)+
		scale_fill_manual(values=col)+
		geom_hline(aes(yintercept=0), linetype="dashed", color="black")+
		theme(
			plot.background = element_blank(),
			panel.border=element_rect(fill='transparent',color="black"),
			panel.background=element_rect(fill='transparent',color="black"),
			panel.grid=element_blank(),
			axis.line=element_blank(),
			axis.ticks.y=element_blank(),
			axis.text.y=element_blank(),
			axis.title.y=element_blank(),
			legend.position = "none"
		)
		
}

get_plot=function(df2,map,x_text_rotate,main_theme,x_theme,y_theme,legend_theme,col,index, add_background)
{

	ngroup=length(levels(df2[,"Group"]))
	nindex=length(index)
	 
	# 先进行差异检验，筛选P值重新给index赋值，过滤不显著的
	q2_result=add_signif1(df2,ngroup,index,x_theme)
	q2 = q2_result[[1]]
	index = q2_result[[2]]
	if(length(index) > 0){
		filter_df2 <- df2[, match(index, colnames(df2)), drop=FALSE]
		filter_df2=t(filter_df2)
		filter_df2=cbind(rownames(filter_df2),filter_df2)
		filter_df2=as.data.frame(filter_df2)
		write.table(filter_df2,paste0(outdir,"/Abundance.tsv"),sep="\t",row.names = F,col.names =T,quote=FALSE)
		q1 = alpha_plot(df2,nindex,ngroup,main_theme,x_theme,y_theme,legend_theme,col,map,index, add_background)
		if(ngroup == 2){
			q3=add_errorbar(df2,map,ngroup,index,x_theme, col, add_background)
			df_blank <- data.frame()
			write.table(df_blank,paste0(outdir,"/posthoc_test-P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)
		}else{
			if(dunn){
				q3=add_signif2(df2,ngroup,index,x_theme)
			}else{
				df_blank <- data.frame()
				write.table(df_blank,paste0(outdir,"/posthoc_test-P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)
			}
		}
		## 绘图
		theme_temp=theme(legend.position = "none",plot.title=element_blank())

		# 图片是水平展示，原来高宽进行调换
		height=(ngroup*box_height/3+0.5)*nindex+1
		height=ifelse(height <= 10, 10, height)

		width=max(nchar(as.character(index))-1)/7/90*x_text_rotate*1+8
		if(ngroup > 2 && !dunn){
			width = width - 1.5
		}
		width=ifelse(width <= 10, 10, width)
		main_title = ifelse(main_title=="", test_method, main_title)
		if(ngroup == 2 || dunn){
			p= q1 + q3 + q2 + plot_layout(ncol = 3, widths  = c(5,1.5,1.5)) + 
			plot_layout(guides = "collect")+
			plot_annotation(title = main_title, theme=main_theme)
		}else{
			p= q1 + q2 + plot_layout(ncol = 2, widths  = c(5,1.5)) + 
			plot_layout(guides = "collect")+
			plot_annotation(title = main_title, theme=main_theme)
		}
		title=plot_grid(get_title(q1))
		return_list = list(p,width,height,title)
	}else{
		df_blank <- data.frame()
		write.table(df_blank,paste0(outdir,"/Abundance.tsv"),sep="\t",row.names = F,col.names =T,quote=FALSE)
		write.table(df_blank,paste0(outdir,"/bartlett_P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)
		write.table(df_blank,paste0(outdir,"/Normality_P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)
		write.table(df_blank,paste0(outdir,"/P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)
		write.table(df_blank,paste0(outdir,"/posthoc_test-P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)
		return_list = list(index)
	}
	return(return_list)
}

color.ls <- list(
    "Classic" = c("#FF6A00","#FFD11A","#99A7FF","#5FC6FF","#ABDF67","#56D6B6","#80B3FF","#C06CEB","#FF6262","#815293","#FFCA80","#FF9933","#FF964A","#FFCC00","#A2AFFF","#6CCBFF","#ACEC58","#49F7CC","#A9CAFB","#DB94FF","#FF7F7F","#AF67CA","#FFB958","#FF9831"),
    "Vitality" = c("#2E8AE6","#36aee4","#17D8E6","#3DCC9C","#8ADD7A","#C0E673","#FEDF66","#FEB266","#FA7D64","#e7688f","#b45fbb","#717cdc","#5CA3EA","#7ACEF3","#79EEF6","#63DDB4","#A6EC99","#D2EC9D","#FCEBA7","#F9CA9B","#FCA08D","#F085A6","#C880CE","#939BDE"),
    "Retro" = c("#187679","#E6382E","#FFD119","#C27F35","#FE685F","#887A81","#CC525C","#FA9A4B","#44BBA4","#A7A09F","#5C9966","#86BF8B","#418E91","#E7716A","#F3D455","#C78E4F","#FE685F","#A8959F","#CC666E","#FAB276","#5BBDAA","#BDB3B1","#6DB97A","#A0D6A5"),
    "Summer" = c("#FF735A","#5BCCB6","#4782B3","#EB7BC0","#FFAA33","#FED840","#AB80F5","#EACC93","#5C9966","#A0C896","#EB8D8E","#CEAFAF","#FF9B89","#89DACB","#6FA1CA","#EC9DCE","#FFC26D","#FEE376","#C3A4F8","#F0DAB2","#77C785","#BBD8B4","#F1ADAE","#DCC6C6"),
    "Funny" = c("#7C54FF","#FED840","#2D928F","#69DFDC","#FFA367","#BDB6E4","#5AA4C9","#98D9AD","#5490DE","#C574D5","#FF5C5C","#9777FF","#FEE067","#58A8A6","#86E3E1","#FFB686","#CAC5E9","#7CB6D4","#8BE5A8","#77A7E5","#D190DD","#FF7D7D","#89D7E1"),
    "Refresh" = c("#F3C625","#439A60","#6DBDF2","#4670C2","#FF9139","#BDB6E4","#F35F25","#ECBA82","#8F6DD0","#C24664","#CF8395","#A69374","#F5D151","#69AE80","#8ACAF5","#6B8DCE","#FFA761","#CAC5E9","#F57F51","#F0C89B","#A58AD9","#CE6B83","#D99CAA","#B8A990"),
    "Mystery" = c("#3369FF","#36D2B8","#7258F3","#DB8D76","#47B359","#8545E6","#E0B18B","#0095FF","#52CC8F","#6675FF","#CC516E","#5CA9E6","#5C87FF","#5EDBC6","#8E79F5","#E2A491","#6CC27A","#9D6AEB","#E6C1A2","#33AAFF","#75D6A5","#8591FF","#D6748B","#7DBAEB"),
    "Spring" = c("#0F8299","#12B3B3","#6CD99F","#ACE673","#FFDB4D","#4C993D","#91B362","#40886A","#84C2B7","#278C78","#35CF80","#99B3CC","#3EA4B8","#57C8C8","#82DBAC","#B9EA8A","#FDDA54","#71B963","#B8D78E","#63B794","#98E2D5","#4EAC99","#78DCA9","#96B8D9"),
	"Science" =c ("#3A4892", "#ED0000", "#008B45", "#631879", "#00817F", "#BA0020", "#5F549A", "#A10055", "#7F807F", "#A93B52", "#9F6D25", "#545A62", "#50557D", "#8A5C4E", "#9A3C5D", "#883A77", "#96546A", "#626489", "#7E4372", "#C85113", "#407354", "#5C3B7B", "#657167", "#AC293F", "#774989", "#9D3860", "#727285", "#CC2C31", "#6E7E35", "#5E3F6D", "#3C6C7E", "#A54137", "#83497B", "#962866", "#8D6B75", "#51568D", "#666D92", "#ED7777", "#468B68", "#6E4979", "#418180", "#BA5D6D", "#7D779A", "#A1517B", "#808080", "#A9727D", "#9F8662", "#5B5E62", "#67697D", "#8A736C"),
    "NPG" = c("#E54B34", "#4CBAD4", "#009F86", "#3B5387", "#F29A7F", "#8491B3", "#91D1C1", "#DC0000", "#7E6047", "#CCCCCC", "#BC8B83", "#33ADAD", "#347988", "#9F7685", "#C1969A", "#8BB0BB", "#CE8662", "#B04929", "#A59487", "#E3907E", "#D46F5B", "#41B4C1", "#278C87", "#726486", "#DA988C", "#88A0B7", "#B9AC91", "#C63517", "#927A66", "#DBAEA4", "#97A4AB", "#21A69A", "#3A6688", "#C98882", "#A593A7", "#8EC0BE", "#D85935", "#985738", "#B9AFA9", "#E67059", "#E5BFB9", "#B2CED4", "#779F99", "#747A87", "#F2DCD5", "#A7ABB3", "#C1D1CD", "#DCA5A5", "#7E7770", "#CCCCCC"),
    "Lancet" = c("#00468B", "#EC0000", "#41B43F", "#0099B3", "#925E9F", "#FDAE91", "#AC002A", "#ACB6B6", "#9E3B4F", "#B2811E", "#40A67D", "#6B7EA9", "#C98599", "#D7675A", "#B86E6B", "#6A7BA1", "#6E4D6D", "#D1793E", "#5EAD73", "#6694AE", "#AE80A1", "#EAA392", "#B46264", "#949CAB", "#C65356", "#8D9D4B", "#4D9F9B", "#8A7CA4", "#E3ABA9", "#C26161", "#B69C9A", "#596D96", "#46688B", "#EC7676", "#7BB47A", "#5AA6B3", "#997F9F", "#FDD6C7", "#AC566B", "#B1B6B6", "#9E6D76", "#B29A68", "#73A692", "#8A93A9", "#C9A7B1", "#D79F99", "#B89392", "#868EA1", "#6E5E6E", "#D1A588"),
    "NEJM" = c("#BB3B28", "#0072B4", "#E08626", "#1F854D", "#7876B1", "#6E99AC", "#DC91FF", "#ED4C97", "#905B6E", "#A07C74", "#928A3C", "#587F7F", "#7487AF", "#A897D5", "#E970C9", "#D6435F", "#A86463", "#7C7E95", "#BA9554", "#52826E", "#858BB0", "#99A2C1", "#E39AE4", "#E26E95", "#747291", "#C2916E", "#6E8856", "#758298", "#8198AE", "#CCAAEA", "#EC82BF", "#C96265", "#BB7B72", "#5A93B4", "#E0B383", "#528569", "#9494B1", "#8DA3AC", "#EEC8FF", "#ED9DC2", "#90767F", "#A08E8A", "#928E67", "#6C7F7F", "#929BAF", "#BFB6D5", "#E9ACD9", "#D68D9B", "#A88686", "#898A95"),
    "JAMA" = c("#374D54", "#DF8E44", "#00A0D4", "#B24645", "#79AE97", "#6A6599", "#7F796B", "#8E6E4F", "#A7988E", "#91778A", "#A07F6C", "#748998", "#776E82", "#5C6360", "#655E52", "#C6946A", "#6F8CAF", "#AB6558", "#779C98", "#71698D", "#6E6D65", "#B67E4A", "#7B9DB1", "#A56167", "#919781", "#6F7798", "#7C7376", "#49585A", "#465154", "#DFB792", "#6ABAD4", "#B27C7B", "#93AEA3", "#827F99", "#7F7C75", "#8E7E6F", "#A7A09B", "#91848E", "#A09086", "#869198", "#7D7882", "#606362", "#65625C", "#C6AD98", "#8F9EAF", "#AB8882", "#8A9C9A", "#7F7B8D", "#6E6E6A", "#B69A80"),
    "Report1" = c("#BF5A17", "#F0017F", "#386CB0", "#FDBF85", "#BEADD3", "#7FC97F", "#FA7F72", "#666666", "#B4B1B1", "#D8434F", "#A95597", "#AA949D", "#E1B6AD", "#A2BCAA", "#C7A978", "#B1746C", "#8C8A8A", "#C28665", "#CC5036", "#CE3F8A", "#7B7FA7", "#EFBB9A", "#B1B5BF", "#A7BA7B", "#D67A6F", "#787878", "#BE9B8A", "#E42F67", "#7D63A3", "#D4A992", "#D0B2C1", "#92C394", "#E29675", "#8C6D69", "#A09D9D", "#C27140", "#BF8C6B", "#F079B7", "#748EB0", "#FDDEC1", "#C9C0D3", "#A4C9A4", "#FABCB6", "#666666", "#B4B3B3", "#D88D93", "#A97FA0", "#AA9FA3", "#E1CCC7", "#AFBCB3"),
    "Report2" = c("#FF6600", "#FFFF66", "#009966", "#FF6666", "#666600", "#CCFFCC", "#669933", "#339966", "#FFB637", "#9ACB68", "#A78A65", "#B26B39", "#9BAF69", "#99CA7E", "#52984E", "#B0893E", "#FFAC57", "#D3E587", "#7A9371", "#D88671", "#83894D", "#BEE4B4", "#6C9857", "#849262", "#FFE47A", "#77B27A", "#D49281", "#8C723C", "#BDD6A8", "#8BB16F", "#5A986A", "#D99353", "#FFB380", "#FFFFB3", "#4D9980", "#FFB3B3", "#666633", "#E6FFE6", "#809966", "#669980", "#FFDB9B", "#B3CB9A", "#A79986", "#B28F76", "#A5AF8C", "#B2CAA4", "#759873", "#B09C77", "#FFD6AB", "#DCE5B6")
)

color.continuous.ls <- list(
    'BrBG' = c('#8C510A', '#F5F1E7', '#01665E'),
    'PiYG' = c('#C51B7D', '#F8F0F4', '#4D9221'),
    'RdBu' = c('#B2182B', '#F8EFE9', '#2166AC'),
    'bl2rd' = c('#0000FF', '#0CE2F2', '#FF0000'),
    'gnbu' = c('#084081', '#6AC1C8', '#F7FCF0'),
    'matlablike' = c('#0000AA', '#A1FFDB', '#AA0000'),
    'Gwr' = c('#0000FF', '#ffffff', '#FF0000'),
    'Bwr' = c('#3399CC', '#ffffff', '#FF6666')
)

###################################################
###################################################

main <- function(conf){
    if (conf != "") {
        source(conf)
    }else{
		mapdir <- indir
	}
	map <- load_map(paste0(mapdir, "/", map))
    df <- load_data(paste0(indir, "/", input),map)
    # index取丰度排序后的物种名
    index <- colnames(df)[2:(ncol(df)-1)]

    if(!file.exists(outdir)){
		dir.create(outdir,recursive=TRUE)
	}

    ## 主题
	add_grid = ifelse(add_background, FALSE, add_grid)
    main_theme <- get_main_theme(grid_color=ifelse(add_grid,"grey","white"),sub_family = "serif",sub_colour = "#000000",sub_size = 12,family=main_title_font,colour=main_title_color,size=main_title_size)
	
    y_theme <- get_y_theme(title_family = y_title_font, title_colour = y_title_color, title_size = y_title_size,
	                       text_family = y_text_font, text_colour = y_text_color, text_size = y_text_size, text_angle = y_text_rotate)
	
    x_theme <- get_x_theme(title_family = x_title_font, title_colour = x_title_color, title_size = x_title_size,
							text_family = x_text_font,text_colour = x_text_color,text_size = x_text_size,text_hjust=ifelse(x_text_rotate==0,0.5,1),text_vjust = ifelse((x_text_rotate==0|x_text_rotate==90),0.5,1),text_angle =  x_text_rotate)
	
    legend_theme <- get_legend_theme(title_family = legend_title_font,title_colour = legend_title_color,title_size =legend_title_size,text_family = legend_text_font,text_colour = legend_text_color,text_size = legend_text_size)

    ## 颜色样式选择
    if (color_name %in% names(color.ls)){
		col=color.ls[[color_name]]
    }else if (color_name %in% names(color.continuous.ls)){
		col=color.continuous.ls[[color_name]]
		col=colorRampPalette(col,space = "Lab")(length(unique(map[,3])))
    }
	
    if (as.numeric(min(summary(map[,3])))<3 | as.numeric(length(unique(map[,3])))<2){  ###有分组小于3个样本	
		p=ggplot(data=data.frame(x=1,y=1))+geom_text(x=0.5,y=0.75,size=12,label=c("分组数、总样本数或分组最小样本数不满足分"))+geom_text(x=0.5,y=0.65,size=12,label=c("析要求; 可调整后再尝试"))+theme_bw()
		df_blank <- data.frame()
		write.table(df_blank,paste0(outdir,"/Abundance.tsv"),sep="\t",row.names = F,col.names =T,quote=FALSE)
		write.table(df_blank,paste0(outdir,"/bartlett_P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)
		write.table(df_blank,paste0(outdir,"/Normality_P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)
		write.table(df_blank,paste0(outdir,"/P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)
		write.table(df_blank,paste0(outdir,"/posthoc_test-P.value.xls"),sep="\t",row.names = F,col.names =T,quote=FALSE)		
		svg(paste0(outdir,"/DiffStat.svg"),
			width=10,
			height=7,
			bg="white", onefile=F)
		print(p)
		dev.off()
		png(paste0(outdir,"/DiffStat.png"),
			width=10,
			height=7,
			units="in",
			res=300,
			bg="white")
		print(p)
		dev.off()

		pdf(paste0(outdir,"/DiffStat.pdf"),
			width=10,
			height=7,
			bg="white", onefile=F)
		print(p)
		dev.off()

    }else{
		if(length(rownames(df)) > 0){
			gg=get_plot(df2=df,map,x_text_rotate,main_theme,x_theme,y_theme,legend_theme,col,index, add_background)
		}else{
			gg=list
		}

		if(length(gg) == 4){
			p=gg[[1]]
			title=gg[[4]]
			vpmain<-viewport(layout =grid.layout(100, 1) ,name = "main")
			vp1=viewport(layout.pos.col=1, layout.pos.row=2)
			vp2=viewport(layout.pos.col=1, layout.pos.row=3:98)
				
			## 输出PDF
			pdf(paste0(outdir,"/DiffStat.pdf"),
					width=gg[[2]],
					height=gg[[3]],
					bg="white", onefile=F)
			grid.newpage()
			pushViewport(vpmain)
			print(title, vp=vp1)
			print(p, vp=vp2)
			dev.off()

			## 输出svg
			svg(paste0(outdir,"/DiffStat.svg"),
					width=gg[[2]],
					height=gg[[3]],
					bg="white", onefile=F)
			grid.newpage()
			pushViewport(vpmain)
			print(title, vp=vp1)
			print(p, vp=vp2)
			dev.off()

			## 输出PNG
			png(paste0(outdir,"/DiffStat.png"),
					width=gg[[2]],
					height=gg[[3]],
					units="in",
					res=100,
					bg="white")
			grid.newpage()
			pushViewport(vpmain)
			print(title, vp=vp1)
			print(p, vp=vp2)
			dev.off()
		}else{
			p=ggplot(data=data.frame(x=1,y=1))+geom_text(x=0.5,y=0.75,size=12,label=c("显著性水平筛选后topN中没有差异显著的功能"))+
			geom_text(x=0.5,y=0.65,size=12,label=c("请适当调整筛选条件!"))+theme_bw()
					
			svg(paste0(outdir,"/DiffStat.svg"),
					width=10,
					height=7,
					bg="white", onefile=F)
			print(p)
			dev.off()

			png(paste0(outdir,"/DiffStat.png"),
				width=10,
				height=7,
				units="in",
					res=300,
					bg="white")
			print(p)
			dev.off()

			pdf(paste0(outdir,"/DiffStat.pdf"),
					width=10,
					height=7,
					bg="white", onefile=F)
			print(p)
			dev.off()
		}
    }
}

## 开始运行
argv <- commandArgs(TRUE)
if (length(argv) == 1) {
    main(argv[1])
}
