# Suppressing library messages to avoid clutter
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(ggalluvial))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))

# Passing arguments and defining variables
argv <- commandArgs(TRUE)
funcAbdDir <- argv[1]
category <- argv[2]
grouped <- argv[3]
mapDir <- argv[4]

# Reading the grouping table
map <- read.table(mapDir, sep="\t", header = F)

# Reading functional abundance data and preprocessing
funcAbd <- read.table(funcAbdDir, skip = 0, sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
colnames(funcAbd)[1] <- 'Name'

# If the last column is a string, remove it
if (!is.numeric(funcAbd[, ncol(funcAbd)])) {
  funcAbd <- funcAbd[-ncol(funcAbd)]
}

# Calculating abundance by group
if (grouped == TRUE) {
  # Reshaping data from wide to long format
  meltfuncAbd <- melt(funcAbd, variable.name = "sampleid", value.name = "abundance", id.vars = "Name")

  # Calculating abundance by group
  colnames(map)[4:5] <- c("Group", "sampleid")
  funcAbd_group <- left_join(meltfuncAbd, map, "sampleid")
  funcAbd_mean_group <- ddply(funcAbd_group, .(Name, Group), summarize, mean = round(mean(abundance), 2))

  # Reshaping data from long to wide format
  df3_group <- dcast(funcAbd_mean_group, Name ~ Group)
  funcAbd <- df3_group
}
rownames(funcAbd) <- funcAbd[, 1]
funcAbd <- funcAbd[, -1]

# Removing rows with specific strings
# filter <- grep("Unclassified|Other|unidentified|uncultured|Incertae|NA", rownames(funcAbd))
# funcAbd <- funcAbd[-filter, , drop = F]

# Sorting the data
funcAbd$sum <- apply(funcAbd, 1, sum)
funcAbd <- funcAbd[order(funcAbd$sum, decreasing = T), ]
funcAbd$sum <- NULL
funcAbd <- funcAbd[match(unique(rownames(funcAbd)), rownames(funcAbd)), , drop = F]

# Converting to relative abundance
funcAbd_rel <- as.data.frame(lapply(funcAbd, function(x) x / sum(x))
row.names(funcAbd_rel) <- row.names(funcAbd)

# Preparing for the plot
#num_Category <- ""
df1 = data.frame(Category = rownames(funcAbd_rel), funcAbd_rel)
df3 = melt(df1, id = c("Category"))
df3 = df3[order(df3$value, df3$variable, decreasing = F), ]
df3$Category = factor(as.character(df3$Category), levels = rownames(funcAbd_rel))
#df3 = df3[which(df3$Category != "Others"), ]
col_source = rep(c(c(brewer.pal(8, "Accent"), brewer.pal(12, "Set3"))[-4][-9][-9][-11], "#B87333", "#8B00FF", "#CCB38C", "#FF2400"), t = 100)
col = c(rev(col_source[1:1000]), "grey30")
num_x_text = max(nchar(colnames(funcAbd_rel))
if (num_x_text > 2) {
  angle = 45
  hjust = 1
  vjust = 1
} else {
  angle = 0
  hjust = 0.5
  vjust = 1
}

# Creating the plot
p <- ggplot(df3, aes(x = variable, stratum = Category, alluvium = Category, y = value * 100, fill = Category)) +
  geom_stratum(alpha = 1, linetype = 0, width = 2/3) +
  scale_color_manual(values = col) +
  ylab("Relative Abundance (%)") +
  scale_fill_manual(values = col, guide = guide_legend(nrow = min(35, nrow(funcAbd_rel)), title = paste("Category"))) +
  scale_x_discrete(expand = c(1/nrow(funcAbd_rel), 1/nrow(funcAbd_rel))) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12, angle = angle, hjust = hjust, vjust = vjust), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        legend.justification = c(0, 0.5), legend.position = "right",
        legend.text = element_text(face = c(rep("italic", t = nrow(funcAbd_rel)), "plain"), size = 8), legend.title = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.key.width = unit(1, 'cm'), legend.key.height = unit(0.5, 'cm'),
        plot.margin = unit(c(0.2, 0.2, 0, 0.2), "cm"))

## With trend
p1 <- p + geom_alluvium(alpha = 0.3, linetype = 0, width = 2/3)

## Without trend
p2 <- p + geom_alluvium(alpha = 0, linetype = 0, width = 2/3)

# Setting output parameters
## Width and height for PDF
if (length(colnames(df3)) <= 30) {
  pdf.width <- 11
} else {
  pdf.width <- (ncol(df3) / 4) + 5
}
pdf.height <- 8

# Width and height for SVG
svg.width <- 11
svg.height <- 8

# Width and height for PNG
png.width <- 11
png.height <- 8

## Output file names
if (grouped == FALSE) {
  pdf.trend.name <- paste0("bar_", category, ".pdf")
  pdf.notrend.name <- paste0("bar_", category, "_without_trend.pdf")
  
  svg.trend.name <- paste0("bar_", category, ".svg")
  svg.notrend.name <- paste0("bar_", category, "_without_trend.svg")
  
  png.trend.name <- paste0("bar_", category, ".png")
  png.notrend.name <- paste0("bar_", category, "_without_trend.png")
} else {
  pdf.trend.name <- paste0("bar_", category, "_grouped.pdf")
  pdf.notrend.name <- paste0("bar_", category, "_grouped_without_trend.pdf")
  
  svg.trend.name <- paste0("bar_", category, "_grouped.svg")
  svg.notrend.name <- paste0("bar_", category, "_grouped_without_trend.svg")
  
  png.trend.name <- paste0("bar_", category, "_grouped.png")
  png.notrend.name <- paste0("bar_", category, "_grouped_without_trend.png")
}

# Output with trend
pdf(pdf.trend.name, width = pdf.width, height = pdf.height, bg = "white")
print(p1)
dev.off()

# Output in SVG format
svg(svg.trend.name, width = svg.width, height = svg.height, bg = "white", onefile = F)
print(p1)
dev.off()

# Output in PNG format
png(png.trend.name, units = "in", res = 900, width = png.width, height = png.height, bg = "white")
print(p1)
dev.off()

# Output without trend
pdf(pdf.notrend.name, width = pdf.width, height = pdf.height, bg = "white")
print(p2)
dev.off()

# Output in SVG format
svg(svg.notrend.name, width = svg.width, height = svg.height, bg = "white", onefile = F)
print(p2)
dev.off()

# Output in PNG format
png(png.notrend.name, units = "in", res =900, width = png.width, height = png.height, bg = "white")
print(p2)
dev.off()
