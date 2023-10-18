# Load necessary libraries while suppressing messages
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(ggalluvial))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))

# Pass parameters and define variables
argv <- commandArgs(TRUE)
taxAbdDir <- argv[1]
taxa <- argv[2]
num_taxa <- as.numeric(argv[3])
grouped <- argv[4]
mapDir <- argv[5]

# Read the mapping table
map <- read.table(mapDir, sep="\t", header = FALSE)

# Read and preprocess abundance table
taxAbd <- read.table(taxAbdDir, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
taxAbd <- taxAbd[, -ncol(taxAbd)]
colnames(taxAbd)[1] <- 'Taxon'

# Calculate abundance by group if grouped is TRUE
if(grouped==TRUE){
  # Reshape the data from wide to long
  melttaxAbd <- melt(taxAbd, variable.name = "sampleid", value.name = "abundance", id.vars = "Taxon")
  
  # Calculate abundance by group
  colnames(map)[4:5] <- c("Group", "sampleid")
  df1_group <- left_join(melttaxAbd, map, "sampleid")
  df2_group <- ddply(df1_group, .(Taxon, Group), summarize, mean = round(mean(abundance), 2))
  
  # Reshape the data from long to wide
  df3_group <- dcast(df2_group, Taxon ~ Group)
  taxAbd <- df3_group
}
rownames(taxAbd) <- taxAbd[,1]
taxAbd <- taxAbd[,-1]

# Remove rows containing specific terms
filter <- grep("Unclassified|Other|unidentified|uncultured|Incertae|NA", rownames(taxAbd))
taxAbd <- taxAbd[-filter, , drop = FALSE]

# Sort the data by total abundance
taxAbd$sum <- apply(taxAbd, 1, sum)
taxAbd <- taxAbd[order(taxAbd$sum, decreasing = TRUE),]  
taxAbd$sum <- NULL
taxAbd <- taxAbd[match(unique(rownames(taxAbd)), rownames(taxAbd)), , drop = FALSE]

# Calculate relative abundance
taxAbd_rel <- apply(taxAbd, 2, function(x) {x / sum(x)})
taxAbd_rel <- as.data.frame(taxAbd_rel)

# Select the top N taxa and group the rest as "Others"
if(nrow(taxAbd_rel) > num_taxa){ 
  taxAbd_rel <- taxAbd_rel[1:num_taxa, , drop = FALSE]
}
taxAbd_rel["Others",] <- 1 - colSums(taxAbd_rel)
taxAbd_rel <- taxAbd_rel[nrow(taxAbd_rel):1, , drop = FALSE]

# Prepare data for plotting
df1 <- data.frame(taxa_name = rownames(taxAbd_rel), taxAbd_rel)
df3 <- melt(df1, id = c("taxa_name"))
df3 <- df3[order(df3$value, df3$variable, decreasing = FALSE),]
df3$taxa_name <- factor(as.character(df3$taxa_name), levels = rownames(taxAbd_rel))

# Define colors for the plot
col_source <- rep(c(c(brewer.pal(8, "Accent"), brewer.pal(12, "Set3"))[-4][-9][-9][-11], "#B87333", "#8B00FF", "#CCB38C", "#FF2400"), t = 100)
col <- c(rev(col_source[1:num_taxa]), "grey30")

# Define text angle and positioning
num_x_text <- max(nchar(colnames(taxAbd_rel))
if (num_x_text > 2) {
  angle <- 45
  hjust <- 1
  vjust <- 1
} else {
  angle <- 0
  hjust <- 0.5
  vjust <- 1
}

# Create the plot
p <- ggplot(df3, aes(x = variable, stratum = taxa_name, alluvium = taxa_name, y = value * 100, fill = taxa_name)) +
  geom_stratum(alpha = 1, linetype = 0, width = 2/3) + 
  scale_color_manual(values = col) +
  ylab("Relative Abundance (%)") +
  scale_fill_manual(values = col, guide = guide_legend(nrow = min(35, nrow(taxAbd_rel)), title = paste("Top ", min(num_taxa, nrow(taxAbd_rel)), " ", taxa, sep="")) + 
  scale_x_discrete(expand = c(1/nrow(taxAbd_rel), 1/nrow(taxAbd_rel))) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12, angle = angle, hjust = hjust, vjust = vjust), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        legend.justification = c(0, 0.5), legend.position = "right",
        legend.text = element_text(face = c(rep("italic", t = nrow(taxAbd_rel)), "plain"), size = 8), legend.title = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.key.width = unit(1, 'cm'), legend.key.height = unit(0.5, 'cm'),
        plot.margin = unit(c(0.2, 0.2, 0, 0.2), "cm"))

# Create plots with and without trend
p1 <- p + geom_alluvium(alpha = 0.3, linetype = 0, width = 2/3)
p2 <- p + geom_alluvium(alpha = 0, linetype = 0, width = 2/3)

# Set output parameters for different file formats
pdf.width <- if(length(colnames(df)) <= 30) { 11 } else { (ncol(df) / 4) + 5 }
pdf.height <- 8
svg.width <- 11
svg.height <- 8
png.width <- 11
png.height <- 8

# Define output file names
if (grouped == FALSE) {
  pdf.trend.name <- paste0("bar_", taxa, ".pdf")
  pdf.notrend.name <- paste0("bar_", taxa, "_without_trend.pdf")
  svg.trend.name <- paste0("bar_", taxa, ".svg")
  svg.notrend.name <- paste0("bar_", taxa, "_without_trend.svg")
  png.trend.name <- paste0("bar_", taxa, ".png")
  png.notrend.name <- paste0("bar_", taxa, "_without_trend.png")
} else {
  pdf.trend.name <- paste0("bar_", taxa, "_grouped.pdf")
  pdf.notrend.name <- paste0("bar_", taxa, "_grouped_without_trend.pdf")
  svg.trend.name <- paste0("bar_", taxa, "_grouped.svg")
  svg.notrend.name <- paste0("bar_", taxa, "_grouped_without_trend.svg")
  png.trend.name <- paste0("bar_", taxa, "_grouped.png")
  png.notrend.name <- paste0("bar_", taxa, "_grouped_without_trend.png")
}

# Save plots to different file formats
pdf(pdf.trend.name, width = pdf.width, height = pdf.height, bg = "white")
print(p1)
dev.off()

svg(svg.trend.name, width = svg.width, height = svg.height, bg = "white", onefile = FALSE)
print(p1)
dev.off()

png(png.trend.name, units = "in", res = 900, width = png.width, height = png.height, bg = "white")
print(p1)
dev.off()

pdf(pdf.notrend.name, width = pdf.width, height = pdf.height, bg = "white")
print(p2)
dev.off()

svg(svg.notrend.name, width = svg.width, height = svg.height, bg = "white", onefile = FALSE)
print(p2)
dev.off()

png(png.notrend.name, units = "in", res = 900, width = png.width, height = png.height, bg = "white")
print(p2)
dev.off()
