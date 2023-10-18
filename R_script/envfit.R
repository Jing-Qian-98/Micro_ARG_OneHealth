# Environmental fitting analysis: After calculating the distance matrix, we perform regression analysis using the 'envfit' function to fit environmental factors (observed variables) to ordination (unconstrained axes).

# Load necessary libraries
library(RColorBrewer)
library(vegan)
library(ggplot2)
library(ggrepel)

# Parse command line arguments and define variables
argv <- commandArgs(TRUE)
in_taxa <- argv[1]
env_data <- argv[2]
map_file <- argv[3]

# Read species abundance table and preprocess
taxaAbd <- read.table(in_taxa, sep = "\t", header = TRUE, row.names = 1)
taxaAbd <- taxaAbd[, -ncol(taxaAbd)]

# Convert to relative abundance and transpose rows and columns
taxaAbd_rel <- apply(taxaAbd, 2, function(x) { x / sum(x) })
phylum <- t(taxaAbd_rel)

# Read group file
map <- read.table(map_file, row.names = 1, header = FALSE, sep = "\t")
group <- as.factor(map$V4)

# Read environmental factors
env <- read.delim(env_data, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

# Perform MDS (Multidimensional Scaling) analysis
mds <- cmdscale(vegdist(phylum, "bray"), k = 2, eig = TRUE)
mds_point <- data.frame(mds$points)   # Get coordinates of each sample
colnames(mds_point) <- c('X1', 'X2')
eig <- mds$eig

# Define colors for plotting
color <- c(brewer.pal(nrow(map), "Set1"))

# Create a scatter plot with ellipses
ggplot(mds_point, aes(x = X1, y = X2, color = group)) +
  geom_point(aes(color = group), size = 4, alpha = 0.6) +
  stat_ellipse(aes(x = X1, y = X2, fill = group), geom = "polygon", alpha = 1/2, levels = 0.95) +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) 

# Fit environmental vectors to the ordination
fit <- envfit(mds, env, permutations = 999)
fit_val <- scores(fit, display = c("vectors"))
fit_val <- fit_val * vegan::ordiArrowMul(fit_val, fill = 1.5)

# Output the results of the 'envfit' analysis
fit$vectors

# Add vectors to the plot using ggplot
p = ggplot(mds_point, aes(x = X1, y = X2, color = group)) +
  geom_point(size = 4, alpha = 0.6) +
  stat_ellipse(aes(x = X1, y = X2, fill = group), geom = "polygon", alpha = 1/2) +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) +
  geom_segment(data = data.frame(fit_val), 
               aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = 'black', alpha = 1)  + 
  geom_label_repel(data = data.frame(fit_val), aes(Dim1, Dim2, label = rownames(fit_val)),
                   color = 'black', alpha = 1,
                   segment.color = 'grey35',
                   point.padding = unit(0.1, "lines")) +
  labs(x = paste("CAP 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)", sep = ""), 
       y = paste("CAP 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)", sep = "")) 

# Apply visual enhancements such as the 'bw' theme and remove styling
p = p + theme_bw() + theme(panel.grid.major = element_line(colour = NA))

# Output the plot as a PDF file
pdf(paste0("envfit.pdf"), width = 12, height = 10, bg = "white", onefile = FALSE)
print(p)
dev.off()
