##packages_needed
# List of R packages needed for this script.
#RColorBrewer; ggplot2; reshape2; ggalluvial; plyr; dplyr; vegan; ggrepel; data.table; magrittr; pheatmap; ape; igraph; rjson; tidytable; tibble; broom; grid; ggpubr; cowplot; PMCMRplus; ggsignif; patchwork; car; purrr; foreach; doMC; iterators; doParallel; parallel; Matrix; bigmemory; biganalytics; gRbase; gplots; gridExtra; phyloseq; qgraph; stringr; rgexf; phangorn; permute; pROC


## Load the R environment
# Activate the specified R environment using Conda.
source /share/home/LiuChang/software/biosoftware/miniconda3/bin/activate /share/home/LiuChang/project/MetaDenovo_V1/envs/miniconda3/envs/R3.6.1
## Script and data paths
# Define paths to the script and data directories.
bin_path=/share/home/LiuChang/MetaTest/script/
data_path=/share/home/LiuChang/MetaTest_bk/1_rawdata

## Taxa and functional composition bar plots
# Create directories for taxa and functional composition bar plots and run R scripts for generating them.
mkdir -p taxa_bar func_bar && cd taxa_bar
Rscript ${bin_path}/bar.Taxa.R ${data_path}/taxa_table_species_RC.tsv species 20 TRUE ${data_path}/map.txt
cd ../func_bar
Rscript ${bin_path}/bar.Fun.R ${data_path}/KEGG_L1_RC.tsv kegg_L1 FALSE ${data_path}/map.txt
cd ..

## Taxa and functional composition PCA
mkdir -p taxa_PCA func_PCA && cd taxa_PCA
Rscript ${bin_path}/plot_taxa_pca_names_graph.R ${data_path}/taxa_table_species_RC.tsv ${data_path}/map.txt species
cd ../func_PCA
Rscript ${bin_path}/plot_function_pca_names_graph.R ${data_path}/KEGG_L1_RC.tsv ${data_path}/map.txt kegg_L1
cd ..

## Taxa and functional composition PCoA
mkdir -p taxa_PCoA func_PCoA && cd taxa_PCoA
Rscript ${bin_path}/plot_PCoA_names.R ${data_path}/taxa_table_species_RC.tsv bray ${data_path}/map.txt
cd ../func_PCoA
Rscript ${bin_path}/plot_function_PCoA_names.R ${data_path}/KEGG_L1_RC.tsv bray ${data_path}/map.txt
cd ..

## Taxa and functional composition NMDS
mkdir -p taxa_NMDS func_NMDS && cd taxa_NMDS
Rscript ${bin_path}/plot_nmds_taxa_names.R ${data_path}/taxa_table_species_RC.tsv ${data_path}/map.txt bray species
cd ../func_NMDS
Rscript ${bin_path}/plot_nmds_function_names.R ${data_path}/KEGG_L1_RC.tsv ${data_path}/map.txt bray kegg_L1
cd ..

## Taxa and functional composition UPGMA clustering
mkdir -p taxa_cluster func_cluster && cd taxa_cluster
Rscript ${bin_path}/plot_cluster_taxa_tree.R ${data_path}/taxa_table_species_RC.tsv bray species
cd ../func_cluster
Rscript ${bin_path}/plot_cluster_function_tree.R ${data_path}/KEGG_L1_RC.tsv bray kegg_L1
cd ..

## Differential analysis box plots
mkdir -p taxa_diff func_diff && cd taxa_diff
Rscript ${bin_path}/taxDiffStat.R ${bin_path}/taxDiffStat_genes.conf
cd ../func_diff
Rscript ${bin_path}/funcDiffStat.R ${bin_path}/funcDiffStat_default.conf
cd ..

## Differential analysis heatmaps
# Abundance.tsv is located in the directory of the previous differential analysis.
mkdir -p taxa_diff_heatmap func_diff_heatmap && cd taxa_diff_heatmap
Rscript ${bin_path}/heatmap.R ../taxa_diff/Abundance.tsv ${data_path}/map.txt 10
mkdir -p taxa_diff_heatmap func_diff_heatmap && cd taxa_diff_heatmap
Rscript ${bin_path}/heatmap.R ../taxa_diff/Abundance.tsv ${data_path}/map.txt 10
cd ../func_diff_heatmap
Rscript ${bin_path}/heatmap.R ../func_diff/Abundance.tsv ${data_path}/map.txt 10
cd ..

## Environmental factor RDA (Redundancy Analysis)
mkdir -p RDA && cd RDA
Rscript ${bin_path}/RDA.R ${data_path}/taxa_table_species_RC.tsv ${data_path}/env.txt ${data_path}/map.txt species
cd ..

## Environmental factor envfit
mkdir -p envfit && cd envfit
Rscript ${bin_path}/envfit.R ${data_path}/taxa_table_species_RC.tsv ${data_path}/env.txt ${data_path}/map.txt
cd ..

## Machine learning (Random Forest)
mkdir -p random_forest && cd random_forest
bash ${bin_path}/random_forest_ROC_genus.sh

## Correlation network analysis
mkdir -p network_sparCC
Rscript ${bin_path}/prepare.sparCC_data_lhy.R  ${data_path}/taxa_table_species_RC.tsv ${data_path}/map.txt network_sparCC/ network_sparCC/ 8
Rscript ${bin_path}/make_network.R ${data_path}/taxa_table_species_RC.tsv ${data_path}/map.txt ${data_path}/taxa_table_species_RC.tsv network_sparCC network_sparCC
