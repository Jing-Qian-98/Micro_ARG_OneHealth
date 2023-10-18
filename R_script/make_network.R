.libPaths("/PERSONALBIO/work/microbio/m13/bin/env/PERSONALBIO/Work/Mb/Mb07/.conda/envs/picrust2/lib/R/library/")
setwd("./")
################################################################################
# Load Packages
suppressMessages(library("foreach"))
suppressMessages(library("doMC"))
suppressMessages(library("iterators"))
suppressMessages(library("doParallel"))
suppressMessages(library("parallel"))
suppressMessages(library("Matrix"))
suppressMessages(library("bigmemory"))
suppressMessages(library("biganalytics"))
suppressMessages(library("gRbase"))
suppressMessages(library("gplots"))
suppressMessages(library("ggplot2"))
suppressMessages(library("grid"))
suppressMessages(library("gridExtra"))
suppressMessages(library("data.table"))
suppressMessages(library("plyr"))
suppressMessages(library("ape"))
suppressMessages(library("phyloseq"))
suppressMessages(library("vegan"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("igraph"))
suppressMessages(library("qgraph"))
suppressMessages(library("stringr"))
################################################################################
################################################################################

argv <- commandArgs(TRUE)

  if (length(argv) !=5 ){
    stop('The script needs  files:otu_table; sample_metadata; result_outputfile; data_inputfile ;core_number, please put them in a dir  ')
  }
 
#argv=c("data/16s.uparse_no_pynast_failures_even40000.csv", "data/sample_data_Treat1.tsv", "data/taxonomy.txt", "result", "data")

################################################################################
################################################################################

#Preallocate global data structures 
PARAM <- list();
PARAM$folder.data <- argv[5]
PARAM$folder.output <- argv[4]              #PUT RESULTS/OUTPUT FOLDER HERE
PARAM$file.otu_data <- argv[3]
PARAM$file.sample_metadata <- argv[2]
PARAM$file.otu_table <- argv[1]
PARAM$use.cores <- 1

###################################
commond=paste("mkdir -p ",PARAM$folder.output ,sep="") 
system(commond)
###################################
###################################
#Minimum OTU size required across all samples
PARAM$thresh.otu_size <- 10;
#Minimum sample size required
PARAM$thresh.sample_size <- 1000;
#Minimum relative number of sequences retained
PARAM$thresh.retained_seqs <- 0.5;
####################################


###################################
#Load data
load(file=paste(PARAM$folder.data, "/All_samples.otu_table.RData", sep=""));
load(file=paste(PARAM$folder.output, "/All_samples.SparCC_global.RData", sep=""));
#Read sample raw metadata
sample.data.raw <- read.table(PARAM$file.sample_metadata, header=T, comment.char = "@", sep="\t", colClasses = "character");
colnames(sample.data.raw)[4] <- "Group"
rownames(sample.data.raw) <-as.character(sample.data.raw[,1]);
idx=rownames(sample.data.raw) %in% colnames(ot.2)
sample.data=sample.data.raw[idx,]
samples=rownames(sample.data)
otu.data=read.table(PARAM$file.otu_data, header=T, sep="\t");

##转换otu数据格式
df <- data.frame(otu.data[,1],otu.data[,ncol(otu.data)])
df_split <- str_split_fixed(df[,2], ";", 7)
df_merge <- cbind(df[,1], df_split)
df_merge[,1] <- df_merge[,8]
colnames(df_merge) <- c("OTU", "domain", "phylum", "class", "order", "family", "genus", "species")
otu.data <- as.data.frame(df_merge)
####
rownames(otu.data) <-as.character(otu.data$OTU);
sample.data$Group=factor(as.character(sample.data$Group),levels=unique(sample.data$Group))


############################
############################
#Plot larger network with context
############################
network.samples=samples
#Get current pruned OTU table
curr.ot <- ot.2[which(rowSums(ot.2[,network.samples]) > 0 ), network.samples];
curr.otus <- rownames(curr.ot);

############################
#Generate network plots
############################
#Get SparCC cooc between current OTUs
curr.sparcc <- RMT.sparcc[curr.otus, curr.otus];
curr.sparcc2 <- curr.sparcc
curr.sparcc2[curr.sparcc2 < 0] <- 0;
#Coerce into igraph object
# sprcc结果转换为igraph画图处理的结果
curr.graph <- graph.adjacency(curr.sparcc, mode="upper", weighted=T, diag=F);
curr.graph2 <- graph.adjacency(curr.sparcc2, mode="upper", weighted=T, diag=F);
remove.vertices <- which(degree(curr.graph) == 0);
keep.vertices <- which(degree(curr.graph) >= 1);
curr.graph <- delete.vertices(curr.graph, remove.vertices);
curr.graph2 <- delete.vertices(curr.graph2, remove.vertices);
E(curr.graph)$sim=E(curr.graph)$weight
E(curr.graph)$weight=abs(E(curr.graph)$weight)
curr.ot.rel.raw=data.frame(dgCMatrix2matrix(ot.rel.raw[curr.otus[keep.vertices],samples]))
curr.ot=data.frame(dgCMatrix2matrix(ot[curr.otus[keep.vertices],samples]))

####################################
####################################
#for ggraph 
####################################
####################################
curr.graph.dat <- list();
curr.graph.dat$size <- rowSums(ot[curr.otus[keep.vertices],samples]);

curr.graph.dat$node <- rownames(ot[curr.otus[keep.vertices],samples]);

#get first 9 phylum
temp=aggregate(curr.ot.rel.raw, by=list(otu.data[rownames(curr.ot.rel.raw),"phylum"]), sum)
rownames(temp)=temp$Group.1
temp$Group.1=NULL
taxa10=names(sort(rowSums(temp),decreasing=T))[1:10]

#Color by taxonomy
curr.graph.dat$color.tax <- rep.int("#bdbdbd", length(curr.otus[keep.vertices]));
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[1]] <- "#8DD3C7";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[2]] <- "#FFFFB3";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[3]] <- "#BEBADA";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[4]] <- "#FB8072";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[5]] <- "#80B1D3"
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[6]] <- "#FDB462";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[7]] <-"#B3DE69"
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[8]] <-"#BC80BD"
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[9]] <-"#CCEBC5"
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[10]] <-"#FFED6F"


col_p=colorRampPalette(brewer.pal(9, "Reds"),bias=0.3,alpha=0.3)(100)
col_n=rev(colorRampPalette(brewer.pal(9, "Greens"),bias=0.3,alpha=0.3)(100))
col_edge=c(col_n,col_p)


############################
############################
#Compute a layout

E(curr.graph)$curved <- 0.25
E(curr.graph)$color <- 0
E(curr.graph)$color=col_edge[100+signif(E(curr.graph)$sim*100,1)]

#community_split using co-occurrence network
com = multilevel.community(curr.graph2, weights = E(curr.graph2)$sim)
subgroup = split(com$names, com$membership)
## subgroup
V(curr.graph)$sg = paste("module",com$membership,sep="_")
modu=names(sort(summary(as.factor(V(curr.graph)$sg),maxsum=10000),decreasing=T))
V(curr.graph)$color="grey"
for (i in 1:10) {
V(curr.graph)$color[which(V(curr.graph)$sg==modu[i])] = c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[c(3,5:7,9,12:13,15,20,19)][i]
}

####
node.id <- c(1:length(V(curr.graph)$name))
names(node.id) <- V(curr.graph)$name

edge <- get.edgelist(curr.graph)
e <- cbind(node.id[edge[,1]],node.id[edge[,2]])
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(curr.graph),area=10*(vcount(curr.graph)^2),repulse.rad=(vcount(curr.graph)^3.1))

############################
############################
# filled by modules
pdf(paste(PARAM$folder.output, "/network_module.pdf",sep=""),width=10,height=10,bg="white")
plot(curr.graph, layout=l, 
    fade=TRUE,
    bg="white",
    trans=TRUE,
    vertex.label=NA,
    vertex.size=log2(curr.graph.dat$size*1000000/length(samples))/6, 
    vertex.color=V(curr.graph)$color,
    vertex.label.cex=0.5,
    edge.arrow.size=0, 
    edge.color=E(curr.graph)$color, 
    edge.width=E(curr.graph)$weight*4-min(E(curr.graph)$weight)*4+0.001, 
    edge.alpha=0.6,
#    vertex.label.family=NA, 
#    vertex.label.color="black",
#    vertex.label.degree=pi/2,
#    vertex.label.dist=rep(0,length(V(curr.graph)$name)),
#    vertex.pie=pie.values, 
#    vertex.pie.border=list(rep("white"),length(colorg)),
#    vertex.pie.color=list(colorg), 
    vertex.frame.color=NA)
with(data.frame(V(curr.graph)$sg,V(curr.graph)$color),legend("topright",title="Modules", c(as.character(modu[1:min(10,length(modu))])), inset=c(-0.05,0),
    title.adj = 0,
    pt.bg=c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[c(3,5:7,9,12:13,15,20,19)][1:min(10,length(modu))],bty = "n",col ="black", pch = 22, cex=1.0))
with(data.frame(V(curr.graph)$sg,V(curr.graph)$color),legend("bottomright",title="Similarity", c("-1.0","-0.9","-0.8","-0.7","0.7","0.8","0.9","1.0"), 
    title.adj = 0,
    col=col_edge[100+c(-0.99,-0.9,-0.8,-0.7,0.7,0.8,0.9,0.99)*100],bty = "n", pch = "--", cex=1.0,pt.cex=2.5,pt.lwd=2))
dev.off();

############################
# filled by Phylum 
pdf(paste(PARAM$folder.output, "/network_phylum.pdf",sep=""),width=10,height=10,bg="white")
plot(curr.graph, layout=l, 
    fade=TRUE,
    bg="white",
    trans=TRUE,
    vertex.label=NA,
    vertex.size=log2(curr.graph.dat$size*1000000/length(samples))/6, 
    vertex.color=curr.graph.dat$color.tax,
    vertex.label.cex=0.5,
    edge.arrow.size=0, 
    edge.color=E(curr.graph)$color, 
    edge.width=E(curr.graph)$weight*4-min(E(curr.graph)$weight)*4+0.001, 
    edge.alpha=0.6,
    vertex.frame.color=NA)
with(data.frame(V(curr.graph)$sg,V(curr.graph)$color),legend("topright",title="Phylum", c(as.character(taxa10[!is.na(taxa10)])), inset=c(-0.05, 0), 
    title.adj = 0,
    pt.bg=c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#BC80BD", "#CCEBC5", "#FFED6F")[1:min(10,length(taxa10[!is.na(taxa10)]))], 
    bty = "n",col ="black", pch = 22, cex=.7))
with(data.frame(V(curr.graph)$sg,V(curr.graph)$color),legend("bottomright",title="Similarity", c("-1.0","-0.9","-0.8","-0.7","0.7","0.8","0.9","1.0"), 
    title.adj = 0,
    col=col_edge[100+c(-0.99,-0.9,-0.8,-0.7,0.7,0.8,0.9,0.99)*100],bty = "n", pch = "--", cex=1.0,pt.cex=2.5,pt.lwd=2))
dev.off();

############################
# filled by Groups
if(length(levels(sample.data$Group)) <= 10){

pie_size <- aggregate(t(curr.ot),list(sample.data$Group), mean) # using mean values as the sizes of the pie slices
Groups <- pie_size[,1]
pie_size <- pie_size[,V(curr.graph)$name]

col_source=c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[c(3,5:7,9,12:13,15,20,19)]
col_pie=rep(list(col_source[1:length(levels(Groups))]),vcount(curr.graph))

pdf(paste(PARAM$folder.output, "/network_group.pdf",sep=""),width=10,height=10,bg="white")
plot(curr.graph, layout=l, 
    fade=TRUE,
    bg="white",
    trans=TRUE,
    vertex.label=NA,
    vertex.size=log2(curr.graph.dat$size*1000000/length(samples))/7, 
    vertex.shape="pie", 
    edge.arrow.size=0, 
    edge.color=E(curr.graph)$color, 
    edge.width=E(curr.graph)$weight*4-min(E(curr.graph)$weight)*4+0.001, 
    edge.alpha=0.6,
    vertex.pie=pie_size, 
    vertex.pie.color=col_pie, 
    vertex.frame.color=NA)
with(data.frame(V(curr.graph)$sg,V(curr.graph)$color),legend("topright",title="Group", 
    title.adj = 0,
    c(as.character(Groups)), pt.bg=col_source[1:length(levels(Groups))],bty = "n",col ="black", pch = 22, cex=1.0))
with(data.frame(V(curr.graph)$sg,V(curr.graph)$color),legend("bottomright",title="Similarity", c("-1.0","-0.9","-0.8","-0.7","0.7","0.8","0.9","1.0"), 
    title.adj = 0,
    col=col_edge[100+c(-0.99,-0.9,-0.8,-0.7,0.7,0.8,0.9,0.99)*100],bty = "n", pch = "--", cex=1.0,pt.cex=2.5,pt.lwd=2))
dev.off();

}

################################################################
################################################################
#for topology
################################################################
################################################################
relations_OTU <- E(curr.graph,directed=F)

nodes_OTU <- V(curr.graph)

ER_OTU=sample_gnm(length(nodes_OTU),length(relations_OTU),directed=F,loops=F )

ER_OTU.list <- list()
for(i in 1:9){
    ER_OTU.list[[i]] <- sample_gnm(length(nodes_OTU),length(relations_OTU),directed=F,loops=F )
}
ER_OTU.list[[10]] <- ER_OTU

#3 Degree
Degree_OTU=degree(curr.graph, loops =F, normalized = F)

Degree_ER_OTU=degree(ER_OTU, loops =F, normalized = F)

#4 Transitivity
Transitivity_OTU= transitivity(curr.graph, type="global")

Transitivity_ER_OTU= transitivity(ER_OTU, type="global")

mean_distance_OTU=mean_distance(curr.graph,directed = F,unconnected=T)

mean_distance_ER_OTU=mean_distance(ER_OTU,directed = F,unconnected=T)



###plot_degree_distribution
count_OTU=rep(1,t=length(Degree_OTU))
count_ER_OTU=rep(1,t=length(Degree_ER_OTU))
data_OTU=aggregate(x=count_OTU,by=list(Degree_OTU),"sum")
data_ER_OTU=aggregate(x=count_ER_OTU,by=list(Degree_ER_OTU),"sum")
data_OTU$Group="Empirical network"
data_ER_OTU$Group="Random network"
data_degree=rbind(data_OTU,data_ER_OTU)
pdf(paste(PARAM$folder.output, "/network_degree_distribution.pdf",sep=""),width=10,height=10,bg="white")
ggplot(data_degree)+geom_point(aes(x=Group.1,y=x,group=Group,col=Group,fill=Group))+
  theme_bw() +theme(axis.ticks=element_line(size = 1,colour = "black"),axis.ticks.length = unit(-2,units="mm"),axis.text.y = element_text(margin=unit(c(0,3,0,0),units="mm"),color="black",size=14),
                    axis.text.x = element_text(margin=unit(c(3,0,0,0),units="mm"),color ="black",size=14),axis.title= element_text(size=14))+
  theme(legend.text = element_text(size = 14),legend.title = element_blank(),legend.key.size = unit(12,"mm"),legend.position = c(0.5,0.7))+
  labs(y ="Degrees",x="node counts")+
  theme(plot.background=element_rect(fill = "transparent",colour = NA),panel.background=element_rect(fill = "transparent",colour = NA),plot.margin = unit(c(0, 0, 0, 3), "lines"),panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA))
dev.off()

######################################for each sample
vids_OTU=graph_OTU=graph_OTU2=list()
Betw_centrality_OTU= Betw_centrality_nd_OTU=list()
Betw_edge_centrality_OTU=Betw_edge_centrality_nd_OTU=Betw_edge_centrality_nt=OTU=list()
Clos_centrality_OTU=Clos_centrality_nd_OTU=list()
Degree_centrality_OTU=Degree_centrality_nd_OTU=list()
Transitivity_nd_OTU=list()
knn_nd_OTU=list()


Betw_centrality_nt_OTU=c()
Clos_centrality_nt_OTU=c()
Degree_centrality_nt_OTU=c()
Transitivity_nt_OTU=c()
knn_nt_OTU=c()
Aver_path_nt_OTU=c()
Degree_assortativity_nt_OTU=c()
Density_nt_OTU=c()
cluster_num_nt_OTU=c()
Diameter_OTU=c()
edge_number_OTU=c()
vertice_number_OTU=c()
modularity_nt_OTU=c()


ot=ot[,sort(colnames(ot))]
for (i in 1:length(samples)){
  vids_OTU[[i]]=ot[names(nodes_OTU),sort(samples)[i]]
  for (j in 1:length(vids_OTU[[i]])) {
    if (vids_OTU[[i]][j]!=0) vids_OTU[[i]][j]=j
  }
  vids_OTU[[i]]=vids_OTU[[i]][vids_OTU[[i]]!=0]
  graph_OTU[[i]]=induced_subgraph(curr.graph, vids_OTU[[i]],impl = c("copy_and_delete"))
  remove.vertices_temp <- which(degree(graph_OTU[[i]]) == 0);
  graph_OTU[[i]] <- delete.vertices(graph_OTU[[i]], remove.vertices_temp);
  
  graph_OTU2[[i]]=induced_subgraph(curr.graph2, vids_OTU[[i]],impl = c("copy_and_delete"))
  remove.vertices_temp <- which(degree(graph_OTU2[[i]]) == 0);
  graph_OTU2[[i]] <- delete.vertices(graph_OTU2[[i]], remove.vertices_temp);
  
  #Topological features 
  #1 Betweenness centrality 
  Betw_centrality_OTU[[i]]=centr_betw(graph_OTU[[i]], directed = FALSE,normalized=F)
  Betw_centrality_nd_OTU[[i]]=Betw_centrality_OTU[[i]]$res
  Betw_centrality_nt_OTU[i]=Betw_centrality_OTU[[i]]$centralization
  Betw_edge_centrality_OTU[[i]]=edge_betweenness(graph_OTU[[i]], directed = FALSE)
  
  #2 Closeness centrality 
  Clos_centrality_OTU[[i]]=centr_clo(graph_OTU[[i]], normalized = F)
  Clos_centrality_nd_OTU[[i]]=Clos_centrality_OTU[[i]]$res
  Clos_centrality_nt_OTU[i]=Clos_centrality_OTU[[i]]$centralization
  
  #3 Degree 
  Degree_centrality_OTU[[i]]=centr_degree(graph_OTU[[i]], loops = F, normalized = F)
  Degree_centrality_nd_OTU[[i]]=Degree_centrality_OTU[[i]]$res
  Degree_centrality_nt_OTU[i]=Degree_centrality_OTU[[i]]$centralization
  
  #4 Transitivity 
  Transitivity_nd_OTU[[i]] = transitivity(graph_OTU[[i]], type="local")
  Transitivity_nt_OTU[i] = transitivity(graph_OTU[[i]], type="global")
  
  #5 Average nearest neighbor degree 
  knn_nd_OTU[[i]]=knn(graph_OTU[[i]])$knn
  knn_nt_OTU[i]=mean(knn_nd_OTU[[i]],na.rm=T)
  
  #6 Average path 
  Aver_path_nt_OTU[i]=mean_distance(graph_OTU[[i]], directed = F, unconnected = T)
  
  #7 Degree assortativity 
  Degree_assortativity_nt_OTU[i]=assortativity_degree(graph_OTU[[i]]) 
  
  #8 density 
  Density_nt_OTU[i]=edge_density(graph_OTU[[i]], loops=FALSE)
  
  #9 Cluster number 
   cluster_num_nt_OTU[i] <- components(graph_OTU[[i]])$no
  
  #10 Diameter
  Diameter_OTU[i]=diameter(graph_OTU[[i]])
  
  #11 edges number 
  edge_number_OTU[i] <- ecount(graph_OTU[[i]])
  
  #12 vertices number 
  vertice_number_OTU[i] <- vcount(graph_OTU[[i]])  
  
  #13 modularity 
  modularity_nt_OTU[i] <- modularity(multilevel.community(graph_OTU2[[i]], weights = E(graph_OTU2[[i]])$weight))
  print(i)
}


network_topo_table=data.frame(samples=sort(samples),
                              Average_nearest_neighbor_degree=knn_nt_OTU,
                              Average_path_length=Aver_path_nt_OTU,
                              Betweenness_centrality=Betw_centrality_nt_OTU,
                              Closeness_centrality=Clos_centrality_nt_OTU,
                              Degree_assortativity=Degree_assortativity_nt_OTU,
                              Degree_centralization=Degree_centrality_nt_OTU,
                              Density=Density_nt_OTU,
                              Cluster_num=cluster_num_nt_OTU,
                              Diameter=Diameter_OTU,
                              Transitivity=Transitivity_nt_OTU,
                              Num_vertice=vertice_number_OTU,
                              Num_edge=edge_number_OTU,
                              Modularity=modularity_nt_OTU)
write.table(network_topo_table,paste(PARAM$folder.output, "/subnetwork_topo_list_per_sample.xls",sep=""),sep="\t",row.names=F,quote = F)
########################################## for all 
  #1 Betweenness centrality 
  Betw_centrality=centr_betw(curr.graph, directed = FALSE,normalized=F)
  Betw_centrality_nd=Betw_centrality$res
  Betw_centrality_nt=Betw_centrality$centralization
  Betw_edge_centrality=edge_betweenness(curr.graph, directed = FALSE)
  BET=betweenness(curr.graph,directed=F,weights=E(curr.graph)$weight)
  
  #2 Closeness centrality 
  Clos_centrality=centr_clo(curr.graph, normalized = F)
  Clos_centrality_nd=Clos_centrality$res
  Clos_centrality_nt=Clos_centrality$centralization
  
  #3 Degree 
  Degree_centrality=centr_degree(curr.graph, loops = F, normalized = F)
  Degree_centrality_nd=Degree_centrality$res
  Degree_centrality_nt=Degree_centrality$centralization
  DEG=degree(curr.graph, loops =F, normalized = F)
  
  #4 Transitivity 
  Transitivity_nd = transitivity(curr.graph, type="local")
  Transitivity_nt = transitivity(curr.graph, type="global")
  
  #5 Average nearest neighbor degree 
  knn_nd=knn(curr.graph)$knn
  knn_nt=mean(knn_nd,na.rm=T)
  
  #6 Average path 
  Aver_path_nt=mean_distance(curr.graph, directed = F, unconnected = T)
  
  #7 Degree assortativity 
  Degree_assortativity_nt=assortativity_degree(curr.graph) 
  
  #8 density 
  Density_nt=edge_density(curr.graph, loops=FALSE)
  
  #9 Cluster number 
  cluster_num_nt <- components(curr.graph)$no
  
  #10 Diameter
  Diameter=diameter(curr.graph)
  
  #11 edges number 
  edge_number <- ecount(curr.graph)
  
  #12 vertices number 
  vertice_number <- vcount(curr.graph)  
  
  #13 modularity 
  modularity_nt <- modularity(multilevel.community(curr.graph2, weights = E(curr.graph2)$weight))
  Module_group <- V(curr.graph)$sg

  ## for random network
  Betw_centrality_nt.s <- function(g){centr_betw(g, directed = FALSE,normalized=F)$centralization}
  Betw_centrality_nt_ER.avg = mean(laply(ER_OTU.list, Betw_centrality_nt.s))
  Betw_centrality_nt_ER.sd = sd(laply(ER_OTU.list, Betw_centrality_nt.s))
  
  Clos_centrality_nt.s <- function(g){centr_clo(g, normalized = F)$centralization}
  Clos_centrality_nt_ER.avg = mean(laply(ER_OTU.list, Clos_centrality_nt.s))
  Clos_centrality_nt_ER.sd = sd(laply(ER_OTU.list, Clos_centrality_nt.s))

  Degree_centrality_nt.s <- function(g){centr_degree(g, loops = F, normalized = F)$centralization}
  Degree_centrality_nt_ER.avg = mean(laply(ER_OTU.list,Degree_centrality_nt.s))
  Degree_centrality_nt_ER.sd = sd(laply(ER_OTU.list, Degree_centrality_nt.s))

  Transitivity_nt.s <- function(g){transitivity(g, type="global")}
  Transitivity_nt_ER.avg = mean(laply(ER_OTU.list, Transitivity_nt.s))
  Transitivity_nt_ER.sd = sd(laply(ER_OTU.list, Transitivity_nt.s))

  knn.nt.s <- function(g){mean(knn(g)$knn, na.rm=T)}
  knn_nt_ER.avg = mean(laply(ER_OTU.list, knn.nt.s))
  knn_nt_ER.sd = sd(laply(ER_OTU.list, knn.nt.s))
  
  Aver_path_nt.s <- function(g){mean_distance(g, directed = F, unconnected = T)}
  Aver_path_nt_ER.avg = mean(laply(ER_OTU.list, Aver_path_nt.s))
  Aver_path_nt_ER.sd = sd(laply(ER_OTU.list, Aver_path_nt.s))

  Degree_assortativity_nt_ER.avg = mean(laply(ER_OTU.list, assortativity_degree))
  Degree_assortativity_nt_ER.sd = sd(laply(ER_OTU.list, assortativity_degree))

  Density_nt_ER.avg = mean(laply(ER_OTU.list, edge_density))
  Density_nt_ER.sd = sd(laply(ER_OTU.list, edge_density))

  cluster_num_nt.s <- function(g){components(g)$no}
  cluster_num_nt_ER.avg = mean(laply(ER_OTU.list, cluster_num_nt.s))
  cluster_num_nt_ER.sd = sd(laply(ER_OTU.list, cluster_num_nt.s))

  Diameter_ER.avg = mean(laply(ER_OTU.list, diameter))
  Diameter_ER.sd = sd(laply(ER_OTU.list, diameter))
  
  edge_number_ER.avg = mean(laply(ER_OTU.list, ecount))
  edge_number_ER.sd = sd(laply(ER_OTU.list, ecount))

  vertice_number_ER.avg = mean(laply(ER_OTU.list, vcount))
  vertice_number_ER.sd = sd(laply(ER_OTU.list, vcount))

  modularity_nt.s <- function(g){modularity(multilevel.community(g))}
  modularity_nt_ER.avg = mean(laply(ER_OTU.list, modularity_nt.s))
  modularity_nt_ER.sd = sd(laply(ER_OTU.list, modularity_nt.s))

  random_network_avg = c("--", knn_nt_ER.avg, Aver_path_nt_ER.avg, Betw_centrality_nt_ER.avg, Clos_centrality_nt_ER.avg, Degree_assortativity_nt_ER.avg, 
                         Degree_centrality_nt_ER.avg, Density_nt_ER.avg, cluster_num_nt_ER.avg, Diameter_ER.avg, Transitivity_nt_ER.avg, 
                         vertice_number_ER.avg, edge_number_ER.avg, modularity_nt_ER.avg)
  random_network_sd = c("--", knn_nt_ER.sd, Aver_path_nt_ER.sd, Betw_centrality_nt_ER.sd, Clos_centrality_nt_ER.sd, Degree_assortativity_nt_ER.sd, 
                         Degree_centrality_nt_ER.sd, Density_nt_ER.sd, cluster_num_nt_ER.sd, Diameter_ER.sd, Transitivity_nt_ER.sd, 
                         vertice_number_ER.sd, edge_number_ER.sd, modularity_nt_ER.sd)

network_topo_list=data.frame(topo=c("RMT_threshold",
                                    "Average_nearest_neighbor_degree",
                                    "Average_path_length",
                                    "Betweenness_centrality",
                                    "Closeness_centrality",
                                    "Degree_assortativity",
                                    "Degree_centralization",
                                    "Density",
                                    "Cluster_num",
                                    "Diameter",
                                    "Transitivity",
                                    "Num_vertice",
                                    "Num_edge",
                                    "Modularity"),
                              empirical_network=as.character(c(rmt,
                                      knn_nt,
                                      Aver_path_nt,
                                      Betw_centrality_nt,
                                      Clos_centrality_nt,
                                      Degree_assortativity_nt,
                                      Degree_centrality_nt,
                                      Density_nt,
                                      cluster_num_nt,
                                      Diameter,
                                      Transitivity_nt,
                                      vertice_number,
                                      edge_number,
                                      modularity_nt)),
                                      random_network_avg=as.character(random_network_avg),
                                      random_network_sd=as.character(random_network_sd))
write.table(network_topo_list, paste(PARAM$folder.output, "/network_topo_list.xls",sep=""),sep="\t",row.names=F,quote = F)

############################################
############################################

####################################################
# Within-module degree and participation coefficient
####################################################

names(com$membership) <- com$names
edge <- as.data.frame(edge)
edge$V1 <- as.character(edge$V1)
edge$V2 <- as.character(edge$V2)

edge$V1.m <- com$membership[edge$V1]
edge$V2.m <- com$membership[edge$V2]

# compute degree: nodes ==> modules
degree.mat <- table(data.frame(c(edge$V1,edge$V2),c(edge$V2.m,edge$V1.m)))

# compute Z, P parameters
Z_P.dat <- matrix(0,nrow=vcount(curr.graph), ncol=2)
rownames(Z_P.dat) <- V(curr.graph)$name
colnames(Z_P.dat) <- c("Z", "P")

for(i in rownames(Z_P.dat)){
  Z_P.dat[i, "Z"] <- (degree.mat[i, as.character(com$membership[i])] -
                      mean(degree.mat[, as.character(com$membership[i])][names(com$membership)[com$membership == com$membership[i]]]))/
                      sd(degree.mat[, as.character(com$membership[i])][names(com$membership)[com$membership == com$membership[i]]])
  Z_P.dat[i, "P"] <- 1 - sum((degree.mat[i,]/DEG[i])^2)
}
Z_P.dat.0 <- as.data.frame(Z_P.dat)
Z_P.dat <- Z_P.dat[!is.na(Z_P.dat[,1]),] # remove nodes in modules which have only 2 members / nodes connect to each other within a module
Z_P.dat <- as.data.frame(Z_P.dat)
Z_P.dat$module <- com$membership[rownames(Z_P.dat)]
Z_P.dat$Phylum <- as.character(otu.data[rownames(Z_P.dat), "phylum"])
Z_P.dat$Phylum[!Z_P.dat$Phylum %in% taxa10] <- "Other"
Z_P.dat$Phylum <- factor(Z_P.dat$Phylum, levels=c(taxa10,"Other"))
Z_P.dat$size <- log2(curr.graph.dat$size[rownames(Z_P.dat)]*1000000/length(samples))


#col_roles <- c("Peripherals"="#ffdfba", "Connectors"="#baffc9", "Module hubs"="#bae1ff", "Network hubs"="#e3c6f0")
region.label <- data.frame(x = c(0.28, 0.84, 0.12, 0.88), 
                           y = c(min(Z_P.dat$Z)*1.05, min(Z_P.dat$Z)*1.05, max(max(Z_P.dat$Z)*1.05, 3*0.99),  max(max(Z_P.dat$Z)*1.05, 3*0.99)),
                           label = c("Peripherals", "Connectors", "Module hubs", "Network hubs")
                          )

pdf(paste(PARAM$folder.output, "/network_ZiPi_plot.pdf",sep=""),width=10,height=7.5,bg="white")
ggplot(Z_P.dat,aes(P,Z, colour = Phylum)) + # P threshold = 0.62, Z threshold = 2.5
  #geom_rect(aes(fill = "Peripherals", xmin = -Inf, xmax = 0.62, ymin = -Inf, ymax = 2.5), colour = NA) + 
  #geom_rect(aes(fill = "Connectors", xmin = 0.62, xmax = Inf, ymin = -Inf, ymax = 2.5), colour = NA) + 
  #geom_rect(aes(fill = "Module hubs", xmin = -Inf, xmax = 0.62, ymin = 2.5, ymax = Inf), colour = NA) + 
  #geom_rect(aes(fill = "Network hubs", xmin = 0.62, xmax = Inf, ymin = 2.5, ymax = Inf), colour = NA) + 
  ylim(min(Z_P.dat$Z)*1.05, max(max(Z_P.dat$Z)*1.05, 3)) +
  xlim(0, 1) +
  geom_hline(yintercept = 2.5, linetype=2) +
  geom_vline(xintercept = 0.62, linetype=2) +
  #geom_point(colour = "black" , aes(size=size)) + geom_point(aes(size = size*0.8)) + 
  geom_text(data = region.label, color = "black", size = 5, aes(x = x, y = y, label = label)) +
  geom_point(aes(size = size*0.8, color = Phylum)) +
  xlab("Among-module connectivity(Pi)") + ylab("Within-module connectivity(Zi)") + 
  #scale_fill_manual(name = "Roles", values=col_roles) + scale_size_continuous(range=c(0.5,4), guide="none") + 
  scale_size_continuous(range=c(0.5,4), guide="none") +
  scale_color_manual(values=c(
    "#8DD3C7",
    "#BF8686",
    "#BEBADA",
    "#FB8072",
    "#80B1D3",
    "#FDB462",
    "#B3DE69",
    "#BC80BD",
    "#CCEBC5",
    "#FFED6F",
    "#bdbdbd")) +
  theme_bw() + 
  theme(axis.text.y = element_text(color = "black", size = 15),
        axis.text.x = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 18, face = "bold"),
        axis.title.y = element_text(color = "black", size = 18, face = "bold"),
        legend.text=element_text(face="plain",size=15),legend.title = element_text(size=15,face="bold"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(face = "bold")) +
  guides(colour = guide_legend(override.aes = list(size = 5)))
dev.off()

Z_P_list <- data.frame(OTU_ID=rownames(Z_P.dat), Zi_score=Z_P.dat$Z, Pi_score=Z_P.dat$P)
Z_P.dat.0$Roles = "Peripherals"
Z_P.dat.0$Roles[Z_P.dat.0$Z > 2.5 & Z_P.dat.0$P < 0.62] = "Module hubs"
Z_P.dat.0$Roles[Z_P.dat.0$Z < 2.5 & Z_P.dat.0$P > 0.62] = "Connectors"
Z_P.dat.0$Roles[Z_P.dat.0$Z > 2.5 & Z_P.dat.0$P > 0.62] = "Network hubs"
write.table(Z_P_list, paste(PARAM$folder.output, "/network_ZiPi_list.xls",sep=""),sep="\t",row.names=F,quote = F)

 

#for gephi
############################################
############################################
suppressMessages(library("rgexf"))

node=data.frame(ID = c(1:vcount(curr.graph)), label= as.character(V(curr.graph)$name))
edge=as.data.frame(get.edges(curr.graph, c(1:ecount(curr.graph))))
dimnames(edge)[[2]]=c("source","target")
WGH=abs(E(curr.graph)$weight)
SIM=E(curr.graph)$sim
SIZ=log2(curr.graph.dat$size*1000000/length(samples))
module.color=V(curr.graph)$color
phylum.color=V(curr.graph)$color.tax
page_rank=page.rank(curr.graph,weights=E(curr.graph)$weight)$vector
Module_group=V(curr.graph)$sg
grafo=write.gexf(node,edge, defaultedgetype = "undirected",keepFactors=F,edgesWeight=WGH, nodesAtt = data.frame(SIZ,DEG,BET,page_rank,Module_group,Degree_centrality_nd,Clos_centrality_nd,Betw_centrality_nd,Transitivity_nd,knn_nd ),edgesAtt = data.frame(SIM,WGH))
if ("a"=="b"){
sink(paste(PARAM$folder.output, "/network_large.gexf",sep=""),append=FALSE,split=TRUE)
print(grafo)
sink()
}
#############################
#############################

## for node topo
node_topo_list <- data.frame(OTU_ID=names(SIZ),Size=SIZ,Degree=DEG,Betweenness=BET,page_rank,Module_group,Clos_centrality_nd,
                             Betw_centrality_nd,Transitivity_nd,knn_nd,Z=Z_P.dat.0$Z,P=Z_P.dat.0$P,Roles=Z_P.dat.0$Roles)
write.table(node_topo_list, paste(PARAM$folder.output, "/node_topo_list.xls",sep=""),sep="\t",row.names=F,,quote = F)

################################################################################################################
#
#         TOP50 OTU SUBGRAPH
#
################################################################################################################

# Get subgraph

curr.otus.top50 <- names(head(sort(rowSums(curr.ot),decreasing = T),50))
curr.graph.top50 <- induced.subgraph(curr.graph, vids=curr.otus.top50)
curr.otus.top50 <- names(degree(curr.graph.top50))[degree(curr.graph.top50) >= 1]

curr.graph.top50 <- induced.subgraph(curr.graph, vids=curr.otus.top50)
curr.graph2.top50 <- induced.subgraph(curr.graph2, vids=curr.otus.top50)

#Compute a layout

E(curr.graph.top50)$curved <- 0.25
E(curr.graph.top50)$color <- 0
E(curr.graph.top50)$color=col_edge[100+signif(E(curr.graph.top50)$sim*100,1)]

#community_split using co-occurrence network
com = multilevel.community(curr.graph2.top50, weights = E(curr.graph2.top50)$sim)
subgroup = split(com$names, com$membership)
## subgroup
V(curr.graph.top50)$sg = paste("module",com$membership,sep="_")
modu=names(sort(summary(as.factor(V(curr.graph.top50)$sg),maxsum=10000),decreasing=T))
V(curr.graph.top50)$color="grey"
for (i in 1:10) {
V(curr.graph.top50)$color[which(V(curr.graph.top50)$sg==modu[i])] = c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[c(3,5:7,9,12:13,15,20,19)][i]
}

####
node.id <- c(1:length(V(curr.graph.top50)$name))
names(node.id) <- V(curr.graph.top50)$name

edge <- get.edgelist(curr.graph.top50)
e <- cbind(node.id[edge[,1]],node.id[edge[,2]])
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(curr.graph.top50),area=20*(vcount(curr.graph.top50)^2),
                                       repulse.rad=(vcount(curr.graph.top50)^3.1))
## vertex label 
V(curr.graph.top50)$genus <- as.character(otu.data[V(curr.graph.top50)$name, "genus"])
#V(curr.graph.top50)$genus[V(curr.graph.top50)$genus == ""] <- as.character(otu.data[V(curr.graph.top50)$name[V(curr.graph.top50)$genus == ""], "family"])
#V(curr.graph.top50)$genus[V(curr.graph.top50)$genus == ""] <- as.character(otu.data[V(curr.graph.top50)$name[V(curr.graph.top50)$genus == ""], "order"])
#V(curr.graph.top50)$genus[V(curr.graph.top50)$genus == ""] <- as.character(otu.data[V(curr.graph.top50)$name[V(curr.graph.top50)$genus == ""], "class"])
#V(curr.graph.top50)$genus[V(curr.graph.top50)$genus == ""] <- as.character(otu.data[V(curr.graph.top50)$name[V(curr.graph.top50)$genus == ""], "phylum"])
#V(curr.graph.top50)$genus[V(curr.graph.top50)$genus == ""] <- as.character(otu.data[V(curr.graph.top50)$name[V(curr.graph.top50)$genus == ""], "domain"])
#V(curr.graph.top50)$genus[V(curr.graph.top50)$genus == ""] <- "Unassigned"

############################
# filled by modules 
pdf(paste(PARAM$folder.output, "/subnetwork_module_top50.pdf",sep=""),width=10,height=10,bg="white")
plot(curr.graph.top50, layout=l, 
    fade=TRUE,
    bg="white",
    trans=TRUE,
    vertex.label=V(curr.graph.top50)$genus,
    vertex.label.color="black",
    vertex.label.family="sans", 
    #vertex.label.cex=1.1, 
    vertex.size=log2(curr.graph.dat$size[curr.graph.dat$node %in% curr.otus.top50]*1000000/length(samples))/2, 
    vertex.color=V(curr.graph.top50)$color,
    vertex.label.cex=0.7,
    edge.arrow.size=0, 
    edge.color=E(curr.graph.top50)$color, 
    edge.width=E(curr.graph.top50)$weight*4-min(E(curr.graph.top50)$weight)*4+0.005, 
    edge.alpha=0.85,
    vertex.frame.color=NA)
with(data.frame(V(curr.graph.top50)$sg,V(curr.graph.top50)$color),legend("topright",title="Modules", c(as.character(modu[1:min(10,length(modu))])), inset=c(-0.05, 0),  
    title.adj = 0, 
    pt.bg=c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[c(3,5:7,9,12:13,15,20,19)][1:min(10,length(modu))],bty = "n",col ="black", pch = 22, cex=1.0))
with(data.frame(V(curr.graph.top50)$sg,V(curr.graph.top50)$color),legend("bottomright",title="Similarity", 
    title.adj = 0,
    c("-1.0","-0.9","-0.8","-0.7","0.7","0.8","0.9","1.0"), col=col_edge[100+c(-0.99,-0.9,-0.8,-0.7,0.7,0.8,0.9,0.99)*100],bty = "n", pch = "--", cex=1.0,pt.cex=2.5,pt.lwd=2))
dev.off();

############################
# filled by Phylum 
pdf(paste(PARAM$folder.output, "/subnetwork_phylum_top50.pdf",sep=""),width=10,height=10,bg="white")
plot(curr.graph.top50, layout=l, 
    fade=TRUE,
    bg="white",
    trans=TRUE,
    vertex.label=V(curr.graph.top50)$genus,
    vertex.label.color="black",
    vertex.label.family="sans", 
    #vertex.label.cex=1.1, 
    vertex.size=log2(curr.graph.dat$size[curr.graph.dat$node %in% curr.otus.top50]*1000000/length(samples))/2, 
    vertex.color=curr.graph.dat$color.tax[curr.graph.dat$node %in% curr.otus.top50],
    vertex.label.cex=0.5,
    edge.arrow.size=0, 
    edge.color=E(curr.graph.top50)$color, 
    edge.width=E(curr.graph.top50)$weight*4-min(E(curr.graph.top50)$weight)*4+0.005, 
    edge.alpha=0.85,
    vertex.frame.color="Black")
with(data.frame(V(curr.graph)$sg,V(curr.graph)$color),legend("topright",title="Phylum", c(as.character(taxa10[!is.na(taxa10)])), inset=c(-0.05, 0), 
    title.adj = 0, 
    pt.bg=c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#BC80BD", "#CCEBC5", "#FFED6F")[1:min(10,length(taxa10[!is.na(taxa10)]))], 
    bty = "n",col ="black", pch = 22, cex=.7))
with(data.frame(V(curr.graph)$sg,V(curr.graph)$color),legend("bottomright",title="Similarity", 
    title.adj = 0, 
    c("-1.0","-0.9","-0.8","-0.7","0.7","0.8","0.9","1.0"), col=col_edge[100+c(-0.99,-0.9,-0.8,-0.7,0.7,0.8,0.9,0.99)*100],bty = "n", pch = "--", cex=1.0,pt.cex=2.5,pt.lwd=2))
dev.off();

############################
# filled by Groups
if(length(levels(sample.data$Group)) <= 10){

pie_size <- pie_size[,V(curr.graph.top50)$name]

col_source=c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[c(3,5:7,9,12:13,15,20,19)]
col_pie=rep(list(col_source[1:length(levels(Groups))]),vcount(curr.graph.top50))

pdf(paste(PARAM$folder.output, "/subnetwork_group_top50.pdf",sep=""),width=10,height=10,bg="white")
plot(curr.graph.top50, layout=l, 
    fade=TRUE,
    bg="white",
    trans=TRUE,
    vertex.label=V(curr.graph.top50)$genus,
    vertex.label.color="black",
    vertex.label.family="sans", 
    #vertex.label.cex=1.1, 
    vertex.size=log2(curr.graph.dat$size[curr.graph.dat$node %in% curr.otus.top50]*1000000/length(samples))/3, 
    vertex.shape="pie", 
    edge.arrow.size=0, 
    edge.color=E(curr.graph.top50)$color, 
    edge.width=E(curr.graph)$weight*4-min(E(curr.graph)$weight)*4+0.005, 
    edge.alpha=0.85,
    vertex.pie=pie_size, 
    vertex.pie.color=col_pie, 
    vertex.frame.color="Black")
with(data.frame(V(curr.graph)$sg,V(curr.graph)$color),legend("topright",title="Group", c(as.character(Groups)), 
    title.adj = 0, 
    pt.bg=col_source[1:length(levels(Groups))],bty = "n",col ="black", pch = 22, cex=1.0))
with(data.frame(V(curr.graph)$sg,V(curr.graph)$color),legend("bottomright",title="Similarity", 
    title.adj = 0, 
    c("-1.0","-0.9","-0.8","-0.7","0.7","0.8","0.9","1.0"), col=col_edge[100+c(-0.99,-0.9,-0.8,-0.7,0.7,0.8,0.9,0.99)*100],bty = "n", pch = "--", cex=1.0,pt.cex=2.5,pt.lwd=2))
dev.off();

}

#####################################################################################
## 转换gml结果文件
node_topo_list$Domain=as.character(otu.data[names(SIZ), "domain"])
node_topo_list$Phylum=as.character(otu.data[names(SIZ), "phylum"])
node_topo_list$Class=as.character(otu.data[names(SIZ), "class"])
node_topo_list$Order=as.character(otu.data[names(SIZ), "order"])
node_topo_list$Family=as.character(otu.data[names(SIZ), "family"])
node_topo_list$Genus=as.character(otu.data[names(SIZ), "genus"])
node_topo_list$Species=as.character(otu.data[names(SIZ), "species"])
node_topo_list$Top50="FALSE"
node_topo_list[names(head(sort(rowSums(curr.ot),decreasing = T),50)), "Top50"] <- "TRUE"
node_topo_list[is.na(node_topo_list)] <- 0

curr.graph <- delete_vertex_attr(curr.graph, "color") %>% 
              delete_vertex_attr("sg") %>% 
              delete_edge_attr("color") %>% 
              delete_edge_attr("curved")
 
node_topo_list$OTU_ID <- as.character(node_topo_list$OTU_ID)
node_topo_list$Module_group <- as.character(node_topo_list$Module_group)
node_topo_list$Roles <- as.character(node_topo_list$Roles)

for(i in colnames(node_topo_list)){
    curr.graph <- set_vertex_attr(curr.graph, i, value = node_topo_list[,i])
}

write_graph(curr.graph, paste(PARAM$folder.output, "/network_large.gml",sep=""), format = "gml")

warnings()
