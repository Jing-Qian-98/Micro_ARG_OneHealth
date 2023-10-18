setwd("./")
################################################################################
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
suppressMessages(library("RColorBrewer"))
suppressMessages(library("RMThreshold"))
################################################################################
################################################################################

argv <- commandArgs(TRUE)

  if (length(argv) !=5 ){
    stop('The script needs  files:otu_table; sample_metadata; result_outputfile; data_inputfile ;core_number, please put them in a dir  ')
  }
 
#Preallocate global data structures 
PARAM <- list();
PARAM$folder.data <- argv[4]
PARAM$folder.output <- argv[3]              #PUT RESULTS/OUTPUT FOLDER HERE
PARAM$file.sample_metadata <- argv[2]
PARAM$file.otu_table <- argv[1]
PARAM$use.cores <- argv[5];
###################################
commond=paste("mkdir -p ",PARAM$folder.output ,sep="") 
system(commond)

###################################
# 根据样本出现频率过滤
remove_rare <- function(table , cutoff ) {
  row2keep <- c()
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}
###################################

#ASV序列总和大于10
PARAM$thresh.otu_size <- 10; ### from 3 to 30
#样品reads数量大于1000
PARAM$thresh.sample_size <- 1000;
#过滤后剩余序列占比大于0.2
PARAM$thresh.retained_seqs <- 0.2;
###################################

#Read & pre-process data
########################
#Read sample raw metadata
sample.data.raw <- read.table(PARAM$file.sample_metadata, header=T, sep="\t",comment.char="@",colClasses = "character");
rownames(sample.data.raw) <-as.character(sample.data.raw[,1]);
#Load raw OTU table
otu_table=read.table(file=PARAM$file.otu_table, header=T, sep="\t", row.names=1,comment.char="@",check.names = F)
otu_table <- otu_table[, -ncol(otu_table)]
idx = colnames(otu_table) %in% rownames(sample.data.raw) 
otu_table_rare_removed <- remove_rare(table=otu_table[,idx],max(5,min(20,length(idx)/5)))     #########OTU至少出现在N个样品中N的为样本数的1/5，但是最小为5最大为20。
ot.raw = Matrix(as.matrix(otu_table_rare_removed), sparse=T)

size.otu.raw1 <- rowSums(otu_table[,idx]);
size.sample.raw1 <- colSums(otu_table[,idx]);
size.sample.raw2 <- colSums(ot.raw);
size.otu.raw2 <- rowSums(ot.raw);

########################
#Filter data
#=> by a minimum of "used" sequences per sample (relative to the raw sequence count per sample)
#=> by minimum sample size
#=> by minimum OTU size (across all samples)
########################
#Minimum relative number of retained sequences & minimum sample size
retained.tmp <- size.sample.raw2 /  size.sample.raw1
samples <- colnames(ot.raw)[size.sample.raw2 > PARAM$thresh.sample_size & retained.tmp > PARAM$thresh.retained_seqs];           #######如果一个样本过多的ASV序列被过滤掉了（只剩20%以下），那么可以考虑这个样本是否真的与其他样本有关联，是否要纳入到网络分析中来。如果大量样本被过滤，此时可以考虑使用属水平数据
ot.tmp <- ot.raw[, samples];
#Minimum OTU size
otus <- rownames(ot.tmp)[rowSums(ot.tmp) > PARAM$thresh.otu_size];
#Prune data accordingly
ot <- ot.raw[otus, samples];
size.sample <- colSums(ot); n.sample <- length(size.sample);
size.otu <- rowSums(ot); n.otu <- length(size.otu);
ot.rel <- t(t(ot) / size.sample);
ot.rel.raw <- t(t(ot) / size.sample.raw1[samples]);
########################


################################################################################
################################################################################
#Calculate pairwise indices on OTUs, for count correction
#=> calculate "raw" pairwise similarities between OTUs (SparCC, cophenetic phylogenetic distance)
#=> for each OTU-wise matrix, get the transformed OTU association matrix S (or "C", "Phi"), based on OTU-OTU-correlation across the raw matrix
#=> in other words, pairwise OTU association is calculated as the similarity of OTUs in the original similarity space
#=> matrix S is a symmetric, transformed correlation matrix: S[i,j] = 0.5 * (1 + rho) => S scales on [0,1]
################################################################################
#SparCC correlation
#if ("SparCC" == "SparCC") {
########################
#SparCC correlation analysis
#=> following Friedman & Alm, PLOS CB, 2012
#=> apparently, the pseudocount has a large influence on the results...
#=> expects an OTU table in which OTUs are rows, samples are columns
########################
#Set parameters
size.thresh <-10;
pseudocount <- 10^-9;
nblocks <- 100;
use.cores <- PARAM$use.cores;
########################

########################
#Load packages for parallel processing
suppressMessages(require("foreach"));
suppressMessages(require("bigmemory"));
suppressMessages(library("doMC", quietly=T));
#Register cluster
registerDoMC(cores=use.cores);
########################

########################
#Filter OTU table by removing all OTUs that are observed less then <size.thresh> times across all samples
#=> their correlations to all other OTUs will be (manually) set to 0 later on
#Add pseudocount to all observed counts (to avoid issues with zeroes in log-space)
cat(paste("SparCC preallocation steps =>", Sys.time()), sep="\n");

remove.otus <- rowSums(ot) < size.thresh; keep.otus <- ! remove.otus;
#o.t <- otu.table[1:1000, ] + pseudocount;   #for test
o.t <- ot[! remove.otus, ] + pseudocount;
otus <- rownames(o.t);
ot.2=ot[otus,]
n.otu <- length(otus);
nblocks <- min(100,floor(n.otu/3))
########################
#Preallocate blocks for parallel processing & Aitchinson's T matrix
#=> based on https://gist.github.com/bobthecat/5024079
size.split <- floor(n.otu / nblocks);
if (size.split < 1) {size.split <- 1}
my.split <- list(); length(my.split) <- nblocks;
my.split[1:(nblocks-1)] <- split(1:(size.split*(nblocks-1)), rep(1:(nblocks-1), each = size.split));
my.split[[nblocks]] <- (size.split*(nblocks-1)):n.otu;
dat.split <- mclapply(my.split, function(g) {o.t[g,]}, mc.cores=use.cores);
#Get combinations of splits
my.combs <- expand.grid(1:length(my.split), 1:length(my.split));
my.combs <- t(apply(my.combs, 1, sort));
my.combs <- unique(my.combs);
#Preallocate Aitchinson's T matrix as big.matrix ("shared" in memory, so accessible from w/in foreach loop)
mat.T <- big.matrix(nrow=n.otu, ncol=n.otu, dimnames=list(otus, otus), shared=T);
mat.T.desc <- describe(mat.T);
cat(paste("Done with preallocations =>", Sys.time()), sep="\n");
########################

########################
#Compute Aitchinson's T matrix
#=> iterate through each block combination, calculate matrix
#between blocks and store them in the preallocated matrix on both
#symmetric sides of the diagonal
cat(paste("Starting parallel T matrix calculation =>", Sys.time()), sep="\n");
results <- foreach(i = 1:nrow(my.combs)) %dopar% {
  #Get current combination
  curr.comb <- my.combs[i, ];
  #Get current data
  g.1 <- my.split[[curr.comb[1]]];
  g.2 <- my.split[[curr.comb[2]]];
  dat.1 <- dat.split[[curr.comb[1]]];
  dat.2 <- dat.split[[curr.comb[2]]];
  #Get current part of Aitchinson's matrix
  curr.T <- apply(dat.1, 1, function(x) {apply(dat.2, 1, function(y) {var(log(x/y))})});
  #Store
  curr.mat.T <- attach.big.matrix(mat.T.desc);
  curr.mat.T[g.2, g.1] <- curr.T;
  curr.mat.T[g.1, g.2] <- t(curr.T);
  #Return
  TRUE
}
cat(paste("Done with parallel T matrix calculation =>", Sys.time()), sep="\n");
########################

########################
#Compute component variations ("t_i")
cat(paste("Computing component variations =>", Sys.time()), sep="\n");
var.t <- colsum(mat.T);
cat(paste("Done with component variation calculations =>", Sys.time()), sep="\n");
#Estimate component variances ("omega_i") from t_i by solving a linear equation system
cat(paste("Estimating component variances from linear equation system =>", Sys.time()), sep="\n");
#mat.a <- diag(rep(n.otu-2,n.otu))+1;
mat.a <- matrix(data=1, nrow=n.otu, ncol=n.otu); diag(mat.a) <- n.otu-1;
omega <- sqrt(solve(a=mat.a, b=var.t));
cat(paste("Done with component variation estimation =>", Sys.time()), sep="\n");
#Estimate pairwise correlations based on these values
cat(paste("Estimating correlations =>", Sys.time()), sep="\n");
global.sparcc=matrix();
global.sparcc <- foreach(i = 1:n.otu, .combine=rbind, .multicombine=T) %do% {(omega[i]^2 + omega^2 - mat.T[i,]) / (2 * omega[i] * omega)}
rownames(global.sparcc) <- colnames(global.sparcc) <- otus;
cat(paste("Done with correlation estimation; returning data matrix =>", Sys.time()), sep="\n");
########################
if (length(which(global.sparcc>1))>0){
system("touch warnings")
system(paste("echo \"produced sparcc simmilarity greater than 1; we change it to 1; this may because of the microbes (i.e. ",nrow(ot),") in the data is too simple after filter!!!  Ths project is not suggested to use network analysis\" > warnings",sep=""))
#print(paste("produced sparcc simmilarity greater than 1; we change it to 1; this may because of the microbes (i.e. ",nrow(ot),") in the data is too simple after filter!!!  Ths project is not suggested to use network analysis",sep=""))
}
if (length(which(global.sparcc< -1))>0){
system("touch warnings")
system(paste("echo \"produced sparcc simmilarity greater than -1; we change it to -1; this may because of the microbes (i.e. ",nrow(ot),") in the data is too simple after filter!!!  Ths project is not suggested to use network analysis\" > warnings",sep=""))
#print(paste("produced sparcc simmilarity smaller than -1; we change it to -1; this may because of the microbes (i.e. ",nrow(ot),") in the data is too simple after filter!!!  Ths project is not suggested to use network analysis",sep=""))
}
global.sparcc[upper.tri(global.sparcc)] <- t(global.sparcc)[upper.tri(global.sparcc)]   ##bug fix
global.sparcc[which(global.sparcc>1)]=1
global.sparcc[which(global.sparcc< -1)]=-1
global.sparcc[which(is.na(global.sparcc))]=0
res=rm.get.threshold(global.sparcc,nr.thresholds=121,interval=c(0.30,min(0.90,max(abs(global.sparcc)[abs(global.sparcc)!=1]))), discard.zeros = T,nr.breaks=100,discard.outliers = T,min.mat.dim = 4)
RMT=data.frame(R=res$tested.thresholds,P=res$p.chi)
rmt=0
for ( i in 1:120){
if (RMT$P[i]>0.05 & RMT$P[i]>=RMT$P[i+1] ){
 rmt=RMT$R[i]
 break
 }
}
if (rmt == 0){
for ( i in 1:120){
if (RMT$P[i]>0.01 & RMT$P[i]>=RMT$P[i+1]){
 rmt=RMT$R[i]
 break
 }
}
}
if (rmt == 0){
for ( i in 1:120){
if (RMT$P[i]>0.001 & RMT$P[i]>=RMT$P[i+1]){
 rmt=RMT$R[i]
 break
 }
}
}
if (rmt == 0){
for ( i in 1:120){
if (RMT$P[i]>0.0001 & RMT$P[i]>=RMT$P[i+1]){
 rmt=RMT$R[i]
 break
 }
}
}
if ( rmt < 0.6 ) {
rmt=0.6
}

RMT.sparcc <- global.sparcc;
RMT.sparcc[RMT.sparcc < rmt & RMT.sparcc > -rmt] <- 0;
# 将过滤后以及网络图构建数据存储至RData
save(ot,ot.2,o.t,ot.rel,ot.rel.raw, otu_table,otu_table_rare_removed, n.otu, file=paste(PARAM$folder.data, "/All_samples.otu_table.RData", sep=""));
save(rmt,global.sparcc,RMT.sparcc, file=paste(PARAM$folder.output, "/All_samples.SparCC_global.RData", sep=""));

########################
