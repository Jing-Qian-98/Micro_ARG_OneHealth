library(randomForest)
library(pROC)
library("plyr")
probabilities <- read.csv("metadata.tsv", sep="\t", header=T, stringsAsFactors=FALSE, comment.char="")
probabilities <- arrange(probabilities, probabilities[,1])

map <- read.table("map.txt", sep="\t", header=F, stringsAsFactors=FALSE)
map <- arrange(map, map$V1)
Group <- as.factor(map$V4)
Group.l <- levels(Group)
levels(Group) <- 0:2

pdf("ROC.pdf")
rocobj <- plot.roc(Group, as.numeric(probabilities[, 2]),
                   main="", 
                   of="se",
                   #  of="sp",
                   percent=TRUE,
                   ci=TRUE,# compute AUC (of AUC by default)
                   print.auc=TRUE,  # print the AUC (will contain the CI)
                   print.auc.x=50,
                   print.auc.y=20,
                   #specificities=seq(0, 100, 5), # on a select set of specificities
                    sensitivity=seq(0, 100, 5), # on a select set of specificities
                   ci.type="shape", 
                   ci.col="grey")  
CI=ci.auc(rocobj)
text(40,10, paste("95% CI: ",round(as.numeric(CI)[1], 2),"% - ",round(as.numeric(CI)[3], 2),"%",sep=""))
dev.off()

ci.auc(rocobj)
