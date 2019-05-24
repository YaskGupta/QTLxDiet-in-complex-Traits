#!/usr/bin/env R
# For functional fungal communities (FFCs) use SparCC algorithm to find co-occurring communities
# this file can be used from preprocess_Mycobiome.sh "otu_table.txt"
# Here is the command for creating co-occurring communities "python SparCC.py otu_table.txt --cor_file=CorMyco.txt"
# To access the co-occurring communities we used functions from R package WGCNA 

library(gdata)
library(WGCNA)
library("DESeq2")
sparcc.mat <- read.table(file="CorMyco.txt",header=T,row.names=1,sep="\t");
adj.matrix <- (sparcc.mat+0.5)**4
adj.matrix.scaled <- apply(adj.matrix, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
lowerTriangle(adj.matrix.scaled) <- upperTriangle(adj.matrix.scaled,byrow=T)
dissTOM <- TOMdist(adj.matrix.scaled)
geneTree <- hclust(as.dist(dissTOM),method="average")
modules.artem <- cutreeDynamic(dendro=geneTree,distM=dissTOM,deepSplit=2,pamRespectsDendro=F,minClusterSize=5)

data.otu <- read.table(file="otu_table.txt",header=T,row.names=1,sep="\t")
data.clr <- varianceStabilizingTransformation(as.matrix(data.otu+1))
data.clr <- data.clr[match(rownames(adj.matrix.scaled),rownames(data.clr)),]
colnames(data.clr) <- gsub("X","",colnames(data.clr))
me.myco <- moduleEigengenes(t(data.clr),labels2colors(modules.artem))
rownames(me.myco$eigengenes) <- colnames(data.clr)

datTraits <- read.table(file="TraitsMyco.tab",sep="\t",row.names=1,header=T)
datTraits <- datTraits[match(colnames(data.clr),datTraits$sample_content),]

datTraits4cor <- data.frame(Sex=ifelse(datTraits$Sex == "m",1,2),Diet=ifelse(datTraits$Diet == "cal",1,ifelse(datTraits$Diet == "con",2,3)),MiceDisease=ifelse(datTraits$Disease == "Yes",2,1),Stage=ifelse(datTraits$Category == "Final",3,ifelse(datTraits$Category == "Trans",2,1)))
datcor <- corAndPvalue(me.myco$eigengenes,datTraits4cor)$cor
datp <- corAndPvalue(me.myco$eigengenes,datTraits4cor)$p
datp[,1] = as.numeric(apply(me.myco$eigengenes,2,function(x){wilcox.test(x ~ datTraits$Sex)$p.value}))
datp[,2] = as.numeric(apply(me.myco$eigengenes,2,function(x){kruskal.test(x ~ datTraits$Diet)$p.value}))
datp[,3] = as.numeric(apply(me.myco$eigengenes,2,function(x){wilcox.test(x ~ datTraits$MiceDisease)$p.value}))
datp[,5] = as.numeric(apply(me.myco$eigengenes,2,function(x){kruskal.test(x ~ datTraits$Category)$p.value}))

datcor[,1] = as.numeric(apply(me.myco$eigengenes,2,function(x){cor.test(x,datTraits4cor$Sex,method="spearman")$estimate}))
datcor[,2] = as.numeric(apply(me.myco$eigengenes,2,function(x){cor.test(x,datTraits4cor$Diet,method="spearman")$estimate}))
datcor[,3] = as.numeric(apply(me.myco$eigengenes,2,function(x){cor.test(x,datTraits4cor$MiceDisease,method="spearman")$estimate}))
datcor[,5] = as.numeric(apply(me.myco$eigengenes,2,function(x){cor.test(x,datTraits4cor$Stage,method="spearman")$estimate}))

textMatrix = paste(signif(datcor, 2), "\n(",signif(datp, 2), ")", sep = "");
names.module <- c(paste("FFC",c(1:6),sep=""))
pdf(file="myco_wgcna.pdf")
labeledHeatmap(datcor,xLabels=colnames(datcor),yLabels=names.module,colorLabels = TRUE,textMatrix=textMatrix,colors=greenWhiteRed(200),setStdMargins = FALSE,zlim = c(-1,1),font.lab.x=2,font.lab.y=2)
dev.off()

# After that we create a file which contain eigenOTU values and RNAseq data reference together
# to identify differentially expressed genes for FMCs. "MicroRNAseqTraits.tab"
# RNAseq_raw.RData contain gene.expr object from DESeq2 analaysis of RNAseq data  
lnames <- load("RNAseq_raw.RData")

datTraits.rnaseq <- read.table(file="MycoRNAseqTraits.tab",header=T,sep="\t")
gene.expr <- gene.expr[,match(datTraits.rnaseq$sampleName,colnames(gene.expr))]
dds.wgcna.myco <- DESeqDataSetFromMatrix(countData = gene.expr,colData = datTraits.rnaseq,design = ~conc + Sex + FFC4)
dds.wgcna.lrt<-DESeq(dds.wgcna.myco,test="LRT",reduced = ~conc + Sex)
res.wgcna.lrt<-results(dds.wgcna.lrt)
res.wgcna.lrt<-res.wgcna.lrt[order(res.wgcna.lrt$padj),]
write.table(subset(res.wgcna.lrt,padj < 0.05),file="Results_FFC4_RNAseq.tab",col.names=T,row.names=F,sep="\t",quote=F)
