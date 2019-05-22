#!/usr/bin/env R
# For functional microbiome communities (FMCs) use SparCC algorithm to find co-occurring communities
# use otutab.txt from script preprocess_Microbiome.sh & remove the otu where 30% of samples have less than 25 counts
# this file can be defined as core microbiota "otuCore.txt"
# Here is the command for creating co-occurring communities "python SparCC.py otuCore.txt --cor_file=CorMicro.txt"
# To access the co-occurring communities we used functions from R package WGCNA 
library(gdata)
library(WGCNA)
library(DESeq2)
sparcc.mat <- read.table(file="CorMicro.txt",header=T,row.names=1,sep="\t");
adj.matrix <- (sparcc.mat+0.5)**4
adj.matrix.scaled <- apply(adj.matrix, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
lowerTriangle(adj.matrix.scaled) <- upperTriangle(adj.matrix.scaled,byrow=T)
dissTOM <- TOMdist(adj.matrix.scaled)
geneTree <- hclust(as.dist(dissTOM),method="average")
modules.micro <- cutreeDynamic(dendro=geneTree,distM=dissTOM,method="hybrid",deepSplit=3,pamRespectsDendro=F,minClusterSize=10)

# transform the data according to McMurdle et al. https://doi.org/10.1371/journal.pcbi.1003531
library(DESeq2)
data.otu <- read.table(file="otuCore.txt",header=T,row.names=1,sep="\t",comment.char="")
data.otu.vst <- varianceStabilizingTransformation(as.matrix(data.otu+1))
data.otu.vst <- data.otu.vst[match(rownames(adj.matrix.scaled),rownames(data.otu.vst)),]
colnames(data.otu.vst) <- gsub("X","",colnames(data.otu.vst))
# Function from WGCNA R package
me.micro <- moduleEigengenes(t(data.otu.vst),labels2colors(modules.micro))
rownames(me.micro$eigengenes) <- colnames(data.otu.vst)
colnames(me.micro) <- paste("FMC",1:ncol(me.micro),sep="")
write.table(me.micro,file="MEs_micro.tab",col.names=T,row.names=T,sep="\t",quote=F)
# After that we create a file which contain eigenOTU values and RNAseq data reference together
# to identify differentially expressed genes for FMCs. "MicroRNAseqTraits.tab"
microdatTraits <- read.table(file="MicroRNAseqTraits.tab",header=T,sep="\t")
# the following r "RNAseq_raw.RData" object gene.expr object which is count data for RNAseq generated from the same mice
lnames <- load("/home/yask/artem_data/Mycobiome/RNAseq_raw.RData")
gene.expr <- gene.expr[,match(microdatTraits$SampleRef,colnames(gene.expr))]

# results for FMC1
dds.FMC1 <- DESeqDataSetFromMatrix(countData = gene.expr,colData = datTraits.rnaseq,design = ~Sex + FMC1)
dds.FMC1 <- DESeq(dds.FMC1,test="LRT",reduced = ~Sex)
res.FMC1 <- results(dds.FMC1)
res.FMC1 <- res.FMC1[order(res.FMC1$padj),]
res.FMC1 <- cbind(Gene=rownames(res.FMC1),as.data.frame(res.FMC1))
write.table(subset(res.FMC1,padj < 0.05),file="Results_FMC1_RNAseq.tab",col.names=T,row.names=F,sep="\t",quote=F)

# results for FMC3
dds.FMC3 <- DESeqDataSetFromMatrix(countData = gene.expr,colData = datTraits.rnaseq,design = ~Sex + FMC3)
dds.FMC3 <- DESeq(dds.FMC3,test="LRT",reduced = ~Sex)
res.FMC3 <- results(dds.FMC3)
res.FMC3 <- res.FMC3[order(res.FMC3$padj),]
res.FMC3 <- cbind(Gene=rownames(res.FMC3),as.data.frame(res.FMC3))
write.table(subset(res.FMC3,padj < 0.05),file="Results_FMC3_RNAseq.tab",col.names=T,row.names=F,sep="\t",quote=F)

# results for FMC5
dds.FMC5 <- DESeqDataSetFromMatrix(countData = gene.expr,colData = datTraits.rnaseq,design = ~Sex + FMC5)
dds.FMC5 <- DESeq(dds.FMC5,test="LRT",reduced = ~Sex)
res.FMC5 <- results(dds.FMC5)
res.FMC5 <- res.FMC5[order(res.FMC5$padj),]
res.FMC5 <- cbind(Gene=rownames(res.FMC5),as.data.frame(res.FMC5))
write.table(subset(res.FMC5,padj < 0.05),file="Results_FMC5_RNAseq.tab",col.names=T,row.names=F,sep="\t",quote=F)

# results for FMC6
dds.FMC6 <- DESeqDataSetFromMatrix(countData = gene.expr,colData = datTraits.rnaseq,design = ~Sex + FMC6)
dds.FMC6 <- DESeq(dds.FMC6,test="LRT",reduced = ~Sex)
res.FMC6 <- results(dds.FMC6)
res.FMC6 <- res.FMC6[order(res.FMC6$padj),]
res.FMC6 <- cbind(Gene=rownames(res.FMC6),as.data.frame(res.FMC6))
write.table(subset(res.FMC6,padj < 0.05),file="Results_FMC6_RNAseq.tab",col.names=T,row.names=F,sep="\t",quote=F)
