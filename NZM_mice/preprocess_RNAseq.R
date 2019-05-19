# The following R script identifies differentialy expressed genes for disease and diet from
# RNAseq data. The file "NZM_summary.txt" is obtained from file "preprocess_RNAseq.sh" Please
# make sure that R package DESeq2 is installed in the system.
# Copyright (c) 2016 by Yask Gupta, Dept. of Dermatology Lubeck, Germany. All Rights Reserved.

#load library DESeq2 
library("DESeq2")

# read the feature count output 
gene.expr <- read.table(file="NZM_summary.txt",header=T,row.name=1,sep="\t")

# load attributes for RNAseq data
datTraits.rnaseq <- read.table(file="datTraits.tab",header=T,sep="\t")

#convert disease as factor rather than numeric
datTraits.rnaseq$Disease <- as.factor(ifelse(datTraits.rnaseq$Disease == 0,"N","Y")) 

# load data in DESeq2 and get differentially expressed genes for disease
dds.rnaseq <- DESeqDataSetFromMatrix(countData = gene.expr,colData = datTraits.rnaseq,design = ~Sex+Disease)
dds.rnaseq.dis <- DESeq(dds.rnaseq,test="LRT",reduced = ~Sex)
res.dis <- results(dds.rnaseq.dis)
res.dis <- res.dis[order(res.dis$padj),]
res.dis = cbind(Gene=rownames(res.dis),as.data.frame(res.dis))
write.table(subset(res.dis,padj < 0.05),file="Results_Disease_RNAseq.tab",col.names=T,row.names=F,sep="\t",quote=F)

# load data in DESeq2 and get differentially expressed genes for diet
dds.rnaseq <- DESeqDataSetFromMatrix(countData = gene.expr,colData = datTraits.rnaseq,design = ~Sex+Diet)
dds.rnaseq.diet <- DESeq(dds.rnaseq,test="LRT",reduced = ~Sex)
res.diet <- results(dds.rnaseq.diet)
res.diet <- res.diet[order(res.diet$padj),]
res.diet = cbind(Gene=rownames(res.diet),as.data.frame(res.diet))
write.table(subset(res.diet,padj < 0.05),file="Results_Diet_RNAseq.tab",col.names=T,row.names=F,sep="\t",quote=F)

# Another analysis performed for Figure 5 is differentially expressed
# genes for anti-nuclear (ANA) and calorie-reduced diet vs all the other
# diets (control and western)
# A seperate trait file for ANA is created "sampleTableANA.tab" for this
# analysis
  
# load data in DESeq2 and get differentially expressed genes for ANA
sampleTable <- read.table(file="sampleTableANA.tab",header=T,sep="\t")
dds.ana <- DESeqDataSetFromMatrix(countData = gene.expr,colData = sampleTable,design = ~Sex+ANA100)
dds.ana <- DESeq(dds.ana,test="LRT", reduced= ~Sex)
res.ana <- results(dds.ana)
res.ana <-res.ana[order(res.ana$padj),]
res.ana = cbind(Gene=rownames(res.ana),as.data.frame(res.ana))
write.table(subset(res.ana,padj < 0.05),file="Results_ANA_RNAseq.tab",col.names=T,row.names=F,sep="\t",quote=F)


# load data in DESeq2 and get differentially expressed genes for Calorie-reduced Diet
sampleTable.cal <- sampleTable 
sampleTable.cal$DietCal = as.factor(ifelse(sampleTable.cal$Diet == "cal","Y","N"))
dds.cal <- DESeqDataSetFromMatrix(countData = gene.expr,colData = sampleTable.cal,design = ~Sex+Diet)
dds.cal <- DESeq(dds.cal,test="LRT",reduced= ~Sex)
res.cal <- results(dds.cal)
res.cal <-res.cal[order(res.cal$padj),]
res.cal = cbind(Gene=rownames(res.cal),as.data.frame(res.cal))
write.table(subset(res.cal,padj < 0.05),file="Results_DietCal_RNAseq.tab",col.names=T,row.names=F,sep="\t",quote=F)


