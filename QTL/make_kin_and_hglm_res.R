# install R packages DOQTL and hglm

library(DOQTL)
library(hglm)

# load happy2array output RData file ("QTLRelInput.RData"),
# box-cox transformed phenotype (NCMB_Sept2018.tab) and 
# SNP annotation (snp.info) in the R workspace
# Copyright (c) 2016 by Yask Gupta, Dept. of Dermatology Lubeck, Germany. All Rights Reserved.

load("QTLRelInput.RData")
ncmb.data <- read.table(file="NCMB_Sept2018.tab",header=T,row.names=1,sep="\t")

### calculate kinship matrix using DOQTL R package
K.G.mice <- kinship.probs(probs=array.G.mice,snps=snp.info)

ncmb.hglm <- ncmb.data

# construct new data-frame with hglm residual for each trait
# first 4 columns contains arbitary information about samples
for(i in 5:ncol(ncmb.data)){
	phen.idx <- i
	phen <- colnames(ncmb.data)[phen.idx]
	data.idx <- which(!is.na(ncmb.data[,phen.idx]));
	ncmb.phen <- ncmb.data[data.idx,c("Sex","Diet",phen)];
	K.G.phen <- K.G.mice[data.idx,data.idx];
	if(length(unique(ncmb.phen[,c(phen)])) == 2){
	 ncmb.phen$presid <- hglm(X = as.matrix(rep(1,nrow(ncmb.phen))), y = ncmb.phen[,c(phen)], Z = K.G.phen, family = binomial(link = logit))$resid
	}else{
	 ncmb.phen$presid <- hglm(X = as.matrix(rep(1,nrow(ncmb.phen))), y = ncmb.phen[,c(phen)], Z = K.G.phen)$resid
	}
	ncmb.hglm[,phen.idx] <- ncmb.phen$presid[match(rownames(ncmb.hglm),rownames(ncmb.phen))]
}
save(array.G.mice,K.G.mice,ncmb.hglm,snp.info, file="QTLRelInputHglm.RData")
