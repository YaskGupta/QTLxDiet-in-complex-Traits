#!/usr/bin/env Rscript
# Copyright (c) 2016 by Yask Gupta, Dept. of Dermatology Lubeck, Germany. All Rights Reserved.

library(DOQTL)
library(QTLRel)
load("QTLRelInputHglm.RData")

args = commandArgs(trailingOnly=TRUE);
if (length(args)==0) {
  stop("At least one phenotype must be supplied", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  phen <- NULL;
  phen <- args[1];
}

filename.res <- paste("Output/",phen,"_addSex_addDiet_hglmResPerm.RData",sep="");
phen.idx <- which(colnames(ncmb.hglm) == phen);
ncmb.hglm[,phen.idx] <- as.numeric(as.matrix(ncmb.hglm[,phen.idx]));
data.idx <- which(!is.na(ncmb.hglm[,phen.idx]));
ncmb.phen <- ncmb.hglm[data.idx,c("Sex","Diet",phen)];
array.G.phen <- array.G.mice[data.idx,,];


prdat.phen <- list(pr = array.G.phen, chr = snp.info[,2], dist = snp.info[,3],snp = snp.info[,1]);
class(prdat.phen) <- c(class(prdat.phen), "addEff");

rm(array.G.mice,snp.info,ncmb.hglm)
gc()

res.phen <- scanOne(y = ncmb.phen[,c(phen)], x = ncmb.phen[c("Sex","Diet")], prdat = prdat.phen, vc = NULL,intcovar = NULL,numGeno = TRUE, test = "None");

ntimes <- 1000
perm.data <- matrix(0,nrow=length(res.phen$snp),ncol=ntimes);
nn <- nrow(ncmb.phen);

for(n in 1:ntimes){
	idx<- sample(1:nn,replace=FALSE)
	perm.data[,n]<- scanOne(y=ncmb.phen[idx,c(phen)], x = ncmb.phen[c("Sex","Diet")],prdat=prdat.phen, vc = NULL,intcovar = NULL,numGeno = TRUE, test = "None")$p;
}

save(res.phen,perm.data,file=filename.res)
