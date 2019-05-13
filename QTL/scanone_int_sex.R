#!/usr/bin/env Rscript
#function for fast interaction model
source("scanoneIntYask.R")

# load objects array.G.mice,K.G.mice,ncmb.hglm,snp.info
load("QTLRelInputHglm.RData")

# give phenotype in command line. This will help is 
# running each phenotype seperately in the cluster

args = commandArgs(trailingOnly=TRUE);
if (length(args)==0) {
  stop("At least one phenotype must be supplied", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  phen <- NULL;
  phen <- args[1];
}


#create array for interaction model

array.G.mice <- array.G.mice[,-1,]
array.G.mice.sex <- array.G.mice
n.snps <- dim(array.G.mice)[3]

for(i in c(1:n.snps)){
	sex.probs <- ncmb.hglm$Sex*array.G.mice[,,i]
	array.G.mice.sex[,,i] <- cbind(array.G.mice.sex[,,i],sex.probs)
	colnames(array.G.mice.sex[,,i])[c(4:6)] <- c("CASTD","MRLD","NZMD") 
}

#reserve name of the output .RDATA for each phenotype  
filename.res <- paste("/data/qtl_mice/artem2015/WGSQTL/QtlRel2018/Output/",phen,"_intSex_hglmResPerm.RData",sep="");
phen.idx <- which(colnames(ncmb.hglm) == phen);
ncmb.hglm[,phen.idx] <- as.numeric(as.matrix(ncmb.hglm[,phen.idx]));
data.idx <- which(!is.na(ncmb.hglm[,phen.idx]));
ncmb.phen <- ncmb.hglm[data.idx,c("Sex","Sex",phen)];

# array corresponding to non-missing samples for a phenotype
array.G.phen <- array.G.mice.sex[data.idx,,];

# free memory
rm(array.G.mice.sex,ncmb.hglm)
gc()

# run the regression model
res.phen <- scanone.int(pheno=ncmb.phen[,c(phen)],prdat=array.G.phen,snp=snp.info,nstr=4)

# permutate data for accessing significance 
perm.data<- scanone.int.perm(pheno=ncmb.phen[,c(phen)],prdat=array.G.phen,snp=snp.info,nstr=4,nperm=1000);

# save output
save(res.phen,perm.data,file=filename.res)



