# The following script summarizes the results of all the data.
# The script uses output from scanone_add_diet.R
# For every phenotype the data is stored phenotype_addSex_addDiet_hglmResPerm.RData
# Copyright (c) 2016 by Yask Gupta, Dept. of Dermatology Lubeck, Germany. All Rights Reserved.

library(QTLRel)
library(hglm)
getOutput <- function(phen=NULL,lc=NULL,phen.res=NULL,chrs=NULL,type=NULL){
	lc.frame <- as.data.frame(lc)
	out.mat <- matrix(NA,nrow=length(chrs),ncol=7)
	colnames(out.mat) <- c("PHEN","CHR","SNP","Pos","CI","LOD","Type")
	for(i in 1:length(chrs)){
	 lc.chr <- subset(lc.frame,chr==chrs[i])
	 ci.start <- round(min(lc.chr$lower),digit=3)
	 ci.end <- round(max(lc.chr$upper),digit=3)
	 ci.tot <- paste(ci.start,"-",ci.end,sep="")
	 Tmp.chrs <- subset(phen.res,chr==chrs[i])
	 peak.snp <- rownames(Tmp.chrs)[which.max(Tmp.chrs$y)]
	 lodscore <- round(max(Tmp.chrs$y),digit=3)
	 dist.chr <- round(Tmp.chrs[which.max(Tmp.chrs$y),]$dist,digit=3)
	 out.mat[i,] <- c(phen,chrs[i],peak.snp,dist.chr,ci.tot,lodscore,type)
	}
	return(out.mat)
}


get.hglm.out <- function(phen=NULL,ncmb.data=NULL,array.G.data=NULL,Kin=NULL,out.data.mat=NULL){
	var.mat <- matrix(0,nrow=nrow(out.data.mat),ncol=5)
	colnames(var.mat) <- c("VarSex","VarDiet","VarGeno","VarGenoxDiet","VarGenoxSex")
	phen.idx <- which(colnames(ncmb.data) == phen);
	ncmb.data[,phen.idx] <- as.numeric(as.matrix(ncmb.data[,phen.idx]));
	data.idx <- which(!is.na(ncmb.data[,phen.idx]));
	ncmb.phen <- ncmb.data[data.idx,c("Sex","Diet",phen)];
	K.G.phen <- Kin[data.idx,data.idx];
	if(length(unique(ncmb.phen[,c(phen)])) == 2){
	 phen.resid <- hglm(X = as.matrix(rep(1,nrow(ncmb.phen))), y = ncmb.phen[,c(phen)], Z = K.G.phen, family = binomial(link = logit))$resid
	}else{
	 phen.resid <- hglm(X = as.matrix(rep(1,nrow(ncmb.phen))), y = ncmb.phen[,c(phen)], Z = K.G.phen)$resid
	}
	array.G.phen <- array.G.data[data.idx,-1,];
	for(i in 1:length(out.data.mat$SNP)){
	 geno <- as.matrix(array.G.phen[,,which(dimnames(array.G.phen)[[3]] == out.data.mat$SNP[i])])
	 geno.diet <- geno*ncmb.phen$Diet
	 geno.sex <- geno*ncmb.phen$Sex
	 af <- anova(lm(phen.resid ~ ncmb.phen$Sex + ncmb.phen$Diet + geno + geno.diet + geno.sex),test="Chisq")
	 afss <- af$"Sum Sq"
	 var.mat[i,] <- round((afss/sum(afss)*100)[1:5],digit=3)
	}
	out.var.data <- as.data.frame(cbind(out.data.mat,var.mat))
	return(out.var.data)

}

load("QTLRelInput.RData")
final.table.add <- as.data.frame(matrix(nrow=1,ncol=12))
colnames(final.table.add) <- c("PHEN","CHR","SNP","Pos","CI","LOD","Type","VarSex","VarDiet","VarGeno","VarGenoxDiet","VarGenoxSex")
for(k in 5:ncol(ncmball)){
	phen=colnames(ncmball)[k]
	filename.res <- paste("Output/",phen,"_addSex_addDiet_hglmResPerm.RData",sep="");
	load(filename.res)
	cat(paste("Reading",phen,"\n",sep=" "))
	out.data <- NULL
	lc.phen.g=NULL
	lc.phen.l=NULL
	chr.g <- NULL
	chr.l <- NULL

	thres.g <- quantile(apply(perm.data,2,max),p=0.95)/(2*log(10));
	thres.l <- quantile(apply(perm.data,2,max),p=0.90)/(2*log(10));
	Tmp.phen<- data.frame(chr=res.phen$chr,dist=res.phen$dist,y=res.phen$p/(2*log(10)));
	Tmp.phen$chr<- reorder(Tmp.phen$chr)
	Tmp.phen<- Tmp.phen[order(Tmp.phen$chr,Tmp.phen$dist),]
	chr.g <- unique(Tmp.phen$chr[which(Tmp.phen$y > thres.g)])
	chr.l <- unique(Tmp.phen$chr[which(Tmp.phen$y > thres.l)])
	if(length(chr.g) == 0) chr.g <- NULL
	if(length(chr.l) == 0) chr.l <- NULL

	if(!is.null(chr.l) & !is.null(chr.g)){
		 lc.phen.g <- lodci(Tmp.phen,cv=thres.g,lod=1.5,drop=1.5)
		 out.data.g <- as.data.frame(getOutput(phen=phen,lc=lc.phen.g,phen.res=Tmp.phen,chrs=chr.g,type="G"))
		 lc.phen.l <- lodci(Tmp.phen,cv=thres.l,lod=1.5,drop=1.5)
		 chr.l <- chr.l[!chr.l %in% chr.g]
		 if(length(chr.l) > 0){
		 	out.data.l <- as.data.frame(getOutput(phen=phen,lc=lc.phen.l,phen.res=Tmp.phen,chrs=chr.l,type="L"))
			out.data <- rbind(out.data.g,out.data.l)
		 }else{
			out.data <- out.data.g
		 }		 
		 out.data.var <- get.hglm.out(phen=phen,ncmb.data=ncmball,array.G.data=array.G.mice,Kin=K.G.mice,out.data.mat=out.data)
	}else if(is.null(chr.l) & !is.null(chr.g)){
		 lc.phen.g <- lodci(Tmp.phen,cv=thres.g,lod=1.5,drop=1.5)
		 out.data <- as.data.frame(getOutput(phen=phen,lc=lc.phen.g,phen.res=Tmp.phen,chrs=chr.g,type="G"))
		 out.data.var <- get.hglm.out(phen=phen,ncmb.data=ncmball,array.G.data=array.G.mice,Kin=K.G.mice,out.data.mat=out.data)
	}else if(!is.null(chr.l) & is.null(chr.g)){
		 lc.phen.l <- lodci(Tmp.phen,cv=thres.l,lod=1.5,drop=1.5)
		 out.data <- as.data.frame(getOutput(phen=phen,lc=lc.phen.l,phen.res=Tmp.phen,chrs=chr.l,type="L"))
		 out.data.var <- get.hglm.out(phen=phen,ncmb.data=ncmball,array.G.data=array.G.mice,Kin=K.G.mice,out.data.mat=out.data)
	}else{
		 out.data.var <- c(phen,rep(NA,11))
	}
	final.table.add <- rbind(final.table.add,out.data.var)
	cat(paste(phen,"has been summarized","\n",sep=" "));

}
save(final.table.add,file="Summarize_result_add.RData")
write.table(final.table.int.diet,file="Summarize_results_int_diet.tab",col.names=T,row.names=F,sep="\t",quote=F)


################################# hglm for color #########################

get.hglm.out.color <- function(phen=NULL,ncmb.data=NULL,array.G.data=NULL,Kin=NULL,out.data.mat=NULL){
	var.mat <- matrix(0,nrow=nrow(out.data.mat),ncol=5)
	colnames(var.mat) <- c("VarSex","VarDiet","VarGeno","VarGenoxDiet","VarGenoxSex")
	phen.idx <- which(colnames(ncmb.data) == phen);
	ncmb.data[,phen.idx] <- as.numeric(as.factor(ncmb.data[,phen.idx]));
	data.idx <- which(!is.na(ncmb.data[,phen.idx]));
	ncmb.phen <- ncmb.data[data.idx,c("Sex","Diet",phen)];
	K.G.phen <- Kin[data.idx,data.idx];
	if(length(unique(ncmb.phen[,c(phen)])) == 2){
	 phen.resid <- hglm(X = as.matrix(rep(1,nrow(ncmb.phen))), y = ncmb.phen[,c(phen)], Z = K.G.phen, family = binomial(link = logit))$resid
	}else{
	 phen.resid <- hglm(X = as.matrix(rep(1,nrow(ncmb.phen))), y = ncmb.phen[,c(phen)], Z = K.G.phen)$resid
	}
	array.G.phen <- array.G.data[data.idx,-1,];
	for(i in 1:length(out.data.mat$SNP)){
	 geno <- as.matrix(array.G.phen[,,which(dimnames(array.G.phen)[[3]] == out.data.mat$SNP[i])])
	 geno.diet <- geno*ncmb.phen$Diet
	 geno.sex <- geno*ncmb.phen$Sex
	 af <- anova(lm(phen.resid ~ ncmb.phen$Sex + ncmb.phen$Diet + geno + geno.diet + geno.sex),test="Chisq")
	 afss <- af$"Sum Sq"
	 var.mat[i,] <- round((afss/sum(afss)*100)[1:5],digit=3)
	}
	out.var.data <- as.data.frame(cbind(out.data.mat,var.mat))
	return(out.var.data)

}

