# Copyright (c) 2016 by Yask Gupta, Dept. of Dermatology Lubeck, Germany. All Rights Reserved.
scanone.int.perm <- function(pheno=NULL,prdat=NULL,snp=NULL,nstr=NULL,nperm=NULL){
 perm.int <- matrix(NA,nrow=nrow(snp),ncol=nperm);
 P <- rep(Inf,nrow(snp));
 n <- length(pheno);
 nstr= nstr-1;
 for(i in 1:nperm){
  pheno.tmp <- pheno[sample(1:n,replace=FALSE)];
  for(k in 1:nrow(snp)){
   qr.null = qr(prdat[,c(1:nstr),k]);
   ss.null = sum(qr.resid(qr.null, pheno.tmp)^2)
   qr.int = qr(prdat[,,k]);
   ss.int = sum(qr.resid(qr.int, pheno.tmp)^2)
   P[k] = -length(pheno) * log(ss.int / ss.null)
  }
 perm.int[,i] <- P/(2*log(10));
 }
 return(perm.int)
}

scanone.int <- function(pheno=NULL,prdat=NULL,snp=NULL,nstr=NULL){
 P <- rep(Inf,nrow(snp));
 n <- length(pheno);
 nstr= nstr-1;
 for(k in 1:nrow(snp)){
  qr.null = qr(prdat[,c(1:nstr),k]);
  ss.null = sum(qr.resid(qr.null, pheno)^2)
  qr.int = qr(prdat[,,k]);
  ss.int = sum(qr.resid(qr.int, pheno)^2)
  P[k] = -length(pheno) * log(ss.int / ss.null)
 }
 lod=P/(2*log(10));
 return(list(snp=snp$SNP_ID,chr=snp$chr,pos=snp$Mb,LOD=lod));
}

