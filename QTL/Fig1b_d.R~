#!/usr/bin/env Rscript
# The following script is for producing manhattan plot Fig1b & Fig1d in the
# Note: Summary results are part of the repository but .RDATA files should be 
# generated using previous R scripts
# Copyright (c) 2016 by Yask Gupta, Dept. of Dermatology Lubeck, Germany. All Rights Reserved.

###### Fig 1b , manhattan plot #########

# Summary results for the additive model
add.phen <- read.table(file="../Data/Summarize_result_add.tab",header=T)
add.phen.g <- unique(c(as.matrix(subset(add.phen,Type == "G",PHEN))))

#Load .Rdata for files for Additive QTL mapping
load(paste("/data/qtl_mice/artem2015/WGSQTL/QtlRel2018/Output/Color_addSex_addDiet_hglmResPerm.RData",sep=""));
Tmp.phen<- data.frame(chr=res.phen$chr,dist=res.phen$dist,y=res.phen$p/(2*(log(10))));
rownames(Tmp.phen) = res.phen$snp
Tmp.phen$chr<- reorder(Tmp.phen$chr)
Tmp.phen<- Tmp.phen[order(Tmp.phen$chr,Tmp.phen$dist),]
res.add.lod <- Tmp.phen
colnames(res.add.lod)[3] = "Color"
for(k in 1:length(add.phen.g)){
	phen <- add.phen.g[k]
	load(paste("/data/qtl_mice/artem2015/WGSQTL/QtlRel2018/Output/",phen,"_addSex_addDiet_hglmResPerm.RData",sep=""));
	cat(paste("Reading",phen,"\n",sep=" "))
	Tmp.phen<- data.frame(chr=res.phen$chr,dist=res.phen$dist,y=res.phen$p/(2*(log(10))));
	rownames(Tmp.phen) = res.phen$snp
	Tmp.phen$chr<- reorder(Tmp.phen$chr)
	Tmp.phen<- Tmp.phen[order(Tmp.phen$chr,Tmp.phen$dist),]
	res.add.lod <- cbind(res.add.lod,Tmp.phen$y)
	colnames(res.add.lod)[3+k] = phen
}

# Summary results for the Int Diet model
int.diet.phen <- read.table(file="../Data/Summarize_result_int_diet.tab",header=T)
int.diet.phen.g <- unique(c(as.matrix(subset(int.diet.phen,Type == "G",PHEN))))

#Load .Rdata for files for Int Diet QTL mapping
load(paste("Output/",int.diet.phen.g[1],"_addSex_intDiet2_hglmResPerm.RData",sep=""));
Tmp.phen<- data.frame(chr=res.phen$chr,dist=res.phen$pos,y=res.phen$LOD);
rownames(Tmp.phen) = res.phen$snp
Tmp.phen$chr<- reorder(Tmp.phen$chr)
Tmp.phen<- Tmp.phen[order(Tmp.phen$chr,Tmp.phen$dist),]
res.int.diet.lod <- Tmp.phen
colnames(res.int.diet.lod)[3] = int.diet.phen.g[1]
for(k in 2:length(int.diet.phen.g)){
	phen <- int.diet.phen.g[k]
	load(paste("Output/",phen,"_addSex_intDiet2_hglmResPerm.RData",sep=""));
	cat(paste("Reading",phen,"\n",sep=" "))
	Tmp.phen<- data.frame(chr=res.phen$chr,dist=res.phen$pos,y=res.phen$LOD);
	rownames(Tmp.phen) = res.phen$snp
	Tmp.phen$chr<- reorder(Tmp.phen$chr)
	Tmp.phen<- Tmp.phen[order(Tmp.phen$chr,Tmp.phen$dist),]
	res.int.diet.lod <- cbind(res.int.diet.lod,Tmp.phen$y)
	colnames(res.int.diet.lod)[2+k] = phen
}

# Summary results for the Int Sex model
int.sex.phen <- read.table(file="../Data/Summarize_result_int_sex.tab",header=T)
int.sex.phen.g <- unique(c(as.matrix(subset(int.sex.phen,Type == "G",PHEN))))

#Load .Rdata for files for Int Sex QTL mapping
load(paste("/data/qtl_mice/artem2015/WGSQTL/QtlRel2018/Output/",int.sex.phen.g[1],"_intSex_hglmResPerm.RData",sep=""));
Tmp.phen<- data.frame(chr=res.phen$chr,dist=res.phen$pos,y=res.phen$LOD);
rownames(Tmp.phen) = res.phen$snp
Tmp.phen$chr<- reorder(Tmp.phen$chr)
Tmp.phen<- Tmp.phen[order(Tmp.phen$chr,Tmp.phen$dist),]
res.int.sex.lod <- Tmp.phen
colnames(res.int.sex.lod)[3] = int.sex.phen.g[1]
for(k in 2:length(int.sex.phen.g)){
	phen <- int.sex.phen.g[k]
	load(paste("/data/qtl_mice/artem2015/WGSQTL/QtlRel2018/Output/",phen,"_intSex_hglmResPerm.RData",sep=""));
	cat(paste("Reading",phen,"\n",sep=" "))
	Tmp.phen<- data.frame(chr=res.phen$chr,dist=res.phen$pos,y=res.phen$LOD);
	rownames(Tmp.phen) = res.phen$snp
	Tmp.phen$chr<- reorder(Tmp.phen$chr)
	Tmp.phen<- Tmp.phen[order(Tmp.phen$chr,Tmp.phen$dist),]
	res.int.sex.lod <- cbind(res.int.sex.lod,Tmp.phen$y)
	colnames(res.int.sex.lod)[2+k] = phen
}


res.comb <- cbind(res.add.lod,res.int.diet.lod[,-c(1,2)],res.int.sex.lod[,-c(1,2)])

Tmp.cumsum.dist = cumsum(Tmp.phen$dist)

add.phen.g = subset(add.phen,Type == "G")
int.diet.phen.g = subset(int.diet.phen,Type == "G")
int.sex.phen.g = subset(int.sex.phen,Type == "G")

write.table(rbind(add.phen.g,rbind(int.diet.phen.g,int.sex.phen.g),file="Summarize_result_genome_wide.tab",col.names=T,row.names=F,sep="\t",quote=F)
res.g <- read.table(file="Output/Summarize_result_genome_wide.tab",header=T)
res.g$cdist = Tmp.cumsum.dist[match(res.g$SNP,rownames(Tmp.phen))]
res.g$Col = ifelse(res.g$MODEL == "ADD","blue",ifelse(res.g$MODEL == "INTD","gold","violetred"))

res.unmat <- c(as.matrix(res.comb[,-c(1,2)]))
unorm.idx = which(res.unmat > 16)
res.unmat[unorm.idx] <- 16+(20-16)*((res.unmat[unorm.idx]-min(res.unmat[unorm.idx]))/(max(res.unmat[unorm.idx])-min(res.unmat[unorm.idx])))
res.norm = matrix(res.unmat,nrow=nrow(res.comb),ncol=ncol(res.comb)-2)
rownames(res.norm) = rownames(res.comb)
colnames(res.norm) = colnames(res.comb)[-c(1,2)]

unorm.idx = which(res.g$LOD > 16)
res.g$LOD[unorm.idx] <- 16+(20-16)*((res.g$LOD[unorm.idx]-min(res.g$LOD[unorm.idx]))/(max(res.g$LOD[unorm.idx])-min(res.g$LOD[unorm.idx])))

chrPos = data.frame(chr=c(1:19,"X"),pos.med = rep(0,length(c(1:19,"X"))));k=0;for(i in c(1:19,"X")){k=k+1;chrPos$pos.med[k] = mean(Tmp.cumsum.dist[which(Tmp.phen$chr == i)])}

res.g$CHR = reorder(res.g$CHR)
res.g<- res.g[order(res.g$CHR,res.g$Pos),]

res.chr <- unique(as.character(res.g$CHR))
res.txt <- res.g
for(i in 1:length(res.chr)){
	res.tmp <- subset(res.g,CHR==res.chr[i])
	if(nrow(res.tmp) > 1){
		res.idx <- which(res.g$CHR==res.chr[i])
		res.tmp <- res.tmp[order(res.tmp$LOD),]
		n.r <- length(res.idx)*0.5
		n.s <- seq(from=0,to=n.r,by=0.5)
		n.s <- n.s[-length(n.s)]
		res.tmp$LOD <- res.tmp$LOD+n.s
		res.txt[res.idx,] = res.tmp
	}
}


for(i in 1:ncol(res.norm)){
	points(Tmp.cumsum.dist,as.numeric(res.norm[,i]),col=ifelse(as.integer(Tmp.phen$chr)%%2 == 0,8,1),pch=".")
}

# make a manhattan plot
svg("Plot_Manhattan_Gmice.svg",width=7,height=3)
plot(Tmp.cumsum.dist,rep(0,nrow(Tmp.phen)),xlab="Chromosomes", ylab="LOD", col=ifelse(as.integer(Tmp.phen$chr)%%2 == 0,8,1),ylim=c(0,22),type="l",pch=".",cex.lab=0.8,axes = FALSE,mgp = c(2, 1, 0))
axis(side=1,at=chrPos$pos.med,labels=chrPos$chr,cex.axis=0.3,cex.lab=0.3,font = 2,mgp = c(1, 1, 0))
axis(side=2,at=c(0,5,10,15,22),labels=c(0,5,10,15,40),cex.axis=0.5,cex.lab=0.5,mgp = c(1, 1, 0))
text(x=res.txt$cdist, y=res.txt$LOD+1.5, labels = res.txt$PHEN,cex = 0.3,adj=1,font=2)

points(res.g$cdist,res.g$LOD,pch=21,col="black",bg=res.g$Col,cex=0.5)
dev.off()

###### Fig 1b , coefficient plot for Weight 2 months #########

phen="wt2Month"

load(paste("/data/qtl_mice/artem2015/WGSQTL/QtlRel2018/Output/",phen,"_addSex_addDiet_hglmResPerm.RData",sep=""));
Tmp.wt.add<- data.frame(chr=res.phen$chr,dist=res.phen$dist,y=res.phen$p/(2*(log(10))));
rownames(Tmp.wt.add) = res.phen$snp
Tmp.wt.add$chr<- reorder(Tmp.wt.add$chr)
Tmp.wt.add<- Tmp.wt.add[order(Tmp.wt.add$chr,Tmp.wt.add$dist),]
thresh.wt.add <- quantile(apply(perm.data,2,max),p=0.95)/(2*(log(10)))

load(paste("/data/qtl_mice/artem2015/WGSQTL/QtlRel2018/Output/",phen,"_addSex_intDiet2_hglmResPerm.RData",sep=""));
Tmp.wt.intd<- data.frame(chr=res.phen$chr,dist=res.phen$pos,y=res.phen$LOD);
rownames(Tmp.wt.intd) = res.phen$snp
Tmp.wt.intd$chr<- reorder(Tmp.wt.intd$chr)
Tmp.wt.intd<- Tmp.wt.intd[order(Tmp.wt.intd$chr,Tmp.wt.intd$dist),]
thresh.wt.intd <- quantile(apply(perm.data,2,max),p=0.95)

Tmp.wt.add.g <- Tmp.wt.add[Tmp.wt.add$chr %in% unique(res.g$CHR[which(res.g$PHEN == "WT(2M)")]),]
Tmp.wt.intd.g <- Tmp.wt.intd[Tmp.wt.intd$chr %in% unique(res.g$CHR[which(res.g$PHEN == "WT(2M)")]),]

Tmp.wt.add.g$chr <- as.numeric(Tmp.wt.add.g$chr)
Tmp.wt.intd.g$chr <- as.numeric(Tmp.wt.intd.g$chr)

Tmp.cumsum.dist.wt = cumsum(Tmp.wt.add.g$dist)

res.g.wt.add <- subset(res.g,PHEN=="WT(2M)" & MODEL=="ADD")
res.g.wt.intd <- subset(res.g,PHEN=="WT(2M)" & MODEL=="INTD")

res.g.wt.add$cdist = Tmp.cumsum.dist[match(res.g.wt.add$SNP,rownames(Tmp.wt.add.g))]
res.g.wt.intd$cdist = Tmp.cumsum.dist[match(res.g.wt.intd$SNP,rownames(Tmp.wt.intd.g))]

chrPos.wt = data.frame(chr=unique(Tmp.wt.add.g$chr),pos.med = rep(0,length(unique(Tmp.wt.add.g$chr))));k=0;for(i in unique(Tmp.wt.add.g$chr)){k=k+1;chrPos.wt$pos.med[k] = mean(Tmp.cumsum.dist.wt[which(Tmp.wt.add.g$chr == i)])}

jpeg("plot_wt2.jpeg",width=7,height=3,units="in",res=400)
par(mfrow=c(1,2))
plot(Tmp.cumsum.dist.wt,Tmp.wt.add.g$y,xlab="Chromosomes", ylab="LOD", col=ifelse(as.integer(as.factor(Tmp.wt.add.g$chr))%%2 == 0,8,1),ylim=c(0,11),type="p",pch=20,cex=0.2,cex.axis=0.5,cex.lab=0.5,xaxt = 'n')
axis(side=1,at=chrPos.wt$pos.med,labels=chrPos.wt$chr,cex.axis=0.5,cex.lab=0.5,mgp = c(1, 1, 0))
abline(h=thresh.wt.add,lty=2,col="red")
plot(Tmp.cumsum.dist.wt,Tmp.wt.intd.g$y,xlab="Chromosomes", ylab="LOD", col=ifelse(as.integer(as.factor(Tmp.wt.intd.g$chr))%%2 == 0,8,1),ylim=c(0,11),type="p",pch=20,cex=0.2,cex.axis=0.5,cex.lab=0.5,xaxt = 'n')
axis(side=1,at=chrPos.wt$pos.med,labels=chrPos.wt$chr,cex.axis=0.5,cex.lab=0.5,mgp = c(1, 1, 0))
abline(h=thresh.wt.intd,lty=2,col="red")
dev.off()





