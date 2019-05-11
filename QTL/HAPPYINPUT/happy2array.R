# The following R scipt takes the output of make_all_happy_files.pl and generates 3D array
# which can be used in DOQTL or QTLRel. The array contains prosterior probabilities for each 
# SNP across 4 founders for all the samples.
# Before running the script please make sure that you are in folder containing the HAPPY input
# files chromosome-wise.
# For this script to run install R packages HAPPY and abind  
# Copyright (c) 2016 by Yask Gupta, Dept. of Dermatology Lubeck, Germany. All Rights Reserved.
library(happy.hbrem)
library(abind)
chrm = c(1:19,"X")
dfile=paste("chr",chrm[1],".AAIL.data",sep="")
afile=paste("chr",chrm[1],".AAIL.alleles",sep="")
h=NULL
h = happy(datafile=dfile,allelesfile=afile,generations=20)
array.G.mice = array(0,c(length(h$subjects),length(h$strains),length(h$markers)-1),dimnames=list(h$subjects,h$strains,h$markers[-length(h$markers)]))
for(i in c(1:(length(h$markers)-1))){model.probs[,,i] = hdesign(h,marker=h$markers[i]);}

for(i in c(2:length(chrm))){
	dfile=paste("chr",chrm[i],".AAIL.data",sep="")
	afile=paste("chr",chrm[i],".AAIL.alleles",sep="")
	h=NULL
	h = happy(datafile=dfile,allelesfile=afile,generations=20)
	model.probs.temp = array(0,c(length(h$subjects),length(h$strains),length(h$markers)-1),dimnames=list(h$subjects,h$strains,h$markers[-length(h$markers)]))
	for(j in c(1:(length(h$markers)-1))){model.probs.temp[,,j] = hdesign(h,marker=h$markers[j]);}
	array.G.mice = abind(model.probs,model.probs.temp, along = 3)
	gc()
}

save(array.G.mice, file="QTLRelInput.RData")
