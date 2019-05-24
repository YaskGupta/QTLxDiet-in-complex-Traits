# This is script for pre-processing the micrbiota data
# The script uses third party PIPITS, RDP classifier and
# QIIME extensively. Please install PIPITS, RDP, QIIME 
# before using the script.
# Copyright (c) 2016 by Yask Gupta, Dept. of Dermatology Lubeck, Germany. All Rights Reserved.

# Server side scripts for PIPITS pipeline
INPATH='RawFastQ/'
mkdir pipitsdata
PIPIN='pipitsdata'

cp -r ${INPATH}/rawData/ ${PIPIN}/

pipits_getreadpairslist -i ${PIPIN}/rawData -o ${PIPIN}/readpairslist.txt

pipits_prep --FASTX-q 10 --FASTX-p 40 -t 48 -i ${PIPIN}/rawData -o ${PIPIN}/pipits_prep -l ${PIPIN}/readpairslist.txt

cp -r ${PIPIN}/pipits_prep ${INPATH}

pipits_process -t 48 -i ${PIPIN}/pipits_prep/prepped.fasta -o ${PIPIN}/pipits_process_RAW --Xmx 120G

cp -r ${PIPIN}/pipits_process_RAW ${INPATH}

pipits_funits -i ${PIPIN}/pipits_prep/prepped.fasta -o ${PIPIN}/pipits_funits -x ITS2 -t 48

cp -r ${PIPIN}/pipits_funits ${INPATH}

pipits_process -t 48 -i ${PIPIN}/pipits_funits/ITS.fasta -o ${PIPIN}/pipits_process_ITS --Xmx 120G

cp -r ${PIPIN}/pipits_process_ITS ${INPATH}

# The output files from the above will be in ${PIPIN}/pipits_process_ITS/ folder (otu_table.biom, otu_table.txt and repseqs.fasta)
# divide meta-data file for mycobiome samples stage-wise (Naive,Trans,Final)
# Please not that in Naive and Trans stage the Disease status represent Final disease state of Mice
#------------------------- In R shell one can do something like this --------------------------------------#
#													   #
## otudata <- read.table(file="otu_table.txt",comment.char="",sep="\t",header=T,row.names=1)		   #
## colnames(otudata) <- gsub("X",colnames(otudata))							   #
## metadata.naive <- read.table(file="Metadata_naive_myco.tab",comment.char="",sep="\t",header=T)	   #
## metadata.trans <- read.table(file="Metadata_trans_myco.tab",comment.char="",sep="\t",header=T)	   #
## metadata.final <- read.table(file="Metadata_final_myco.tab",comment.char="",sep="\t",header=T)	   #
# match the samples in metadata file and otu table							   #
## otudata.naive <- otudata[,match(as.chararcter(metadata.naive[,1]),colnames(otudata))]		   #
## otudata.trans <- otudata[,match(as.chararcter(metadata.trans[,1]),colnames(otudata))]		   #
## otudata.final <- otudata[,match(as.chararcter(metadata.final[,1]),colnames(otudata))]		   #
# remove singletons from the data									   #
## otudata.naive <- otudata.naive[which(apply(otudata.naive,1,function(x){length(which(x > 0))}) > 1),]	   #
## otudata.trans <- otudata.trans[which(apply(otudata.trans,1,function(x){length(which(x > 0))}) > 1),]	   #
## otudata.final <- otudata.final[which(apply(otudata.final,1,function(x){length(which(x > 0))}) > 1),]	   #
# output a text file for filter_fasta.py to retain OTU fasta sequences for each stage			   #
## write.table(rownames(otudata.naive),file="otus_naive.txt",sep="",quote=F,row.names=F,col.names=F)	   #
## write.table(rownames(otudata.trans),file="otus_trans.txt",sep="",quote=F,row.names=F,col.names=F)	   #
## write.table(rownames(otudata.final),file="otus_final.txt",sep="",quote=F,row.names=F,col.names=F)	   #
# output a tab delimited OTU abundance matrix for each stage to be used in QIIME			   #
## write.table(otudata.naive,file="otus_naive.tab",sep="\t",quote=F,row.names=T,col.names=T)	   	   #
## write.table(otudata.trans,file="otus_trans.tab",sep="\t",quote=F,row.names=T,col.names=T)		   #
## write.table(otudata.final,file="otus_final.tab",sep="\t",quote=F,row.names=T,col.names=T)		   #
#													   #
# In output tab file replace X with #OTUID manually for QIIME preprocessing				   #
#----------------------------------------------------------------------------------------------------------#


# convert tab file to biom for QIIME usage
biom convert -i otus_naive.tab -o otus_naive.biom --table-type="OTU table" --to-json
biom convert -i otus_trans.tab -o otus_trans.biom --table-type="OTU table" --to-json
biom convert -i otus_final.tab -o otus_final.biom --table-type="OTU table" --to-json

# QIIME workflow for processing the data for one (naive) of the stage

# calculate alpha-diversity 
## create alpha parameter text file "alpha_param.txt" with text "alpha_diversity:metrics shannon"
## and then run the following commands in qiime
alpha_rarefaction.py -i otus_naive.biom -m Metadata_naive_myco.tab -o arare -p alpha_params.txt -e 10000
## Test significance for a trait of interest (e.g. diet)
compare_alpha_diversity.py -i arare/alpha_div_collated/shannon.txt -m Metadata_naive_myco.tab -c Diet -d 10000 -o Diet_naive

# calculate beta-diversity
## create beta parameter text file "beta_param.txt" with text "beta_diversity:metrics jaccard"
## run the following commands in qiime
beta_diversity.py -i otus_naive.biom -m beta_param.txt -o bdiv_even10k -e 10000
## Test significance for a trait of interest (e.g. diet)
compare_categories.py --method adonis -i bdiv_even10k/jaccard_dm.txt -m Metadata_naive.tab -c Diet -o adonis_out -n 999





