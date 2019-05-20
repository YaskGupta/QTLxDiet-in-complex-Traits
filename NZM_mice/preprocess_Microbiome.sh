# This is script for pre-processing the micrbiota data
# The script uses third party usearch package, RDP 
# classifier and QIIME extensively. Please install 
# usearch, RDP, QIIME before using the script.
# Copyright (c) 2016 by Yask Gupta, Dept. of Dermatology Lubeck, Germany. All Rights Reserved.

# Example for one sample

# Merge fastq paired end reads
usearch -fastq_mergepairs 100_GTCAAGAG-CGTTAAGC_L001_R1_001.fastq -reverse 100_GTCAAGAG-CGTTAAGC_L001_R2_001.fastq  -fastq_minovlen 100 -fastq_minlen 200 -fastq_minmergelen 270 -fastq_maxmergelen 330 -fastqout 100_GTCAAGAG-CGTTAAGC.merged.fastq -relabel 100_GTCAAGAG-CGTTAAGC. -threads 40

# Filter low bad quality reads and convert fastq to fasta
usearch -fastq_filter 100_GTCAAGAG-CGTTAAGC.merged.fastq -fastq_maxee 1.0 -relabel 100_GTCAAGAG-CGTTAAGC. -fastaout 100_GTCAAGAG-CGTTAAGC.filtered.fa -threads 40

# remove chimeric sequenced from the fasta file using rdp gold database
usearch -uchime_ref 100_GTCAAGAG-CGTTAAGC.filtered.fa -db rdp_gold.fa -nonchimeras 100_GTCAAGAG-CGTTAAGC.nonchimera.fa -uchime_ref_minpctid 0.99 -strand plus -threads 24

########### run the same commands for all the samples and store all the sample files in nonchimera_fasta folder ####################

# For phylum to genus classification run rdp classifier. This will create nzm_hier.txt which can be
# used for downstream analysis.
java -Xmx100g -jar rdp_classifier_2.11/dist/classifier.jar -g 16srrna -f fixrank -o nzm_classified.txt -h nzm_hier.txt nonchimera_fasta/*

# For OTU level clustering run uclust which can be used for downstream analysis.
cat nonchimera_fasta/* > filtered.fa
usearch -derep_fulllength filtered.fa -relabel Uniq -sizeout -fastaout uniques.fa -threads 40
usearch -cluster_otus uniques.fa -minsize 2 -otus otus.fa -relabel Otu -threads 40
usearch -usearch_global filtered.fa -db otus.fa -strand plus -id 0.97 -otutabout otutab.txt -biomout otutab.biom -threads 40

# divide meta-data file for microbiome samples stage-wise (Naive,Trans,Final)

#------------------------- In R shell one can do something like this --------------------------------------#
#													   #
## otudata <- read.table(file="otutab.txt",comment.char="",sep="\t",header=T,row.names=1)		   #
## colnames(otudata) <- gsub("X",colnames(otudata))							   #
## metadata.naive <- read.table(file="Metadata_naive.tab",comment.char="",sep="\t",header=T)		   #
## metadata.trans <- read.table(file="Metadata_trans.tab",comment.char="",sep="\t",header=T)		   #
## metadata.final <- read.table(file="Metadata_final.tab",comment.char="",sep="\t",header=T)		   #
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
## write.table(otudata.naive,file="otus_naive.tab",sep="\t",quote=F,row.names=T,col.names=T)		   #
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

# filter fasta file for each stage
filter_fasta.py -f otus.fa -o otus_naive.fa -s otus_naive.txt

# assign taxanomy with greengene database (default)
assign_taxonomy.py -i otus_naive.fa

# generate phylogenetic tree using fasttree
align_seqs.py -i otus_naive.fa -o pynast_naive
filter_alignment.py -i pynast_naive/otus_naive_aligned.fasta -o pynast_naive
make_phylogeny.py -i pynast_naive/otus_naive_aligned_pfiltered.fasta -o pynast_naive/otu_naive.tre

# calculate alpha-diversity 
## create alpha parameter text file "alpha_param.txt" with text "alpha_diversity:metrics shannon,chao1"
## and then run the following commands in qiime
alpha_rarefaction.py -i otus_naive.biom -m Metadata_naive.tab -o arare -p alpha_params.txt -t pynast_naive/otu_naive.tre -e 20000
## Test significance for a trait of interest (e.g. diet)
compare_alpha_diversity.py -i arare/alpha_div_collated/chao1.txt -m Metadata_naive.tab -c Diet -d 20000 -o Diet_naive

# calculate beta-diversity 
## run the following commands in qiime
beta_diversity_through_plots.py -i otus_naive.biom -m Metadata_naive.tab -o bdiv_even20k -t pynast_naive/otu_naive.tre -e 20000
## Test significance for a trait of interest (e.g. diet)
compare_categories.py --method adonis -i bdiv_even20k/unweighted_unifrac_dm.txt -m Metadata_naive.tab -c Diet -o adonis_out -n 999

