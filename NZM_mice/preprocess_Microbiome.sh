# This is script for pre-processing the micrbiota data
# The script uses third party usearch package and RDP 
# classifier extensively. Please install usearch and RDP
# before using the script.
# Copyright (c) 2016 by Yask Gupta, Dept. of Dermatology Lubeck, Germany. All Rights Reserved.

# Example for one sample

# Merge fastq paired end reads
usearch8.1.1861_i86linux32 -fastq_mergepairs 100_GTCAAGAG-CGTTAAGC_L001_R1_001.fastq -reverse 100_GTCAAGAG-CGTTAAGC_L001_R2_001.fastq  -fastq_minovlen 100 -fastq_minlen 200 -fastq_minmergelen 270 -fastq_maxmergelen 330 -fastqout 100_GTCAAGAG-CGTTAAGC.merged.fastq -relabel 100_GTCAAGAG-CGTTAAGC. -threads 40

# Filter low bad quality reads and convert fastq to fasta
usearch8.1.1861_i86linux32 -fastq_filter 100_GTCAAGAG-CGTTAAGC.merged.fastq -fastq_maxee 1.0 -relabel 100_GTCAAGAG-CGTTAAGC. -fastaout 100_GTCAAGAG-CGTTAAGC.filtered.fa -threads 40

# remove chimeric sequenced from the fasta file using rdp gold database
usearch8.1.1861_i86linux32 -uchime_ref 100_GTCAAGAG-CGTTAAGC.filtered.fa -db rdp_gold.fa -nonchimeras 100_GTCAAGAG-CGTTAAGC.nonchimera.fa -uchime_ref_minpctid 0.99 -strand plus -threads 24

########### run the same commands for all the samples and store all the sample files in nonchimera_fasta folder ####################

# For phylum to genus classification run rdp classifier. This will create nzm_hier.txt which can be
# used for downstream analysis.
java -Xmx100g -jar rdp_classifier_2.11/dist/classifier.jar -g 16srrna -f fixrank -o nzm_classified.txt -h nzm_hier.txt nonchimera_fasta/*

# For OTU level clustering run uclust which can be used for downstream analysis.
cat nonchimera_fasta/* > filtered.fa
usearch8.1.1861_i86linux32 -derep_fulllength filtered.fa -relabel Uniq -sizeout -fastaout uniques.fa -threads 40
usearch8.1.1861_i86linux32 -cluster_otus uniques.fa -minsize 2 -otus otus.fa -relabel Otu -threads 40
usearch8.1.1861_i86linux32 -usearch_global filtered.fa -db otus.fa -strand plus -id 0.97 -otutabout otutab.txt -biomout otutab.json -threads 40
