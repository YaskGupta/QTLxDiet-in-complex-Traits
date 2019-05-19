# The following script describes the pre-processing RNAseq data. Please make sure that tophat2
# cufflinks, cuffmerge and subread (featureCounts) is install in the system.
# Copyright (c) 2016 by Yask Gupta, Dept. of Dermatology Lubeck, Germany. All Rights Reserved.

# This is an example for one FastQ file.
# Align paired-end data with tophat2
tophat2 --GTF Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genes.gtf --library-type fr-firststrand --num-threads 48 --output-dir Tophat_out/S10_S10 Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genome FastQ/S10_S10_R1_001.fastq.gz FastQ/S10_S10_R2_001.fastq.gz

# Cufflink command
cufflinks --library-type fr-firststrand --GTF-guide Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genes.gtf --library-type fr-firststrand --num-threads 48 --max-bundle-length 6000000 -o Tophat_out/S10_S10 Tophat_out/S10_S10/accepted_hits.bam

# create a .txt file assemeble.txt containing transcripts.gtf
# from all the sample and run cuffmerge to create merged.gtf

cuffmerge -p 8 -g Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genes.gtf -s Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genome.fa assemblies.txt

# use featureCounts from subread package to count the fragments
# for merged.gtf file.

featureCounts -p -t exon -g gene_id -a merged.gtf -o NZM_summary.txt *.bam

# Rest of the data processing is done by DESeq2 R package.

 
