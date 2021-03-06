<h1>Diet unmasks hidden susceptibility in complex traits</h1>
<p></p>
<p>This repository hosts scripts used for investigating Gene-Diet interation for mouse advance intercross line (AIL) in complex traits.The advance intercross line consists of four-way cross between BxD2/TyJ, CAST/EiJ, NZM2410/J and MRL/MpJ mouse strain for 15-20 generations. It also hosts script for analyzing multi-omics data (RNAseq, microbiome and mycobiome) from NZM2410/J which plays critical role in these gene-diet interactions</p>
<br><b> RAW DATA DEPOSITS </b>
<p><br>The raw data used in the analysis is deposited in:
<br>ENA: Whole genome sequencing FastQ file for 3 out of 4 strain. Note CAST/EIJ whole genome sequencing data can be obtained from https://www.sanger.ac.uk/science/data/mouse-genomes-project
<br>DRYAD database: Raw plink formated genotype data (BED and BIM files), input files for HAPPY R package and VCF file for whole genome sequencing.https://datadryad.org/review?doi=doi:10.5061/dryad.c8gc64n
<br>NCBI SRA: Data for the analysis of NYM2410/J mice that includes raw FastQ files for RNAseq, Microbiome [V1-V2] and Mycobiome [ITS2].https://www.ncbi.nlm.nih.gov/sra/PRJNA543200</p>
<br><b> SCRIPTS TO PROCESS RAW DATA </b>
<p><br>The repository is broken down in 4 parts:
<br>WGS: Pipline for processing whole genome sequencing data. 
<br>QTL: Pipline for processing genotype data and performing QTL mapping in complex traits which includes additive, QTLxDiet and QTLxSex interaction models.
<br>NZM_mice: Pipelines for processing and analysis of RNAseq, microbiome and mycobiome data.
<br>Data: Folder containing the trait values for QTL mapping, RNAseq, microbiome and mycobiome data.</p>
