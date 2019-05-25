<h1>Diet unmasks hidden susceptibility in complex traits</h1>
<p></p>
<p>This repository hosts scripts used for investigating Gene-Diet interation for advance intercross line in complex traits.</p>
<p>The advance intercross line consists of four-way cross between BxD2/TyJ, CAST/EiJ, NZM2410/J and MRL/MpJ mouse strain for 15-20 generations.</p>
<br>The raw data used in the analysis is deposited in:</br>
<br>ENA => Whole genome sequencing FastQ file for 3 out of 4 strain. Note CAST/EIJ whole genome sequencing data can be obtained from (https://www.sanger.ac.uk/science/data/mouse-genomes-project)</br>
<br>DRYAD database => Raw plink formated genotype data (Bed and Bim files), input files for HAPPY R package and VCF file for whole genome sequencing. (https://datadryad.org/review?doi=doi:10.5061/dryad.c8gc64n) </br> 
<br>NCBI SRA => Data for the analysis of NYM2410/J mice that includes raw FastQ files for RNAseq, Microbiome [V1-V2] and Mycobiome [ITS2]. (https://www.ncbi.nlm.nih.gov/sra/PRJNA543200) </br>
<p></p><br>The repository is broken down in 4 parts:</br>
WGS: Pipline for processing whole genome sequencing data.
QTL: Pipline for processing genotype data and performing QTL mapping in complex traits which includes additive, QTLxDiet and QTLxSex interaction models.
NZM_mice: Pipelines for processing and analysis of RNAseq, microbiome and mycobiome data.
