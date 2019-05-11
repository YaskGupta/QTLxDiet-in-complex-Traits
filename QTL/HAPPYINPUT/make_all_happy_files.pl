use strict;
use warnings;
# The following perl script creats HAPPY input file from founder alleles and sample geneotype which is splitted chromosome-wise.
# The HAPPY R program is used for creating 3D array containing prosterior probabilitlies and imputation for each SNP across the samples.
# The script can be used as perl make_all_happy.pl 1 
# Copyright (c) 2016 by Yask Gupta, Dept. of Dermatology Lubeck, Germany. All Rights Reserved.

my %founder;
my $i=0;
my @happyM;
# get chromsome number from commandline.
my $chr = $ARGV[0];
my $allelefile="founder_allele_".$chr.".txt";
#make marker file
open FOUNDER,$allelefile or die "file for $chr donot exists!!\n";
my @f_data = <FOUNDER>;
close FOUNDER;
my $num_snps = scalar(@f_data)-1;
push(@happyM,"markers"."\t".$num_snps."\t"."strains"."\t"."4");
push(@happyM,"strain_names"."\t"."Bxd2"."\t"."CAST"."\t"."MRL"."\t"."NZM");
for(my $i=1; $i<scalar(@f_data);$i++){
	chomp $f_data[$i];
	my ($marker,$chr,$pos,$cM,$A1,$A2,$bxd2,$cast,$mrl,$nzm) = split/\t/,$f_data[$i];
	push(@happyM,"marker\t$marker\t3\t$chr\t$cM");
	my @parents = ($bxd2,$cast,$mrl,$nzm);
	my $genoTmp = $bxd2.$cast.$mrl.$nzm;
	my $A1_t = () = $genoTmp =~ /$A1/g;
	my $A2_t = () = $genoTmp =~ /$A2/g;
	my @A1_ar;
	my @A2_ar;
# calculate probability of SNPs comming from each founder
	foreach my $fGeno (@parents){
		my $A1_f = () = $fGeno =~ /$A1/g;
		my $A2_f = () = $fGeno =~ /$A2/g;
		my $A1_ff = $A1_f/$A1_t;
		my $A2_ff = $A2_f/$A2_t;
		if(($A1_ff > 0) && ($A1_ff < 1)){$A1_ff = sprintf("%.2f",$A1_ff);}
		if(($A2_ff > 0) && ($A2_ff < 1)){$A2_ff = sprintf("%.2f",$A2_ff);}
		push(@A1_ar,$A1_ff);
		push(@A2_ar,$A2_ff);
	}
	my $A1_out = "allele"."\t".$A1."\t".join("\t",@A1_ar);
	my $A2_out = "allele"."\t".$A2."\t".join("\t",@A2_ar);
	my $na_out = "allele"."\t"."NA"."\t"."0.25"."\t"."0.25"."\t"."0.25"."\t"."0.25";
	push(@happyM,$na_out); push(@happyM,$A1_out); push(@happyM,$A2_out);
}
#make ped file for for all the samples for each chromsome 
my $markerout = join("\n",@happyM);
my $allelefileO = "chr".$chr.".AAIL.alleles";
open ALLELE,">",$allelefileO or die "Problem with permission file cannot be created\n";
print ALLELE $markerout."\n";
close ALLELE;

my $datafile ="sample_genotype_".$chr.".txt" or die "file for $chr donot exists!!\n";
open SAMPLE, $datafile or die "$datafile donot exists!!\n";
my $k=0;
my @ped;
while(<SAMPLE>){
	chomp $_;
	my @tmp;
	my @s = split("\t",$_);
	for(my $i=1;$i<scalar(@s);$i++){
		if($s[$i] eq "NA"){
			push(@tmp,"NA"." "."NA");
		}else{
			my($a1,$a2) = split("",$s[$i]);
			push(@tmp,$a1." ".$a2);
		}
	}
	my $pedO = $s[0]." "."1"." ".join(" ",@tmp);
	push(@ped,$pedO);
}
close SAMPLE;
my $pedout = join("\n",@ped);
my $datafileO = "chr".$chr.".AAIL.data";
open DATAOUT,">",$datafileO or die "Problem with permission file cannot be created\n";
print DATAOUT $pedout."\n";
close DATAOUT;
