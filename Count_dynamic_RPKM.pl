#!/usr/bin/perl
use strict;
use warnings;


my $list = $ARGV[0]; # "/projects/ps-renlab/chz272/02.REACH-seq/11.final_processing/02.ENCODE_COMBINE/03.all_combine/test1/0502_nor_final/enhancer/2019_05_01_All_Merge_Cluster_final.xls";
my $ann = $ARGV[1];
my $promoter = "/projects/ps-renlab/chz272/genome_ref/mm10_tss_all.xls";
#my $bam = "/projects/ps-renlab/chz272/02.REACH-seq/08.merge_DNA_all/CZ195_mm10_rmdup.bam_merge.bam.clean.bam";
my $bam = "/projects/ps-renlab/chz272/02.REACH-seq/08.merge_DNA_all/CZ195_mm10_rmdup.bam_merge.bam";

if($ann eq "TSS"){
	$promoter = "/projects/ps-renlab/chz272/genome_ref/mm10_tss_all.xls";
}
elsif($ann eq "TTS"){
	$promoter = "/projects/ps-renlab/chz272/genome_ref/mm10_utr3_all.xls";
	$bam = "/projects/ps-renlab/chz272/02.REACH-seq/09.merge_RNA_all/CZ196_mm10Aligned.out.sam_sorted.bam_clean_rmdup.bam_merge.bam.clean.bam"
}
elsif($ann eq "REG"){
	$promoter= "/projects/ps-renlab/chz272/17.0425_reach/01.peaks/all_DNA_reg.bed";
}

my %tss;
open IN, $promoter or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $gene_name = $tmp[0];
	my $chr = $tmp[1];
	my $pss = int($tmp[2]/100);
	my $pse = int($tmp[3]/100);
	foreach my $pos ($pss .. $pse){
		$tss{$chr}{$pos} = $gene_name;
	}
}
close IN;


open IN, $list or die $!;
<IN>;
my %clist;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	#my $ident = $tmp[3];
	my $ident = int($tmp[0]/10) + 1;
	$ident = 10 if $ident > 10;
	#$clist{$tmp[0]} = 0 if $ARGV[0] eq "0";
	#next if $ident ne $ARGV[0];
	$clist{$tmp[1]} = $ident;
}
close IN;

my %total;
my %data;
open IN, "samtools view $bam|" or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $chr = $tmp[2];
	my $pos = int($tmp[3]/100);
	next if not exists $tss{$chr}{$pos};
	my $cid = substr($tmp[0], -22, 11);
	next if not exists $clist{$cid};
	my $cluster = $clist{$cid};
#	my $id = substr($cid, -2, 2);
#	if($id < 40){
#		$cluster = 3;
#	}
#	else{
#		my $idd = substr($id, -1, 1);
#		$cluster = 1 if $idd < 3;
#		$cluster = 2 if $idd > 2;
#	}
	$total{$cluster} = 0 if not exists $total{$cluster};
	$total{$cluster}++; #= 0 if not exists $total{$cluster};
	my $gene_name = $tss{$chr}{$pos};
	$data{$gene_name}{$cluster} = 0 if not exists $data{$gene_name}{$cluster};
	$data{$gene_name}{$cluster} ++;
}
close IN;

open OUT, ">C$ARGV[0]\_$ARGV[1]\_count.xls" or die $!;
foreach my $gene_name(sort keys %data){
	my $output = "$gene_name";
	foreach my $cluster (1 .. 10){
		my $value = 0;
		if(exists $data{$gene_name}{$cluster}){
			$value = $data{$gene_name}{$cluster} / $total{$cluster} * 1000000;
		}
		$output .= "\t$value";
	}	
	$output .= "\n";
	print OUT $output;
}
close OUT;
