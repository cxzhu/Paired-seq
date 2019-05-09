#!/usr/bin/perl
use strict;
use warnings;

open IN, "/projects/ps-renlab/chz272/02.REACH-seq/11.final_processing/03.link/04.stage_specific_link/01.output/merge_link.xls" or die $!;
my %enh;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $gene = $tmp[1];
	my $enha = $tmp[2];
	$enh{$enha}{$gene} = 1;
}
close IN;

open IN, "/projects/ps-renlab/chz272/02.REACH-seq/10.merge_Combined_Analysis/02.all_combine/Combined_matrix/bw_rna_seq_seperate/TTS_enrichment.xls" or die $!;
my %expression;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $gene = $tmp[0];
	my $E12 = ($tmp[1] + $tmp[3] + $tmp[5] + $tmp[9] + $tmp[11])/5;
	my $E16 = ($tmp[2] + $tmp[4] + $tmp[6] + $tmp[10] + $tmp[12])/5;
	$expression{$gene}{"E12"} = $E12;
	$expression{$gene}{"E16"} = $E16;
}
close IN;

my %accessibility;
open IN, "/projects/ps-renlab/chz272/02.REACH-seq/10.merge_Combined_Analysis/02.all_combine/Combined_matrix/bw_seperate/TSS_enrichment.xls" or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $gene = $tmp[0];
	my $E12 = ($tmp[1] + $tmp[3] + $tmp[5] + $tmp[9] + $tmp[11])/5;
	my $E16 = ($tmp[2] + $tmp[4] + $tmp[6] + $tmp[10] + $tmp[12])/5;
	$accessibility{$gene}{"E12"} = $E12;
	$accessibility{$gene}{"E16"} = $E16;

}
close IN;

open OUT, ">Gene_Enhancers_Dynamics_2.xls" or die $!;
open IN, "/projects/ps-renlab/chz272/02.REACH-seq/10.merge_Combined_Analysis/02.all_combine/Combined_matrix/bw_seperate/REG_enrichment.xls" or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $enha = $tmp[0];
	next if not exists $enh{$enha};
	my $E12 = ($tmp[1] + $tmp[3] + $tmp[5] + $tmp[9] + $tmp[11])/5;
	my $E16 = ($tmp[2] + $tmp[4] + $tmp[6] + $tmp[10] + $tmp[12])/5;
	foreach my $gene (sort keys %{$enh{$enha}}){
		next if not exists $expression{$gene};
		next if not exists $accessibility{$gene};
		my $expr12 = $expression{$gene}{"E12"};
		my $expr16 = $expression{$gene}{"E16"};
		my $acce12 = $accessibility{$gene}{"E12"};
		my $acce16 = $accessibility{$gene}{"E16"};
		print OUT "$gene\t$enha\t$expr12\t$expr16\t$acce12\t$acce16\t$E12\t$E16\n";
	}
}
close IN;
close OUT;
