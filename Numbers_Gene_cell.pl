#!/usr/bin/perl
use strict;
use warnings;
my %anno;
open IN, "/projects/ps-renlab/chz272/annotations/merge.txt" or die $!;
#open IN, "/projects/ps-renlab/chz272/annotations/human_hg19_genomic_region_annotation/hg19.gencode.v19.txt" or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $chr = $tmp[0];
	my $pss = int($tmp[1]/1000);
	my $pse = int($tmp[2]/1000);
	foreach my $pos ($pss .. $pse+1){
		$anno{$chr}{$pos} = $tmp[3];
	}
}
close IN;

open IN, "samtools view $ARGV[0]|" or die $!;
my %hash;
while(<IN>){
	chomp;
	my @tmp = split/\t+/, $_;
	my $cell_id = substr($tmp[0], -22, 11);
	my $umi = substr($tmp[0], -10, 10);
	my $chr = $tmp[2];
	my $pos = int($tmp[3]/1000);
	next if not exists $anno{$chr}{$pos};
	my $gene = $anno{$chr}{$pos};
	$hash{$cell_id}{$gene}{$umi} = 1;
}
close IN;

open OUT, ">$ARGV[0].countGene.xls" or die $!;
foreach my $cell_id (sort keys %hash){
	my $n = 0;
	foreach my $gene (sort keys %{$hash{$cell_id}}){
		$n++;
	}
	print OUT "$cell_id\t$n\n";
}
close OUT;
