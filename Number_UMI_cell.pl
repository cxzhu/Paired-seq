#!/usr/bin/perl
use strict;
use warnings;

my %hash;
open IN, "samtools view $ARGV[0]|" or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\t+/, $_;
	my $cell_id = substr($tmp[0], -22, 11);
	my $umi = substr($tmp[0], -10, 10);
	my $chr = $tmp[2];
	my $pos = $tmp[3];
	next if $cell_id =~ m/\*/;
	$hash{$cell_id}{$chr}{$pos}{$umi} = 1;# if not exists $hash{$cell_id}{$chr}{$pos};
}
close IN;

open OUT, ">$ARGV[0]\_count.xls";
foreach my $cell_id (sort keys %hash){
	my $total = 0;
	foreach my $chr (sort keys %{$hash{$cell_id}}){
		foreach my $pos (sort keys %{$hash{$cell_id}{$chr}}){
			my $n = keys %{$hash{$cell_id}{$chr}{$pos}};
			$total+=$n;
		}
	}
	print OUT "$cell_id\t$total\n";
}
close OUT;

