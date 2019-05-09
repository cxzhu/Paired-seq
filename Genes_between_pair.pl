#!/usr/bin/perl
use strict;
use warnings;

my %data;
my %hash;
open IN, "/projects/ps-renlab/chz272/genome_ref/mm10_tss_all.xls" or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $chr = $tmp[1];
	next if $_ =~ m/chrM/;
	my $pss = int($tmp[2]/100);
	my $pse = int($tmp[3]/100);
	foreach my $pos ($pss .. $pse){
		$hash{$chr}{$pos} = $tmp[0];
	}
}
close IN;

open IN, "$ARGV[0]" or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my @tss = split/[\:\-]/, $tmp[0];
	my @reg = split/[\:\-]/, $tmp[2];
	my $chr = $tss[0];
	my $pss = $tss[2];
	my $pse = $reg[1];
	if($tss[2]>$reg[2]){
		$pss = $reg[2];
		$pse = $tss[1];
	}
	my $link = "$tmp[1]_$tmp[2]";
	$data{$link}{"empty"} = 1;
	foreach my $pos (int($pss / 100) .. int($pse / 100)){
		next if not exists $hash{$chr}{$pos};
		my $ele = $hash{$chr}{$pos};
		next if $ele eq "$tmp[1]";
		$data{$link}{$ele} = 1;
	}
}
close IN;

open OUT, ">$ARGV[0]_NumGenesBetweenLink.xls" or die $!;
foreach my $link (sort keys %data){
	my $num = keys %{$data{$link}};
	$num=$num-1;
	print OUT "$link\t$num\n";
}
close OUT;
