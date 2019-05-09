#!/usr/bin/perl
use strict;
use warnings;

open IN1, "dna_jar.all.xls" or die $!;
open IN2, "rna_jar_full.xls" or die $!;
open OU1, ">had_full.xls" or die $!;
open OU2, ">had_fourth.xls" or die $!;
while(<IN1>){
	chomp;
	my @dna = split/\s+/, $_;
	my $l2 = <IN2>;
	chomp($l2);
	my @rna = split/\s+/, $_;
	my $output1;
	my $output2;
	my $c = 0;
	foreach my $i ( 0 .. $#dna){
		$c++;
		my $v = ($dna[$i]+6)*($rna[$i]+6);
		$output1 .= "$v\t";
		next if $c < 3;
		$output2 .= "$v\t";
		$c = 0;
	}
	chomp($output1);
	chomp($output2);
	print OU1 "$output1\n";
	print OU2 "$output2\n";
}
close IN1;
close IN2;
close OU1;
close OU2;
