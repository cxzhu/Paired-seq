#!/usr/bin/perl
use strict;
use warnings;

open IN, $ARGV[0] or die $!;
open OUT, ">$ARGV[0]_best.xls" or die $!;
my %hash;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $cur = $tmp[1];
	my $cor = $tmp[2];
	next if $cor eq "NA";
	next if $cor < 0;
	$hash{$cur}{"cor"}=$cor if not exists $hash{$cur}{"cor"};
	$hash{$cur}{"line"} = $_ if not exists $hash{$cur}{"line"};

	next if $cor <= $hash{$cur}{"cor"};
	$hash{$cur}{"cor"} = $cor;
	$hash{$cur}{"line"} = $_;

}
close IN;
foreach my $cur (sort keys %hash){
	print OUT $hash{$cur}{"line"}."\n";
}

close OUT;
