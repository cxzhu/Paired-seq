#!/usr/bin/perl
use strict;
use warnings;

open IN, "$ARGV[0]" or die $!;
my $name = substr($ARGV[0],0,length($ARGV[0])-4);
open OUT, "|gzip - >  $name\_cov.fq.gz";

while(<IN>){
	next if m/^\@/;
	my @tmp = split/\s+/, $_;
	next if $tmp[2] eq '*';
	my $cell_id = $tmp[2];
	my @sp = split/\:/, $tmp[0];
	my $readname = "@".$sp[0].":".$sp[1].":".$sp[2].":".$sp[3].":".$sp[4].":".$sp[5].":".$sp[6].":".$cell_id.":".$sp[7];
	my $read = $sp[8];
	my $l = length($read);
	my $mark = "+";
	my $qual = substr($tmp[0], -$l, $l);
	print OUT "$readname\n$read\n$mark\n$qual\n";
}
close IN;
close OUT;
