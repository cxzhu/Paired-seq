#!/usr/bin/perl
use strict;
use warnings;

my %black;
open IN, "zcat /projects/ps-renlab/mandyjiang/Blacklist/lists/mm10-blacklist.v2.bed.gz|" or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $chr = $tmp[0];
	my $pss = int($tmp[1]/100);
	my $pse = int($tmp[2]/100);
	foreach my $pos ($pss .. $pse){
		$black{$chr}{$pos} = 1;
	}
}
close IN;

open IN, "./DNA/peaks.bed" or die $!;
my $i = 0;
my $j = 0;
my %old_peak_list;
open OUT, ">./DNA_filtered/peaks.bed" or die $!;
while(<IN>){
	chomp;
	$i++;
	my @tmp = split/\t+/, $_;
	my $chr = $tmp[0];
	my $pss = int($tmp[1]/100);
	my $pse = int($tmp[2]/100);
	my $valid = 1;
	foreach my $pos ($pss .. $pse){
		$valid = 0 if exists $black{$chr}{$pos};
	}
	next if $valid == 0;
	$j++;
	$old_peak_list{$i} = $j;
	print OUT $_."\n";
}
close IN;
close OUT;
system("cp ./DNA/barcodes.tsv ./DNA_filtered/");

open IN, "./DNA/matrix.mtx" or die $!;
my %hash;
my $total = 0;
<IN>;
my $t = <IN>;
chomp($t);
my @sp = split/\s+/, $t;
my $n_cells = $sp[1];
my $n_peaks = $j;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $old_peak_id = $tmp[0];
	my $cell_id = $tmp[1];
	next if not exists $old_peak_list{$old_peak_id};
	$total++;
	my $peak_id = $old_peak_list{$old_peak_id};
	$hash{$peak_id}{$cell_id} = 1;
}
close IN;

open OUT, ">./DNA_filtered/matrix.mtx" or die $!;
print OUT "\%\%MatrixMarket matrix coordinate real general\n\%\n$n_peaks $n_cells $total\n";
foreach my $peak_id (sort {$a<=>$b} keys %hash){
	foreach my $cell_id (sort {$a<=>$b} keys %{$hash{$peak_id}}){
		print OUT "$peak_id $cell_id 1\n";
	}
}
close OUT;
