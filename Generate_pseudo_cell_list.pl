#!/usr/bin/perl
use strict;
use warnings;

my %value;
my %clust;
my %n_each_clust;

if(0){
	open IN, "$ARGV[0]" or die $!;
	my $fist = <IN>;

	chomp($fist);
	my @cell_list = split/\s+/, $fist;
	my $m = -1;
	while(<IN>){
		chomp;
		$m++;
		my $cell_a = $cell_list[$m];
		my @sp = split/\s+/, $_;
		my $n = -1;
		foreach my $v (@sp){
			$n++;
			my $cell_b = $cell_list[$n];
			#$value{$v}{$cell_a}{$cell_b} = 1;
			print "$v\t$cell_a\t$cell_b\n";
		}
	}
	close IN;
}


if(1){
	open IN, "$ARGV[0]" or die $!;
	open OUT, ">$ARGV[0]_pseudo_cell_E12.xls" or die $!;
	#open OUT, ">$ARGV[0]_pseudo_cell_E16.xls" or die $!;
	my $full = 0;
	my $cid = 0;
	my %report;
	my $line = 0;
	while(<IN>){
		$line ++;
		if(int($line / 1000000) == $line / 1000000){
			print STDERR "E12".$line."\n";
		}
		chomp;
		my @tmp = split/\s+/, $_;
		my $v = $tmp[0];
		next if $v == 0;
		my $cell_a = $tmp[1];
		my $cell_b = $tmp[2];
		next if (exists $clust{$cell_a} && exists $clust{$cell_b});
		my $t1 = substr($cell_a, -1, 1);
		my $t2 = substr($cell_b, -1, 1);
		next if $t1>=3;
		next if $t2>=3;
		## if one assigned
		if(exists $clust{$cell_a} || exists $clust{$cell_b}){
			my $c = 0;
			$c = $clust{$cell_a} if exists $clust{$cell_a};
			$c = $clust{$cell_b} if exists $clust{$cell_b};
			### check if full
			if($n_each_clust{$c}>50){
				if(not exists $report{$c}){
					$full++;
					print STDERR "$full\n";
					$report{$c} = 1;
				}
				next;
			}
			$clust{$cell_a} = $c;
			$clust{$cell_b} = $c;
			$n_each_clust{$c}++;
		}
		else{
		### nobody is assigned
			$cid++;
			next if $cid > 100;
			$n_each_clust{$cid} = 1;
			print STDERR "E12 clust $cid\n";	
			$clust{$cell_a} = $cid;
			$clust{$cell_b} = $cid;
		}
	}



	foreach my $cell (sort keys %clust){
		my $cid = $clust{$cell};
		print OUT "$cell\t$cid\n";
	}
	close OUT;
}
