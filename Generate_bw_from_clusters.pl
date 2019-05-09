#!/usr/bin/perl
use strict;
use warnings;
use FileHandle;

#my $bam = "/projects/ps-renlab/chz272/02.REACH-seq/05.salk_tissue/05.batch/01.bam/DNA_merge.bam";
my $bam = "/projects/ps-renlab/chz272/02.REACH-seq/08.merge_DNA_all/CZ195_mm10_rmdup.bam_merge.bam.clean.bam";

my %fhn;
my %fh;
my %cid;


my $in = "louvain_tSNE_cluster.xls";
$in = $ARGV[0];

my %clusters;
my $total_cells = 0;
my %cluster_cells;


open IN, $in or die $!;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $id = $tmp[0];
	my $cluster = $tmp[1];
	$fhn{$cluster} = 1;
	$clusters{$id} = $cluster;
	$total_cells ++;
	$cluster_cells{$cluster} = 0 if not exists $cluster_cells{$cluster};
	$cluster_cells{$cluster}++;
}
close IN;

foreach my $i (sort {$a<=>$b} keys %fhn){
	open $fh{$i}, "|samtools view -b - > ATAC_cluster_$i.bam" or die $!;
}


open IN, "samtools view -h $bam|" or die $!;
while(<IN>){
	my $line = $_;
	if($line =~ m/^\@/){
		foreach my $i (sort {$a<=>$b} keys %fhn){
			$fh{$i}->print($line);
		}
	}
	else{
		chomp($line);
		my @tmp = split/\t+/, $line;
		my $cid = substr($tmp[0], -22, 11);
		next if not exists $clusters{$cid};
		my $i = $clusters{$cid};
		$fh{$i}->print("$line\n");
	}
}

close IN;

foreach my $i (sort {$a<=>$b} keys %fhn){
	close $fh{$i};
}


$bam = "/projects/ps-renlab/chz272/02.REACH-seq/09.merge_RNA_all/CZ196_mm10Aligned.out.sam_sorted.bam_clean_rmdup.bam_merge.bam.clean.bam";

foreach my $i (sort {$a<=>$b} keys %fhn){
	open $fh{$i}, "|samtools view -b - > RNA_cluster_$i.bam" or die $!;
}


open IN, "samtools view -h $bam|" or die $!;
while(<IN>){
	my $line = $_;
	if($line =~ m/^\@/){
		foreach my $i (sort {$a<=>$b} keys %fhn){
			$fh{$i}->print($line);
		}
	}
	else{
		chomp($line);
		my @tmp = split/\t+/, $line;
		my $cid = substr($tmp[0], -22, 11);
		next if not exists $clusters{$cid};
		my $i = $clusters{$cid};
		$fh{$i}->print("$line\n");
	}
}

close IN;
open OUT, ">proc.sh" or die $!;
foreach my $i (sort {$a<=>$b} keys %fhn){
	close $fh{$i};
	my $cur_clust_number = $cluster_cells{$i};
	my $sf = $total_cells/$cur_clust_number;
	print OUT "samtools index ATAC_cluster_$i.bam\n";
	print OUT "samtools index RNA_cluster_$i.bam\n";
#	print OUT "cxrun 1 4 bamCoverage -p 4 -e 100 -b ATAC_cluster_$i.bam --normalizeUsing None --scaleFactor $sf -o PLSA_DNA_sca_C$i.bw\n";
	print OUT "cxrun 1 4 bamCoverage -p 4 -e 100 -b ATAC_cluster_$i.bam --normalizeUsing RPKM -o DNA_C$i.bw\n";
	print OUT "sleep 1\n";
#	print OUT "cxrun 1 4 bamCoverage -p 4 -e 100 -b RNA_cluster_$i.bam --normalizeUsing None --scaleFactor $sf -o PLSA_RNA_sca_C$i.bw\n";
	print OUT "cxrun 1 4 bamCoverage -p 4 -e 100 -b RNA_cluster_$i.bam --normalizeUsing RPKM -o RNA_C$i.bw\n";
	print OUT "sleep 1\n";
}
close OUT;

system("cxrun 1 4 sh proc.sh")

#foreach my $i (sort {$a<=>$b} keys %fhn){
#	system("cxrun 1 1 sh proc.sh ATAC_cluster_$i");
#	system("sleep 1");
#}