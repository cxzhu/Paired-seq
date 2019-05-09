#!/usr/bin/perl
use strict;
use warnings;

## cut off
my $peak_cutoff = 10;
my $gene_cutoff = 10;

## DNA matrix
my %dna_data;
my %dna_barcodes_list;
my %dna_peaks_list;

## RNA matrix
my %rna_data;
my %rna_barcodes_list;
my %rna_genes_list;

## QC matrix
my %dna_nPeaks;
my %rna_nGenes;
my %val_cells;

my %data;


### DNA

open IN, "./DNA_matrix/barcodes.tsv" or die $!;
my $i=0;
while(<IN>){
	chomp;
	$i++;
	$dna_barcodes_list{$i} = $_;
}
close IN;

open IN, "./DNA_matrix/genes.tsv" or die $!;
$i = 0;
while(<IN>){
	chomp;
	$i++;
	$dna_peaks_list{$i} = $_;
}
close IN;

open IN, "./DNA_matrix/matrix.mtx" or die $!;
<IN>;
<IN>;
<IN>;

while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $barcode_id = $tmp[1];
	my $peak_id = $tmp[0];
	my $value = $tmp[2];
	my $barcode = $dna_barcodes_list{$barcode_id};
	my $peak = $dna_peaks_list{$peak_id};
	$dna_data{$peak}{$barcode} = 1;
	$dna_nPeaks{$barcode} = 0 if not exists $dna_nPeaks{$barcode};
	$dna_nPeaks{$barcode}++;
	$data{$barcode}{$peak} = $value;
}
close IN;

### RNA

open IN, "./RNA_matrix/barcodes.tsv" or die $!;
$i = 0;
while(<IN>){
	chomp;
	$i++;
	$rna_barcodes_list{$i} = $_;
}
close IN;

open IN, "./RNA_matrix/genes.tsv" or die $!;
$i = 0;
while(<IN>){
	chomp;
	$i++;
	$rna_genes_list{$i} = $_;
}
close IN;

open IN, "./RNA_matrix/matrix.mtx" or die $!;
<IN>;
<IN>;
<IN>;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $barcode_id = $tmp[1];
	my $gene_id = $tmp[0];
	my $value = $tmp[2];
	my $barcode = $rna_barcodes_list{$barcode_id};
	my $gene = $rna_genes_list{$gene_id};
	$rna_data{$gene}{$barcode} = 1;
	$rna_nGenes{$barcode} = 0 if not exists $rna_nGenes{$barcode};
	$rna_nGenes{$barcode}++;
	$data{$barcode}{$gene} = $value;
}
close IN;


foreach my $cell (sort keys %dna_nPeaks){
	next if not exists $rna_nGenes{$cell};
	next if $rna_nGenes{$cell} < $gene_cutoff;
	next if $dna_nPeaks{$cell} < $peak_cutoff;
	$val_cells{$cell} = 1;
}

### output matrix;
system("mkdir Combined_matrix");
system("mkdir Combined_matrix/DNA");
system("mkdir Combined_matrix/RNA");

my %feature_list_dna;
my $total_dna = 0;
## DNA first
foreach my $peak (sort keys %dna_data){
	foreach my $cell (sort keys %{$dna_data{$peak}}){
		next if not exists $val_cells{$cell};
		$total_dna++;
		$feature_list_dna{$peak} = 1;
	}
}
## RNA second
my %feature_list_rna;
my $total_rna = 0;
foreach my $gene (sort keys %rna_data){
	foreach my $cell (sort keys %{$rna_data{$gene}}){
		next if not exists $val_cells{$cell};
		$total_rna++;
		$feature_list_rna{$gene} = 1;
	}
}

$i=0;
open OUT, ">./Combined_matrix/DNA/genes.tsv" or die $!;
foreach my $peak (sort keys %feature_list_dna){
	$i++;
	$feature_list_dna{$peak} = $i;
	print OUT $peak;
	print OUT "\n";
}
close OUT;
my $num_features_DNA = $i;


$i=0;
open OUT, ">./Combined_matrix/RNA/genes.tsv" or die $!;
foreach my $gene (sort keys %feature_list_rna){
	$i++;
	$feature_list_rna{$gene} = $i;
	print OUT $gene;
	print OUT "\n";
}
my $num_features_RNA = $i;
close OUT;

my $num_cells = keys %val_cells;
open OU1, ">./Combined_matrix/DNA/barcodes.tsv" or die $!;
open OU2, ">./Combined_matrix/RNA/barcodes.tsv" or die $!;

$i = 0;
foreach my $barcode (sort keys %val_cells){
	$i++;
	print OU1 $barcode;
	print OU1 "\n";
	print OU2 $barcode;
	print OU2 "\n";
	$val_cells{$barcode} = $i;
}
close OU1;
close OU2;


open OUT, ">./Combined_matrix/DNA/matrix.mtx" or die $!;
print OUT "\%\%MatrixMarket matrix coordinate real general\n\%\n$num_features_DNA $num_cells $total_dna\n";
foreach my $cell (sort keys %val_cells){
	my $cell_id = $val_cells{$cell};
	foreach my $feature (sort keys %{$data{$cell}}){
		next if ( (not exists $feature_list_dna{$feature})  );
		my $feature_id = 0;
		$feature_id = $feature_list_dna{$feature} if exists $feature_list_dna{$feature};
		my $value = $data{$cell}{$feature};
		print OUT "$feature_id $cell_id $value\n";
	}
}
close OUT;


open OUT, ">./Combined_matrix/RNA/matrix.mtx" or die $!;
print OUT "\%\%MatrixMarket matrix coordinate real general\n\%\n$num_features_RNA $num_cells $total_rna\n";
foreach my $cell (sort keys %val_cells){
	my $cell_id = $val_cells{$cell};
	foreach my $feature (sort keys %{$data{$cell}}){
		next if ( (not exists $feature_list_rna{$feature}) );
		my $feature_id = 0;
		$feature_id = $feature_list_rna{$feature} if exists $feature_list_rna{$feature};
		my $value = $data{$cell}{$feature};
		print OUT "$feature_id $cell_id $value\n";
	}
}
close OUT;




















