s=$1

reacthools combine CZ${s}_R1.fq.gz CZ${s}_R2.fq.gz
zcat CZ${s}_combined.fq.gz | bowtie /projects/ps-renlab/chz272/genome_ref/cell_id/cell_id -p 8 -v 1 -m 1 --norc - -S CZ${s}.sam

reachtools convert CZ${s}.sam
trim_galore -a AAAAAAAAAAAAAAAACCTGCAGGNNNNNNNNNN CZ${s}_cov.fq.gz
trim_galore -a GGGGGGNNNNNNNNNNNNNNNN CZ${s}_cov_trimmed.fq.gz

STAR --runThreadN 6 --genomeDir /projects/ps-renlab/chz272/genome_ref/refdata-cellranger-mm10-3.0.0/star/ --readFilesIn CZ${s}_cov_trimmed_trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix CZ${s}_mm10 
samtools sort CZ${s}_mm10Aligned.sam CZ${s}_mm10
reachtools rmdup CZ${s}_mm10.bam

reachtools bam2Mtx CZ${s}_mm10.bam mm10
