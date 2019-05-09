s=$1

reacthools combine CZ${s}_R1.fq.gz CZ${s}_R2.fq.gz
zcat CZ${s}_combined.fq.gz | bowtie /projects/ps-renlab/chz272/genome_ref/cell_id/cell_id -p 8 -v 3 --norc - -S CZ${s}.sam

reachtools convert CZ${s}.sam
trim_galore CZ${s}_cov.fq.gz

bowtie2 -x /projects/ps-renlab/chz272/genome_ref/mm10/mm10 -U CZ${s}_cov_trimmed.fq.gz --no-unal -p 8 -S CZ${s}_mm10.sam
samtools sort CZ${s}_mm10.sam CZ${s}_mm10
reachtools rmdup CZ${s}_mm10.bam

reachtools bam2Mtx CZ${s}_mm10.bam mm10_1kb
