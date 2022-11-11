#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 20

# input
base='/scratch3/rogerlab_databases/Fastq_files/Ergobibamus/Illumina/gDNA/qc_seq/trimmomatic-0.39'
R1=${base}/eb_dna_fw_prd.fastq.gz
R2=${base}/eb_dna_rv_prd.fastq.gz
U1=${base}/eb_dna_fw_unp.fastq.gz
U2=${base}/eb_dna_rv_unp.fastq.gz
ASSEMBLY='canu_ergo_contigs_clean.fasta'

# output
PREFIX=dnaseq_vs_canu_ergo_assembly
INDEX=canu_ergo_assembly

# load hisat2 software
source activate hisat2

# align rnaseq reads
hisat2-build -f $ASSEMBLY $INDEX
## q: quiet, k: how many a single read can map
## --no-spliced-alignment to specify dna read mapping
hisat2 \
    -q --threads 20 --phred33 -k 2 --no-spliced-alignment \
    -x $INDEX -1 $R1 -2 $R2 -U $U1,$U2 -S $PREFIX.sam
## -X/--maxins <int> is not set, (expected library fragment size)
## so default value of 500 bp is used

conda deactivate

# convert to bam
samtools view --threads 20 -b $PREFIX.sam > $PREFIX.bam
samtools sort --threads 20 $PREFIX.bam > $PREFIX.sorted.bam
samtools index -@ 20 $PREFIX.sorted.bam $PREFIX.sorted.bam.bai

rm $PREFIX.sam $PREFIX.bam

# notes on HISAT2 #

## it seems that HISAT2 reports only three possible MAPQ values: 0, 1 and 60
## it is not explained on the website or manual what it means, but a biostars website suggests the following
## 0  - a non-uniquely mapped read, i.e. this read is also mapped somewhere else. there are many indels/mismatches
## 1  - a non-uniquely mapped read, i.e. this read is also mapped somewhere else. there are few or no indels/mismatches
## 60 - a uniquely mapped read, i.e. this read is not mapped anywhere else. indels/mismatches are not taken into account in calculating MAPQ

## for HISAT2, any read mapping that is "secondary" (i.e. there is a better mapping for this read somewhere else),
## the 256 bitwise FLAG is set (i.e. "not primary alignment")
