#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 20

# input
base='/scratch3/rogerlab_databases/Fastq_files/Ergobibamus/Illumina/RNAseq/trimmomatic-0.39'
R1=${base}/eb_rna_fw.prd.fastq.gz
R2=${base}/eb_rna_rv.prd.fastq.gz
U1=${base}/eb_rna_fw.unp.fastq.gz
U2=${base}/eb_rna_rv.unp.fastq.gz
ASSEMBLY='ergo_cyp_genome.fasta.masked'

# output
PREFIX=rnaseq_vs_masked_ergo_genome
INDEX=masked_ergo_genome

# load hisat2 software
source activate hisat2

# align rnaseq reads
hisat2-build -f $ASSEMBLY $INDEX
## q: quiet, k: how many a single read can map
hisat2 \
    -q --threads 20 --phred33 \
    --max-intronlen 1000 \
    --rna-strandness RF -k 2 \
    -x $INDEX -1 $R1 -2 $R2 -U $U1,$U2 | \
	samtools sort --threads 20 -o $PREFIX.sort.bam -

samtools index -@ 20 $PREFIX.sort.bam $PREFIX.sort.bam.bai

conda deactivate

