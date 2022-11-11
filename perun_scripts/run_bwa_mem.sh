#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 20
#$ -q 256G-batch

# input
base='/scratch3/rogerlab_databases/Fastq_files/Ergobibamus/Illumina/gDNA/qc_seq/trimmomatic-0.39'
R1=${base}/eb_dna_fw_prd.fastq.gz
R2=${base}/eb_dna_rv_prd.fastq.gz
U1=${base}/eb_dna_fw_unp.fastq.gz
U2=${base}/eb_dna_rv_unp.fastq.gz
ASSEMBLY='canu_medaka_polish.fasta'

# output
PREFIX=dnaseq_vs_canu_medaka_polish

# load bwa
source activate bwa

# align dnaseq reads
## first index
bwa index $ASSEMBLY
## then align paired reads
bwa mem -t 20 $ASSEMBLY $R1 $R2 > ${PREFIX}-pe.sam
## and unpaired reads
bwa mem -t 20 $ASSEMBLY $U1 > ${PREFIX}-se1.sam
bwa mem -t 20 $ASSEMBLY $U2 > ${PREFIX}-se2.sam

conda deactivate

# convert to bam
for SAM in *.sam; do
    samtools view --threads 20 -b $SAM > ${SAM/sam/bam}
    samtools sort --threads 20 ${SAM/sam/bam} > ${SAM/sam/sorted.bam}
    samtools index -@ 20 ${SAM/sam/sorted.bam} ${SAM/sam/sorted.bam.bai}
    rm $SAM ${SAM/sam/bam}
done

