#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 40
#$ -q 256G-batch

# input
base='/scratch3/rogerlab_databases/Fastq_files/Ergobibamus/Illumina/gDNA/qc_seq/trimmomatic-0.39'
R1=${base}/eb_dna_fw_prd.fastq.gz
R2=${base}/eb_dna_rv_prd.fastq.gz
U1=${base}/eb_dna_fw_unp.fastq.gz
U2=${base}/eb_dna_rv_unp.fastq.gz
ASSEMBLY='canu_ergo_contigs.fasta'

# output
PREFIX=dnaseq_vs_canu_medaka_pilon_polish

# load bwa-mem2 software
source activate bwa-mem2

# align dnaseq reads
## first index
bwa-mem2 index -p ${ASSEMBLY%.fasta} $ASSEMBLY
## then align paired reads
bwa-mem2 mem -t 40 ${ASSEMBLY%.fasta} $R1 $R2 | samtools sort --threads 40 -o ${PREFIX}-pe.sort.bam -
## and unpaired reads
bwa-mem2 mem -t 40 ${ASSEMBLY%.fasta} $U1 | samtools sort --threads 40 -o ${PREFIX}-se1.sort.bam -
bwa-mem2 mem -t 40 ${ASSEMBLY%.fasta} $U2 | samtools sort --threads 40 -o ${PREFIX}-se2.sort.bam -

conda deactivate

# convert to bam
for BAM in *.bam; do
    samtools index -@ 40 $BAM $BAM.bai
done

