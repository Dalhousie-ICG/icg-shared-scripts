#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -m bea
#$ -M joran.martijn@dal.ca
#$ -pe threaded 20
#$ -q 768G-batch

base='/scratch3/rogerlab_databases/Fastq_files/Ergobibamus/Illumina/RNAseq/trimmomatic-0.39'
R1=${base}/eb_rna_fw.prd.fastq.gz
R2=${base}/eb_rna_rv.prd.fastq.gz
# U1=${base}/eb_rna_fw.unp.fastq.gz
# U2=${base}/eb_rna_rv.unp.fastq.gz
RNASEQ_VS_REF_BAM=rnaseq_vs_masked_ergo_cyp_genome.sort.bam

# This is the latest version of Trinity we have, v2.13.2(released Sep 4th 2021)

# library was prepped with dUTP method,
# so  all forward reads (read1) map to antisense
# and all reverse reads (read2) map to sense
# strands of the reference genome

# Trinity considers this type RF

source activate trinity

Trinity \
        --seqType fq \
        --max_memory 768G \
        --genome_guided_bam $RNASEQ_VS_REF_BAM \
        --genome_guided_max_intron 1000 \
        --left $R1 \
        --right $R2 \
        --SS_lib_type RF \
        --CPU 20
