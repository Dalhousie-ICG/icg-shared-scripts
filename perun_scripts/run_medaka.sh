#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q 144G-batch,256G-batch
#$ -pe threaded 8

source activate medaka

READS='eb_flongle_reads_pooled.fastq.gz'
REFERENCE='canu_ergo_contigs.fasta'
OUTDIR='medaka_out'
THREADS=8
MODEL='r941_min_sup_g507'

medaka_consensus \
	-i $READS \
	-d $REFERENCE \
	-o $OUTDIR \
	-t 8 \
	-m $MODEL

conda deactivate
