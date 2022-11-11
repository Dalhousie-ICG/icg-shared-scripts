#!/bin/bash
#$ -S /bin/bash
# . /etc/profile
#$ -cwd
#$ -m bea
#$ -M joran.martijn@dal.ca
#$ -pe threaded 20

# input
READS='eb_flongle_reads_pooled.fastq.gz'
REF='canu_all_contigs.fasta'

# output
SAM='pooled_vs_all_contigs.sam'

source activate ngmlr

# align reads
ngmlr \
	--query $READS \
	--reference $REF \
	--threads 20 \
	--output $SAM \
	--presets ont

BASE=${SAM%.sam}
# convert sam to bam
samtools view --threads 20 -b $BASE.sam > $BASE.bam
samtools sort --threads 20 $BASE.bam > $BASE.sort.bam
samtools index $BASE.sort.bam $BASE.sort.bam.bai
