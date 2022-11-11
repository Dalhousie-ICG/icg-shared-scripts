#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe threaded 20
#$ -q 144G-batch,256G-batch

source activate pilon

# Be sure to allocate at least 8GB for pilon
# Pilon is a java program, so you can use the -Xmx<memory> option;
# 16384m equals 16GB 

ASSEMBLY='canu_medaka_pilon.fasta'
PAIRED_BAM='dnaseq_vs_canu_medaka_pilon_polish-pe.sort.bam'
UNP_FW_BAM='dnaseq_vs_canu_medaka_pilon_polish-se1.sort.bam'
UNP_RV_BAM='dnaseq_vs_canu_medaka_pilon_polish-se2.sort.bam'
OUTDIR='pilon_out'
OUTPREFIX='canu_medaka_pilon_x2_polish'
THREADS=20

pilon -Xmx16384m \
	--genome $ASSEMBLY \
	--frags $PAIRED_BAM \
	--unpaired $UNP_FW_BAM \
	--unpaired $UNP_RV_BAM \
	--output $OUTPREFIX \
	--outdir $OUTDIR \
	--changes \
	--vcf \
	--threads $THREADS \
	--tracks

conda deactivate
