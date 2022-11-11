#!/bin/bash                                                                                                                         
#$ -S /bin/bash
#$ -cwd
#$ -pe threaded 20

# input
CONFIG='pasa.config'
GENOME='ergo_cyp_genome.fasta.masked'
TRANSCRIPTOME='Trinity-GG.fasta'
THREADS=20

source activate pasa-2.5.2

# run pasa
Launch_PASA_pipeline.pl \
        --create --run \
        -c $CONFIG \
        -g $GENOME \
        -t $TRANSCRIPTOME \
        --transcribed_is_aligned_orient \
        --ALIGNERS blat,gmap,minimap2 \
        --CPU $THREADS

conda deactivate
