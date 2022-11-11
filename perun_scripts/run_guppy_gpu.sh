#!/bin/bash
#SBATCH --account=def-roger
#SBATCH --gres=gpu:p100:1
#SBATCH --time=1-00:00:00
#SBATCH --mem=2000M

INPUT_DIR='fast5'
OUTPUT_DIR='01_guppy-5.0.16/fastq'
#FLOWCELL=FLO-FLG001
#KIT=SQK-LSK109

# config file for super high accuracy model
CFG='/home/jmartijn/tools/ont-guppy-5.0.16/data/dna_r9.4.1_450bps_sup.cfg'

/home/jmartijn/tools/ont-guppy-5.0.16/bin/guppy_basecaller \
    -i $INPUT_DIR \
    --recursive \
    -s $OUTPUT_DIR \
    --config $CFG \
    -x "cuda:0" \
    --records_per_fastq 0 \
    --calib_detect \
    --trim_strategy none


#    --flowcell $FLOWCELL \
#    --kit $KIT \
