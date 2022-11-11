#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q 768G-batch
#$ -m bea
#$ -M joran.martijn@dal.ca
#$ -pe threaded 20

source activate braker2

# setting relevant environment variables, if necessary
#
## braker requires gmes_petap.pl, a GeneMark tool. 
## braker will look for it under $PATH or $GENEMARK_PATH
## GENEMARK_PATH seems to be already set to /scratch2/software/gmes_linux_4.69/ by default
## which contains gmes_petap.pl
## Another option:
## export GENEMARK_PATH=/scratch2/software/gmes_linux_64-aug-2020/
#
## as braker requires GeneMark, you will also need a valid .gm_key
## it will expire after 200 days
## download from http://topaz.gatech.edu/genemark/license_download.cgi
## and unzip it in your $HOME
#
## braker needs writing access to the config/ directory of AUGUSTUS
## this directory is located at /scratch2/software/anaconda/envs/braker2/config/
## when using the braker2 environment. Regular users do not have writing access here
## you will therefore need to copy this directory to a place where you do have writing access:
##
## cp -r /scratch2/software/anaconda/envs/braker2/config/ ~/tools/augustus/
##
## braker will look for this directory using $AUGUSTUS_CONFIG_PATH
## so set this environment variable to where you copied the config directory on your home
export AUGUSTUS_CONFIG_PATH='/home/jmartijn/tools/augustus/config'
#
## braker will look for the augustus binaries and scripts under $PATH
## if it can't find them, it will try to guess them from your $AUGUSTUS_CONFIG_PATH
## under (braker2) environment they are in your $PATH so no need to configure anything extra here
## alternatively, set
## export AUGUSTUS_BIN_PATH=/scratch2/software/anaconda/envs/braker2-2.1.6/bin
## export AUGUSTUS_SCRIPTS_PATH=/scratch2/software/anaconda/envs/braker2-2.1.6/bin
#
## braker also need python3, and will look for it under $PATH
## under (braker2) environment, python3 is under $PATH
#
## braker also need bamtools, and will look for it under $PATH
## under (braker2) environment, bamtools is under $PATH
#
## braker also needs cdbfasta
## under (braker2) environment, cdbfasta is under $PATH
#
## braker also needs ncbiblast or diamond
## under (braker2), diamond is already under $PATH
#
## braker also optionally uses samtools
## under (braker2), samtools is already under $PATH


# input
ORIGINAL_GENOME='ergo_cyp_genome.fasta.masked'
RNASEQ='rnaseq_vs_masked_ergo_cyp_genome.sort.bam'
SPECIES='ergobibamus'

# run braker
braker.pl \
        --species=$SPECIES \
        --useexisting \
        --gff3 \
        --softmasking \
        --cores=20 \
        --genome=$ORIGINAL_GENOME \
        --bam=$RNASEQ 

# options explained
##
## --useexisting        use existing genemark training set if it exists
##                      the training set is created with --species='species'
##                      so essentially if a run fails at a later stage, and you rerun it,
##                      it won't have to recreate the training set
##
## --gff3               output in GFF3 instead of GTF (=GFF2) format
##
## --softmasking        input genome is softmasked (repetitive regions are in lower case)
##
## --cores=20           use 20 cores


# relevant information from the online manual
##
## Repeat Masking
## In order to predict genes accurately in a novel genome, the genome should be masked for repeats.
## This will avoid the prediction of false positive gene structures in repetitive and low complexitiy regions.
## Repeat masking is also essential for mapping RNA-Seq data to a genome with some tools (other RNA-Seq mappers, such as HISAT2, ignore masking information).
## In case of GeneMark-EX and AUGUSTUS, softmasking (i.e. putting repeat regions into lower case letters and all other regions into upper case letters) 
## leads to better results than hardmasking (i.e. replacing letters in repetitive regions by the letter N for unknown nucleotide).
## If the genome is masked, use the --softmasking flag of braker.pl.

