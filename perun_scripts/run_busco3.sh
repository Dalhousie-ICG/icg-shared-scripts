#!/bin/bash
#$ -S /bin/bash
# . /etc/profile
#$ -cwd
#$ -pe threaded 8
#$ -q 144G-batch,256G-batch

source activate busco-3

INPUT='canu_ergo_contigs_clean.fasta'
OUTDIR='toxoplasma_busco3_odb9_long_out'
MODE='genome'
TRAINING_SPECIES='toxoplasma'

# busco v3 only works with odb9 databases
LINEAGEDB='/home/dsalas/Shared/BUSCO/eukaryota_odb9'

# in the busco-3 environment, AUGUSTUS_CONFIG_PATH is set to
# /scratch2/software/anaconda/envs/busco-3/config/
# but we don't have writing permissions there
# not sure why we need writing permissions but it doesnt work anyway
# but we copied that dir to a place where we do have writing permissions:
export AUGUSTUS_CONFIG_PATH="/home/jasons/tools/busco-3/config/"

# run busco
## do not specify output dir with a trailing slash, it will lead to a fatal error
## modes are genome, proteins, transcriptome
## the below command will use Augustus as gene predictor,
## using a specified species gene finding parameters $TRAINING_SPECIES
## if you want to do self-training, specify --long instead
run_BUSCO.py \
    --in $INPUT \
    --out $OUTDIR \
    --mode $MODE \
    --lineage_path $LINEAGEDB \
    --cpu 8 \
    --species $TRAINING_SPECIES \
#    --long

conda deactivate
