#!/bin/bash
#$ -S /bin/bash
# . /etc/profile
#$ -cwd
#$ -m bea
#$ -pe threaded 8
#$ -q 144G-batch,256G-batch

# I specified the high memory nodes, because I had a memory error
# with the makeblastdb step when it happened to run on a 16G node.
# It may have been a fluke but to keep it safe I just specified
# the high memory nodes

source activate busco-5

# in the busco-5 environment, AUGUSTUS_CONFIG_PATH is set to
# /scratch2/software/anaconda/envs/busco-5/config/
# but we don't have writing permissions there
# not sure why we need writing permissions but it doesnt work anyway
# but we copied that dir to a place where we do have writing permissions:
# you may want to copy it to your own home
export AUGUSTUS_CONFIG_PATH="/home/jasons/tools/busco-5/config/"

INPUT='canu_ergo_contigs_clean.fasta'
OUTDIR='busco5_augustus_leishmania_out'
MODE='genome'
AUGUSTUS_SPECIES='leishmania_tarentolae'

# setting the lineage db
## the latest version (as of writing) is odb10
LINEAGEDB='/scratch3/rogerlab_databases/other_dbs/BUSCO/lineages/eukaryota_odb10/'
## busco v5 only works with odb10
## it will not work with odb9

# run busco
## do not specify output dir with a trailing slash, it will lead to a fatal error
## modes are genome, proteins, transcriptome
## the below command will use Augustus as gene predictor, Leishmania tarentolae as species
## use --long instead of --augustus_species if you want to do self-training
busco \
    --in $INPUT \
    --out $OUTDIR \
    --mode $MODE \
    --lineage_dataset $LINEAGEDB \
    --cpu 8 \
    --augustus \
    --augustus_species $AUGUSTUS_SPECIES \
#    --long

conda deactivate
