#!/bin/bash
#$ -S /bin/bash
# . /etc/profile
#$ -cwd
#$ -m bea
#$ -M joran.martijn@dal.ca
#$ -pe threaded 8

source activate busco-5

INPUT='canu_ergo_contigs_clean.fasta'
OUTDIR='busco5_out'
MODE='genome'


# setting the lineage db
## the latest busco db for eukaryota is odb10
LINEAGEDB='/scratch3/rogerlab_databases/other_dbs/BUSCO/lineages/eukaryota_odb10/'
## busco v5 only works with odb10
## it will not work with odb9


# run busco
## do not specify output dir with a trailing slash, it will lead to a fatal error
## modes are genome, proteins, transcriptome
## the below command will use Metaeuk as gene predictor
busco \
    --in $INPUT \
    --out $OUTDIR \
    --mode $MODE \
    --lineage_dataset $LINEAGEDB \
    --cpu 8

conda deactivate
