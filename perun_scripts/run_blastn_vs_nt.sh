#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe threaded 10

QUERY=Ergo_cyp_seqs_list.fasta
DB=/db1/blastdb-sep1-2021/nt
THREADS=10
OUT=Ergo_cyp_seqs_vs_nt.blastn

source activate blast

blastn \
    -query $QUERY \
    -db $DB \
    -evalue 1e-5 \
    -num_threads $THREADS \
    -outfmt '6 std qcovhsp sskingdoms stitle' \
    -out $OUT

sort -k1,1 -k12,12gr -k11,11g -k3,3gr $OUT | sort -u -k1,1 > ${OUT/.blastn/.tophits.blastn}

conda deactivate
