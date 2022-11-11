#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe threaded 24

# load repeatmodeler
## the repeatmodeler environment actually contains 
## repeatmodeler AND repeatmasker
## so no need to seperately activate repeatmasker environment
source activate repeatmodeler

# build database
ASSEMBLY="ergo_cyp_genome.fasta"
NAME=${ASSEMBLY%.fasta}
BuildDatabase -name $NAME $ASSEMBLY

# run repeatmodeler
JOBS=6
RepeatModeler -pa $JOBS -database $NAME > repeatmodeler.out

# run repeatmasker
OUTDIR='repeatmasker_out'
RepeatMasker \
    -engine rmblastn \
    -parallel $JOBS \
    -xsmall -s -nolow -a -inv \
    -gff -html -ace -source \
    -dir $OUTDIR \
    -lib ${NAME}-families.fa \
    $ASSEMBLY

conda deactivate

##########################
## RepeatMasker options ##
##########################

# -parallel defines the number of parallel jobs to be run.
## note that 1 job requires 4 cores!

# -xsmall returns repetitive regions in lower case letters (actg) instead of a sequence of N's (i.e. masked)
## lowercase repetitive regions = soft-masking
## NNNNNNNNN repetitive regions = hard-masking
## tools like hisat2 will treat soft-masked regions as regular unmasked regions!

# -s invokes slow search, the most sensitive setting

# -nolow does not mask low complexity DNA regions or simple repeats

# -a invokes writing alignments in .align output file

# -inv, in conjunction with -a, ensures that alignments are presented in the orientation of the repeat

# -dir defines output directory

# -lib defines the custom library. In our case, the output of RepeatModeler

# -gff, -html, -ace and create respective output files


#########################
## RepeatMasker output ##
#########################

# The main output is the ${ASSEMBLY}.out file
## It contains more or less the same information as the GFF output file, but additionally has the repeat class/family annotated
## To view this in IGV we would need to conver this .out file to a GFF3 file

# The gff output is in GFF2 format
## this means that, amongst other things, the second column is not 'type' (as in GFF3), but 'method'
## terms in this field are not Sequence Ontology compliant (as in GFF3)

# The html output seems to be a funny interactive html page where you can see alignments for each inferred repeat feature
