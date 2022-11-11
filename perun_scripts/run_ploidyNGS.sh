#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -q 256G-batch

# ploidyNGS uses quite a bit of memory, so that is why
# I specified a large memory node above

source activate ploidyNGS-dependencies
# add ploidyNGS install location to PATH so it can find the R histogram script as well
export PATH="/scratch2/software/ploidyNGS:$PATH"

## ploidyNGS will create an indexed BAM file if it is not in the same directory as the BAM file
## but in the source code it says #TODO: check that indexing worked OK
## so maybe its best to create and provide the indexed BAM file yourself with `samtools index`
BAMFILE='dnaseq_vs_canu_ergo_assembly.sorted.bam'
if [[ ! -e "${BAMFILE}.bai" ]]; then
    echo "No indexed BAM file found, exiting"
    exit
fi

## ploidyNGS will also consider regions with very low coverage (i.e. 1x, 2x, etc)
## random sequence errors in reads of low coverage regions could explain 50/50, 33/67 etc allele-frequencies
## it may be a good idea to only use higher coverage regions when estimating ploidy
MINCOV=10

## read each alignment column until a maxdepth is reached
## default is 100, but I found setting it much higher gives
## much more accurate results, at little cost to memory usage
MAXDEPTH=2000

# basename for outfiles
OUTBASE='dnaseq_vs_canu_ergo_assembly'

# ploidyNGS_minCov.py is an edit of ploidyNGS.py

## that allows one to set a minimum coverage threshold
## for a alignment position to be considered

## it also makes the reported positions in the .tab file 1-indexed
## (the original script reported 0-indexed, which can be confusing)

## it will also print out to the STDOUT the total number of
## heteromorphic positions
ploidyNGS_minCov.py \
    --out $OUTBASE \
    --bam $BAMFILE \
    --min_cov $MINCOV \
    --max_depth $MAXDEPTH

conda deactivate

# some more comments about ploidyNGS

## the python script relies heavily on the 'pysam', a python module for working with SAM/BAM files etc
## it calls the .pileup() function, which by default ignores read mappings with bitwise FLAGs
## 4    0x4   BAM_FUNMAP       - unmapped reads
## 256  0x100 BAM_FSECONDARY   - not primary read mapping (the read has a better mapping, aka the primary read mapping)
## 512  0x200 BAM_FQCFAIL      - read mapping fails some kind of quality control
## 1024 0x400 BAM_FDUP         - read mapping is a PCR or optical duplicate
## so many read mappings in the SAM/BAM file that you can see in Tablet, may not actually be taken into account by ploidyNGS
