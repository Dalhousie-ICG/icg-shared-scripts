#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe threaded 40

source activate evidencemodeler

# add to evm perl scripts used below to $PATH
export PATH="$EVM_HOME/EvmUtils:$PATH"

# input
# note that the braker gff3 file is not the final gff3 output
# of braker, but rather the augustus.hints file
GENOME='ergo_cyp_genome.fasta.masked'
BRAKER_GFF3='augustus.hints.gff3'
PASA_GFF3='ergobibamus.sqlite.pasa_assemblies.gff3'

# validate GFF3s
gff3_gene_prediction_file_validator.pl $BRAKER_GFF3
gff3_gene_prediction_file_validator.pl $PASA_GFF3
# NOTE: if there are no issues with the GFF3 files, 
# these scripts will not yield any STDERR or STDOUT output
# so no output is good: the GFF3 files have been validated


# partition the genome and the GFF3 files
partition_EVM_inputs.pl \
	--genome $GENOME \
	--gene_predictions $BRAKER_GFF3 \
	--transcript_alignments $PASA_GFF3 \
	--segmentSize 1000000 \
	--overlapSize 200000 \
	--partition_listing partitions_list.out
# the genome sequences and gff3 files are partitioned into individual contigs
# the output is a <partitions_list.out> and a bunch of contig directoriescp
# so each output directory has one FASTA file (the contig), and the GFF3 files
# (the GFF3 entries for that contig).
# large contigs (> 1Mbp) are segmented into smaller overlapping chunks
# this would allow you to process the chunks in parallel

# be sure to specify the weights file in full path, otherwise it complains
# the verbose option -S does not get recognized :/

# input
WEIGHT_FILE=$(realpath weights.txt)
# the weight file is something like this
#
# ABINITIO_PREDICTION AUGUSTUS    1
# TRANSCRIPT  assembler-ergobibamus.sqlite    10
#
# 'AUGUSTUS' and 'assembler-ergobibamus.sqlite'
# refer to the 2nd column of the GFF3 files
#
# 1 and 10 are the weights assigned to each prediction


# write the evm commands used for parallelization 
# in the next step
write_EVM_commands.pl \
	--genome $GENOME \
	--gene_predictions $BRAKER_GFF3 \
	--transcript_alignments $PASA_GFF3 \
	--weights $WEIGHT_FILE \
	--output_file_name evm.out \
	--partitions partitions_list.out \
	> commands.list


# execute commands in parallel
parallel -j 40 < commands.list


# Combine <evm.out> for each contig chunk
# only for contigs that were partitioned
recombine_EVM_partial_outputs.pl \
	--partitions partitions_list.out \
	--output_file_name evm.out

# convert 'evm.out' to 'evm.out.gff3'
convert_EVM_outputs_to_GFF3.pl \
	--partitions partitions_list.out \
 	--output evm.out \
 	--genome $GENOME

conda deactivate

# combine all gff3 files into a single one, for all contigs
cat */evm.out.gff3 > EVM.all.gff3

