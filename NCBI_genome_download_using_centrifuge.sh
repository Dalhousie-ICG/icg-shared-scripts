#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -o centrifuge-download_Protozoa.log
#$ -pe threaded 10
cd $PWD
source activate centrifuge

echo "download genbank protozoa genomes @ Any assembly level"
centrifuge-download -o protozoa_genbank -P 10 -a 'Any' -d "protozoa" genbank > seqid2taxid_protozoa.map

echo "download genbank protozoa genomes @ complete assembly level"
centrifuge-download -o protozoa_genbank_complete -P 10  -d "protozoa" genbank > seqid2taxid_complete_protozoa.map

echo "download refseq protozoa genomes @ complete assembly level"
centrifuge-download -o protozoa_refseq_complete -P 10  -d "protozoa" refseq > seqid2taxid_refseq_protozoa.map

echo "download NCBI taxonomy"
centrifuge-download -o taxonomy taxonomy

conda deactivate

# >source activate centrifuge
# >centrifuge-download
#
# centrifuge-download [<options>] <database>
# 
# ARGUMENT
#  <database>        One of refseq, genbank, contaminants or taxonomy:
#                      - use refseq or genbank for genomic sequences,
#                      - contaminants gets contaminant sequences from UniVec and EmVec,
#                      - taxonomy for taxonomy mappings.
# 
# COMMON OPTIONS
#  -o <directory>         Folder to which the files are downloaded. Default: '.'.
#  -P <num of threads>      Number of processes when downloading (uses xargs). Default: '1'
# 
# WHEN USING database refseq OR genbank:
#  -d <domain>            What domain to download. One or more of bacteria, viral, archaea, fungi, protozoa, invertebrate, plant, vertebrate_mammalian, vertebrate_other (comma separated).
#  -a <assembly level>    Only download genomes with the specified assembly level. Default: 'Complete Genome'. Use 'Any' for any assembly level.
#  -c <refseq category>   Only download genomes in the specified refseq category. Default: any.
#  -t <taxids>            Only download the specified taxonomy IDs, comma separated. Default: any.
#  -g <program>           Download using program. Options: rsync, curl, wget. Default curl (auto-detected).
#  -r                     Download RNA sequences, too.
#  -u                     Filter unplaced sequences.
#  -m                     Mask low-complexity regions using dustmasker. Default: off.
#  -l                     Modify header to include taxonomy ID. Default: off.
#  -g                     Download GI map.
#  -v                     Verbose mode
