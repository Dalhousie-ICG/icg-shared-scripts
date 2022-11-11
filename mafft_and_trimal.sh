#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#-pe threaded 2

mafft --auto --thread 2 name.fas > name.aligned.fas 
trimal -in name.aligned.fas -out name.trimal -gappyout
