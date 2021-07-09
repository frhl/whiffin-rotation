#!/bin/bash

#$ -wd /well/lindgren/flassen/projects/utrs/whiffin-rotation
#$ -N GO
#$ -o logs/GO.log
#$ -e logs/GO.errors.log
#$ -q short.qe
#$ -P lindgren.prjc
#$ -pe shmem 1

module load Anaconda3/2020.07
module load java/1.8.0_latest
module load R

Rscript workflows/210707_UTR_GO_analysis.R 

