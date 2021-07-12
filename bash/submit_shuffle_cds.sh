#!/bin/bash

#$ -wd /well/lindgren/flassen/projects/utrs/whiffin-rotation
#$ -N CDS
#$ -o logs/CDS.log
#$ -e logs/CDS.errors.log
#$ -q short.qe
#$ -P lindgren.prjc
#$ -pe shmem 1

module load Anaconda3/2020.07
module load java/1.8.0_latest
module load R

Rscript workflows/210703_shuffle_CDS.R 

