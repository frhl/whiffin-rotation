#!/bin/bash

#$ -wd /well/lindgren/flassen/projects/utrs/whiffin-rotation
#$ -N submit_shuffle_5utr
#$ -o logs/submit_shuffle_5utr.log
#$ -e logs/submit_shuffle_5utr.errors.log
#$ -q short.qe
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -t 1-10


set +eu
module load Anaconda3/2020.07
module load java/1.8.0_latest
source "/apps/eb/skylake/software/Anaconda3/2020.07/etc/profile.d/conda.sh"
conda activate reticulate
set -eu

readonly array_replicate="${SGE_TASK_ID}"

readonly rscript="workflows/shuffle_utrs.R"

readonly path_mane="/well/lindgren/flassen/projects/210629_MANE.GRCh38.v0.95.combined-table.txt"
readonly path_features="/well/lindgren/flassen/projects/utrs/whiffin-rotation/derived/tables/210629_MANE.v0.95.UTR_features.txt"
readonly subset="five_prime_UTR"
readonly iterations=1000


readonly out_dir="/well/lindgren/flassen/projects/utrs/whiffin-rotation/derived/run001"
readonly out_prefix="${out_dir}/test5"

mkdir -p ${out_dir}

Rscript ${rscript} \
  --path_mane "${path_mane}" \
  --path_features "${path_features}" \
  --subset "${subset}" \
  --array_replicate "${array_replicate}" \
  --iterations "${iterations}" \
  --out_prefix "${out_prefix}"

