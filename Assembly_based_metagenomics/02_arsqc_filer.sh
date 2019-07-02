#!/bin/bash
#SBATCH --job-name=rqc_filter_Ametagenomics
#SBATCH --output=rqc_Ametagenomics_%A_%a.out
#SBATCH --error=rqc_Ametagenomics_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --array=1-156
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 10

module load miniconda
source activate Ametagenomics
module load pigz


RUN=${SLURM_ARRAY_TASK_ID}
echo "My Slurm RUN_ID: '${RUN}'"
echo "My TMPDIR IS: " $TMPDIR

rootdir="/project/"
echo "$rootdir"


itl_file=$(ls data/*_interleave.fastq.gz | sed -n ${RUN}p)
echo "$itl_file"
time rqcfilter.py --fastq $itl_file --output "$rootdir"rqc_data -p
