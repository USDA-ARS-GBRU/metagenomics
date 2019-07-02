#!/bin/bash
#SBATCH --job-name=rqc_interleave_Ametagenomics
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
module load bbtools/gcc/64/36.86

RUN=${SLURM_ARRAY_TASK_ID}
echo "My Slurm RUN_ID: '${RUN}'"
echo "My TMPDIR IS: " $TMPDIR

infile1=$(ls data/*_R1.fastq.gz | sed -n ${RUN}p)
echo "$infile1"

infile2=$(echo $infile1 | sed 's/R1/R2/')
echo "$infile2"

outbase=$(basename $infile1 _R1.fastq.gz)
outfile=data/${outbase}_interleave.fastq.gz
echo "$outfile"

reformat.sh in=$infile1 in2=$infile2 out=$outfile
