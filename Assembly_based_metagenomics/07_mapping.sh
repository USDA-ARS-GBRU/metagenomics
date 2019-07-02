#!/bin/bash

#SBATCH --job-name=bowtie_map
#SBATCH --output=bowtie_map_%A_%a.out
#SBATCH --error=bowtie_map_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --array=1-12
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 40

module load bowtie2/2.3.4
module load samtools

rootdir="/project/"
echo "$rootdir"

RUN=${SLURM_ARRAY_TASK_ID}
echo "My Slurm RUN_ID: '${RUN}'"

iarray=("$rootdir"rqc_data/merged_fastqByRelicates/*_interleaved_merge.fq.gz)
echo "$iarray"

infile="${iarray[$SLURM_ARRAY_TASK_ID - 1]}"
echo "$infile"

bname=`basename $infile _interleaved_merge.fq.gz`
echo "$bname"

time bowtie2 -p 40 -x megahit_output/contigs --interleaved $infile -S "$rootdir"mapping/"$bname".sam
time samtools sort -o "$rootdir"mapping/"$bname".bam "$rootdir"mapping/"$bname".sam
