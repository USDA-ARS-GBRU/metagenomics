#!/bin/bash

#SBATCH --job-name=index_contig
#SBATCH --output=index_contig_%A_%a.out
#SBATCH --error=index_contig_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 40

module load bowtie2/2.3.4
module load samtools


bowtie2-build megahit_output/final.contigs.fa megahit_output/contigs
