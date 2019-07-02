#!/bin/bash
#SBATCH --job-name=metabat
#SBATCH --output=anvi-metabat_%A_%a.out
#SBATCH --error=anvi-metabat_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH -p mem
#SBATCH -N 1
#SBATCH -n 40

module load miniconda
source activate metabat


runMetaBat.sh -m 5000 /project/megahit_output/final.contigs.fa /project/mapping/*.bam
