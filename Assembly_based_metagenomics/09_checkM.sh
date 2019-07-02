#!/bin/bash
#SBATCH --job-name=CHECKM
#SBATCH --output=anvi-metabat_%A_%a.out
#SBATCH --error=anvi-metabat_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 40

module load miniconda
#source activate metabat
module load checkm


checkm lineage_wf -t 40 -x fa final.contigs.fa.metabat-bins5000/ /projectq/checkM
