#!/bin/bash
#SBATCH --job-name=PORKA
#SBATCH --output=porka_bybins_%A_%a.out
#SBATCH --error=porka_bybins_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --array=1-19
#SBATCH -p mem
#SBATCH -N 1
#SBATCH -n 40


module load miniconda
source activate megan


RUN=${SLURM_ARRAY_TASK_ID}
echo "My Slurm RUN_ID: '${RUN}'"
echo "My TMPDIR IS: " $TMPDIR

infile=$(ls filter_bins_metabat_checkM/*.fa | sed -n ${RUN}p)
echo "$infile"

outbase=$(basename $infile .fa)
outdirect=porka_bybins/${outbase}
echo "$outdirect"

prokka $infile --outdir $outdirect --prefix $outbase --addgenes --addmrna --centre X --compliant --genus --species --strain --kingdom --cpus 40
