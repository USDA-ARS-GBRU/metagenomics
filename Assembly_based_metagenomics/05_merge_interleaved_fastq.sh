#!/bin/bas
#SBATCH --job-name=merge
#SBATCH --output=merge_%A_%a.out
#SBATCH --error=merge_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 40

cat reads/A1*.fq.gz > merged_fastqByRelicates/A1_interleaved_merge.fq.gz
cat reads/A2*.fq.gz > merged_fastqByRelicates/A2_interleaved_merge.fq.gz
cat reads/A3*.fq.gz > merged_fastqByRelicates/A3_interleaved_merge.fq.gz
cat reads/M1*.fq.gz > merged_fastqByRelicates/M1_interleaved_merge.fq.gz
cat reads/M2*.fq.gz > merged_fastqByRelicates/M2_interleaved_merge.fq.gz
cat reads/M3*.fq.gz > merged_fastqByRelicates/M3_interleaved_merge.fq.gz
cat reads/O1S*.fq.gz > merged_fastqByRelicates/O1S_interleaved_merge.fq.gz
cat reads/O2S*.fq.gz > merged_fastqByRelicates/O2S_interleaved_merge.fq.gz
cat reads/O3S*.fq.gz > merged_fastqByRelicates/O3S_interleaved_merge.fq.gz
cat reads/C1W*.fq.gz > merged_fastqByRelicates/C1W_interleaved_merge.fq.gz
cat reads/C2W*.fq.gz > merged_fastqByRelicates/C2W_interleaved_merge.fq.gz
