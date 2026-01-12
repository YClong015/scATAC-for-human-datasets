#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --job-name=Sobject
#SBATCH --time=8:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-22

export PEAKPATH="/QRISdata/Q8448/Human_disease_MACS/consensus_regions_SCPM1.bed" # universal peak filtered by SCPM >= 1
export SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [Samples.txt]) # one sample name per line for the text file
export OUTDIR="/QRISdata/Q8448/Human_disease_Signac/Heart/SObjects/universal_SCPM1/"

module load r/4.4.2

cd $TMPDIR
Rscript /home/s4799891/scripts/SeuratObject.R

mv $TMPDIR/* ${OUTDIR}
