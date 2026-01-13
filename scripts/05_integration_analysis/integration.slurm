#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=512G
#SBATCH --job-name=integration
#SBATCH --time=16:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-n

#export IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /home/s4799891/dataLists/Muscle_inds.txt)
#export DIR="/QRISdata/Q8448/Human_disease_Signac/Muscle/SObjects/universal_SCPM1/integration/"

module load r/4.4.2

cd $TMPDIR
Rscript /home/s4799891/scripts/integration.R

#mv $TMPDIR/* $DIR
