#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --job-name=Dobject
#SBATCH --time=8:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-n

export IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [inds.txt]) # one donor ID per line
export DIR="/QRISdata/Q8448/Human_disease_Signac/Muscle/Sobject/universal_SCPM1"

module load r/4.4.2

cd $TMPDIR
Rscript /home/s4799891/scripts/DonorObject.R

mv $TMPDIR/* ${DIR}
