#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --job-name=doublet_filtering
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --account=a_nefzger
#SBATCH -o /home/s4799891/outputs/%x_%A_%a.out
#SBATCH -e /home/s4799891/outputs/%x_%A_%a.err
#SBATCH --array=1-n

export DIR="/QRISdata/Q8448/Human_disease_Signac/Muscle/SObjects/universal_SCPM1/doublet_removed/"
export SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /home/s4799891/dataLists/Muscle_Samples.txt)
export IND=$(echo "${SAMPLE}" | cut -d "_" -f 2)
export FCOUNT=$(grep "${IND}" /home/s4799891/dataLists/Muscle_fragCounts.txt | cut -d " " -f 2)

module load r/4.4.2

cd $TMPDIR

Rscript /home/s4799891/scripts/doublet_filtering.R

mv $TMPDIR/* ${DIR}
