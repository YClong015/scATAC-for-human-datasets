#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=Header_modify
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-n

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [Liver_Samples.txt])

cd $TMPDIR
mkdir ${SAMPLE}
cd ${SAMPLE}

zcat /QRISdata/Q8448/Human_disease_data/Liver/${SAMPLE}/${SAMPLE}_S1_L001_R1_001.fastq.gz | sed 's/^@\([^ ]*\) \([ACTG]\+\):/@\1_\2:/' | gzip > ${TMPDIR}/${SAMPLE}/${SAMPLE}_R1_modified.fastq.gz
zcat /QRISdata/Q8448/Human_disease_data/Liver/${SAMPLE}/${SAMPLE}_S1_L001_R2_001.fastq.gz | sed 's/^@\([^ ]*\) \([ACTG]\+\):/@\1_\2:/' | gzip > ${TMPDIR}/${SAMPLE}/${SAMPLE}_R2_modified.fastq.gz

cp ${SAMPLE}_R1_modified.fastq.gz /QRISdata/Q8448/Human_disease_data/Liver/${SAMPLE}
cp ${SAMPLE}_R2_modified.fastq.gz /QRISdata/Q8448/Human_disease_data/Liver/${SAMPLE}
