#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=Rename
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-n

NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [SampleNames.txt] | cut -f 1) # first field of SampleNames.txt per line is the sample ID
SRAID1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [SampleNames.txt] | cut -f 2) #second field of SampleNames.txt per line is the SRAID


cd /QRISdata/Q8448/Human_disease_data/Brain_Neuron/

tar -xvf ${SRAID1}_Fastq.tar.gz
mkdir ${NAME}
mv ${SRAID1}_Fastq/${SRAID1}/${SRAID1}_1.fastq.gz ${NAME}/${NAME}_S1_L001_I1_001.fastq.gz
mv ${SRAID1}_Fastq/${SRAID1}/${SRAID1}_2.fastq.gz ${NAME}/${NAME}_S1_L001_R1_001.fastq.gz
mv ${SRAID1}_Fastq/${SRAID1}/${SRAID1}_3.fastq.gz ${NAME}/${NAME}_S1_L001_R2_001.fastq.gz
mv ${SRAID1}_Fastq/${SRAID1}/${SRAID1}_4.fastq.gz ${NAME}/${NAME}_S1_L001_R3_001.fastq.gz
rm -r ${SRAID1}_Fastq

#cd ${SRAID1}
#mv *_1.fq.gz ../${NAME}/${NAME}_S1_L001_R1_001.fastq.gz
#mv *_2.fq.gz ../${NAME}/${NAME}_S1_L001_R2_001.fastq.gz
#cd ../
#rm -r ${SRAID1}
