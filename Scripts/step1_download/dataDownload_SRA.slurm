#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=DWLD
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4799891/outputs/%x_%A_%a.out
#SBATCH -e /home/s4799891/outputs/%x_%A_%a.err
#SBATCH --array=1-n

DATA=Eye
SRAID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [DataIDs.txt]) # one SRA ID per line in the text file
DATADIR=/QRISdata/Q8448/Human_disease_data/${DATA}

module load sra-toolkit

cd $TMPDIR
mkdir ${SRAID}_Fastq
cd ${SRAID}_Fastq

OUTDIR=${TMPDIR}/${SRAID}_Fastq/${SRAID}
#prefetch ${SRAID} -O ${TMPDIR}/${SRAID}_Fastq/${SRAID}
fasterq-dump -e 8 -O ${OUTDIR} --progress --split-files --include-technical $SRAID

cd $SRAID
gzip *

cd ../../
tar -zcvf ${SRAID}_Fastq.tar.gz ${SRAID}_Fastq

cp ${SRAID}_Fastq.tar.gz ${DATADIR}/${SRAID}_Fastq.tar.gz
