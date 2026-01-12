#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=DWLD_CNGB
#SBATCH --time=24:00:00
#SBATCH -o /home/s4896245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4896245/outputs/%x_%A_%a.err
#SBATCH --array=1-n

LINK=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [DataLinks.txt]) # one CNGB download link per line in the text file
CNGBID=$(echo "$LINK" | sed -n 's#.*/\([^/]*\)/[^/]*$#\1#p') # extract the ID from the full link
DATADIR=/QRISdata/Q8448/Human_disease_data/Muscle

cd $TMPDIR
mkdir ${CNGBID}_Fastq
cd ${CNGBID}_Fastq

${LINK}

gzip *

mkdir ${DATADIR}/${CNGBID}
cp * ${DATADIR}/${CNGBID}
