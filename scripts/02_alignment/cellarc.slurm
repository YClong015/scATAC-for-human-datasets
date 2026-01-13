#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --job-name=cellarc
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-n

DATA=Eye
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [sampleIDs.txt]) # one sample ID per line in the text file

export PATH="/scratch/user/s4799891/cellranger-arc-2.0.2":$PATH # path to cellranger-arc-2.0.2

## Copy fastq files to scratch
cd $TMPDIR

mkdir ${SAMPLE}_atac
cp -r /QRISdata/Q8448/Human_disease_data/Eye/${SAMPLE} $TMPDIR/${SAMPLE}_atac
mkdir ${SAMPLE}_gex
cp -r /QRISdata/Q8448/Human_disease_data/Eye_gex/${SAMPLE} $TMPDIR/${SAMPLE}_gex

cd $TMPDIR

echo "fastqs,sample,library_type" > Eye_libraries.csv
echo "${TMPDIR}/${SAMPLE}_atac,${SAMPLE},Chromatin Accessibility" >> Eye_libraries.csv
echo "${TMPDIR}/${SAMPLE}_gex,${SAMPLE},Gene Expression" >> Eye_libraries.csv

cellranger-arc count --id=${DATA}_${SAMPLE} \
                       --reference=/scratch/user/s4799891/refdata-cellranger-arc-GRCh38-2024-A \ # path to the reference directory
                       --libraries=/home/s4799891/dataLists/Eye_libraries.csv \
                       --localcores=32 \
                       --localmem=256 \
                       --no-bam

tar -zcvf ${DATA}_${SAMPLE}.tar.gz ${DATA}_${SAMPLE}
cp ${DATA}_${SAMPLE}.tar.gz /QRISdata/Q8448/Human_disease_cellatac/${DATA}/${DATA}_${SAMPLE}.tar.gz

