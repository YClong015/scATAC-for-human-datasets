#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --job-name=cellatac
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-n

DATA=Eye
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [sampleIDs.txt]) # one sample ID per line in the text file

export PATH=/scratch/user/s4799891/cellranger-atac-2.2.0:$PATH # path to cellranger-atac-2.2.0

cd $TMPDIR
mkdir ${SAMPLE}

## Copy fastq files to scratch
cp -r /QRISdata/Q8448/Human_disease_data/${DATADIR}/${SAMPLE} $TMPDIR/${SAMPLE}

cellranger-atac count --id=${DATA}_${SAMPLE} \
                       --reference=/scratch/user/s4799891/refdata-cellranger-arc-GRCh38-2024-A \ # path to the reference directory
                       --fastqs=${TMPDIR}/${SAMPLE} \
                       --sample=${SAMPLE} \
                       --localcores=32 \
                       --localmem=256 \
                       --nosecondary

tar -zcvf ${DATA}_${SAMPLE}.tar.gz ${DATA}_${SAMPLE}
cp ${DATA}_${SAMPLE}.tar.gz /QRISdata/Q8448/Human_disease_cellatac/${DATADIR}/${DATA}_${SAMPLE}.tar.gz
