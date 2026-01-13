#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --job-name=Precellar_barcode
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-n

module load anaconda3/2023.09-0
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate precellar # conda environment with Python 3.10.18

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [sampleIDs.txt]) # one sample ID [er line for the text file


cd $TMPDIR
cp -r /QRISdata/Q8448/Human_disease_data/Liver/${SAMPLE} $TMPDIR/
mv ${SAMPLE}/${SAMPLE}_R1_modified.fastq.gz R1.fq.gz
mv ${SAMPLE}/${SAMPLE}_R2_modified.fastq.gz R2.fq.gz

python /home/s4799891/scripts/precellar_stripBarcode.py

cp ${TMPDIR}/R1_processed.fq.gz /QRISdata/Q8448/Human_disease_data/Liver/${SAMPLE}_processed/${SAMPLE}_S1_L001_R1_001.fastq.gz
cp ${TMPDIR}/barcode.fq.gz /QRISdata/Q8448/Human_disease_data/Liver/${SAMPLE}_processed/${SAMPLE}_S1_L001_R2_001.fastq.gz
cp ${TMPDIR}/R2_processed.fq.gz /QRISdata/Q8448/Human_disease_data/Liver/${SAMPLE}_processed/${SAMPLE}_S1_L001_R3_001.fastq.gz

conda deactivate
