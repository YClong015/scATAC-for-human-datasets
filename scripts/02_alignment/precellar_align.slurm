#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --job-name=Precellar_align
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-n

module load anaconda3/2023.09-0
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate precellar # conde environment with Python 3.10.18

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [sampleIDs.txt]) # one sample ID per line for the text file

cd $TMPDIR

# copy fastq files to the tmp dir
cp /QRISdata/Q8448/Human_disease_data/Liver/${SAMPLE}_processed/${SAMPLE}_S1_L001_R1_001.fastq.gz ${TMPDIR}/R1_processed.fq.gz
cp /QRISdata/Q8448/Human_disease_data/Liver/${SAMPLE}_processed/${SAMPLE}_S1_L001_R2_001.fastq.gz ${TMPDIR}/I1.fq.gz
cp /QRISdata/Q8448/Human_disease_data/Liver/${SAMPLE}_processed/${SAMPLE}_S1_L001_R3_001.fastq.gz ${TMPDIR}/R2_processed.fq.gz

mkdir ${TMPDIR}/${SAMPLE}
export OUTDIR="${TMPDIR}/${SAMPLE}"

python /home/s4799891/scripts/precellar_align.py

mkdir /QRISdata/Q8448/Human_disease_cellatac/Liver_Precellar/Liver_${SAMPLE}
cp ${TMPDIR}/${SAMPLE}/* /QRISdata/Q8448/Human_disease_cellatac/Liver_Precellar/Liver_${SAMPLE}

conda deactivate
