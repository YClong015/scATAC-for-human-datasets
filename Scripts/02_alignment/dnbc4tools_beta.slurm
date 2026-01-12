#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --job-name=dnbc4tools_beta
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-n

DATA=Muscle
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [sampleIds.txt]) # one sample per line for the text file

cd $TMPDIR
## Copy fastq files to scratch
cp /QRISdata/Q8448/Human_disease_data/${DATA}/${SAMPLE}/* ${TMPDIR}/
mv ${TMPDIR}/${SAMPLE}_S1_L001_R1_001.fastq.gz ${TMPDIR}/${SAMPLE}_R1.fastq.gz
mv ${TMPDIR}/${SAMPLE}_S1_L001_R2_001.fastq.gz ${TMPDIR}/${SAMPLE}_R2.fastq.gz

/scratch/user/s4799891/dnbc4tools3.0beta_v3/dnbc4tools atac run \ # path to dnbc4tools3.0beta_v3/dnbc4tools
		--name ${SAMPLE} \
		--fastq1 ${TMPDIR}/${SAMPLE}_R1.fastq.gz \
		--fastq2 ${TMPDIR}/${SAMPLE}_R2.fastq.gz \
                --customize "cb,R2:51-60;cb,R2:61-70;R1,R1:1-50;R2,R2:1-50" \ # specify the locations of cell barcode, R1 and R2
		--genomeDir /scratch/user/s4799891/dnbc4tools_ref \ # path to the reference directory
		--threads 32

mkdir /QRISdata/Q8448/Human_disease_cellatac/${DATA}_dnbc4tools/${DATA}_${SAMPLE}
cp -r ${TMPDIR}/${SAMPLE}/logs /QRISdata/Q8448/Human_disease_cellatac/${DATA}_dnbc4tools/${DATA}_${SAMPLE}/
cp -r ${TMPDIR}/${SAMPLE}/outs /QRISdata/Q8448/Human_disease_cellatac/${DATA}_dnbc4tools/${DATA}_${SAMPLE}/  
