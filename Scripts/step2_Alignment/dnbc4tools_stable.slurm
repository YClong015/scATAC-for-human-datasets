#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --job-name=dnbc4tools_stable
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-n

DATA=Muscle
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [sampleIDs.txt]) # one sample ID per line in the text file

cd $TMPDIR
## Copy fastq files to scratch
cp /QRISdata/Q8448/Human_disease_data/${DATA}/${SAMPLE}/* ${TMPDIR}/
mv ${TMPDIR}/${SAMPLE}_S1_L001_R1_001.fastq.gz ${TMPDIR}/${SAMPLE}_R1.fastq.gz
mv ${TMPDIR}/${SAMPLE}_S1_L001_R2_001.fastq.gz ${TMPDIR}/${SAMPLE}_R2.fastq.gz

/scratch/user/s4799891/dnbc4tools2.1.3/dnbc4tools atac run \ # path to dnbc4tools2.1.3/dnbc4tools
		--name ${SAMPLE} \
		--fastq1 ${TMPDIR}/${SAMPLE}_R1.fastq.gz \
		--fastq2 ${TMPDIR}/${SAMPLE}_R2.fastq.gz \
		--genomeDir /scratch/user/s4799891/dnbc4tools_ref \ # path to the reference directory generated from the command below
		--threads 32

tar -zcvf $TMPDIR/${SAMPLE}.tar.gz ${TMPDIR}/${SAMPLE}
cp ${TMPDIR}/${SAMPLE}.tar.gz /QRISdata/Q8448/Human_disease_cellatac/${DATA}_dnbc4tools/${DATA}_${SAMPLE}.tar.gz


## Make refeerence for the tool##
# /scratch/user/s4799891/dnbc4tools2.1.3/dnbc4tools atac mkref --fasta [genome.fa] --ingtf [genes.gtf] --species Homo_sapiens --threads 32
