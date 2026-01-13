#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --job-name=MACS
#SBATCH --time=20:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-n

module load anaconda3/2023.09-0
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate pycistopic # conda environment with Python 3.11.13

DATA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [DataList.txt]) # one data name per line for the text file

mkdir ${TMPDIR}/data

# copy fragment files to tmp dir
cd /QRISdata/Q8448/Human_disease_cellatac/${DATA}
for f in $(ls /QRISdata/Q8448/Human_disease_cellatac/${DATA} | grep -v ".tar.gz"); do cp ${f}/outs/fragments.tsv.gz ${TMPDIR}/data/${f}_fragments.tsv.gz; done

# Generate a list of paths to samples' fragment files
cd $TMPDIR
for f in $(ls data); do echo "${f%%_fragments.tsv.gz}: ${TMPDIR}/data/${f}" >> ${TMPDIR}/${DATA}_bedPaths.txt; done

mkdir ${TMPDIR}/outs
export TMPDIR=${TMPDIR}
export OUTDIR="${TMPDIR}/outs"
export BEDPATH="${TMPDIR}/${DATA}_bedPaths.txt"

python /home/s4799891/scripts/MACS.py

cp ${TMPDIR}/outs/* /QRISdata/Q8448/Human_disease_MACS/${DATA}/sample_peak

conda deactivate
