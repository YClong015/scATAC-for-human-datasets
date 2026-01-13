#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --job-name=Consensus_peak_ind
#SBATCH --time=20:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-n

module load anaconda3/2023.09-0
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate pycistopic # conda environment with Python 3.11.13

DATA=Muscle
IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [inds.txt]) # one donor ID per line for the text file

# copy sample peak sets to tmp dir
mkdir ${TMPDIR}/data
cd /QRISdata/Q8448/Human_disease_MACS/${DATA}/sample_peak/
for i in $(ls | grep ${IND} ); do f=$(ls $i | grep "narrowPeak"); cp ${i}/${f} ${TMPDIR}/data/${f}; done

# Generate a list of paths to sample peak sets 
cd $TMPDIR
for f in $(ls data); do echo "${f%%_peaks.narrowPeak}: ${TMPDIR}/data/${f}" >> ${TMPDIR}/narrow_peak_paths.txt; done

mkdir ${TMPDIR}/outs
export TMPDIR=${TMPDIR}
export OUTDIR="${TMPDIR}/outs"
export PEAKPATH="${TMPDIR}/narrow_peak_paths.txt"

python /home/s4799891/scripts/consensus_peak_1st.py

cp ${TMPDIR}/outs/* /QRISdata/Q8448/Human_disease_MACS/${DATA}/individual_peak/${IND}/

conda deactivate
