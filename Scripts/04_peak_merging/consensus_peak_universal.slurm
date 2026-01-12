#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --job-name=Consensus_peak_universal
#SBATCH --time=20:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%j.out
#SBATCH -e /home/s4869245/outputs/%x_%j.err

module load anaconda3/2023.09-0
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate pycistopic 

# copy tissue/ disease peak sets to tmp dir
mkdir ${TMPDIR}/data
cd /QRISdata/Q8448/Human_disease_MACS/
for f in $(ls); do cp ${f}/individual_peak/${f}_nonunique_regions.bed ${TMPDIR}/data/${f}_consensus_regions.bed; done

# Generate a list of paths to tissue/disease peak sets 
cd $TMPDIR
for f in $(ls data); do 
sed -i 's/$/\t250/' ${TMPDIR}/data/${f}; 
echo "${f%%_consensus_regions.bed}: ${TMPDIR}/data/${f}" >> ${TMPDIR}/peak_paths.txt; done

mkdir ${TMPDIR}/outs
export TMPDIR=${TMPDIR}
export OUTDIR="${TMPDIR}/outs"
export PEAKPATH="${TMPDIR}/peak_paths.txt"

python /home/s4799891/scripts/consensus_peak_2nd3rdy.py

cp ${TMPDIR}/outs/* /QRISdata/Q8448/Human_disease_MACS/

conda deactivate
