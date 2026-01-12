#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --job-name=peak_uniq
#SBATCH --time=8:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%A_%a.out
#SBATCH -e /home/s4869245/outputs/%x_%A_%a.err
#SBATCH --array=1-9

module load bedtools

DATA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" [DataList.txt]) # one data name per line for the text file

cd /QRISdata/Q8448/Human_disease_MACS/${DATA}/individual_peak

# Merge all peaks from individual peak sets
for IND in $(cat /home/s4799891/${DATA}_inds.txt); do
cat ${IND}/consensus_regions.bed >> all.bed;
done

# Count overlapped peaks
bedtools intersect -a consensus_regions.bed\
                   -b all.bed\
                   -f 0.5 -c > consensus_regions_cnts.bed

# extract unique/ non-unique peaks
awk '$NF == 1' consensus_regions_cnts.bed | cut -f 1-3 > ${DATA}_unique_regions.txt
awk '$NF > 1' consensus_regions_cnts.bed | cut -f 1-3 > ${DATA}_nonunique_regions.txt

grep -wf ${DATA}_unique_regions.txt consensus_regions.bed > ${DATA}_unique_regions.bed
grep -wf ${DATA}_nonunique_regions.txt consensus_regions.bed > ${DATA}_nonunique_regions.bed

# Summary
TOTAL=$(wc -l consensus_regions.bed)
UNIQUE=$(wc -l ${DATA}_unique_regions.bed)
NONUNIQUE=$(wc -l ${DATA}_nonunique_regions.bed)
echo "Total number of peaks in ${DATA}: ${TOTAL}"
echo "Number of unique peaks in ${DATA}: ${UNIQUE}"
echo "Number of filtered peaks in ${DATA}: ${NONUNIQUE}"
