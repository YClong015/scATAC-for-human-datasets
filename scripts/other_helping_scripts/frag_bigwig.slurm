#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=frag_bw
#SBATCH --time=8:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%j.out
#SBATCH -e /home/s4869245/outputs/%x_%j.err

module load anaconda3/2023.09-0
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate scatac_fragment_tools 

scatac_fragment_tools bigwig -i [path to fragments.tsv.gz] -c [path to chrom.sizes] -o [output path]

conda deactivate
