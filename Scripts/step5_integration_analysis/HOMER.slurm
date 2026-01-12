#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=HOMER
#SBATCH --time=5:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%j.out
#SBATCH -e /home/s4869245/outputs/%x_%j.err

PATH=$PATH:/scratch/user/s4799891/homer/.//bin/

findMotifsGenome.pl /QRISdata/Q8448/Human_disease_Signac/Muscle/SObjects/universal_SCPM1/DAR/Type1_myonuclei_AgedNDvYoungND_005pct_closing.bed hg38 motif_res/Type1_myonuclei_AgedNDvYoungND_close -p 8 -nomotif
findMotifsGenome.pl /QRISdata/Q8448/Human_disease_Signac/Muscle/SObjects/universal_SCPM1/DAR/Type1_myonuclei_AgedNDvYoungND_005pct_closing.bed hg38 motif_res/Type1_myonuclei_AgedNDvYoungND_close -p 8 -nomotif
findMotifsGenome.pl /QRISdata/Q8448/Human_disease_Signac/Muscle/SObjects/universal_SCPM1/DAR/Type2_myonuclei_AgedNDvYoungND_005pct_opening.bed hg38 motif_res/Type2_myonuclei_AgedNDvYoungND_open -p 8 -nomotif
findMotifsGenome.pl /QRISdata/Q8448/Human_disease_Signac/Muscle/SObjects/universal_SCPM1/DAR/Type2_myonuclei_AgedNDvYoungND_005pct_closing.bed hg38 motif_res/Type2_myonuclei_AgedNDvYoungND_close -p 8 -nomotif
