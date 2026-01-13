#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=tabix
#SBATCH --time=5:00:00
#SBATCH --partition=general
#SBATCH -o /home/s4869245/outputs/%x_%j.out
#SBATCH -e /home/s4869245/outputs/%x_%j.err

module load tabixpp

DIR=/QRISdata/Q8448/Human_disease_cellatac/Liver_Precellar

cd $TMPDIR

for f in $(ls ${DIR}); do
mkdir $TMPDIR/outs;
cp ${DIR}/${f}/fragments.tsv.gz $TMPDIR/outs/; 
gzip -d ${TMPDIR}/outs/fragments.tsv.gz;
(grep "^#" ${TMPDIR}/outs/fragments.tsv; grep -v "^#" ${TMPDIR}/outs/fragments.tsv | sort -t "`printf '\t'`" -k 1,1 -k 2,2n) | bgzip > ${TMPDIR}/outs/fragments.tsv.gz;
tabix -p bed ${TMPDIR}/outs/fragments.tsv.gz;
mv ${TMPDIR}/outs ${DIR}/${f}/; done
