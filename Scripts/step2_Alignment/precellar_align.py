import os
import precellar

output_file = os.environ['OUTDIR'] + "/fragments.tsv.gz"

assay = precellar.Assay("/scratch/user/s4799891/generic_atac.yaml")

assay.update_read("R1", fastq="R1_processed.fq.gz")
assay.update_read("I1", fastq="I1.fq.gz")
assay.update_read("R2", fastq="R2_processed.fq.gz")

bwa = precellar.aligners.BWAMEM2("/scratch/user/s4799891/bwa_index/GRCh38")
qc = precellar.align(
    assay,
    aligner=bwa,
    modality="atac",
    output=output_file,
    output_type='fragment',
    num_threads=32
)

qc_file = os.environ['OUTDIR'] + "/qc.txt"
with open(qc_file, 'w') as f:
    print(qc, file=f)
