# scATAC-seq Workflow (SLURM)

This repository contains a SLURM-based scATAC-seq workflow:
data download → alignment/preprocessing → peak calling → peak merging →
Seurat/Signac integration & annotation → DAR & motif enrichment,
plus helper utilities (renaming, bigWig, tabix).

## Requirements
- Linux + SLURM
- Tools: samtools, bedtools, htslib/tabix, MACS2, (bigWig tools e.g. bedGraphToBigWig or deepTools)
- R: Seurat, Signac (+ other packages used by the R scripts)
- HOMER (for motif enrichment)
- Optional: synapseclient (Synapse CLI) for Brain_CB

## Inputs
You need:
1) A sample sheet containing at least:
   - sample_id, donor_id, tissue, condition, data_source, fastq_1, fastq_2
2) Reference files:
   - genome fasta
   - annotation gtf
   - chromosome sizes / blacklist

## How to run (overview)
Each step is submitted as a SLURM job.
files before running.

### Step 1. Data download
**SRA** (Brain_N, Brain_SA, Eye, Heart, Kidney_NC, Liver, Lung)
- `dataDownload_SRA.slurm`
```bash
sbatch dataDownload_SRA.slurm
