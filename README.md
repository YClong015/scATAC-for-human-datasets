# scATAC-seq Workflow (SLURM)

This repository contains a SLURM-based scATAC-seq workflow covering the entire analysis pipeline:
**Data Download** → **Alignment/Preprocessing** → **Peak Calling** → **Peak Merging** → **Seurat/Signac Integration & Annotation** → **DAR & Motif Enrichment**.

It also includes helper utilities for file renaming, bigWig conversion, and tabix indexing.

---

## Requirements

* **System**: Linux environment with SLURM scheduler
* **Package Management**: All dependencies are defined in `environment.yml`.
    * Install via: `conda env create -f environment.yml`

* **Language Support**:
    * R (>= 4.3)
    * Python (>= 3.9)
    * Perl (for HOMER)

* **Bioinformatics CLI Tools** (Installed via Conda):
    * **Alignment**: `BWA`, `Samtools`
    * **Data Download**: `SRA Toolkit`, `synapseclient`
    * **Peak Calling**: `MACS2`
    * **Motif Analysis**: `HOMER`
    * **Utilities**: `Bedtools`, `htslib` (tabix), `bedGraphToBigWig`, `pigz`

* **Python Libraries**:
    * `pandas`, `numpy`, `pysam`

* **R Packages**:
    * **Core**: `Seurat`, `Signac`, `dplyr`, `ggplot2`, `patchwork`, `Matrix`
    * **Genomics**: `GenomicRanges`, `EnsDb.Hsapiens.v86`
    * **Integration**: `harmony`
    * **QC & Filtering**: `scDblFinder`
    * **Motif Utilities**: `BiocGenerics`
---

## Inputs

To run this pipeline, you need:

1.  **Sample Sheet**: A CSV/TSV file containing at least the following columns:
    * `sample_id`, `donor_id`, `tissue`, `condition`, `data_source`, `fastq_1`, `fastq_2`
2.  **Reference Files**:
    * Genome FASTA
    * Annotation GTF
    * Chromosome sizes
    * Blacklist regions (optional)

---

## How to Run

Each step is designed to be submitted as a SLURM job.
**Note**: Please edit the paths and parameters inside the `.slurm` files (and related `.py/.R` scripts) to match your environment before running.

> **Directory Note**: The commands below assume scripts are located in the repository root. If your scripts are in a `slurm/` folder, update paths accordingly (e.g., `sbatch slurm/dataDownload_SRA.slurm`).

### Step 1. Data Download

#### SRA (Brain_N, Brain_SA, Eye, Heart, Kidney_NC, Liver, Lung)
```bash
sbatch dataDownload_SRA.slurm
```

#### CNGB (Muscle)
```bash
sbatch dataDownload_CNGB.slurm
```

#### Synapse (Brain_CB)
```bash
synapse login
synapse get-download-list --downloadLocation <DOWNLOAD_DIR>
```

---

### Step 2. Alignment / Preprocessing

#### Standard Tissues (Brain_N, Brain_SA, Brain_CB, Heart, Kidney_NC, Lung)
```bash
sbatch cellatac.slurm
```

#### Eye
```bash
sbatch cellarc.slurm
```

#### Muscle
```bash
sbatch dnbc4tools_stable.slurm
# OR
sbatch dnbc4tools_beta.slurm
```

#### Liver (Precellar Pipeline)
This tissue requires a specific sequence of scripts:
1.  `modify_header.slurm`
2.  `precellar_stripBarcode.py` (called via wrapper if applicable)
3.  `precellar_barcode.slurm`
4.  `precellar_align.py` (called via wrapper if applicable)
5.  `precellar_align.slurm`

**Run order:**
```bash
sbatch modify_header.slurm
sbatch precellar_barcode.slurm
sbatch precellar_align.slurm
```

---

### Step 3. Peak Calling

Uses MACS2 via Python wrapper.

```bash
sbatch MACS.slurm
```

---

### Step 4. Peak Merging

The merging strategy follows this hierarchy:
`Sample` → `Donor` → `Tissue/Disease` → `Universal` (with unique peak filtering).

#### 1. Sample → Donor peak set
```bash
sbatch consensus_peak_ind.slurm
```

#### 2. Donor → Tissue/Disease peak set
```bash
sbatch consensus_peak_tissue.slurm
```

#### 3. Unique Peaks Filtering
```bash
sbatch peak_uniq.slurm
```

#### 4. Tissue/Disease → Universal peak set
```bash
sbatch consensus_peak_universal.slurm
```

---

### Step 5. Integration Analysis (Seurat/Signac)

#### Create Seurat Objects (Per-sample)
```bash
sbatch SeuratObject.slurm
```

#### Merge into Donor Objects
```bash
sbatch DonorObject.slurm
```

#### Doublet Filtering
```bash
sbatch doublet_filtering.slurm
```

#### Integration + Cell Type Annotation
```bash
sbatch integration.slurm
```

#### DAR Calling
If SLURM scripts are available, submit via `sbatch`. Otherwise, run directly in R:
```bash
Rscript cell_types.R
Rscript DAR.R
```

#### Motif Enrichment
```bash
sbatch HOMER.slurm
Rscript motif_visualisation.R
```

---

### Step 6. Helper Scripts

#### Rename Downloaded Sample Names
```bash
sbatch renaming.slurm
```

#### Convert fragments.tsv.gz to bigWig
```bash
sbatch frag_bigwig.slurm
```

#### Generate Tabix Index for Fragments
```bash
sbatch tabix.slurm
```

---

## Notes

* **Version Control**: Do **NOT** commit large files (FASTQ, BAM, fragments, bigWig) to git. Ensure they are listed in your `.gitignore`.
* **Configuration**: It is recommended to keep paths and parameters in a separate `config/` folder or file to keep scripts clean.
* **Logging**: Suggested SLURM log directory structure: `logs/<jobname>_%j.out`.
