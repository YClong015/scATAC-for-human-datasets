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

To run this pipeline, you need to prepare two main components:

### 1. Sample Sheet (`sample_sheet.csv`)
A CSV file containing metadata and file paths for each sample. 
> **Tip**: An example is provided in `examples/sample_sheet.csv`.

**Required Columns:**
* `sample_id`: Unique identifier for the sample.
* `donor_id`: Identifier for the biological donor (for merging/integration).
* `tissue`: Tissue type (e.g., Lung, Kidney).
* `condition`: Experimental condition (e.g., Control, Treated).
* `data_source`: Origin of data (e.g., SRA, In-house).
* `fastq_1`: Absolute path to Read 1 file (`.fastq.gz`).
* `fastq_2`: Absolute path to Read 2 file (`.fastq.gz`).

### 2. Reference Files
* **Genome FASTA**: (`.fa` or `.fasta`) required for alignment (BWA).
* **Annotation GTF**: (`.gtf`) required for peak annotation and TSS enrichment QC.
* **Chromosome Sizes**: (`chrom.sizes`) required for bigWig conversion.
* **Blacklist Regions**: (`.bed`) optional but recommended for filtering artifact peaks (e.g., ENCODE blacklist).

---

## How to Run

Each step is designed to be submitted as a SLURM job.
**Note**: Please edit the paths and parameters inside the `.slurm` files (and related `.py/.R` scripts) to match your environment before running.

> **Directory Note**: The commands below assume scripts are located in the repository root. If your scripts are in a `slurm/` folder, update paths accordingly (e.g., `sbatch slurm/dataDownload_SRA.slurm`).

### Step 1. Data Download

#### SRA (Brain_N, Brain_SA, Eye, Heart, Kidney_NC, Liver, Lung)
```bash
sbatch scripts/01_download/dataDownload_SRA.slurm
```

#### CNGB (Muscle)
```bash
sbatch scripts/01_download/dataDownload_CNGB.slurm
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
sbatch scripts/02_alignment/cellatac.slurm
```

#### Eye
```bash
sbatch scripts/02_alignment/cellarc.slurm
```

#### Muscle
```bash
sbatch scripts/02_alignment/dnbc4tools_stable.slurm
# OR
sbatch scripts/02_alignment/dnbc4tools_beta.slurm
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
sbatch scripts/02_alignment/modify_header.slurm
sbatch scripts/02_alignment/precellar_barcode.slurm
sbatch scripts/02_alignment/precellar_align.slurm
```

---

### Step 3. Peak Calling

Uses MACS2 via Python wrapper.

```bash
sbatch scripts/03_peak_calling/MACS.slurm
```

---

### Step 4. Peak Merging

The merging strategy follows this hierarchy:
`Sample` → `Donor` → `Tissue/Disease` → `Universal` (with unique peak filtering).

#### 1. Sample → Donor peak set
```bash
sbatch scripts/04_peak_merging/consensus_peak_ind.slurm
```

#### 2. Donor → Tissue/Disease peak set
```bash
sbatch scripts/04_peak_merging/consensus_peak_tissue.slurm
```

#### 3. Unique Peaks Filtering
```bash
sbatch scripts/04_peak_merging/peak_uniq.slurm
```

#### 4. Tissue/Disease → Universal peak set
```bash
sbatch scripts/04_peak_merging/consensus_peak_universal.slurm
```

---

### Step 5. Integration Analysis (Seurat/Signac)

#### Create Seurat Objects (Per-sample)
```bash
sbatch scripts/05_integration_analysis/SeuratObject.slurm
```

#### Merge into Donor Objects
```bash
sbatch scripts/05_integration_analysis/DonorObject.slurm
```

#### Doublet Filtering
```bash
sbatch scripts/05_integration_analysis/doublet_filtering.slurm
```

#### Integration + Cell Type Annotation
```bash
sbatch scripts/05_integration_analysis/integration.slurm
```

#### DAR Calling
If SLURM scripts are available, submit via `sbatch`. Otherwise, run directly in R:
```bash
Rscript scripts/05_integration_analysis/cell_types.R
Rscript scripts/05_integration_analysis/DAR.R
```

#### Motif Enrichment
```bash
sbatch scripts/05_integration_analysis/HOMER.slurm
Rscript scripts/05_integration_analysis/motif_visualisation.R
```

---

### Step 6. Helper Scripts

#### Rename Downloaded Sample Names
```bash
sbatch scripts/other_helping_scripts/renaming.slurm
```

#### Convert fragments.tsv.gz to bigWig
```bash
sbatch scripts/other_helping_scripts/frag_bigwig.slurm
```

#### Generate Tabix Index for Fragments
```bash
sbatch scripts/other_helping_scripts/tabix.slurm
```

---

## Notes

* ** Data Management**:
    Large intermediate files (FASTQ, BAM, fragments, bigWig) are strictly excluded. **Do NOT** force commit these files.

* **Ep Logging**:
    SLURM logs are directed to the `logs/` directory. The standard naming convention used in `.slurm` scripts is `logs/<job_name>_%j.out` (where `%j` is the job ID) to facilitate debugging.
