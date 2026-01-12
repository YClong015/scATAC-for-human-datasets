library(Signac)
library(Seurat)
library(GenomicRanges)
library(rtracklayer)
library(rlist)
library(AnnotationHub)
library(future)

plan("multicore", workers = 4)

peak <- rtracklayer::import(Sys.getenv("PEAKPATH"))
s <- paste0(Sys.getenv("SAMPLE"))

# sample metadata
md <- read.csv(
  file = "/QRISdata/Q8448/Human_disease_cellatac/Muscle_dnbc4tools/Muscle_cellMD.csv", # cellMD with the donor info appended manually
  header = TRUE,
  sep = ",",
  row.names = 1,
  na.strings = "-",
  stringsAsFactors = FALSE)[-1,]

# Create Seurat Object
cell.md <- read.table(
  file = paste0("/QRISdata/Q8448/Human_disease_cellatac/Muscle_dnbc4tools/", s, "/output/singlecell.csv"),
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
) [-1, ]

# Subset filtered cells
filtered <- read.delim(paste0("/QRISdata/Q8448/Human_disease_cellatac/Muscle_dnbc4tools/", s, "/output/filter_peak_matrix/barcodes.tsv.gz"),
                       header = FALSE, sep = "\t")

cell.md <- cell.md[filtered$V1,]

fragment <- CreateFragmentObject(
  path = paste0("/QRISdata/Q8448/Human_disease_cellatac/Muscle_dnbc4tools/", s, "/output/fragments.tsv.gz"),
  cells = rownames(cell.md))

counts <- FeatureMatrix(
  fragments = fragment,
  features = peak
)

assay <- CreateChromatinAssay(
  counts = counts,
  fragments = fragment
)

sobject <-  CreateSeuratObject(counts = assay, assay = "ATAC")
sobject$Sample <- s
sobject$Donor <- md[s,"Donor"]
sobject$Age <- md[s,"Age"]
sobject$Muscle <- md[s,"Muscle"]
sobject$CI <- md[s,"Charlson.Index"]
sobject$BI <- md[s,"Barthel Index"]
sobject$Condition <- md[s,"Condition"]

# save.image(paste0(Sys.getenv("OUTDIR"), s, ".RData"))
save(sobject, file = paste0(Sys.getenv("OUTDIR"), s, ".RData"))
