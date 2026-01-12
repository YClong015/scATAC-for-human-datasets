library(Signac)
library(Seurat)
library(scDblFinder)
library(GenomicRanges)
library(AnnotationHub)
library(ggplot2)
library(gridExtra)

# Data setup
donor <- Sys.getenv("IND")
sample.id <- Sys.getenv("SAMPLE")
load(paste0("/QRISdata/Q8448/Human_disease_Signac/Muscle/SObjects/universal_SCPM1/", donor, "_filtered.RData"))

# QC
QCobject <- subset(
  x = iObject,
  subset = peak_region_fragments > 100 &
    peak_region_fragments < as.integer(Sys.getenv("FCOUNT")) &
    TSS.enrichment > 1
)

# Get sample object set
objects <- SplitObject(QCobject, split.by = "Sample")
dblt <- objects[[sample.id]]
sce <- SingleCellExperiment(assays = SimpleList(counts = dblt[["ATAC"]]$counts))

# scDblFinder method
set.seed(123)
sce <- scDblFinder(sce, 
                   aggregateFeatures = TRUE, nfeatures = 25, 
                   processing = "normFeatures") 
# artificialDoublets=NULL ->  the maximum of the number of cells or 5*nbClusters^2 (with a minimum of 1500)

# Amulet method
fragfile <- paste0("/QRISdata/Q8448/Human_disease_cellatac/Muscle_dnbc4tools/", sample.id, "/output/fragments.tsv.gz")
toExclude <- GRanges(c("chrX", "chrY", "chrM"), IRanges(1L, width = 10^8))
res <- amulet(fragfile, regionsToExclude = toExclude)

# p-value combination
row.names(res) <- paste0(sample.id, "_", row.names(res))
res$scDblFinder.p <- 1-colData(sce)[rownames(res), "scDblFinder.score"]
res$combined <- apply(res[,c("scDblFinder.p", "p.value")], 1, FUN = function(x){
  x[x<0.001] <- 0.001
  suppressWarnings(aggregation::fisher(x))
})

# Filter out doublets
doublets <- row.names(res[res$combined<0.05,]) # 0.05 as the threshold
res$doublet <- ifelse(res$combined<0.05, "doublet", "singlet")
dblt <- AddMetaData(object = dblt, metadata = res$doublet, col.name = "doublet")
sglt <- subset(x = dblt, subset = doublet == "singlet")

# # Doublet - normalisation, dimensional reduction and clustering
# dblt <- RunTFIDF(dblt)
# dblt <- FindTopFeatures(dblt, min.cutoff = 'q0')
# dblt <- RunSVD(object = dblt)
# DepthCor(dblt) 
# dblt <- RunUMAP(
#   object = dblt,
#   reduction = "lsi",
#   dims = 2:25           
# )
# dblt <- FindNeighbors(
#   object = dblt,
#   reduction = "lsi",
#   dims = 2:25           
# )
# dblt <- FindClusters(
#   object = dblt,
#   algorithm = 3,
#   resolution = 1.2,
#   verbose = FALSE
# )
# 
# png(paste0(sample.id, "_doublet.png"))
# DimPlot(object = dblt, label = FALSE, group.by = "doublet") + ggtitle(sample.id, subtitle = "Dim: 25")
# dev.off()

# save.image(paste0(Sys.getenv("DIR"), sample.id, "_dbltFiltered.RData"))
save(sglt, file = paste0(Sys.getenv("DIR"), sample.id, "_dbltFiltered.RData"))

