library(Signac)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(future)
plan("multicore", workers = 4)

out.dir <- Sys.getenv("DIR")

## run a for loop over sample IDs and generate a list of objects that can be input to FindIntegrationAnchors.
all.objects <- list()
samples <- readLines("/home/s4799891/dataLists/Muscle_Samples.txt")
for (s in samples) {
  load(paste0("/QRISdata/Q8448/Human_disease_Signac/Muscle/SObjects/universal_SCPM1/doublet_removed/", s, "_dbltFiltered.RData"))
  sglt <- subset(sglt, nCount_ATAC > 1000)
  all.objects[[s]] <- sglt
}

# Merge into one object
if (length(all.objects) == 1) {
  mergedObject <- all.objects[[1]]
} else if (length(all.objects) > 1) {
  mergedObject <- merge(
    x = all.objects[[1]],
    y = all.objects[-1]
  )
}

# load("/QRISdata/Q8448/Human_disease_Signac/Muscle/SObjects/universal_SCPM1/integration/Muscle_singlets_merged.RData")

DefaultAssay(mergedObject) <- "ATAC"
dim(mergedObject)

mergedObject <- RunTFIDF(mergedObject)
mergedObject <- FindTopFeatures(mergedObject, min.cutoff = 'q75')
mergedObject <- RunSVD(object = mergedObject)

all.objects <- SplitObject(mergedObject, split.by = "Donor")

# find integration anchors
features <- list()
for (n in 1:length(all.objects)) {
  all.objects[[n]] <- RunTFIDF(all.objects[[n]])
  all.objects[[n]] <- FindTopFeatures(all.objects[[n]], min.cutoff = 'q75')
  all.objects[[n]] <- RunSVD(object = all.objects[[n]])
  features[[n]] <- VariableFeatures(all.objects[[n]])
}

#anchor.features <- unique(unlist(features))
#length(anchor.features)
anchor.features.table <- table(unlist(features))
sum(anchor.features.table > 1)
anchor.features <- names(anchor.features.table)[which(anchor.features.table > 1)]
length(anchor.features)
head(anchor.features)

# # Union of features
# union.f <- function(features) {
#   if (length(features) == 1) {
#     return(features)
#   } else {
#     return(union(features[[1]], union.f(features[-1])))
#   }
# }
# anchor.features <- union.f(features)

options(future.globals.maxSize = 150000 * 1024^2)
integration.anchors <- FindIntegrationAnchors(
  object.list = all.objects,
  anchor.features = anchor.features, 
  reduction = "rlsi",
  dims = 2:25
)
# integrate LSI embedding
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = mergedObject[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:25
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:25)

integrated <- FindNeighbors(
  object = integrated,
  reduction = "integrated_lsi",
  dims = 2:25           
)
integrated_0.2 <- FindClusters(
  object = integrated,
  algorithm = 3,  # SLM algorithm
  resolution = 0.2,
  verbose = FALSE
)

DimPlot(object = integrated_0.2) + ggtitle("Muscle", subtitle = "Dim: 25; Res: 0.2")
DimPlot(object = integrated_0.2, split.by = "Donor") + ggtitle("Muscle", subtitle = "Dim: 25; Res: 0.2")
# save.image(paste0(Sys.getenv("DIR"), "Muscle_integration.RData"))

# png(paste0(out.dir, "Muscle_UMAP.png"))
# DimPlot(object = integrated_0.2) + ggtitle("Muscle", subtitle = "Dim: 25; Res: 0.2")
# dev.off()

# gene activity
# compute gene activities
gene.activities <- GeneActivity(integrated_0.2)

# add the gene activity matrix to the Seurat object as a new assay
integrated_0.2[['RNA']] <- CreateAssayObject(counts = gene.activities)
integrated_0.2 <- NormalizeData(
  object = integrated_0.2,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated_0.2$nCount_RNA)
)
DefaultAssay(integrated_0.2) <- 'RNA'

# VlnPlot(integrated_0.2, 
#         features = c("AQP7", "CDH5", "SCN7A", "CD3D", "PAX7", "F13A1", 
#                                  "KCNJ8", "MYH11", "SCG2", "TNNT1", "TNNT3"),
#         y.max = 5,
#         pt.size = 0.01,
#         alpha = 0.05,
#         ncol = 3
#         )

# png(paste0("Muscle_GeneActivity.png"))
gene.plots <- FeaturePlot(
  object = integrated_0.2,
  features = c("AQP7", "CDH5", "SCN7A", "CD3D", "PAX7", "F13A1",
               "KCNJ8", "MYH11", "SCG2", "TNNT1", "TNNT3"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  label.size = 4
)
gene.plots <- lapply(X = gene.plots, 
                     FUN = function(x) x + theme(
                       plot.title = element_text(size = 18), 
                       axis.text = element_text(size = 14),
                       axis.title = element_text(size = 14),
                       legend.text = element_text(size = 12),
                       legend.title = element_text(size = 12)))
# dev.off()

cell.types <- as.data.frame(integrated_0.2$ATAC_snn_res.0.2)
colnames(cell.types) <- c("cluster")
cell.types$type <- ifelse(cell.types$cluster == 0, "Type I myonuclei", 
                          ifelse(cell.types$cluster == 1, "Type II myonuclei",
                                 ifelse(cell.types$cluster == 7, "Muscle stem cells",
                                        ifelse(cell.types$cluster == 4, "Endothelial cells",
                                               ifelse(cell.types$cluster == 6, "Smooth muscle cells", "unclassified")))))
integrated_0.2 <- AddMetaData(integrated_0.2, metadata = cell.types$type, col.name = "type")
DimPlot(object = integrated_0.2, group.by = "type") + 
  ggtitle("Muscle") +
  theme(plot.title = element_text(size = 18), 
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10))
ggarrange(plotlist = gene.plots, ncol = 3, nrow = 4)

## Activity score for marker genes in clusters
markers <- c("AQP7", "CDH5", "SCN7A", "CD3D", "PAX7", "F13A1",
             "KCNJ8", "MYH11", "SCG2", "TNNT1", "TNNT3")
Idents(integrated_0.2) <- "seurat_clusters"
FindMarkers(integrated_0.2, features = markers, ident.1 = 2)

# save.image(paste0(out.dir, "Muscle_repSample_integration.RData"))
save(integrated_0.2, file=paste0(out.dir, "Muscle_repSample_integration.RData"))

## Track Visualisation ##
roi = "chr12-6530589-6542299"
png(paste0(Sys.getenv("DIR"), "Universal_Gapdh_Muscle.png"))
CoveragePlot(
  object = integrated_0.2[rownames(integrated_0.2) %in% DAR.peaks,],
  group.by = 'Donor',
  region = roi,
  annotation = TRUE,
  peaks = TRUE
)
dev.off()
