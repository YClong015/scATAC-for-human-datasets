library(Signac)
library(Seurat)
library(GenomicRanges)
library(rtracklayer)
library(rlist)
library(AnnotationHub)
library(future)

plan("multicore", workers = 4)
i <- Sys.getenv("IND")

# object list
objects <- list()
names <- list()

samples <- readLines("[Samples.txt]")
for (s in samples) {
  if (grepl(pattern = paste0(i, "_"), s)) {
  load(paste0(Sys.getenv("DIR"), s, ".RData"))

  objects[[s]] <- sobject
  names[[s]] <- sobject$Sample[[1]]
  }
}

# Merge into patient level objects
if (length(objects) == 1) {
  iObject <- objects[[1]]
  RenameCells(iObject, new.names = paste0(names[[1]], "_", Cells(iObject)))
} else if (length(objects) > 1) {
  iObject <- merge(
    x = objects[[1]],
    y = objects[-1],
    add.cell.ids = names
  )
}
iObject$Donor <- i

# save.image(paste0(Sys.getenv("DIR"), i, "_filtered.RData"))
save(iobject, file = paste0(Sys.getenv("OUTDIR"), i, "_filtered.RData"))

