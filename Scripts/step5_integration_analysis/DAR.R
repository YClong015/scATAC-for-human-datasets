library(dplyr)
library(ggplot2)
library(ggpubr)

### Function ###
apply_DESeq2_test_seurat <- function(seurat.object, 
                                     population.1 = NULL, 
                                     population.2 = NULL, 
                                     exp.thresh = 0.075,
                                     fc.thresh=0.1, 
                                     adj.pval.thresh = 0.05, 
                                     num.splits = 6, 
                                     seed.use = 1,
                                     replicates.1 = NULL,
                                     replicates.2 = NULL,
                                     verbose = TRUE, 
                                     return.deseq2.res = FALSE, 
                                     assay.use = "ATAC",
                                     ncores = 1) {
  
  if (!'DESeq2' %in% rownames(x = installed.packages())) {
    stop("Please install DESeq2 before using this function
         (http://bioconductor.org/packages/release/bioc/html/DESeq2.html)")
  }
  
  if (is.null(population.1) & is.null(replicates.1)) {
    stop("Both population.1 and replicates.1 cannot be NULL")
  }
  
  if (!is.null(population.1) & !is.null(replicates.1)) {
    stop("Values for either population.1 OR replicates.1 should be provided - not both")
  }
  
  ## reduce counts in a cluster to num.splits cells for genes with > 1 peak
  if (is.null(replicates.1)) {
    high.expressed.peaks <- GetExpressedPeaks(seurat.object, population.1, population.2, threshold = exp.thresh, assay.use = assay.use)
    length(high.expressed.peaks)
  } else {
    high.expressed.peaks <- GetExpressedPeaks(seurat.object, unlist(replicates.1), unlist(replicates.2), threshold = exp.thresh, assay.use = assay.use)
  }

  if (verbose) print(paste(length(high.expressed.peaks), "expressed peaks passing threshold ", toString(exp.thresh)))
  
  peaks.use <- high.expressed.peaks
  if (verbose) print(paste(length(peaks.use), "individual peak sites to test"))
  
  ## make pseudo-bulk profiles out of cells
  ## set a seed to allow replication of results
  set.seed(seed.use)
  if (is.null(replicates.1)) {
    
    if (length(population.1) == 1) {
      cells.1 <- names(Seurat::Idents(seurat.object))[which(Seurat::Idents(seurat.object) == population.1)]
    } else{
      cells.1 <- population.1
    }
    
    cells.1 = sample(cells.1)
    cell.sets1 <- split(cells.1, sort(1:length(cells.1)%%num.splits))
  } else{
    ## user has provided cells for replicates - use these instead
    cell.sets1 <- replicates.1
  }
  
  ## create a profile set for first cluster
  profile.set1 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets1))
  for (i in 1:length(cell.sets1)) {
    this.set <- cell.sets1[[i]]
    sub.matrix <- Seurat::GetAssayData(seurat.object, layer = "counts", assay = assay.use)[peaks.use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x) sum(x)))
      profile.set1[, i] <- this.profile
    } else {
      profile.set1[, i] <- sub.matrix
    }
  }
  rownames(profile.set1) <- peaks.use
  colnames(profile.set1) <- paste0("Population1_", 1:length(cell.sets1))
  
  ## create a profile set for second cluster
  if (is.null(replicates.2)) {
    if (is.null(population.2)) {
      cells.2 <- setdiff(colnames(seurat.object), cells.1)
    } else {
      if (length(population.2) == 1) {
        cells.2 <- names(Seurat::Idents(seurat.object))[which(Seurat::Idents(seurat.object) == population.2)]
      } else {
        cells.2 <- population.2
      }
    }
    
    cells.2 = sample(cells.2)
    cell.sets2 <- split(cells.2, sort(1:length(cells.2)%%num.splits))
  } else{
    ## user has provided cells for replicates - use these instead
    cell.sets2 <- replicates.2
  }
  
  
  profile.set2 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets2))
  for (i in 1:length(cell.sets2)) {
    this.set <- cell.sets2[[i]]
    sub.matrix <- Seurat::GetAssayData(seurat.object, layer = "counts", assay = assay.use)[peaks.use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x) sum(x)))
      profile.set2[, i] <- this.profile
    } else {
      profile.set2[, i] <- sub.matrix
    }
  }
  rownames(profile.set2) <- peaks.use
  colnames(profile.set2) <- paste0("Population2_", 1:length(cell.sets2))
  
  ## merge the count matrices together
  peak.matrix <- cbind(profile.set1, profile.set2)
  
  ## Create the DEXSeq sample table
  sampleTable <- data.frame(row.names = c(colnames(profile.set1), colnames(profile.set2)),
                            condition = c(rep("target", ncol(profile.set1)),
                                          rep("comparison", ncol(profile.set2))))
  
  #dexseq.feature.table <- Tool(apa.seurat.object, "Sierra")[, c("Gene_name", "Gene_part", "Peak_number")]
  #dexseq.feature.table$Peak <- rownames(dexseq.feature.table)
  #dexseq.feature.table <- dexseq.feature.table[rownames(peak.matrix), ]
  
  # removing potential colons and spaces from gene names to match output of DEXSeq                                     
  #pid_gene_names <- gsub('[: ]', '', dexseq.feature.table$Gene_name)
  #rownames(dexseq.feature.table) <- paste0(pid_gene_names, ":", dexseq.feature.table$Peak_number) 
  #rownames(peak.matrix) <- rownames(dexseq.feature.table)
  
  #peak_ID_set = dexseq.feature.table[rownames(peak.matrix), "Peak_number"]
  #gene_names = dexseq.feature.table[rownames(peak.matrix), "Gene_name"]
  
  ## Build the DESeq2 object
  #dxd = DEXSeq::DEXSeqDataSet(peak.matrix, sampleData=sampleTable, groupID = gene_names,
  #                            featureID = peak_ID_set, design= ~sample+exon+condition:exon)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=peak.matrix, 
                                        colData=sampleTable, 
                                        design=~condition)
  
  ## Run DESeq2
  if (verbose) print("Running DESeq2 test...")
  
  dds <- DESeq2::DESeq(dds)
  res.specific <- DESeq2::results(dds, contrast=c( "condition", "target", "comparison" ), test="Wald")
  
  if (return.deseq2.res) {
    return(dds)
  } else {
    return(res.specific)
  } 

}



GetExpressedPeaks <- function(seurat.object, population.1, population.2=NULL, threshold=0.05, assay.use = "ATAC") {
  
  if (length(population.1) == 1){ # cluster identity used as input
    foreground.set = names(Seurat::Idents(seurat.object)[Seurat::Idents(seurat.object)==population.1])
  } else { # cell identity used as input
    foreground.set = population.1
  }
  if (is.null(population.2)) {
    remainder.set = names(Seurat::Idents(seurat.object)[Seurat::Idents(seurat.object)!=population.1])
  } else {
    if (length(population.2) == 1) { # cluster identity used as input
      remainder.set = names(Seurat::Idents(seurat.object)[Seurat::Idents(seurat.object)==population.2])
    } else { # cell identity used as input
      remainder.set = population.2
    }
  }
  
  peak.names = rownames(seurat.object)
  
  # Get the peaks expressed in the foreground set based on proportion of non-zeros
  this.data <- Seurat::GetAssayData(seurat.object, layer = "counts", assay=assay.use)
  nz.row.foreground = tabulate(this.data[, foreground.set]@i + 1, nbins = nrow(seurat.object))
  nz.prop.foreground = nz.row.foreground/length(foreground.set)
  peaks.foreground = peak.names[which(nz.prop.foreground > threshold)]
  
  # Now identify the peaks expressed in the background set
  nz.row.background = tabulate(this.data[, remainder.set]@i + 1, nbins = nrow(seurat.object))
  nz.prop.background = nz.row.background/length(remainder.set)
  peaks.background = peak.names[which(nz.prop.background > threshold)]
  
  return(union(peaks.foreground, peaks.background))
}

###
out.dir <- "/QRISdata/Q8448/Human_disease_Signac/Muscle/SObjects/universal_SCPM1/DAR/"
DefaultAssay(integrated_0.2) <- 'ATAC'
###
# DAR Calling (Type I myonuclei - age matched diseased vs control)
output.label <- "Type1_myonuclei"
res.table <- apply_DESeq2_test_seurat( integrated_0.2, population.1 = type1.diseased.cells, population.2 = type1.control.cells,
                                       exp.thresh = 0.05, num.splits = 10 )
hist(res.table$pvalue)
res.table <- as.data.frame(res.table)
res.table.sig <- subset(res.table, padj < 0.05)
nrow(res.table.sig)
sum(res.table.sig$log2FoldChange > 0)
sum(res.table.sig$log2FoldChange < 0)

## Write out all the tested peaks
write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
            sep = "\t", quote = FALSE)
## Write out the significant DARs
write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
            sep = "\t", quote = FALSE)

### Write out opening DARs
peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
length(peaks.up)
peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
                         Start = sub(".*-(.*)-.*", "\\1", peaks.up),
                         End = sub(".*-.*-(.*)", "\\1", peaks.up))
write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
### Write out closing DARs
peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
length(peaks.down)
peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
                         Start = sub(".*-(.*)-.*", "\\1", peaks.down),
                         End = sub(".*-.*-(.*)", "\\1", peaks.down))
write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

###
# DAR Calling (Type II myonuclei - age matched diseased vs control)
output.label <- "Type2_myonuclei"
res.table <- apply_DESeq2_test_seurat( integrated_0.2, population.1 = type2.diseased.cells, population.2 = type2.control.cells,
                                       exp.thresh = 0.05, num.splits = 10 )
hist(res.table$pvalue)
res.table <- as.data.frame(res.table)
res.table.sig <- subset(res.table, padj < 0.05)
nrow(res.table.sig)
sum(res.table.sig$log2FoldChange > 0)
sum(res.table.sig$log2FoldChange < 0)

## Write out all the tested peaks
write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
            sep = "\t", quote = FALSE)
## Write out the significant DARs
write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
            sep = "\t", quote = FALSE)

### Write out opening DARs
peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
length(peaks.up)
peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
                         Start = sub(".*-(.*)-.*", "\\1", peaks.up),
                         End = sub(".*-.*-(.*)", "\\1", peaks.up))
write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
### Write out closing DARs
peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
length(peaks.down)
peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
                         Start = sub(".*-(.*)-.*", "\\1", peaks.down),
                         End = sub(".*-.*-(.*)", "\\1", peaks.down))
write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# # DAR Calling (MuSC)
# output.label <- "MuSC"
# res.table <- apply_DESeq2_test_seurat( integrated_0.2, population.1 = MuSC.diseased.cells, population.2 = MuSC.control.cells,
#                                        exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# # DAR Calling (EC)
# output.label <- "EC"
# res.table <- apply_DESeq2_test_seurat( integrated_0.2, population.1 = EC.diseased.cells, population.2 = EC.control.cells,
#                                        exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# # DAR Calling (SMC)
# output.label <- "SMC"
# res.table <- apply_DESeq2_test_seurat( integrated_0.2, population.1 = SMC.diseased.cells, population.2 = SMC.control.cells,
#                                        exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# # DAR Calling (Type I myonuclei - aged vs young)
# output.label <- "Type1_myonuclei_AvY"
# res.table <- apply_DESeq2_test_seurat( integrated_0.2, population.1 = type1.aged.cells, population.2 = type1.young.cells,
#                                        exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# # DAR Calling (Type II myonuclei - aged vs young)
# output.label <- "Type2_myonuclei_AvY"
# res.table <- apply_DESeq2_test_seurat( integrated_0.2, population.1 = type2.aged.cells, population.2 = type2.young.cells,
#                                        exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# # DAR Calling (MuSC - aged vs young)
# output.label <- "MuSC_AvY"
# res.table <- apply_DESeq2_test_seurat( integrated_0.2, population.1 = MuSC.aged.cells, population.2 = MuSC.young.cells,
#                                        exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# # DAR Calling (EC - aged vs young)
# output.label <- "EC_AvY"
# res.table <- apply_DESeq2_test_seurat( integrated_0.2, population.1 = EC.aged.cells, population.2 = EC.young.cells,
#                                        exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# # DAR Calling (SMC - aged vs young)
# output.label <- "SMC_AvY"
# res.table <- apply_DESeq2_test_seurat( integrated_0.2, population.1 = SMC.aged.cells, population.2 = SMC.young.cells,
#                                        exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# # DAR Calling (Type I myonuclei - Aged healthy vs Young healthy)
# output.label <- "Type1_myonuclei_AgedHvYoungH"
# res.table <- apply_DESeq2_test_seurat(integrated_0.2, population.1 = type1.control.cells, population.2 = type1.young.healthy.cells,
#                                        exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# # DAR Calling (Type II myonuclei - Aged healthy vs Young healthy)
# output.label <- "Type2_myonuclei_AgedHvYoungH"
# res.table <- apply_DESeq2_test_seurat(integrated_0.2, population.1 = type2.control.cells, population.2 = type2.young.healthy.cells,
#                                       exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# # DAR Calling (Type I myonuclei - Aged non-diseased vs Young non-diseased)
# output.label <- "Type1_myonuclei_AgedNDvYoungND"
# res.table <- apply_DESeq2_test_seurat(integrated_0.2, population.1 = type1.aged.nondiseased.cells, population.2 = type1.young.cells,
#                                       exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# # DAR Calling (Type II myonuclei - Aged non-diseased vs Young non-diseased)
# output.label <- "Type2_myonuclei_AgedNDvYoungND"
# res.table <- apply_DESeq2_test_seurat(integrated_0.2, population.1 = type2.aged.nondiseased.cells, population.2 = type2.young.cells,
#                                       exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# # DAR Calling (Type I myonuclei - Aged healthy vs Young)
# output.label <- "Type1_myonuclei_AgedHvYoung"
# res.table <- apply_DESeq2_test_seurat(integrated_0.2, population.1 = type1.control.cells, population.2 = type1.young.cells,
#                                       exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# # DAR Calling (Type II myonuclei - Aged healthy vs Young)
# output.label <- "Type2_myonuclei_AgedHvYoung"
# res.table <- apply_DESeq2_test_seurat(integrated_0.2, population.1 = type2.control.cells, population.2 = type2.young.cells,
#                                       exp.thresh = 0.05, num.splits = 10 )
# hist(res.table$pvalue)
# res.table <- as.data.frame(res.table)
# res.table.sig <- subset(res.table, padj < 0.05)
# nrow(res.table.sig)
# sum(res.table.sig$log2FoldChange > 0)
# sum(res.table.sig$log2FoldChange < 0)
# 
# ## Write out all the tested peaks
# write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
#             sep = "\t", quote = FALSE)
# ## Write out the significant DARs
# write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
#             sep = "\t", quote = FALSE)
# 
# ### Write out opening DARs
# peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
# length(peaks.up)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.up),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.up))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# ### Write out closing DARs
# peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
# length(peaks.down)
# peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
#                          Start = sub(".*-(.*)-.*", "\\1", peaks.down),
#                          End = sub(".*-.*-(.*)", "\\1", peaks.down))
# write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)

# DAR Calling (Type I myonuclei - Aged unknown vs Young unknown)
output.label <- "Type1_myonuclei_AgedNoBIvYoungNoBI"
res.table <- apply_DESeq2_test_seurat(integrated_0.2, population.1 = type1.aged.unknown.cells, population.2 = type1.young.unknown.cells,
                                      exp.thresh = 0.05, num.splits = 10 )
hist(res.table$pvalue)
res.table <- as.data.frame(res.table)
res.table.sig <- subset(res.table, padj < 0.05)
nrow(res.table.sig)
sum(res.table.sig$log2FoldChange > 0)
sum(res.table.sig$log2FoldChange < 0)

## Write out all the tested peaks
write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
            sep = "\t", quote = FALSE)
## Write out the significant DARs
write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
            sep = "\t", quote = FALSE)

### Write out opening DARs
peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
length(peaks.up)
peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
                         Start = sub(".*-(.*)-.*", "\\1", peaks.up),
                         End = sub(".*-.*-(.*)", "\\1", peaks.up))
write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
### Write out closing DARs
peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
length(peaks.down)
peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
                         Start = sub(".*-(.*)-.*", "\\1", peaks.down),
                         End = sub(".*-.*-(.*)", "\\1", peaks.down))
write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# DAR Calling (Type II myonuclei - Aged unknown vs Young unknown)
output.label <- "Type2_myonuclei_AgedNoBIvYoungNoBI"
res.table <- apply_DESeq2_test_seurat(integrated_0.2, population.1 = type2.aged.unknown.cells, population.2 = type2.young.unknown.cells,
                                      exp.thresh = 0.05, num.splits = 10 )
hist(res.table$pvalue)
res.table <- as.data.frame(res.table)
res.table.sig <- subset(res.table, padj < 0.05)
nrow(res.table.sig)
sum(res.table.sig$log2FoldChange > 0)
sum(res.table.sig$log2FoldChange < 0)

## Write out all the tested peaks
write.table(res.table, file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
            sep = "\t", quote = FALSE)
## Write out the significant DARs
write.table(res.table.sig, file = paste0(out.dir, output.label, "_005pct_DESeq2_padj005.csv"),
            sep = "\t", quote = FALSE)

### Write out opening DARs
peaks.up <- rownames(subset(res.table.sig, log2FoldChange > 0))
length(peaks.up)
peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.up),
                         Start = sub(".*-(.*)-.*", "\\1", peaks.up),
                         End = sub(".*-.*-(.*)", "\\1", peaks.up))
write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_opening.bed"), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
### Write out closing DARs
peaks.down <- rownames(subset(res.table.sig, log2FoldChange < 0))
length(peaks.down)
peak.table <- data.frame(Chr = sub("(.*)-.*-.*", "\\1", peaks.down),
                         Start = sub(".*-(.*)-.*", "\\1", peaks.down),
                         End = sub(".*-.*-(.*)", "\\1", peaks.down))
write.table(peak.table, file = paste0(out.dir, output.label, "_005pct_closing.bed"), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


## Correlation analysis ##
### Read in the diff comparison
diff.label <- "Type2_myonuclei_AgedNoBIvYoungNoBI"
diff.table <- read.table(paste0(out.dir, diff.label, "_005pct_DESeq2_all.csv"),
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
head(diff.table)
nrow(diff.table)
### Read in the MEF/reprogramming comparison
output.label <- "Type2_myonuclei"
res.table <- read.table(file = paste0(out.dir, output.label, "_005pct_DESeq2_all.csv"),
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
nrow(res.table)
overlapping.peaks <- intersect(rownames(diff.table), rownames(res.table))
length(overlapping.peaks)
### Generate a correlation plot
ggData <- data.frame(Diff_LFC = diff.table[overlapping.peaks, "log2FoldChange"],
                     MEF_LFC = res.table[overlapping.peaks, "log2FoldChange"],
                     Peak = overlapping.peaks)
type2.cor <- ggplot(ggData, aes(x = Diff_LFC, y = MEF_LFC)) +
  geom_point(size = 0.5) + theme_classic(base_size = 15) +
  ylab(paste0("Diseased vs Control", " LFC")) + xlab(paste0("Aged vs Young", " LFC")) +
  geom_smooth(method = "lm", se = TRUE) + ggtitle( "B" ) +
  geom_hline(yintercept = 0, colour = "grey") + geom_vline(xintercept = 0, colour = "grey") +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~','~"))) +
  theme(plot.title = element_text(size = 18), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))
cor.test(ggData$Diff_LFC, ggData$MEF_LFC, method = "spearman")
cor.test(ggData$Diff_LFC, ggData$MEF_LFC, method = "spearman")$p.value
ggarrange(type1.cor, type2.cor, ncol =2)


## DAR peaks tracks ##
DAR.peaks <- GetExpressedPeaks(integrated_0.2, population.1 = type2.aged.cells, population.2 = type2.young.cells, threshold = 0.05)

roi = "chr12-6530500-6540000"
t2_dvc <- CoveragePlot(
  object = integrated_0.2[rownames(integrated_0.2) %in% DAR.peaks,],
  group.by = 'Donor',
  region = roi,
  annotation = TRUE,
  peaks = TRUE
) & theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12))

ggarrange(t1_avy, t2_avy, t1_dvc, t2_dvc, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))

## Peak plot ##
PeakPlot(object = integrated_0.2[rownames(integrated_0.2) %in% DAR.peaks,], region = roi) & ggtitle("Type II Myonuclei - Aged vs Young") & theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))
