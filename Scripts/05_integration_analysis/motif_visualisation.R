library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

#######################################################
### Generate a bubble plot for the motif enrichment ###
#######################################################
out.dir <- "/QRISdata/Q8448/Human_disease_Signac/Muscle/SObjects/universal_SCPM1/DAR/"

comparison.set <- c("Type1_myonuclei_open", "Type1_myonuclei_AgedNoBIvYoungNoBI_open", 
                    "Type2_myonuclei_open", "Type2_myonuclei_AgedNoBIvYoungNoBI_open")
names(comparison.set) <- c("Type1.disease", "Type1.age", "Type2.disease", "Type2.age")


comparisons.combined <- c(comparison.set)
# dar.direction <- "closing_vs_NS"
motifs.set <- c()
for (this.comparison in comparisons.combined) {
  this.file <- paste0(out.dir, "motif_res/", this.comparison, "/knownResults.txt")
  homer.table <- read.csv(this.file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  rownames(homer.table) <- make.unique(homer.table$Motif.Name)
  homer.table <- homer.table[1:10, ]
  motifs.set <- union(motifs.set, rownames(homer.table))
}
length(motifs.set)


# dar.direction <- "closing_vs_NS"
ggData <- c()
for (this.comparison in comparisons.combined) {
  
  this.file <- paste0(out.dir, "motif_res/", this.comparison, "/knownResults.txt")

  homer.table <- read.csv(this.file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  rownames(homer.table) <- make.unique(homer.table$Motif.Name)
  
  table.subset <- homer.table[motifs.set, ]
  
  pct.target.mat.aged <- as.numeric(gsub("%","",table.subset$X..of.Target.Sequences.with.Motif))
  pct.background.mat.aged <- as.numeric(gsub("%","",table.subset$X..of.Background.Sequences.with.Motif))
  
  fc.diff <- ((pct.target.mat.aged+0.001)/(pct.background.mat.aged+0.001))
  
  dar.direction <- ifelse(grepl("Aged", this.comparison), "A vs Y", "D vs C")
  ggData <- rbind(ggData, 
                  data.frame(
                    Comparison = rep(this.comparison, nrow(table.subset)),
                    DAR = rep(dar.direction, nrow(table.subset)),
                    Motif = table.subset$Motif.Name,
                    Target_pct = pct.target.mat.aged,
                    FC = log2(fc.diff+1),
                    Num_hits = table.subset[, 6],
                    Log10_pval = -log10(table.subset$P.value)))
  
}
ggData$Log10_pval[which(ggData$Log10_pval == Inf)] <- max(ggData$Log10_pval[-which(ggData$Log10_pval == Inf)])
ggData$DAR <- factor(ggData$DAR, levels = c("A vs Y", "D vs C"))
# # Repeat motifs in closing regions
# ggData$Motif[ggData$Motif == "GRE(NR),IR3/A549-GR-ChIP-Seq(GSE32465)/Homer"] <- "GRE(NR),IR3 (2)/A549-GR-ChIP-Seq(GSE32465)/Homer"
# motifs.set[motifs.set == "GRE(NR),IR3/A549-GR-ChIP-Seq(GSE32465)/Homer"] <- "GRE(NR),IR3 (2)/A549-GR-ChIP-Seq(GSE32465)/Homer"

ggData$Motif <- factor(ggData$Motif, levels = rev(motifs.set))
ggData$Comparison <- plyr::mapvalues(ggData$Comparison, from = as.character(comparisons.combined), to = names(comparisons.combined))
ggData$Comparison <- factor(ggData$Comparison, levels = c("Type1.age", "Type2.age", "Type1.disease", "Type2.disease"))


pl1 <- ggplot(ggData, aes(x = Comparison, y = reorder(Motif, Log10_pval), size = Num_hits, colour = Log10_pval)) +
  geom_point(alpha = 0.75) +
  facet_wrap(~DAR, nrow = 1) +
  theme_classic(base_size = 15) + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
  scale_size(range = c(1, 15), name="Num", limits = c(1, max(ggData$Num_hits))) +
  ylab("") + 
  scale_colour_gradientn(colours = c("#e3e1e1", "#7a7ef0", "blue"), values = scales::rescale(c(0, 10, 200), c(0,1))) 
pl1


pl1 <- ggplot(ggData, aes(x = Comparison, y = reorder(Motif, Log10_pval), size = FC, colour = Log10_pval)) +
  geom_point(alpha = 0.75) +
  facet_wrap(~DAR, nrow = 1) +
  theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_size(range = c(1, 15), name="LFC", limits = c(0,  max(ggData$FC))) +
  ylab("") + 
  scale_colour_gradientn(colours = c("#e3e1e1", "#7a7ef0", "blue"), values = scales::rescale(c(0, 10, 200), c(0,1))) 
pl1


pl2 <- ggplot(ggData, aes(x = Comparison, y = reorder(Motif, Log10_pval), size = Target_pct, colour = Log10_pval)) +
  geom_point(alpha = 0.75) + 
  labs(x = NULL, y = NULL) +
  ggtitle("B") +
  facet_wrap(~DAR, nrow = 1, scales = "free_x") +
  scale_size(range = c(1, 15), name="Target_pct", limits = c(0,  max(ggData$Target_pct))) +
  scale_y_discrete(labels = function(x) sapply(strsplit(x, "/"), '[', 1)) +
  theme_classic(base_size = 16) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 16)) +
  scale_colour_gradientn(colours = c("#e3e1e1", "#7a7ef0", "blue"), values = scales::rescale(c(0, 10, 200), c(0,1))) 
pl2

pl1 <- ggplot(ggData, aes(x = Comparison, y = reorder(Motif, Log10_pval), size = Target_pct, colour = Log10_pval)) +
  geom_point(alpha = 0.75) +
  ggtitle("A") +
  labs(x = NULL, y = NULL) +
  facet_wrap(~DAR, nrow = 1, scales = "free_x") +
  scale_size(range = c(1, 15), name="Target_pct", limits = c(0,  max(ggData$Target_pct))) +
  scale_y_discrete(labels = function(x) sapply(strsplit(x, "/"), '[', 1)) +
  theme_classic(base_size = 16) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 16)) +
  scale_colour_gradientn(colours = c("#e3e1e1", "#f57f7f","red"), values = scales::rescale(c(0, 10, 200), c(0,1))) 
pl1

ggarrange(pl1, pl2, ncol = 2)

