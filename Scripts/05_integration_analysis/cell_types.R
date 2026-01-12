load("/QRISdata/Q8448/Human_disease_Signac/Muscle/SObjects/universal_SCPM1/integration/Muscle_repSample_integration.RData")

DefaultAssay(integrated_0.2) <- 'ATAC'

diseased.md <- subset(integrated_0.2@meta.data, BI < 100)
table(diseased.md$Donor)
control.md <- subset(integrated_0.2@meta.data, BI == 100 & Age > 75) # aged healthy
table(control.md$Donor)
aged.md <- subset(integrated_0.2@meta.data, Age > 75)
table(aged.md)
young.md <- subset(integrated_0.2@meta.data, Age < 50)
table(young.md)
young.healthy.md <- subset(integrated_0.2@meta.data, BI == 100 & Age < 50)
table(young.healthy.md)
aged.nondiseased.md <- subset(aged.md, Donor != "P17" & Donor != "P21")
table(young.healthy.md)
aged.unknown.md <- subset(integrated_0.2@meta.data, grepl("OM", Donor))
table(aged.unknown.md)
young.unknown.md <- subset(integrated_0.2@meta.data, grepl("YM", Donor))
table(aged.unknown.md)

# Typr I myonuclei (cluster 0)
type1.diseased.cells <- rownames(subset(diseased.md, ATAC_snn_res.0.2 == "0"))
length(type1.diseased.cells)

type1.control.cells <- rownames(subset(control.md, ATAC_snn_res.0.2 == "0"))
length(type1.control.cells)

type1.aged.cells <- rownames(subset(aged.md, ATAC_snn_res.0.2 == "0"))
length(type1.aged.cells)

type1.young.cells <- rownames(subset(young.md, ATAC_snn_res.0.2 == "0"))
length(type1.young.cells)

type1.young.healthy.cells <- rownames(subset(young.healthy.md, ATAC_snn_res.0.2 == "0"))
length(type1.young.healthy.cells)

type1.aged.nondiseased.cells <- rownames(subset(aged.nondiseased.md, ATAC_snn_res.0.2 == "0"))
length(type1.aged.nondiseased.cells)

type1.aged.unknown.cells <- rownames(subset(aged.unknown.md, ATAC_snn_res.0.2 == "0"))
length(type1.aged.unknown.cells)

type1.young.unknown.cells <- rownames(subset(young.unknown.md, ATAC_snn_res.0.2 == "0"))
length(type1.young.unknown.cells)

# Type II Myonuclei (cluster 1)
type2.diseased.cells <- rownames(subset(diseased.md, ATAC_snn_res.0.2 == "1"))
length(type2.diseased.cells)

type2.control.cells <- rownames(subset(control.md, ATAC_snn_res.0.2 == "1"))
length(type2.control.cells)

type2.aged.cells <- rownames(subset(aged.md, ATAC_snn_res.0.2 == "1"))
length(type2.aged.cells)

type2.young.cells <- rownames(subset(young.md, ATAC_snn_res.0.2 == "1"))
length(type2.young.cells)

type2.young.healthy.cells <- rownames(subset(young.healthy.md, ATAC_snn_res.0.2 == "1"))
length(type2.young.healthy.cells)

type2.aged.nondiseased.cells <- rownames(subset(aged.nondiseased.md, ATAC_snn_res.0.2 == "1"))
length(type2.aged.nondiseased.cells)

type2.aged.unknown.cells <- rownames(subset(aged.unknown.md, ATAC_snn_res.0.2 == "1"))
length(type2.aged.unknown.cells)

type2.young.unknown.cells <- rownames(subset(young.unknown.md, ATAC_snn_res.0.2 == "1"))
length(type2.young.unknown.cells)

# MuSC (cluster 7)
MuSC.diseased.cells <- rownames(subset(diseased.md, ATAC_snn_res.0.2 == "7"))
length(MuSC.diseased.cells)

MuSC.control.cells <- rownames(subset(control.md, ATAC_snn_res.0.2 == "7"))
length(MuSC.control.cells)

MuSC.aged.cells <- rownames(subset(aged.md, ATAC_snn_res.0.2 == "7"))
length(MuSC.aged.cells)

MuSC.young.cells <- rownames(subset(young.md, ATAC_snn_res.0.2 == "7"))
length(MuSC.young.cells)

MuSC.aged.unknown.cells <- rownames(subset(aged.unknown.md, ATAC_snn_res.0.2 == "7"))
length(MuSC.aged.unknown.cells)

MuSC.young.unknown.cells <- rownames(subset(young.unknown.md, ATAC_snn_res.0.2 == "7"))
length(MuSC.young.unknown.cells)

# EC (cluster 4)
EC.diseased.cells <- rownames(subset(diseased.md, ATAC_snn_res.0.2 == "4"))
length(EC.diseased.cells)

EC.control.cells <- rownames(subset(control.md, ATAC_snn_res.0.2 == "4"))
length(EC.control.cells)

EC.aged.cells <- rownames(subset(aged.md, ATAC_snn_res.0.2 == "4"))
length(EC.aged.cells)

EC.young.cells <- rownames(subset(young.md, ATAC_snn_res.0.2 == "4"))
length(EC.young.cells)

EC.aged.unknown.cells <- rownames(subset(aged.unknown.md, ATAC_snn_res.0.2 == "4"))
length(EC.aged.unknown.cells)

EC.young.unknown.cells <- rownames(subset(young.unknown.md, ATAC_snn_res.0.2 == "4"))
length(EC.young.unknown.cells)

# SMC (cluster 6)
SMC.diseased.cells <- rownames(subset(diseased.md, ATAC_snn_res.0.2 == "6"))
length(SMC.diseased.cells)

SMC.control.cells <- rownames(subset(control.md, ATAC_snn_res.0.2 == "6"))
length(SMC.control.cells)

SMC.aged.cells <- rownames(subset(aged.md, ATAC_snn_res.0.2 == "6"))
length(SMC.aged.cells)

SMC.young.cells <- rownames(subset(young.md, ATAC_snn_res.0.2 == "6"))
length(SMC.young.cells)

SMC.aged.unknown.cells <- rownames(subset(aged.unknown.md, ATAC_snn_res.0.2 == "6"))
length(SMC.aged.unknown.cells)

SMC.young.unknown.cells <- rownames(subset(young.unknown.md, ATAC_snn_res.0.2 == "6"))
length(SMC.young.unknown.cells)

