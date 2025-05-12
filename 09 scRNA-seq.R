#09 scRNA-seq Figure 6 and Supplementary Figure 5 
getwd()
#整合
setwd("/home/gejiachen/scrnaseq/CRA001160")
DimPlot(data, reduction = "umap", label = TRUE,raster=FALSE) 
Idents(data) <- data@meta.data$celltype
table(Idents(data))
Macrophage <- subset(x = data,idents=c("Mono/Macro"))
Macrophage$dataset <- "CRA001160"
Ductal <- subset(x = data,idents=c("Ductal"))
Ductal$dataset <- "CRA001160"
saveRDS(Macrophage,"Macrophage_CRA001160.rds")
saveRDS(Ductal,"Ductal_CRA001160.rds")
DimPlot(data, reduction = "umap",label = T)
list_genes=list(Immu=c("PTPRC"),
                Tcell=c("CD3D","CD3E"),
                CD8T=c("CD8A","CD8B"),
                CD8Tex=c("PDCD1","LAG3"),
                B=c("CD79A","MS4A1"),
                Plasma=c("SLAMF7"),
                MonoMacro=c("CD14","CSF1R","CD68","MMP9"),
                Endo=c("VWF","CHGB"),
                Fibro=c("COL1A2"),
                Malignant=c("EGLN3"),
                Peric=c("RGS5"),
                Acinar=c("PRSS1"),
                Ductal=c("AMBP","KRT19"))
DotPlot(data,features=list_genes,cols = c("grey", "red"),cluster.idents = T) + 
  RotatedAxis() + theme(strip.text = element_text(margin=margin(b=3, unit="mm"), size = 8))
VlnPlot(data, features = c("A2ML1","FAM83A"), group.by = "celltype", stack = TRUE, sort = TRUE) +
  theme(legend.position = "none") + ggtitle("CRA001160 cohort")

setwd("/home/gejiachen/scrnaseq/GSE111672")
sc.data <- Read10X_h5("PAAD_GSE111672_expression.h5")
metadata <- read.delim(file.choose(),row.names = 1)
data <- CreateSeuratObject(counts = sc.data, project = "GSE111672",meta.data = metadata)  #如读取了metadata运行这个代码，不需要运行上一个代码
data <- NormalizeData(object = data)
data <- FindVariableFeatures(object = data)
data <- ScaleData(object = data)
data <- RunPCA(object = data)
data <- RunUMAP(object = data, dims = 1:10)
data <- RunTSNE(object = data, dims = 1:10)
Idents(data) <- data@meta.data$Celltype..minor.lineage.
table(data@meta.data$Celltype..minor.lineage.) 
save(data, file = "GSE111672.Rdata")
Macrophage <- subset(x = data,idents=c("M1"))
Macrophage$dataset <- "GSE111672"
Ductal <- subset(x = data,idents=c("Ductal"))
Ductal$dataset <- "GSE111672"
saveRDS(Macrophage,"Macrophage_GSE111672.rds")
saveRDS(Ductal,"Ductal_GSE111672.rds")
DimPlot(data, reduction = "umap",label = T)
DimPlot(data, reduction = "tsne",label = T)
list_genes=list(Tcell=c("CD3D","CD3E"),
                CD8T=c("CD8A","CD8B"),
                Tproli=c("MKI67"),
                Mast=c("CD79A","MS4A1"),
                Neutr=c("S100A9"),
                MonoMacro=c("CD14","CD68"),
                Endo=c("VWF","PECAM1"),
                Fibro=c("COL1A1"),
                Acinar=c("PRSS1"),
                Ductal=c("AMBP","KRT19"))
DotPlot(data,features=list_genes,cols = c("grey", "red"),cluster.idents = T) + 
  RotatedAxis() + theme(strip.text = element_text(margin=margin(b=3, unit="mm"), size = 8))
VlnPlot(data, features = c("A2ML1","FAM83A"), group.by = "Celltype..major.lineage.", stack = TRUE, sort = TRUE) +
  theme(legend.position = "none") + ggtitle("GSE111672 cohort")

setwd("/home/gejiachen/scrnaseq/GSE154778")
sc.data <- Read10X_h5("PAAD_GSE154778_expression.h5")
metadata <- read.delim(file.choose(),row.names = 1)
data <- CreateSeuratObject(counts = sc.data, project = "GSE154778",meta.data = metadata)  #如读取了metadata运行这个代码，不需要运行上一个代码
data <- NormalizeData(object = data)
data <- FindVariableFeatures(object = data)
data <- ScaleData(object = data)
data <- RunPCA(object = data)
data <- RunUMAP(object = data, dims = 1:10)
data <- RunTSNE(object = data, dims = 1:10)
Idents(data) <- data@meta.data$Celltype..minor.lineage.
table(data@meta.data$Celltype..minor.lineage.) 
save(data, file = "GSE154778.Rdata")
Macrophage <- subset(x = data,idents=c("Mono/Macro"))
Macrophage$dataset <- "GSE154778"
Ductal <- subset(x = data,idents=c("Epithelial","Malignant"))
Ductal$dataset <- "GSE154778"
saveRDS(Macrophage,"Macrophage_GSE154778.rds")
saveRDS(Ductal,"Ductal_GSE154778.rds")
DimPlot(data, reduction = "umap",label = T)
DimPlot(data, reduction = "tsne",label = T)
list_genes=list(Immu=c("PTPRC"),
                Tcell=c("CD3D","CD3E"),
                CD8T=c("CD8A","CD8B"),
                Plasma=c("JCHAIN"),
                MonoMacro=c("CD14","CD68","CST3","CSF1R"),
                Fibroblasts=c("COL6A1","MYL9","COL1A1","MMP2","COL1A2","MYLK","COL3A1"),
                Epith=c("EPCAM"),
                Malignant=c("CDH1","CCND1","MET","SOX4"))
DotPlot(data,features=list_genes,cols = c("grey", "red"),cluster.idents = T) + 
  RotatedAxis() + theme(strip.text = element_text(margin=margin(b=3, unit="mm"), size = 8))
VlnPlot(data, features = c("A2ML1","FAM83A"), group.by = "Celltype..major.lineage.", stack = TRUE, sort = TRUE) +
  theme(legend.position = "none") + ggtitle("GSE154778 cohort")

setwd("/home/gejiachen/scrnaseq/GSE155698")
Idents(TotalTissue.combined.harmony) <- TotalTissue.combined.harmony@meta.data$Expanded_labels
table(TotalTissue.combined.harmony@meta.data$Expanded_labels) 
Macrophage <- subset(x = TotalTissue.combined.harmony,idents=c("Mono&Macro"))
Macrophage$dataset <- "GSE155698"
Ductal <- subset(x = TotalTissue.combined.harmony,idents=c("Epithelial"))
Ductal$dataset <- "GSE155698"
saveRDS(Macrophage,"Macrophage_GSE155698.rds")
saveRDS(Ductal,"Ductal_GSE155698.rds")
DimPlot(TotalTissue.combined.harmony, reduction = "umap",label = T)
DimPlot(TotalTissue.combined.harmony, reduction = "tsne",label = T)
list_genes=list(Acinar=c('PRSS1','CTRB2','REG1A','SPINK1'),
                Fibroblasts=c('ACTA2','CDH11','PDGFRB','COL1A1','COL3A1','RGS5','IGFBP7','PDPN','DCN',
                              'MCAM','PDGFA'),
                Endothelial=c('VWF','CDH5'),
                Mast=c('TPSAB1','CPA3'),
                Myeloid=c('ITGAM','FCGR3A','FCGR3B','APOE','C1QA','MARCO','LYZ','HLA-DRA'), 
                DC=c('BATF3','IRF8','CD1C','CD101'),
                TandNK=c('CD3D','CD3E','NKG7','CD8A','PRF1','GZMB','CD69','TIGIT'),
                MonoMacro=c("CD14","CD68","CSF1R"),
                B=c('CD79A','MS4A1','CD20','CD138','IGJ','IGLL5','CXCR4','CD117','CD27'),
                Epithelial=c("EPCAM","CDH1","CCND1","MET","SOX4"))
DotPlot(TotalTissue.combined.harmony,features=list_genes,cols = c("grey", "red"),cluster.idents = T) + 
  RotatedAxis() + theme(strip.text = element_text(margin=margin(b=3, unit="mm"), size = 8))
VlnPlot(TotalTissue.combined.harmony, features = c("A2ML1","FAM83A"), group.by = "Expanded_labels", stack = TRUE, sort = TRUE) +
  theme(legend.position = "none") + ggtitle("GSE155698 cohort")

setwd("/home/gejiachen/scrnaseq/GSE158356")
sc.data <- Read10X_h5("PAAD_GSE158356_expression.h5")
metadata <- read.delim(file.choose(),row.names = 1)
data <- CreateSeuratObject(counts = sc.data, project = "GSE158356",meta.data = metadata)  #如读取了metadata运行这个代码，不需要运行上一个代码
data <- NormalizeData(object = data)
data <- FindVariableFeatures(object = data)
data <- ScaleData(object = data)
data <- RunPCA(object = data)
data <- RunUMAP(object = data, dims = 1:10)
data <- RunTSNE(object = data, dims = 1:10)
Idents(data) <- data@meta.data$Celltype..minor.lineage.
table(data@meta.data$Celltype..minor.lineage.) 
save(data, file = "GSE158356.Rdata")
Macrophage <- subset(x = data,idents=c("Mono/Macro"))
Macrophage$dataset <- "GSE158356"
Ductal <- subset(x = data,idents=c("Epithelial"))
Ductal$dataset <- "GSE158356"
saveRDS(Macrophage,"Macrophage_GSE158356.rds")
saveRDS(Ductal,"Ductal_GSE158356.rds")
DimPlot(data, reduction = "umap",label = T)
DimPlot(data, reduction = "tsne",label = T)
list_genes=list(Immu=c("PTPRC"),
                Tcell=c("CD3D","CD3E"),
                CD8T=c("CD8A","CD8B"),
                B=c("CD79A","CD19","MS4A1"),
                MonoMacro=c("CD14","CD68","CST3","CSF1R"),
                Fibroblasts=c("COL6A1","MYL9","COL1A1","MMP2","COL1A2","PDGFA","MYLK","COL3A1"),
                Epith=c("EPCAM"),
                Acinar=c("PRSS1"))
DotPlot(data,features=list_genes,cols = c("grey", "red"),cluster.idents = T) + 
  RotatedAxis() + theme(strip.text = element_text(margin=margin(b=3, unit="mm"), size = 8))
VlnPlot(data, features = c("A2ML1","FAM83A"), group.by = "Celltype..major.lineage.", stack = TRUE, sort = TRUE) +
  theme(legend.position = "none") + ggtitle("GSE158356 cohort")

setwd("/home/gejiachen/scrnaseq/GSE212966")
pdac <- readRDS("~/scrnaseq/GSE212966/pdac.rds")
Idents(pdac) <- pdac@meta.data$celltype
table(pdac@meta.data$celltype) 
Macrophage <- subset(x = pdac,idents=c("Mono/Macro"))
Macrophage$dataset <- "GSE212966"
Ductal <- subset(x = pdac,idents=c("Ductal"))
Ductal$dataset <- "GSE212966"
saveRDS(Macrophage,"Macrophage_GSE212966.rds")
saveRDS(Ductal,"Ductal_GSE212966.rds")
Idents(data) <- data@meta.data$celltype
DimPlot(data, reduction = "umap",label = T)
DimPlot(data, reduction = "tsne",label = T)
list_genes=list(Immu=c("PTPRC"),
                Tcell=c("CD3D","CD3E","GZMA","CCL5"),
                Stellelate=c("C11orf96","ADIRF","MYH11"),
                Plasma=c("CD79A","CD19"),
                Neutrophil=c("CXCL8","G0S2","IL1B"),
                MonoMacro=c("C1QB","CD68","CD14","SPP1"),
                Fibroblasts=c("COL1A1","COL1A2","DUM","LUM"),
                Ductal=c("EPCAM","KRT17","TFF1","TFF2","APOE"),
                Endothelial=c("CDH5","PLVAP","PECAM1","VWF"),
                Schwann=c("ITGB8"),
                Mast=c("CPA3","MS4A2"))
DotPlot(data,features=list_genes,cols = c("grey", "red"),cluster.idents = T) + 
  RotatedAxis() + theme(strip.text = element_text(margin=margin(b=3, unit="mm"), size = 8))
VlnPlot(data, features = c("A2ML1","FAM83A"), group.by = "celltype", stack = TRUE, sort = TRUE) +
  theme(legend.position = "none") + ggtitle("GSE212966 cohort")

setwd("/home/gejiachen/scrnaseq/GSE214295")
data <- readRDS("~/scrnaseq/GSE214295/data.rds")
Idents(data) <- data@meta.data$celltype
table(data@meta.data$celltype) 
Macrophage <- subset(x = data,idents=c("Mono/Macro"))
Macrophage$dataset <- "GSE214295"
Ductal <- subset(x = data,idents=c("Ductal"))
Ductal$dataset <- "GSE214295"
saveRDS(Macrophage,"Macrophage_GSE214295.rds")
saveRDS(Ductal,"Ductal_GSE214295.rds")
Idents(data) <- data@meta.data$celltype
DimPlot(data, reduction = "umap",label = T)
DimPlot(data, reduction = "tsne",label = T)
list_genes=list(Immu=c("PTPRC"),
                B=c("CD79A","MS4A1"),
                Tcell=c("CD3D","CD3E","GZMA","CCL5"),
                Stellelate=c("C11orf96","ADIRF","MYH11"),
                Pericyte=c("RGS5"),
                Acinar=c("PRSS1"),
                MonoMacro=c("C1QB","CD68","CD14","SPP1"),
                Fibroblasts=c("COL1A1","COL1A2","DUM","LUM"),
                Ductal=c("EPCAM","KRT17","TFF1","TFF2","APOE"),
                Endothelial=c("CDH5","PLVAP","PECAM1","VWF"),
                Schwann=c("ITGB8"),
                Mast=c("CPA3","MS4A2"))
DotPlot(data,features=list_genes,cols = c("grey", "red"),cluster.idents = T) + 
  RotatedAxis() + theme(strip.text = element_text(margin=margin(b=3, unit="mm"), size = 8))
VlnPlot(data, features = c("A2ML1","FAM83A"), group.by = "celltype", stack = TRUE, sort = TRUE) +
  theme(legend.position = "none") + ggtitle("GSE214295 cohort")

setwd("/home/gejiachen/scrnaseq/GSE217845")
data <- readRDS("~/scrnaseq/GSE217845/data.rds")
Idents(data) <- data@meta.data$celltype
table(data@meta.data$celltype) 
Macrophage <- subset(x = data,idents=c("Mono/Macro"))
Macrophage$dataset <- "GSE217845"
Ductal <- subset(x = data,idents=c("Ductal","Malignant"))
Ductal$dataset <- "GSE217845"
saveRDS(Macrophage,"Macrophage_GSE217845.rds")
saveRDS(Ductal,"Ductal_GSE217845.rds")
Idents(data) <- data@meta.data$celltype
DimPlot(data, reduction = "umap",label = T)
DimPlot(data, reduction = "tsne",label = T)
list_genes=list(Immu=c("PTPRC"),
                B=c("MS4A1","CD79A","CD79B","JCHAIN","MZB1","IGHA2","IGHG2"),
                Tcell=c("CD3E","GZMA","KLRD1","PRF1"),
                MonoMacro=c("C1QB","CD68","CD14"),
                Mast=c("CPA3","KIT"),
                Acinar=c("CELA3A","CLPS"),
                Ductal=c("CFTR"),
                Endocrine=c("INS","GCG"),
                Malignant=c("KRT17","CEACAM6","LAMA3","SOX9","SPINK1"),
                Endothelial=c("CDH5","PLVAP","PECAM1","VWF"),
                Fibroblasts=c("COL1A1","COL1A2","PDGFRA","PDPN","ACTA2"))
DotPlot(data,features=list_genes,cols = c("grey", "red"),cluster.idents = T) + 
  RotatedAxis() + theme(strip.text = element_text(margin=margin(b=3, unit="mm"), size = 8))
VlnPlot(data, features = c("A2ML1","FAM83A"), group.by = "celltype", stack = TRUE, sort = TRUE) +
  theme(legend.position = "none") + ggtitle("GSE217845 cohort")
