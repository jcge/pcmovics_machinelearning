#空间转录组
rm(list=ls())
options(stringsAsFactors = F)
library(ggsci)
library(dplyr) 
library(future)
library(Seurat)
library(clustree)
library(cowplot)
library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(qs)
library(Matrix)
setwd("/home/gejiachen/MOVICS/xiuhui/03/sp")
Spatial_data <- readRDS("PDAC_Updated.rds")
Spatial_data
table(Spatial_data$patient)
table(Spatial_data$orig.ident) 
metadata <- Spatial_data@meta.data
write.csv(metadata,"metadata.csv")
metadata <- read.csv("metadata.csv",header = T,row.names = 1)

exp_data <- FetchData(Spatial_data, vars = c("A2ML1"))
table(exp_data[,1])
Spatial_data@meta.data$A2ML1 <- NULL

Spatial_data$A2ML1 <- metadata$A2ML1
table(Spatial_data$Histology,Spatial_data$A2ML1)
table(Spatial_data$orig.ident,Spatial_data$A2ML1)
Idents(Spatial_data) <- Spatial_data@meta.data$second_type
#原发癌
SpatialDimPlot(Spatial_data, 
               images = "IU_PDA_T4", 
               alpha = 0,
               pt.size.factor = 0)
SpatialDimPlot(Spatial_data, label = TRUE, label.size = 3, 
               images = "IU_PDA_T4", group.by = "first_type", 
               pt.size.factor = 1.2)
SpatialDimPlot(Spatial_data, label = TRUE, label.size = 3, 
               images = "IU_PDA_T4", group.by = "A2ML1", 
               pt.size.factor = 1.2)
?SpatialDimPlot
# SpatialFeaturePlot(Spatial_data, 
#                    images = "IU_PDA_T4",alpha = 0, 
#                    features = NULL)
# ?SpatialFeaturePlot
#肝转移癌
SpatialDimPlot(Spatial_data, 
               images = "IU_PDA_HM4", 
               alpha = 0,
               pt.size.factor = 0)
SpatialDimPlot(Spatial_data, label = TRUE, label.size = 3, 
               images = "IU_PDA_HM4", group.by = "first_type", 
               pt.size.factor = 1.2)
SpatialDimPlot(Spatial_data, label = TRUE, label.size = 3, 
               images = "IU_PDA_HM4", group.by = "A2ML1", 
               pt.size.factor = 1.2)
#淋巴结转移
SpatialDimPlot(Spatial_data, 
               images = "IU_PDA_LNM10", 
               alpha = 0,
               pt.size.factor = 0)
SpatialDimPlot(Spatial_data, label = TRUE, label.size = 3, 
               images = "IU_PDA_LNM10", group.by = "first_type", 
               pt.size.factor = 1.2)
SpatialDimPlot(Spatial_data, label = TRUE, label.size = 3, 
               images = "IU_PDA_LNM10", group.by = "A2ML1", 
               pt.size.factor = 1.2)
#正常胰腺
SpatialDimPlot(Spatial_data, 
               images = "IU_PDA_NP2", 
               alpha = 0,
               pt.size.factor = 0)
SpatialDimPlot(Spatial_data, label = TRUE, label.size = 3, 
               images = "IU_PDA_NP2", group.by = "first_type", 
               pt.size.factor = 1.2)
SpatialDimPlot(Spatial_data, label = TRUE, label.size = 3, 
               images = "IU_PDA_NP2", group.by = "A2ML1", 
               pt.size.factor = 1.2)