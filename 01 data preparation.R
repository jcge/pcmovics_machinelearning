#01 data preparation
library(MOVICS)
library(purrr)
library(openxlsx)
library(tidyverse)
library(limma)
library(readr)
getwd()
setwd("/home/gejiachen/MOVICS")
# count
count <- read_tsv("TCGA-PAAD.star_counts.tsv.gz")
probeMap <- read.table("gencode.v36.annotation.gtf.gene.probemap",sep = "\t" , header = T)
count <- count %>%
  inner_join(probeMap, by = c("Ensembl_ID" = "id")) %>%
  dplyr::select(gene, starts_with("TCGA"))
count <- as.data.frame(avereps(count[,-1],ID = count$gene))
# colnames(count) <- substring(colnames(count),1,15) %>% gsub("-",".",.)
TCGA_group_list <- ifelse(as.numeric(substring(colnames(count),14,15)) < 10,
                          "Tumor","Normal") %>% 
  factor(.,levels = c("Normal","Tumor"))
table(TCGA_group_list)
# http://vega.archive.ensembl.org/info/about/gene_and_transcript_types.html http://asia.ensembl.org/
library(GenomicFeatures)
library(rtracklayer)
gc_data = import('gencode.v31.annotation.gtf')
gtf <- as.data.frame(gc_data)
gtf$seqnames %>% table()
table(gtf$transcript_type)
gtf <- subset(gtf, transcript_type == c("protein_coding","miRNA","lncRNA"))
mRNA_info <- subset(gtf, transcript_type == "protein_coding")
lncRNA_info <- subset(gtf, transcript_type == "lncRNA")
miRNA_info <- subset(gtf, transcript_type == "miRNA")
count <- count[rownames(count) %in% gtf$gene_name,]
count <- count[,sample]
saveRDS(count, "count.rds")

# fpkm
fpkm <- read_tsv("TCGA-PAAD.star_fpkm.tsv.gz")
fpkm <- fpkm %>%
  inner_join(probeMap, by = c("Ensembl_ID" = "id")) %>%
  dplyr::select(gene, starts_with("TCGA"))
fpkm = as.data.frame(avereps(fpkm[,-1],ID = fpkm$gene))
# colnames(fpkm) <- substring(colnames(fpkm),1,15) %>% gsub("-",".",.)
TCGA_group_list <- ifelse(as.numeric(substring(colnames(fpkm),14,15)) < 10,
                          "Tumor","Normal") %>% 
  factor(.,levels = c("Normal","Tumor"))
table(TCGA_group_list)
fpkm <- fpkm[rownames(fpkm) %in% gtf$gene_name,]
fpkm <- fpkm[,sample]
saveRDS(fpkm, "fpkm.rds")

#maf 
# maf <- read_tsv("TCGA-PAAD.somaticmutation_wxs.tsv.gz")
devtools::install_local("/home/gejiachen/rpackage/TCGAmutations-master.zip")
library(maftools) 
library(TCGAmutations)
paad.maf <- TCGAmutations::tcga_load(study = "PAAD")
maf <- paad.maf@maf.silent
maf <- maf[,c(23,1:22)]
maf <- maf[,c(1:9)]
saveRDS(maf, "maf.rds")

#segment
segment <- read_tsv("TCGA-PAAD.masked_cnv_DNAcopy.tsv.gz")
saveRDS(segment, "segment.rds")

sample <- colnames(count)
write.csv(sample, "sample.csv")
sample <- read.csv("sample.csv")
sample <- sample[,2]
# sample <- as.matrix(sample)
count <- count[,sample]
fpkm <- fpkm[,sample]
maf <- maf[maf$Tumor_Sample_Barcode_min == "sample",]
write.csv(maf, "maf.csv")
maf <- read.csv("maf.csv")
maf <- maf[which(maf$Tumor_Sample_Barcode_min == c(sample)),]
filter(maf,Tumor_Sample_Barcode_min == "sample")

#methylation
methylation <- read_tsv("TCGA-PAAD.methylation450.tsv.gz") 
methylation <- column_to_rownames(methylation,"Composite Element REF") 
# segment$sample = str_sub(segment$sample,1,15)
# colnames(mut.status) <- gsub("\\.","-",colnames(mut.status))
methylation <- methylation[,sample]
rownames(methylation) <- methylation[,1]
methylation <- na.omit(methylation)
methylation <- methylation[,sample]
# methylation <- mo.data[["meth.beta"]]
saveRDS(methylation, "methylation.rds")

#clinical
clinical <- read.table("clinical.txt",header = T)
rownames(clinical) <- clinical[,1]
clinical <- clinical[,-1]
clinical <- clinical[sample,]
saveRDS(clinical, "clinical.rds")

count <- count[,rownames(clin)]
#构建list
paad.tcga <- list()
paad.tcga[["mRNA.expr"]] <- mRNA.expr
paad.tcga[["miRNA.expr"]] <- miRNA.expr
paad.tcga[["lncRNA.expr"]] <- lncRNA.expr
paad.tcga[["meth.beta"]] <- methylation
paad.tcga[["mut.status"]] <- mut.status
paad.tcga[["count"]] <- count
paad.tcga[["fpkm"]] <- fpkm
paad.tcga[["maf"]] <- maf
paad.tcga[["segment"]] <- segment
paad.tcga[["clin.info"]] <- clinical
saveRDS(paad.tcga,"paad.tcga.rds")
#cox按照p值筛选
if(F){
  names(paad.tcga)
  surv.info <- paad.tcga$clin.info
  mo.data = lapply(paad.tcga[c(1:4)], function(x){
    x %>% 
      getElites(elite.num = 500) %>% 
      pluck('elite.dat') %>% 
      getElites(method= "cox",surv.info = surv.info,p.cutoff = 0.05) %>%
      pluck('elite.dat')
  })
  mo.data$mut.status = getElites(dat = paad.tcga$mut.status,
                                 method = "freq",
                                 elite.pct = 0.02) %>% 
    pluck('elite.dat') 
  sapply(mo.data, dim)
}
saveRDS(mo.data,"mo.data.rds")
# 整理预后相关mRNA,miRNA,lncRNA
expr <- mo.data[["fpkm"]]
mRNA <- mRNA_info$gene_name
mRNA <- unique(mRNA) 
mRNA <- intersect(mRNA,rownames(count))
mRNA.expr <- count[mRNA,]
mRNA.expr <- mRNA.expr[1:1500,]
mRNA.expr <- mRNA.expr[,sample]
saveRDS(mRNA.expr,"mRNA.expr.rds")

miRNA <- miRNA_info$gene_name
miRNA <- unique(miRNA) 
miRNA <- intersect(miRNA,rownames(count))
miRNA.expr <- count[miRNA,]
miRNA.expr <- miRNA.expr[1:1500,]
miRNA.expr <- miRNA.expr[,sample]
saveRDS(miRNA.expr,"miRNA.expr.rds")

lncRNA <- lncRNA_info$gene_name
lncRNA <- unique(lncRNA) 
lncRNA <- intersect(lncRNA,rownames(count))
lncRNA.expr <- count[lncRNA,]
lncRNA.expr <- lncRNA.expr[1:1500,]
lncRNA.expr <- lncRNA.expr[,sample]
saveRDS(lncRNA.expr,"lncRNA.expr.rds")
# MAD筛选分子 
elite.mRNA <- getElites(dat = mRNA.expr,
                        method = "sd",
                        elite.num = 800,
                        elite.pct = 0.1) # 保留MAD前10%的基因
elite.miRNA <- getElites(dat = miRNA.expr,
                         method = "sd",
                         elite.num = 800,
                         elite.pct = 0.1) # 保留MAD前10%的基因
elite.lncRNA <- getElites(dat = lncRNA.expr,
                          method = "sd",
                          elite.num = 800,
                          elite.pct = 0.1) # 保留MAD前10%的基因
elite.methylation <- getElites(dat = methylation,
                               method = "sd",
                               elite.num = 800,
                               elite.pct = 0.1) # 保留MAD前10%的基因
mRNA.expr <- elite.mRNA[["elite.dat"]]
miRNA.expr <- elite.miRNA[["elite.dat"]]
lncRNA.expr <- elite.lncRNA[["elite.dat"]]
methylation <- elite.methylation[["elite.dat"]]

# 突变组：宽数据转变为长数据
paad.maf <- maf[,1:2]
## 1.计算频次：group_by用于分类
dd1 <- paad.maf %>% 
  group_by(Tumor_Sample_Barcode_min, Hugo_Symbol)%>%  ##用
  summarise(num=n()) 
head(dd1)
## 2.数据变宽
library(tidyr)
dd2 <- dd1%>% 
  pivot_wider(names_from = Hugo_Symbol,values_from = num)
## 3.na 变成0
#用is.na(dd2)产生逻辑值，用于定位
dd2[is.na(dd2)]=0
dd2 <- t(dd2)
colnames(dd2) <- dd2[1,]
dd2 <- dd2[-1,]
dd2 <- as.data.frame(dd2)
saveRDS(dd2,"mut.status.rds")
sample <- colnames(mut.status)
mut.status[mut.status >= 1] = 1
mut.status[mut.status == " 0"] = 0
mut.status[mut.status == " 1"] = 1
table(mut.status)
mut.status <- gsub(" ","",mut.status)
mut.status <- mut.status[,sample]

write.csv(mut.status,"mut.status.csv",sep = ",")
write.csv(a,"a.csv",sep = ",")
miRNA.expr <- t(miRNA.expr)
miRNA.expr <- na.omit(miRNA.expr)
mut.status <- read.csv("mut.status.csv",row.names = 1,check.names = F)
colnames(mut.status) <- gsub("\\.","-",colnames(mut.status))
sample <- intersect(sample,colnames(mut.status))