#03 NTP预测外部数据 external validation Figure 3 and Supplementary Figure 1    
# 整合GSE28735
df <- GSE28735[!duplicated(GSE28735[,c("Tag")]),]
GSE28735 <- df
rownames(GSE28735) <- GSE28735$Tag
GSE28735 <- GSE28735[,-1]
GSE28735 <- as.data.frame(GSE28735)
sample <- read.table("sample.txt",header = T)
sample <- subset(sample, Cohort == "GSE28735")
GSE28735 <- GSE28735[,sample$Sample]
df <- list()
df[["mRNA.expr"]] <- GSE28735 
GSE28735 <- df
# run NTP in Yau cohort by using up-regulated biomarkers
GSE28735.ntp.pred <- runNTP(expr  = GSE28735$mRNA.expr,
                            templates  = marker.up$templates, # the template has been already prepared in runMarker()
                            scaleFlag  = TRUE, # scale input data (by default)
                            centerFlag = TRUE, # center input data (by default)
                            doPlot     = TRUE, # to generate heatmap
                            fig.name   = "NTP HEATMAP FOR GSE28735") 
write.csv(GSE28735.ntp.pred[["ntp.res"]],"GSE28735group.csv")
# runKappa(subt1     = cmoic.paad$clust.res$clust,
#          subt2     = GSE28735.ntp.pred$clust.res$clust,
#          subt1.lab = "CMOIC",
#          subt2.lab = "NTP",
#          fig.name  = "CONSISTENCY HEATMAP FOR GSE28735 between CMOIC and NTP")
# 整合GSE57495
df <- GSE57495[!duplicated(GSE57495[,c("Tag")]),]
GSE57495 <- df
rownames(GSE57495) <- GSE57495$Tag
GSE57495 <- GSE57495[,-1]
GSE57495 <- as.data.frame(GSE57495)
setwd("/home/gejiachen/MOVICS")
sample <- read.table("sample.txt",header = T)
sample <- subset(sample, Cohort == "GSE57495")
GSE57495 <- GSE57495[,sample$Sample]
df <- list()
df[["mRNA.expr"]] <- GSE57495 
GSE57495 <- df
# run NTP in Yau cohort by using up-regulated biomarkers
GSE57495.ntp.pred <- runNTP(expr  = GSE57495$mRNA.expr,
                            templates  = marker.up$templates, # the template has been already prepared in runMarker()
                            scaleFlag  = TRUE, # scale input data (by default)
                            centerFlag = TRUE, # center input data (by default)
                            doPlot     = TRUE, # to generate heatmap
                            fig.name   = "NTP HEATMAP FOR GSE57495") 
write.csv(GSE57495.ntp.pred[["ntp.res"]],"GSE57495group.csv")
# 整合GSE62452
df <- GSE62452[!duplicated(GSE62452[,c("Tag")]),]
GSE62452 <- df
rownames(GSE62452) <- GSE62452$Tag
GSE62452 <- GSE62452[,-1]
GSE62452 <- as.data.frame(GSE62452)
sample <- read.table("sample.txt",header = T)
sample <- subset(sample, Cohort == "GSE62452")
GSE62452 <- GSE62452[,sample$Sample]
df <- list()
df[["mRNA.expr"]] <- GSE62452 
GSE62452 <- df
# run NTP in Yau cohort by using up-regulated biomarkers
GSE62452.ntp.pred <- runNTP(expr  = GSE62452$mRNA.expr,
                            templates  = marker.up$templates, # the template has been already prepared in runMarker()
                            scaleFlag  = TRUE, # scale input data (by default)
                            centerFlag = TRUE, # center input data (by default)
                            doPlot     = TRUE, # to generate heatmap
                            fig.name   = "NTP HEATMAP FOR GSE62452") 
write.csv(GSE62452.ntp.pred[["ntp.res"]],"GSE62452group.csv")
# 整合GSE21501
df <- GSE21501[!duplicated(GSE21501[,c("Tag")]),]
GSE21501 <- df
rownames(GSE21501) <- GSE21501$Tag
GSE21501 <- GSE21501[,-1]
GSE21501 <- as.data.frame(GSE21501)
sample <- read.table("sample.txt",header = T)
sample <- subset(sample, Cohort == "GSE21501")
GSE21501 <- GSE21501[,sample$Sample]
df <- list()
df[["mRNA.expr"]] <- GSE21501 
GSE21501 <- df
# run NTP in Yau cohort by using up-regulated biomarkers
GSE21501.ntp.pred <- runNTP(expr  = GSE21501$mRNA.expr,
                            templates  = marker.up$templates, # the template has been already prepared in runMarker()
                            scaleFlag  = TRUE, # scale input data (by default)
                            centerFlag = TRUE, # center input data (by default)
                            doPlot     = TRUE, # to generate heatmap
                            fig.name   = "NTP HEATMAP FOR GSE21501") 
write.csv(GSE21501.ntp.pred[["ntp.res"]],"GSE21501group.csv")
# 整合GSE71729
df <- GSE71729[!duplicated(GSE71729[,c("Tag")]),]
GSE71729 <- df
rownames(GSE71729) <- GSE71729$Tag
GSE71729 <- GSE71729[,-1]
GSE71729 <- as.data.frame(GSE71729)
sample <- read.table("sample.txt",header = T)
sample <- subset(sample, Cohort == "GSE71729")
GSE71729 <- GSE71729[,sample$Sample]
df <- list()
df[["mRNA.expr"]] <- GSE71729 
GSE71729 <- df
# run NTP in Yau cohort by using up-regulated biomarkers
GSE71729.ntp.pred <- runNTP(expr  = GSE71729$mRNA.expr,
                            templates  = marker.up$templates, # the template has been already prepared in runMarker()
                            scaleFlag  = TRUE, # scale input data (by default)
                            centerFlag = TRUE, # center input data (by default)
                            doPlot     = TRUE, # to generate heatmap
                            fig.name   = "NTP HEATMAP FOR GSE71729") 
write.csv(GSE71729.ntp.pred[["ntp.res"]],"GSE71729group.csv")
# 整合GSE78229
df <- GSE78229[!duplicated(GSE78229[,c("Tag")]),]
GSE78229 <- df
rownames(GSE78229) <- GSE78229$Tag
GSE78229 <- GSE78229[,-1]
GSE78229 <- as.data.frame(GSE78229)
sample <- read.table("sample.txt",header = T)
sample <- subset(sample, Cohort == "GSE78229")
GSE78229 <- GSE78229[,sample$Sample]
df <- list()
df[["mRNA.expr"]] <- GSE78229 
GSE78229 <- df
# run NTP in Yau cohort by using up-regulated biomarkers
GSE78229.ntp.pred <- runNTP(expr  = GSE78229$mRNA.expr,
                            templates  = marker.up$templates, # the template has been already prepared in runMarker()
                            scaleFlag  = TRUE, # scale input data (by default)
                            centerFlag = TRUE, # center input data (by default)
                            doPlot     = TRUE, # to generate heatmap
                            fig.name   = "NTP HEATMAP FOR GSE78229") 
write.csv(GSE78229.ntp.pred[["ntp.res"]],"GSE78229group.csv")
# 整合GSE79668
df <- GSE79668[!duplicated(GSE79668[,c("Tag")]),]
GSE79668 <- df
rownames(GSE79668) <- GSE79668$Tag
GSE79668 <- GSE79668[,-1]
GSE79668 <- as.data.frame(GSE79668)
sample <- read.table("sample.txt",header = T)
sample <- subset(sample, Cohort == "GSE79668")
GSE79668 <- GSE79668[,sample$Sample]
df <- list()
df[["mRNA.expr"]] <- GSE79668 
GSE79668 <- df
# run NTP in Yau cohort by using up-regulated biomarkers
GSE79668.ntp.pred <- runNTP(expr  = GSE79668$mRNA.expr,
                            templates  = marker.up$templates, # the template has been already prepared in runMarker()
                            scaleFlag  = TRUE, # scale input data (by default)
                            centerFlag = TRUE, # center input data (by default)
                            doPlot     = TRUE, # to generate heatmap
                            fig.name   = "NTP HEATMAP FOR GSE79668") 
write.csv(GSE79668.ntp.pred[["ntp.res"]],"GSE79668group.csv")
# 整合GSE85916
df <- GSE85916[!duplicated(GSE85916[,c("Tag")]),]
GSE85916 <- df
rownames(GSE85916) <- GSE85916$Tag
GSE85916 <- GSE85916[,-1]
GSE85916 <- as.data.frame(GSE85916)
sample <- read.table("sample.txt",header = T)
sample <- subset(sample, Cohort == "GSE85916")
GSE85916 <- GSE85916[,sample$Sample]
df <- list()
df[["mRNA.expr"]] <- GSE85916 
GSE85916 <- df
# run NTP in Yau cohort by using up-regulated biomarkers
GSE85916.ntp.pred <- runNTP(expr  = GSE85916$mRNA.expr,
                            templates  = marker.up$templates, # the template has been already prepared in runMarker()
                            scaleFlag  = TRUE, # scale input data (by default)
                            centerFlag = TRUE, # center input data (by default)
                            doPlot     = TRUE, # to generate heatmap
                            fig.name   = "NTP HEATMAP FOR GSE85916") 
write.csv(GSE85916.ntp.pred[["ntp.res"]],"GSE85916group.csv")
# 整合E_MTAB_6134
df <- EMTAB6134[!duplicated(EMTAB6134[,c("Tag")]),]
E_MTAB_6134 <- df
rownames(E_MTAB_6134) <- df$Tag
E_MTAB_6134 <- E_MTAB_6134[,-1]
E_MTAB_6134 <- as.data.frame(E_MTAB_6134)
sample <- read.table("sample.txt",header = T)
sample <- subset(sample, Cohort == "E_MTAB_6134")
E_MTAB_6134 <- E_MTAB_6134[,sample$Sample]
df <- list()
df[["mRNA.expr"]] <- E_MTAB_6134 
E_MTAB_6134 <- df
# run NTP in Yau cohort by using up-regulated biomarkers
E_MTAB_6134.ntp.pred <- runNTP(expr  = E_MTAB_6134$mRNA.expr,
                               templates  = marker.up$templates, # the template has been already prepared in runMarker()
                               scaleFlag  = TRUE, # scale input data (by default)
                               centerFlag = TRUE, # center input data (by default)
                               doPlot     = TRUE, # to generate heatmap
                               fig.name   = "NTP HEATMAP FOR E_MTAB_6134") 
write.csv(E_MTAB_6134.ntp.pred[["ntp.res"]],"E_MTAB_6134group.csv")
# 整合PACA_AU_arraySeq
df <- PACA_AU_arraySeq[!duplicated(PACA_AU_arraySeq[,c("gene_id")]),]
PACA_AU_arraySeq <- df
rownames(PACA_AU_arraySeq) <- PACA_AU_arraySeq$gene_id
PACA_AU_arraySeq <- PACA_AU_arraySeq[,-1]
PACA_AU_arraySeq <- as.data.frame(PACA_AU_arraySeq)
sample <- read.table("sample.txt",header = T)
sample <- subset(sample, Cohort == "PACA_AU_arraySeq")
PACA_AU_arraySeq <- PACA_AU_arraySeq[,sample$Sample]
df <- list()
df[["mRNA.expr"]] <- PACA_AU_arraySeq 
PACA_AU_arraySeq <- df
# run NTP in Yau cohort by using up-regulated biomarkers
PACA_AU_arraySeq.ntp.pred <- runNTP(expr  = PACA_AU_arraySeq$mRNA.expr,
                                    templates  = marker.up$templates, # the template has been already prepared in runMarker()
                                    scaleFlag  = TRUE, # scale input data (by default)
                                    centerFlag = TRUE, # center input data (by default)
                                    doPlot     = TRUE, # to generate heatmap
                                    fig.name   = "NTP HEATMAP FOR PACA_AU_arraySeq") 
write.csv(PACA_AU_arraySeq.ntp.pred[["ntp.res"]],"PACA_AU_arraySeqgroup.csv")
# 整合PACA_AU_RNA_seq
df <- PACA_AU_RNA_seq[!duplicated(PACA_AU_RNA_seq[,c("Tag")]),]
PACA_AU_RNA_seq <- df
rownames(PACA_AU_RNA_seq) <- PACA_AU_RNA_seq$Tag
PACA_AU_RNA_seq <- PACA_AU_RNA_seq[,-1]
PACA_AU_RNA_seq <- as.data.frame(PACA_AU_RNA_seq)
sample <- read.table("sample.txt",header = T)
table(sample$Cohort)
sample <- subset(sample, Cohort == "PACA_AU_RNA_seq")
PACA_AU_RNA_seq <- PACA_AU_RNA_seq[,sample$Sample]
df <- list()
df[["mRNA.expr"]] <- PACA_AU_RNA_seq 
PACA_AU_RNA_seq <- df
# run NTP in Yau cohort by using up-regulated biomarkers
PACA_AU_RNA_seq.ntp.pred <- runNTP(expr  = PACA_AU_RNA_seq$mRNA.expr,
                                   templates  = marker.up$templates, # the template has been already prepared in runMarker()
                                   scaleFlag  = TRUE, # scale input data (by default)
                                   centerFlag = TRUE, # center input data (by default)
                                   doPlot     = TRUE, # to generate heatmap
                                   fig.name   = "NTP HEATMAP FOR PACA_AU_RNA_seq") 
write.csv(PACA_AU_RNA_seq.ntp.pred[["ntp.res"]],"PACA_AU_RNA_seqgroup.csv")
# 整合PACA_CA_seq
df <- PACA_CA_seq[!duplicated(PACA_CA_seq[,c("gene_id")]),]
PACA_CA_seq <- df
rownames(PACA_CA_seq) <- PACA_CA_seq$gene_id
PACA_CA_seq <- PACA_CA_seq[,-1]
PACA_CA_seq <- as.data.frame(PACA_CA_seq)
setwd("/home/gejiachen/MOVICS")
sample <- read.table("sample.txt",header = T)
sample <- subset(sample, Cohort == "PACA_CA_seq")
PACA_CA_seq <- PACA_CA_seq[,sample$Sample]
df <- list()
df[["mRNA.expr"]] <- PACA_CA_seq 
PACA_CA_seq <- df
# run NTP in Yau cohort by using up-regulated biomarkers
PACA_CA_seq.ntp.pred <- runNTP(expr  = PACA_CA_seq$mRNA.expr,
                               templates  = marker.up$templates, # the template has been already prepared in runMarker()
                               scaleFlag  = TRUE, # scale input data (by default)
                               centerFlag = TRUE, # center input data (by default)
                               doPlot     = TRUE, # to generate heatmap
                               fig.name   = "NTP HEATMAP FOR PACA_CA_seq") 
write.csv(PACA_CA_seq.ntp.pred[["ntp.res"]],"PACA_CA_seqgroup.csv")