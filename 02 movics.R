#02 MOVICS Figure 1 and Figure 2
library(MOVICS)
library(purrr)
optk.paad <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,F,T), #二元数据写T
                         try.N.clust = 2:8, 
                         fig.name    = "CLUSTER NUMBER OF TCGA-PAAD")
optk.paad[["N.clust"]]
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster","iClusterBayes"),
                         N.clust     = 2,
                         type        = c("gaussian","gaussian", "gaussian", "gaussian", "binomial"))
save(moic.res.list, file = "moic.res.list.rda") # 保存下结果
# 整合多种分型结果
cmoic.paad <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")
clust.res <- cmoic.paad[["clust.res"]]
write.csv(clust.res,"clust.res.csv")
# 使用Silhouette准则判断分型质量
getSilhouette(sil      = cmoic.paad$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)

# 多组学分型热图 
# β值矩阵转换为M值矩阵
indata <- mo.data
indata$meth.beta <- log2(indata$meth.beta / (1 - indata$meth.beta))
# 对数据进行标准化
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,T,F)) # no scale for mutation
#标注前10个分子
#根据贝叶斯方法进行分型
iClusterBayes.res <- getMOIC(data        = mo.data,
                             methodslist = "iClusterBayes",
                             N.clust     = 4,
                             type        = c("gaussian","gaussian", "gaussian", "gaussian", "binomial"))
feat   <- iClusterBayes.res$feat.res 
feat1  <- feat[which(feat$dataset == "mRNA.expr"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "miRNA.expr"),][1:10,"feature"] 
feat3  <- feat[which(feat$dataset == "lncRNA.expr"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "meth.beta"),][1:10,"feature"]
feat5  <- feat[which(feat$dataset == "mut.status"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4, feat5)
# 为每个组学的热图自定义颜色，不定义也可
mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
miRNA.col  <- c("#74ADD1", "#FFFFBF", "#FDAE61", "#D73027")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col,miRNA.col,lncRNA.col, meth.col, mut.col)
# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","miRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","miRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = cmoic.paad$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")

#分型热图
#总subtype
Subtype <- cmoic.paad$clust.res
# Subtype <- Subtype$clust
#SNF
SNF <- moic.res.list[["SNF"]]$clust.res
Subtype$SNF <- SNF$clust
#PINSPlus
PINSPlus <- moic.res.list[["PINSPlus"]]$clust.res
Subtype$PINSPlus <- PINSPlus$clust
#NEMO
NEMO <- moic.res.list[["NEMO"]]$clust.res
Subtype$NEMO <- NEMO$clust
#COCA
COCA <- moic.res.list[["COCA"]]$clust.res
Subtype$COCA <- COCA$clust
#LRAcluster
LRAcluster <- moic.res.list[["LRAcluster"]]$clust.res
Subtype$LRAcluster <- LRAcluster$clust
#ConsensusClustering
ConsensusClustering <- moic.res.list[["ConsensusClustering"]]$clust.res
Subtype$ConsensusClustering <- ConsensusClustering$clust
#IntNMF
IntNMF <- moic.res.list[["IntNMF"]]$clust.res
Subtype$IntNMF <- IntNMF$clust
#CIMLR
CIMLR <- moic.res.list[["CIMLR"]]$clust.res
Subtype$CIMLR <- CIMLR$clust
#MoCluster
MoCluster <- moic.res.list[["MoCluster"]]$clust.res
Subtype$MoCluster <- MoCluster$clust
#iClusterBayes
iClusterBayes <- moic.res.list[["iClusterBayes"]]$clust.res
Subtype$iClusterBayes <- iClusterBayes$clust
saveRDS(Subtype,"subtype.rds")
write.csv(Subtype,"subtype.csv")

#subtype热图 
subtype <- read.csv("subtype.csv",sep = ",",row.names = 1)
# subtype <- subtype[,-1]
subtype <- subtype[order(subtype$clust),]
clin <- subtype[,1:10]
exp <- subtype[,11:21]
#计算显著性
table(subtype$clust)
#画图
pre_heatdata <- exp
annColors <- list()
annColors[['Sex']] <- c('Male'='steelblue','Female'='red')
annColors[['Tstage']] <- c('T1'='cornsilk','T2'='paleturquoise','T3'='goldenrod','T4'='firebrick','TX'='white')
annColors[['Nstage']] <- c('N0'='slategray','N1'='skyblue','NX'='white')
annColors[['Mstage']] <- c('M0'='pink','M1'='cornsilk','MX'='white')
# annColors[['AJCC']] <- c('I'='yellowgreen','II'='lightseagreen','III'='blueviolet','IV'='red','Unknown'='white')
annColors[['Grade']] <- c('G1'='lightseagreen','G2'='plum','G3'='purple','G4'='black','GX'='white')
library(pheatmap)
#pdf(file="pheatmap.pdf",width = 10,height = 12)
pre_heatdata <- t(pre_heatdata)
p <- pheatmap(pre_heatdata,
              
              color = colorRampPalette(c("#5bc0eb",'black',"#ECE700"))(1000),
              
              annotation_col = clin,
              
              annotation_colors = annColors,
              
              treeheight_row = 50,
              
              # gaps_col = 92,
              
              show_rownames = T,
              
              show_colnames = F,
              
              cluster_rows = F,
              
              cluster_cols = F)
p

# 生存差异 
surv.paad <- compSurv(moic.res         = cmoic.paad,
                      surv.info        = surv.info,
                      convt.time       = "d", # 把天变成月
                      surv.median.line = "h", 
                      xyrs.est         = c(1,3,5), # 计算5年和10年生存率
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")
# ?compSurv
print(surv.paad) 
# 比较临床特征 
clin.paad <- compClinvar(moic.res      = cmoic.paad,
                         var2comp      = surv.info, #需要比较的临床信息，行名须是样本名
                         strata        = "Subtype", # 分层变量，这里肯定是分型了
                         factorVars    = NULL, #分类变量名字
                         nonnormalVars = NULL, #非正态的连续性变量
                         exactVars     = NULL, #需要使用精确概率法的变量
                         doWord        = TRUE, # 自动生成Word文档
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES")
?compClinvar
# 突变全景图
mut.brca <- compMut(moic.res     = cmoic.paad,
                    mut.matrix   = paad.tcga$mut.status, # 0/1矩阵
                    doWord       = TRUE, # 生成Word文档
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.05, # 保留在至少5%的样本中突变的基因
                    p.adj.cutoff = 0.05, # 保留padj<0.05的基因
                    innerclust   = TRUE, # 在每个亚型中进行聚类
                    # annCol       = annCol, # same annotation for heatmap
                    # annColors    = annColors, # same annotation color for heatmap
                    width        = 6, 
                    height       = 2,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")
table(mut.matrix$`TCGA-IB-A6UG-01A`)
mut.matrix <- mut.matrix[1:10,]

# DEGs
# run DEA with edgeR
runDEA(dea.method = "edger",
       expr       = count, # raw count data
       moic.res   = cmoic.paad,
       prefix     = "TCGA-PAAD") # prefix of figure name
# run DEA with DESeq2
runDEA(dea.method = "deseq2",
       expr       = count, # deseq2也需要count
       moic.res   = cmoic.paad,
       prefix     = "TCGA-PAAD")
# run DEA with limma
runDEA(dea.method = "limma",
       expr       = count, # normalized expression data
       moic.res   = cmoic.paad,
       prefix     = "TCGA-PAAD")

# 识别特定分子
# choose edgeR result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = cmoic.paad,
                       dea.method    = "edger", # name of DEA method
                       prefix        = "TCGA-PAAD", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 5000, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = count, # use normalized expression as heatmap input
                       #annCol        = annCol, # sample annotation in heatmap
                       #annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")
marker.down<- runMarker(moic.res      = cmoic.paad,
                        dea.method    = "edger", # name of DEA method
                        prefix        = "TCGA-PAAD", # MUST be the same of argument in runDEA()
                        dat.path      = getwd(), # path of DEA files
                        res.path      = getwd(), # path to save marker files
                        p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                        p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                        dirct         = "down", # direction of dysregulation in expression
                        n.marker      = 5000, # number of biomarkers for each subtype
                        doplot        = TRUE, # generate diagonal heatmap
                        norm.expr     = count, # use normalized expression as heatmap input
                        #annCol        = annCol, # sample annotation in heatmap
                        #annColors     = annColors, # colors for sample annotation
                        show_rownames = FALSE, # show no rownames (biomarker name)
                        fig.name      = "DOWNREGULATED BIOMARKER HEATMAP")

# GSEA
# MUST locate ABSOLUTE path of msigdb file
MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)
# run GSEA to identify up-regulated GO pathways using results from edgeR
gsea.up <- runGSEA(moic.res     = cmoic.paad,
                   dea.method   = "edger", # name of DEA method
                   prefix       = "TCGA-PAAD", # MUST be the same of argument in runDEA()
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = count, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "UPREGULATED PATHWAY HEATMAP")
print(gsea.up$gsea.list$CS1[1:6,3:6])
head(round(gsea.up$grouped.es,3))
gsea.down <- runGSEA(moic.res     = cmoic.paad,
                     dea.method   = "edger", # name of DEA method
                     prefix       = "TCGA-PAAD", # MUST be the same of argument in runDEA()
                     dat.path     = getwd(), # path of DEA files
                     res.path     = getwd(), # path to save GSEA files
                     msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                     norm.expr    = count, # use normalized expression to calculate enrichment score
                     dirct        = "down", # direction of dysregulation in pathway
                     p.cutoff     = 0.05, # p cutoff to identify significant pathways
                     p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                     gsva.method  = "gsva", # method to calculate single sample enrichment score
                     norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                     fig.name     = "DOWNREGULATED PATHWAY HEATMAP")

# GSVA
# MUST locate ABSOLUTE path of gene set file
GSET.FILE <- system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)
# run GSVA to estimate single sample enrichment score based on given gene set of interest
gsva.res <- runGSVA(moic.res      = cmoic.paad,
                    norm.expr     = count,
                    gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
                    gsva.method   = "gsva", # method to calculate single sample enrichment score
                    # annCol        = annCol,
                    # annColors     = annColors,
                    fig.path      = getwd(),
                    fig.name      = "GENE SETS OF INTEREST HEATMAP",
                    height        = 5,
                    width         = 8)
print(gsva.res$raw.es[1:3,1:3])