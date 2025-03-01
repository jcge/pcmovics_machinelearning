#07 distinct characteristics of risk groups  Figure 5A-F  
# gsea
library(edgeR) 
sample <- subtype[,1]
count <- count[,sample]
# count <- cbind(rownames(count),count)
# colnames(count)[1] <- "gene"
write.table(count,"count.txt",sep = "\t",row.names = F)
count <- count[,-1]
library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例

# degs
group <- read.table("group.txt",header = T)
count <- count[,group$sample]
data_int <- 2^(count) - 1
data_int <- apply(data_int, 2, as.integer)
rownames(data_int) <- rownames(count)
exprSet <- data_int
# 创建设计矩阵，指定组别信息
group <- group[,2]
group <- factor(group)
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(exprSet)
# 创建 DGEList 对象
dge <- DGEList(counts = exprSet, group = group)
# 这里我们使用上面提到的 filterByExpr() 进行自动过滤，去除低表达基因
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
# 归一化，得到的归一化系数被用作文库大小的缩放系数
dge <- calcNormFactors(dge)
# 使用 voom 方法进行标准化
v <- voom(dge, design, plot = TRUE, normalize = "quantile")
# 如果是芯片数据、TPM数据或已标准化的数据，不需要再进行标准化，可直接从这里开始进行差异分析
# 使用线性模型进行拟合
fit <- lmFit(v, design)
# 和上面两个包一样，需要说明是谁比谁
con <- paste(rev(levels(group)), collapse = "-")
con
# 创建对比矩阵
cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
# 获取差异表达基因结果
tempOutput <- topTable(fit2, coef = con, n = Inf)
DEG_limma_voom <- na.omit(tempOutput)
head(DEG_limma_voom)
write.table(DEG_limma_voom,file = "DEG_limma_voom.txt")
info <- read.table("DEG_limma_voom.txt",row.names = 1)
#gene ID转换
GO_database <- 'org.Hs.eg.db'
gene <- bitr(rownames(info),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
colnames(info) 
info <- cbind(rownames(info),info)
names(info) <- c('SYMBOL','Log2FoldChange',"AveExpr","t",'pvalue','padj','B')
rownames(info) <- NULL
info_merge <- merge(info,gene,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange
GSEA_input <- info_merge$Log2FoldChange
names(GSEA_input) = info_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html
GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 0.05)#GSEA富集分析
GSEA_KEGG@result <- GSEA_KEGG@result[order(GSEA_KEGG@result[["enrichmentScore"]]),]
GSEA_KEGG@result <- GSEA_KEGG@result[order(-GSEA_KEGG@result[["enrichmentScore"]]),]
gseaplot2(GSEA_KEGG,1:5)#30是根据ridgeplot中有30个富集通路得到的

# 药物敏感性
install.packages("/home/gejiachen/rpackage/resource/pRRophetic_Guozi", repos = NULL, type = "source")
library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)
allDrugs=c("AICAR", "AKT.inhibitor.VIII", "ATRA", "Axitinib", "Bexarotene", "Bicalutamide","Bleomycin",
           "Bortezomib", "Bosutinib","Camptothecin", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin",  "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", "Etoposide","Gefitinib", "Gemcitabine", "GSK269962A", "Imatinib", "JNK.Inhibitor.VIII", "Lapatinib", "Lenalidomide", "Metformin", "Methotrexate",  "Midostaurin", "Mitomycin.C", "Nilotinib", "Obatoclax.Mesylate", "Paclitaxel", "Parthenolide", "Pazopanib", "Pyrimethamine",  "Rapamycin", "Roscovitine", "Salubrinal", "Shikonin",  "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "Vinblastine", "Vinorelbine", "Vorinostat" )
exp=as.matrix(count)
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]
drug <- NULL
i = 3
for(drug in allDrugs){
  #预测药物敏感性
  senstivity=pRRopheticPredict(data, drug, selection=1)
  senstivity=senstivity[senstivity!="NaN"]
  senstivity=as.data.frame(senstivity)
  group[,i] = senstivity
  colnames(group)[i] <- drug
  i = i + 1
  #senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
}
write.csv(group,"drug.csv")
a<-boxplot[["plot_env"]][["data"]]
print(a)

# 岭回归
cv.fit = cv.glmnet(x = Train_set,
                   y = Surv(Train_surv[,"OS.time"], Train_surv[,"OS"]),
                   family = "cox", alpha = 0, nfolds = 10)
plot(cv.fit)
fit = glmnet(x = Train_set,
             y = Surv(Train_surv[,"OS.time"], Train_surv[,"OS"]),
             family = "cox", alpha = 0, lambda = cv.fit$lambda.min)
plot(cv.fit$glmnet.fit,"lambda",label = T)

ridge <- model[["Ridge"]]
ridge <- cbind(ridge[["beta"]]@Dimnames[[1]],ridge[["beta"]]@x)
ridge <- as.matrix(ridge)

# 药物敏感性 
drug <- read.csv("drug.csv",sep = ",",header = T,row.names = 1)
x=colnames(drug)[1]
# 转换分类矩阵数据
data=melt(drug,id.vars=c("group"))
colnames(data)=c("group","drug","sensitivity")
p=ggviolin(data, x="drug", y="sensitivity", color = "group", 
           ylab="Drug Sensitivity",
           xlab="Drug",
           legend.title=x,
           add.params = list(fill="white"),
           palette = c("#ECE700","#5bc0eb"),    
           width=1, add = "boxplot")
p=p+rotate_x_text(60)
p
p=p+stat_compare_means(aes(group=group),
                       method="wilcox.test", 
                       # 可换其他统计方法
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                       label = "p.signif")
print(p)

# 高低分组对比
setwd("/home/gejiachen/MOVICS/highorlow")
# grade
library("ggpubr")
grade <- read.table("grade.txt",header = T,sep = "\t")
mycom <- list(c("G1","G2"),c("G2","G3or4"))
ggviolin(grade, x = "Grade",y = "Riskscore",color = "Grade", size = 1.2,
         palette = c("#B0E854","#DEA7CC","#79ACD5","#D99264"),
         add = "boxplot") +
  stat_compare_means(comparisons = mycom, label.y = c(1.8,2),
                     label = "p.signif",size = 6)+
  stat_compare_means(label.y = 3)+
  ylim(-1.5,3) 
# margin
margin <- read.table("margin.txt",header = T,sep = "\t")
mycom <- list(c("R0","R1"))
ggviolin(margin, x = "margin",y = "Riskscore",color = "margin", size = 1.2,
         palette = c("#B0E854","#DEA7CC","#79ACD5","#D99264"),
         add = "boxplot") +
  stat_compare_means(comparisons = mycom, label.y = c(1.8),
                     label = "p.signif",size = 6)+
  stat_compare_means(label.y = 2.5)+
  ylim(-1.5,2.5) 
# tumor_cell_classification
tumor_cell_classification <- read.table("tumor_cell_classification.txt",header = T,sep = "\t")
mycom <- list(c("basal-like","classical"))
ggviolin(tumor_cell_classification, x = "tumor_cell_classification",y = "Riskscore",color = "tumor_cell_classification", size = 1.2,
         palette = c("#B0E854","#DEA7CC","#79ACD5","#D99264"),
         add = "boxplot") +
  stat_compare_means(comparisons = mycom, label.y = c(2),
                     label = "p.signif",size = 6)+
  stat_compare_means(label.y = 2.5)+
  ylim(-1,2.5) 
# tmb突变负荷 
#突变组
library(maftools) 
library(TCGAmutations)
a <- TCGAmutations::tcga_available() 
paad.maf <- TCGAmutations::tcga_load(study = "PAAD")
paad <- paad.maf@data 
# write.table(paad,"paad.txt",sep = "\t") 
total = read.maf(maf = paad,
                 vc_nonSyn=names(tail(sort(table(paad$Variant_Classification)))))
plotmafSummary(maf = total,
               rmOutlier = TRUE,
               addStat = "median",
               dashboard = TRUE,
               titvRaw = FALSE)
##瀑布图(Oncoplots)：
oncoplot(maf = total,
         #draw_titv = FALSE,
         #anno_height = 1,#样本注释区高度
         #legend_height = 4, #图例绘图区高度
         #drawRowBar = T, #是否显示右侧条形图
         #drawColBar = T, #是否显示顶部条形图
         top = 10) #高频突变的Top10基因
#titv函数将SNP分类为Transitions_vs_Transversions，
#并以各种方式返回汇总表的列表。汇总数据也可以显示为一个箱线图，
#显示六种不同转换的总体分布，并作为堆积条形图显示每个样本中的转换比例
laml.titv = titv(maf = total, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)
#rainfall plots，展示超突变的基因组区域。
#detectChangePoints设置为TRUE，rainfall plots可以突出显示潜在变化的区域
#rainfallPlot(maf = total, detectChangePoints = TRUE, pointSize = 0.6)
#somaticInteractions函数使用配对Fisher 's精确检验来分析突变基因之间
#的的co-occurring 或者exclusiveness
#exclusive/co-occurance event analysis on top 10 mutated genes.
Interact <- somaticInteractions(maf = total, top = 15, pvalue = c(0.05, 0.1))
#提取P值结果
Interact$pair
#mafComapre参数比较两个不同队列的差异突变基因，检验方式为fisher检验。
tmb = tmb(maf = total)   #默认以log10转化的TMB绘图
write.csv(tmb,"tmb.csv")
