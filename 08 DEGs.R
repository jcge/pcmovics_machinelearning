#08 DEGs  
# 组间均有差异的基因筛选 
# tcga-paad
setwd("/home/gejiachen/MOVICS/zujianchayi")
# tcga-paad
TCGA_PAAD <- read.table("TCGA_PAAD.txt",sep = "\t",header = T,row.names = 1)
x=colnames(TCGA_PAAD)[1]
# 转换分类矩阵数据
data=melt(TCGA_PAAD,id.vars=c("group"))
colnames(data)=c("group","gene","relative expression")
p=ggviolin(data, x="gene", y="relative expression", color = "group", 
           ylab="relative expression",
           xlab="TCGA-PAAD cohort",
           legend.title=x,
           add.params = list(fill="white"),
           palette = c("#ECE700","#5bc0eb"),    
           width=1, add = "boxplot")
p=p+rotate_x_text(60)
p
p=p+stat_compare_means(aes(group=group),
                       method="t.test", 
                       # 可换其他统计方法
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                       label = "p.signif")
print(p)

# GSE28735
GSE28735 <- read.table("GSE28735.txt",sep = "\t",header = T,row.names = 1)
x=colnames(GSE28735)[1]
# 转换分类矩阵数据
data=melt(GSE28735,id.vars=c("group"))
colnames(data)=c("group","gene","relative expression")
p=ggviolin(data, x="gene", y="relative expression", color = "group", 
           ylab="relative expression",
           xlab="GSE28735 cohort",
           legend.title=x,
           add.params = list(fill="white"),
           palette = c("#ECE700","#5bc0eb"),    
           width=1, add = "boxplot")
p=p+rotate_x_text(60)
p
p=p+stat_compare_means(aes(group=group),
                       method="t.test", 
                       # 可换其他统计方法
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                       label = "p.signif")
print(p)

# GSE62452
GSE62452 <- read.table("GSE62452.txt",sep = "\t",header = T,row.names = 1)
x=colnames(GSE62452)[1]
# 转换分类矩阵数据
data=melt(GSE62452,id.vars=c("group"))
colnames(data)=c("group","gene","relative expression")
p=ggviolin(data, x="gene", y="relative expression", color = "group", 
           ylab="relative expression",
           xlab="GSE62452 cohort",
           legend.title=x,
           add.params = list(fill="white"),
           palette = c("#ECE700","#5bc0eb"),    
           width=1, add = "boxplot")
p=p+rotate_x_text(60)
p
p=p+stat_compare_means(aes(group=group),
                       method="t.test", 
                       # 可换其他统计方法
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                       label = "p.signif")
print(p)