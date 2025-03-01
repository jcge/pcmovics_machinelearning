#04 KM plot  survival analysis 
# 多组生存曲线绘制 
getwd()
setwd("/home/gejiachen/MOVICS/ML/Figures")
bmt <- read.csv("survivalplot.csv", header = T) 
# E_MATB_6134
E_MATB_6134 <- subset(bmt, Cohort == 'E_MATB_6134')
res.cut <- surv_cutpoint(E_MATB_6134,time = "OS.time",
                         event = "OS", variables = c("Riskscore"))
summary(res.cut) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS==1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="E_MATB_6134",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

# GSE21501
GSE21501 <- subset(bmt, Cohort == 'GSE21501')
res.cut <- surv_cutpoint(GSE21501,time = "OS.time",
                         event = "OS", variables = c("Riskscore"))
summary(res.cut) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS==1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="GSE21501",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

# GSE28735
GSE28735 <- subset(bmt, Cohort == 'GSE28735')
res.cut <- surv_cutpoint(GSE28735,time = "OS.time",
                         event = "OS", variables = c("Riskscore"))
summary(res.cut) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS==1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="GSE28735",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

# GSE57495
GSE57495 <- subset(bmt, Cohort == 'GSE57495')
res.cut <- surv_cutpoint(GSE57495,time = "OS.time",
                         event = "OS", variables = c("Riskscore"))
summary(res.cut) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS==1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="GSE57495",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

# GSE62452
GSE62452 <- subset(bmt, Cohort == 'GSE62452')
res.cut <- surv_cutpoint(GSE62452,time = "OS.time",
                         event = "OS", variables = c("Riskscore"))
summary(res.cut) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS==1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="GSE62452",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

# GSE71729
GSE71729 <- subset(bmt, Cohort == 'GSE71729')
res.cut <- surv_cutpoint(GSE71729,time = "OS.time",
                         event = "OS", variables = c("Riskscore"))
summary(res.cut) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS==1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="GSE71729",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

# GSE78229
GSE78229 <- subset(bmt, Cohort == 'GSE78229')
res.cut <- surv_cutpoint(GSE78229,time = "OS.time",
                         event = "OS", variables = c("Riskscore"))
summary(res.cut) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS==1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="GSE78229",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

# GSE79668
GSE79668 <- subset(bmt, Cohort == 'GSE79668')
res.cut <- surv_cutpoint(GSE79668,time = "OS.time",
                         event = "OS", variables = c("Riskscore"))
summary(res.cut) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS==1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="GSE79668",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

# GSE85916
GSE85916 <- subset(bmt, Cohort == 'GSE85916')
res.cut <- surv_cutpoint(GSE85916,time = "OS.time",
                         event = "OS", variables = c("Riskscore"))
summary(res.cut) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS==1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="GSE85916",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

# PACA_AU_arraySeq
PACA_AU_arraySeq <- subset(bmt, Cohort == 'PACA_AU_arraySeq')
res.cut <- surv_cutpoint(PACA_AU_arraySeq,time = "OS.time",
                         event = "OS", variables = c("Riskscore"))
summary(res.cut) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS==1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="PACA_AU_arraySeq",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

# PACA_AU_RNA_seq
PACA_AU_RNA_seq <- subset(bmt, Cohort == 'PACA_AU_RNA_seq')
res.cut <- surv_cutpoint(PACA_AU_RNA_seq,time = "OS.time",
                         event = "OS", variables = c("Riskscore"))
summary(res.cut) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS==1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="PACA_AU_RNA_seq",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

# PACA_CA_seq
PACA_CA_seq <- subset(bmt, Cohort == 'PACA_CA_seq')
res.cut <- surv_cutpoint(PACA_CA_seq,time = "OS.time",
                         event = "OS", variables = c("Riskscore"))
summary(res.cut) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS==1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="PACA_CA_seq",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

# TCGA_PAAD
TCGA_PAAD <- subset(bmt, Cohort == 'TCGA_PAAD')
res.cut <- surv_cutpoint(TCGA_PAAD,time = "OS.time",
                         event = "OS", variables = c("Riskscore"))
summary(res.cut) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS==1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="TCGA_PAAD",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)
