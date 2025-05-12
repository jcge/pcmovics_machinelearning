#06 forest plot  Figure 4C and Supplementary Figure 2
library(randomcoloR)
palette <- distinctColorPalette(13)  #差异明显的60种
palette
# 森林图绘制P值和HR 
library(forestplot)
setwd("/home/gejiachen/MOVICS/senlin") 
# GSE71729
rs_forest <- read.table('GSE71729.txt',header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,c(1:2,5)]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V2, #设置均值
           
           lower = rs_forest$V3, #设置均值的lowlimits限
           
           upper = rs_forest$V4, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1.0,1.5),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,1.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box="#D3E0CB",summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
# GSE21501
rs_forest <- read.table("GSE21501.txt",header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,c(1:2,5)]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V2, #设置均值
           
           lower = rs_forest$V3, #设置均值的lowlimits限
           
           upper = rs_forest$V4, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1.0,1.5),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,1.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box=palette[2],summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
# GSE28735
rs_forest <- read.table("GSE28735.txt",header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,c(1:2,5)]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V2, #设置均值
           
           lower = rs_forest$V3, #设置均值的lowlimits限
           
           upper = rs_forest$V4, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1.0,1.5),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,1.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box=palette[3],summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
# GSE78229
rs_forest <- read.table("GSE78229.txt",header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,c(1:2,5)]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V2, #设置均值
           
           lower = rs_forest$V3, #设置均值的lowlimits限
           
           upper = rs_forest$V4, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1.0,1.5),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,1.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box=palette[4],summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
# GSE85916
rs_forest <- read.table("GSE85916.txt",header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,c(1:2,5)]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V2, #设置均值
           
           lower = rs_forest$V3, #设置均值的lowlimits限
           
           upper = rs_forest$V4, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1.0,1.5),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,1.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box=palette[5],summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
# PACA_AU_arraySeq
rs_forest <- read.table("PACA_AU_arraySeq.txt",header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,c(1:2,5)]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V2, #设置均值
           
           lower = rs_forest$V3, #设置均值的lowlimits限
           
           upper = rs_forest$V4, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1.0,1.5),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,1.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box=palette[6],summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
# TCGA_PAAD
rs_forest <- read.table("TCGA_PAAD.txt",header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,c(1:2,5)]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V2, #设置均值
           
           lower = rs_forest$V3, #设置均值的lowlimits限
           
           upper = rs_forest$V4, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1.0,1.5),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,1.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box=palette[7],summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
# EMTAB6134
rs_forest <- read.table("EMTAB6134.txt",header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,c(1:2,5)]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V2, #设置均值
           
           lower = rs_forest$V3, #设置均值的lowlimits限
           
           upper = rs_forest$V4, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1.0,1.5),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,1.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box=palette[8],summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
# GSE79668
rs_forest <- read.table("GSE79668.txt",header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,c(1:2,5)]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V2, #设置均值
           
           lower = rs_forest$V3, #设置均值的lowlimits限
           
           upper = rs_forest$V4, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1.0,1.5),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,1.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box=palette[9],summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
# PACA_AU_RNA_seq
rs_forest <- read.table("PACA_AU_RNA_seq.txt",header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,c(1:2,5)]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V2, #设置均值
           
           lower = rs_forest$V3, #设置均值的lowlimits限
           
           upper = rs_forest$V4, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1.0,1.5),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,1.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box=palette[10],summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
# PACA_CA_seq
rs_forest <- read.table("PACA_CA_seq.txt",header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,c(1:2,5)]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V2, #设置均值
           
           lower = rs_forest$V3, #设置均值的lowlimits限
           
           upper = rs_forest$V4, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1.0,1.5),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,1.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box=palette[11],summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
# GSE57495
rs_forest <- read.table("GSE57495.txt",header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,c(1:2,5)]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V2, #设置均值
           
           lower = rs_forest$V3, #设置均值的lowlimits限
           
           upper = rs_forest$V4, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1.0,1.5),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,1.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box=palette[12],summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
# GSE62452
rs_forest <- read.table("GSE62452.txt",header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,c(1:2,5)]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V2, #设置均值
           
           lower = rs_forest$V3, #设置均值的lowlimits限
           
           upper = rs_forest$V4, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1.0,1.5),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,1.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box=palette[13],summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 4)#设置森林图的位置，此处设置为4，则出现在第四列
