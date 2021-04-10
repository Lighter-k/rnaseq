##时间序列分析
library(maSigPro)
data.abiotic<- read.csv(file = "time_analysis_data.csv",header = T,stringsAsFactors = F) #加载样本COUNTS数据，必须使用归一化之后的表达量
edesign.abiotic<- read.csv(file = "sample_singal.csv",header = T,stringsAsFactors = F) #加载分组及相关信息

##行名变换处理
rownames(data.abiotic)<- data.abiotic[,2]
data.abiotic<- data.abiotic[,-(1:2)]
rownames(edesign.abiotic)<- edesign.abiotic[,1]
edesign.abiotic<- edesign.abiotic[,-1]

##定义回归模型
design <- make.design.matrix(edesign.abiotic, degree = 2)
design$groups.vector#查看回归参数的分配

##差异基因分析
fit <- p.vector(data.abiotic, design, Q = 1, MT.adjust = "BH", min.obs = 6) ##min.obs低于这个数值的真实基因将被排除在分析之外，6最小；若太大会报错(下标出界),输入的
是差异基因时，可以将Q调整为1
#返回参数
#fit$i # returns the number of significant genes
#fit$alfa # gives p-value at the Q false discovery control level
#fit$SELEC # is a matrix with the significant genes and their expression values

##差异分析
tstep <- T.fit(fit, step.method = "forward", alfa = 0.05) ##当只指定ear和sam时，step.method可以选取"forward"参数，或了解其他的

##提取差异基因列表
sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")

##查看差异基因情况
sigs$summary
suma2Venn(sigs$summary[, c(2:3)])
suma2Venn(sigs$summary[, c(1:3)])
sigs$sig.genes$SAMvsControl$g #查看差异基因结果

##将表达模式相似的基因进行聚类可视化其中k参数为cluster的数量，cluster.method为聚类的方式（hclust，kmeans，Mclust)
sam<- see.genes(sigs$sig.genes$SAMvsControl, show.fit = T, dis =design$dis, cluster.method="hclust" ,cluster.data = 1, k = 9)#当只指定sam和ear时，cluster.method改为kmeans
sam_control<- sam$cut ##提取出差异的基因
sam_control<- as.data.frame(sam_control)

#参考https://zhuanlan.zhihu.com/p/104357304


##根据不同的cluster导出基因
for (i in 1:9) {
sam_control_cluster1<- subset(sam_control,sam_control==i) ##提取出cluster
name<- paste("sam_control_cluster",sub("","",i)) ##
write.csv(sam_control_cluster1,paste0(name,".csv"))
}

##使用PlotGroups()函数可以查看某一特定基因的表达情况
Zm00001d011605<- data.abiotic[rownames(data.abiotic)=="Zm00001d011605", ]
PlotGroups (Zm00001d011605, edesign = edesign.abiotic)