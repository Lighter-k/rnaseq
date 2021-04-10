##样本数据均一化处理
##读取数据
counts_data<- read.csv("C:/Users/在干嘛/Desktop/数据分析/pheatmap/All_data.csv",header=T,stringsAsFactors=F)
rownames(counts_data)<- counts_data[,1]
counts_data<- counts_data[,-1]
##构建condition
##colist<- colnames(counts_data)
condition<- factor(colnames(counts_data))
dds <- DESeqDataSetFromMatrix(countData=counts_data,DataFrame(condition),design=~condition)
rld <- rlogTransformation(dds) ##归一化处理
exprSet_new=assay(rld) ##提取DEseq标准化后的数据
write.csv(exprSet_new,file = "C:/Users/在干嘛/Desktop/数据分析/pheatmap/归一化All_data.csv")
exprSet_new<- exprSet_new[rowSums(exprSet_new)!=0,]##去除在所有样本中都不表达的行
write.csv(exprSet_new,file = "C:/Users/在干嘛/Desktop/数据分析/pheatmap/归一化All_data_去除不表达.csv")
plotPCA(rld, intgroup=c('condition'))
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)##绘图
plot(pcaData[,1:2],pch=15,col=c(rep("red",15),rep("blue",15)),font=2,font.lab=2)+text(pcaData[,1],pcaData[,2]+1.5,row.names(pcaData),cex=1,font = 2)+legend("right",inset = 0.04,title = "Rep",c("Rep1","Rep2"),pch = 15,col = c("red","blue"))
