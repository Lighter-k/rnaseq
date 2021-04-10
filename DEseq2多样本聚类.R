差异基因提取、合并
##读取control
sub_data<- read.csv(file="C:/Users/在干嘛/Desktop/数据分析/pheatmap/DEseq2归一化处理后/EAR_All_normalized_control.csv",header=T,stringsAsFactors=F)
All_data<- read.csv(file="C:/Users/在干嘛/Desktop/数据分析/pheatmap/DEseq2归一化处理后/EAR_All_normalized_input.csv",header=T,stringsAsFactors=F)
Geneid<- All_data$gene_id
Geneid<- data.frame(Geneid)
colnames(Geneid)<- c("gene_id")
colist<- colnames(All_data) ##储存进行差异比对的列名
mylist<- list()##定义一个空list，后面用

for (i in 2:length(colist))
{
  ##构建差异组
  data1<- data.frame(All_data[,1],All_data[,i])
  colnames(data1)<- c("gene_id",colnames(All_data[i]))
  merge.count<- merge(sub_data,data1,by="gene_id")
  data2<- merge.count[,-1]
  data2<- as.matrix(data2)
  rownames(data2)<- merge.count$gene_id
  
  ##差异分析
  condition<- factor(c(rep("control",2),rep(colnames(All_data[i]),2)))
  coldata<- data.frame(row.names = colnames(data2),condition)
  data2<- round(data2) ##取整数
  dds<- DESeqDataSetFromMatrix(countData = data2,colData = coldata,design = ~condition)
  dds<- DESeq(dds)
  
  ##差异基因提取
  res<- results(dds,contrast = c("condition","control",colnames(All_data[i])))
  diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)##获取差异表达的基因
  data3<- rownames(diff_gene_deseq2)##将差异表达的geneid储存在data3
  data3<- list(data3)
  names(data3)<- c(colnames(All_data[i]))##修改list的名
  mylist<- c(data3,mylist)##将list一行一行结合
  ##data3<- data.frame(data3)
  ##colnames(data3)<- c(colnames(All_data[i]))

  ##
  ##write.csv(data3,paste0(control vs colnames(All_data[i]),".csv"))##paste函数
  
  ##差异基因合并> Geneid<- All_data$gene_id
  ##Geneid<- merge(Geneid,data3,by=rownames(data3))
}

##把所有的差异gene整合到一列
#aa<- colnames(All_data)
#aa<- aa[18:27]
#for (i in 1:16) {bb<- c(bb,mylist[[aa[i]]])}
#bb<- as.data.frame(bb)

mylist<- data.frame(sapply(mylist, "[", i = 1:max(sapply(mylist, length))))##把list格式的结果改为数据框格式
##参考https://www.zhihu.com/question/357969098/answer/911931396
write.csv(mylist,file="C:/Users/在干嘛/Desktop/数据分析/pheatmap/DEseq2归一化处理后/EAR_DEG.csv")
