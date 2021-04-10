##多样本归一化处理
counts_data<- read.csv("C:/Users/在干嘛/Desktop/数据分析/pheatmap/EAR_All_data.csv",header=T,stringsAsFactors=F)
rownames(counts_data)<- counts_data[,1]
counts_data<- counts_data[,-1]
head(counts_data)
condition<- factor(c(rep("E4_V5C",2),rep("E4_V7",2),rep("E7_V5A",2),rep("E7_V5B",2),rep("E7_V5C",2),rep("E7_v6",2),rep("E7_V7",2)),levels=c("E4_V5C","E4_V7","E7_V5A","E7_V5B","E7_V5C","E7_v6","E7_V7"))#构建condition
##DEseq2均一化
coldata<- data.frame(row.names = colnames(counts_data), condition)
dds<- DESeqDataSetFromMatrix(countData = counts_data, colData = coldata, design = ~ condition)
dds<- DESeq(dds)
sizeFactors(dds)
head(dds)
res<- results(dds)
resdata<- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
head(resdata)
merge_list<- data.frame(counts_data,resdata)
head(merge_list)
resdata<- merge_list
head(resdata)
nor<- select(resdata,c(?:?))
write.csv(nor,file = "C:/Users/在干嘛/Desktop/数据分析/pheatmap/EAR_All_normalized.csv")