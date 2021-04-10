##1对多比较，得到差异基因
setwd("C:/Users/在干嘛/Desktop/数据分析/pheatmap/所有2个重复样本/差异基因")
filelist<- list.files() #建立工作文件目录
mylist<- list()##定义一个空list
file1<- read.csv(file = "SAM_V5A.csv",header = T,stringsAsFactors = F)
file1<- file1[,-1]#删除第一列的数字序列
for (i in filelist)
{
      if(i!=filelist[9]){
      file2<- read.csv(i,header = T,stringsAsFactors = F)
      file2<- file2[,-1]#删除第一列的数字序列
      data1<- merge(file1,file2,by="gene_id")
      rownames(data1)<- data1[,1]
      data1<- data1[,-1]
      data1<- round(data1)#取整数
      condition<- factor(c(rep("control",2),rep(i,2)))
      coldata<- data.frame(row.names = colnames(data1),condition)
      dds<- DESeqDataSetFromMatrix(countData = data1,colData = coldata,design = ~condition)
      dds<- DESeq(dds)
      res<- results(dds,contrast = c("condition","control",i))
      diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
      head(diff_gene_deseq2)
      name<- paste(sub(".csv","",i),"vs","control")##设置文件名字
      write.csv(diff_gene_deseq2,paste0(name,".csv"))#paste0,调用变量命名文件名
      
    }
}