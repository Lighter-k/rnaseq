#两两比较差异基因
library(DESeq2)

##设定工作目录
setwd("C:/Users/在干嘛/Desktop/数据分析/pheatmap/所有2个重复样本")
filelist<- list.files() #建立工作文件目录
mylist<- list()##定义一个空list

for (i in filelist)
{
  for (j in filelist) {
    if(i!=j){
      file1<- read.csv(i,header = T,stringsAsFactors = F)
      file2<- read.csv(j,header = T,stringsAsFactors = F)
      data1<- merge(file1,file2,by="gene_id")
      rownames(data1)<- data1[,1]
      data1<- data1[,-1]
      condition<- factor(c(rep(i,2),rep(j,2)))
      coldata<- data.frame(row.names = colnames(data1),condition)
      data1<- round(data1)
      dds<- DESeqDataSetFromMatrix(countData = data1,colData = coldata,design = ~condition)
      dds<- DESeq(dds)
      res<- results(dds,contrast = c("condition",i,j))
      diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
      data2<- rownames(diff_gene_deseq2)
      data2<- list(data2)
      name<- colnames(data1)
      names(data2)<- paste(sub(".csv","",i),"vs",sub(".csv","",j))#sub字符串替换，将i中的".csv"替换成"".paste连接不同类型的变量及常量
      mylist<- c(data2,mylist)
    }
  }
}

##将所有的list合并成1列
all_DEG<- mylist[[1]]
for (l in 2:length(mylist)) {
  all_DEG<- c(all_DEG,mylist[[l]])
  
}
all_DEG<- as.data.frame(all_DEG)

##去除重复的gene_id
aa<- all_DEG[complete.cases(all_DEG),]##complete.cases函数,先去NA
aa<- as.data.frame(aa)
bb<- aa[!duplicated(aa$aa),]
bb<- as.data.frame(bb)

write.csv(bb,file = "all_DEG.csv")#保存到当前工作目录下


##取出所有样本差异表达gene的counts
library(dplyr)
all_DEG_counts<- semi_join(所有基因的counts,差异表达的gene,by='gene_id')


##pheatmap绘制热图
library(pheatmap)
 p1<- pheatmap::pheatmap(average_counts,show_rownames = F,scale = "row",fontsize = 12)
 A_cluster <- average_counts[p1$tree_row$order,]



