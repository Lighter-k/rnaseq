# RNAseq批量处理

### Trimmomatic过滤

ls *.fastq.gz|while read id; do java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 12 -phred33 ${id%.R1.fastq.gz}.R1.fastq.gz ${id%.R1.fastq.gz}.R2.fastq.gz ${id%.R1.fastq.gz}.R1.paired.fastq.gz ${id%.R1.fastq.gz}.R1.unpaired.fastq.gz ${id%.R1.fastq.gz}.R2.paired.fastq.gz ${id%.R1.fastq.gz}.R2.unpaired.fastq.gz ILLUMINACLIP:/opt/Trimmomatic-0.36/adapters/Smart.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25; done

##{id%.R1.fastq.gz}指取id字符串从.R1.fastq.gz开始的左边字符串，提取指定位置字符串参考 https://blog.csdn.net/qq_26442553/article/details/79916527

##<u>/opt/Trimmomatic-0.36/adapters/Smart.fa</u>指定接头文件的位置

##根据测序的接头和数据质量设定接头文件和最小丢弃序列长度

##过滤后数据会变小，因为过滤的无用数据

### Fastqc质检

ls *gz | while read id; do fastqc -t 4 $id; done

### multiqc整合

multiqc . #扫描当前文件夹

### hista2比对

ls *R1.paired.fastq.gz|while read id; do /opt/hisat2/hisat2 -p 6 -x /lustre2/gmzheng/maize_data/hisat2-index/Zea_mays -1 /lustre2/gmzheng/epi-50/clean_data/${id%.R1.paired.fastq.gz}.R1.paired.fastq.gz -2 /lustre2/gmzheng/epi-50/clean_data/${id%.R1.paired.fastq.gz}.R2.paired.fastq.gz -S /lustre2/gmzheng/epi-50/alig_data/${id%.R1.paired.fastq.gz}.sam; done

建议同时生成比对的日志文件，查看比对率。

ls *.sam|while read id; do samtools flagstat -@ 20 $id >$id.log; done ##未生成日志文件，使用flagstat

### Samtools格式转化

#### sam文件转bam文件

ls *.sam|while read id; do samtools view -S $id -b > ${id%.sam}.bam; done

#### bam文件sort

ls *.bam|while read id; do samtools sort $id -o ${id%.bam}_sort.bam ; done 

#### sort后index

ls *sort.bam|while read id; do samtools index $id; done ##index之后文件变小，大约有1-2M

### htseq-count定量

### featureCounts定量

以sort.bam为输入文件

```
~/biosoft/subread/subread-1.6.0-Linux-x86_64/bin/featureCounts -T 6 -p -t exon -g gene_id -a ~/annotation/mm10/gencode.vM13.annotation.gtf -o SRR3589959_featureCounts222.txt SRR3589959.bam
-a 输入GTF/GFF基因组注释文件
-p 这个参数是针对paired-end数据
-F 指定-a注释文件的格式，默认是GTF
-g 从注释文件中提取Meta-features信息用于read count，默认是gene_id
-t 跟-g一样的意思，其是默认将exon作为一个feature
-o 输出文件
-T 多线程数

ls *bam|while read id; do featureCounts -T 8 -p -t exon -g gene_id -a /lustre2/gmzheng/maize_data/Zea_mays.B73_RefGen_v4.48.gtf -o ${id%_sort.bam}_featureCounts.txt $id; done
multiqc *.txt.summary ##对featureCounts结果整合及可视化
```

ftp://ftp.gramene.org/pub/gramene/release-63/gtf/ ##gtf文件下载

只保留文件的gene_id和样本counts

ls *featureCounts.txt|while read id; do cat $id| cut -f1,7- > ${id%featureCounts.txt}counts.txt; done

### PCA

R中操作

setwd() ##设置工作路径

list.files() ##查看文件

v3<- read.table("epi50-1_L3_Y09D74_counts.txt",sep = "\t",col.names = c("gene_id","v3")) ##读入表达数据

count.merge<- merge(v3,v4,by="gene_id") ##合并表达数据

批量读取合并文件

```
setwd() ##指定工作环境，即文件的储存位置
filelist<- list.files() ##获取文件夹的文件名
merge.counts<- read.table("V3-SAM-2.txt",sep = "\t",col.names = c("gene_id",filelist[1])) ##先指定一个merge.counts文件
merge.counts<- subset(merge.counts,select = -V3.SAM.2.txt) ##删除merge.counts的表达量行，留下进行合并的“by行”(merge.conut<- v3[-1,]),即最开始的merge.counts
for (i in 1:length(filelist)) 
{
new_data<- read.table(filelist[i],sep = "\t",col.names = c("gene_id",filelist[i]))
merge.counts<- merge(new_data,merge.counts,by="gene_id")
} ##批量读取并根据gene_id合并,合并后删除Geneid的行否则非数值的字符会影响后续的分析
```

#### DEseq2 PCA plot

读入整合后的count数据

```
counts_data<- read.csv("C:/Users/在干嘛/Desktop/数据分析/merge_counts.csv",header=T,stringsAsFactors=F)
```

删除gene_id

```<_ 
count_datas2<- counts_data[,-1]
count_datas2<- as.matrix(count_datas2)
rownames(count_datas2)<-counts_data$gene_id ##给count_datas2加上id，以gene_id作为行名
##rownames(data)返回行名 colnames(data)返回列名

rownames(counts_data)<- counts_data[,1]
counts_data<- counts_data[,-1] ##删除带x的行
```

构建condition

```
condition<- factor(c(rep("c1f",15),rep("cf2",15))) ##指定两次重复的样品数量都是15
condition ##查看condition构建的是否正确
```

构建dds对象

```
dds <- DESeqDataSetFromMatrix(countData=count_datas2, DataFrame(condition), design= ~ condition )
or
dds <- DESeqDataSetFromMatrix(countData=count_datas2, colData=coldata,design=~condition)
结果好像没啥变化
```

数据标准化处理

```
rld <- rlogTransformation(dds)##均一化处理，小于30用rlog，大于30用vst
exprSet_new=assay(rld) ##提取归一化处理后的表达数据
write.table(exprSet_new, file="FZH.DESeq2.normalization.txt", sep="\t",quote=F) #保存归一化处理的数据
```

DEseq自带的PCA plot

```
plotPCA(rld, intgroup=c('condition'))
```

精确到具体哪个样本

```
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
plot(pcaData[,1:2],pch=15,col=c(rep("red",15),rep("blue",15)),font=2,font.lab=2)+text(pcaData[,1],pcaData[,2]+1.5,row.names(pcaData),cex=1,font = 2)+legend("right",inset = 0.04,title = "Rep",c("Rep1","Rep2"),pch = 15,col = c("red","blue")) ##加上legend参数时不知道为什么会报错（二进列运算符中有非数值参数），但是依然能出图
```

参考https://www.jianshu.com/p/53258de21108?clicktime=1578979788

### 差异基因提取

前面构建的condition可直接使用或重新构建

```
 condition<- factor(c(rep("c1",15),rep("c2",15)))
 coldata<- data.frame(row.names = colnames(count_datas2),condition)
```

标准化

```
dds<- DESeqDataSetFromMatrix(countData = count_datas2,colData = coldata,design = ~condition)
dds<- dds[rowSums(counts(dds))>1,] ##去除在所有样本中都不表达的基因，也可在构建dds之前就过滤 
##All_data2<- All_data2[rowSums(All_data2)>1,] 去掉行不大于1的行
dds<- DESeq(dds)
```

差异基因提取

```
res <- results(dds,contrast=c("condition","c","t"))
res<- res[order(res$pvalue),] #根据p值排序
summary(res)#对res矩阵总结，统计有多少个genes上调和下调
write.csv(res,file = "C:/Users/在干嘛/Desktop/数据分析/all_result.csv")
```



#### 多样本的差异基因构建脚本



```
##声明变量
sub_data<- read.csv(file="C:/Users/在干嘛/Desktop/数据分析/controls.csv",header=T,stringsAsFactors=F)
Geneid<- All_data$gene_id ##预设进行merge的"gene_id"
Geneid<- data.frame(Geneid)
colnames(Geneid)<- c("gene_id")
colist<- colnames(All_data) ##储存进行差异比对的列名

for (i in 2:length(colist))
{
  ##构建差异组
  data1<- data.frame(All_data[,1],All_data[,i])##根据列号读取列，分别是"gene_id"列和“样本列”,保存在data1
  colnames(data1)<- c("gene_id",colnames(All_data[i]))##对data1命名行名
  merge.count<- merge(sub_data,data1,by="gene_id")##合并对照组和data1，储存在merge.count
  data2<- merge.count[,-1]##对merge.count重新整理，主要是把"gene_id"行变成行名，储存在data2
  data2<- as.matrix(data2)
  rownames(data2)<- merge.count$gene_id
  
  ##差异分析
  condition<- factor(c(rep("wt",2),rep(colnames(All_data[i]),1)))##wt即对照是两个重复，两个分析的样本中至少有一个样本是两次重复或者更多
  coldata<- data.frame(row.names = colnames(data2),condition)
  dds<- DESeqDataSetFromMatrix(countData = data2,colData = coldata,design = ~condition)
  dds<- DESeq(dds)
  
  ##差异基因提取
  res<- results(dds,contrast = c("condition","wt",colnames(All_data[i])))
  data3<- cbind(rownames(res),res$baseMean)##取出res的行名"gene_id"和"baseMean"合并为data3，即输入样本[i]和对照的差异基因
  colnames(data3)<- c("gene_id",colnames(All_data[i]))##给data3的列命名，分别为"gene_id"和输入样本名
  
  ##差异基因合并
  Geneid<- merge(Geneid,data3,by="gene_id")
}
```

```
Geneid2<- read.csv2("C:/Users/在干嘛/Desktop/数据分析/Geneid2.csv",header=T,stringsAsFactors=F,sep=",") ##以","做分隔符

Geneid2<- as.data.frame(lapply(Geneid2,as.numeric)) ##将data.frame中character转换为数值型

Geneid2<- Geneid2[rowSums(Geneid2)>1,] ##去除0的行
```

bam文件转bw文件

```
bedtools bamtobed -i ${name}.clean.bam > ${name}.bed
bedtools bamtobed -i ${name}.sorted.bam > ${name}.sorted.bed
bedtools genomecov -i ${name}.bed -split -bg -g $INDEX.genome.sizes > ${name}.bg
wigToBigWig ${name}.bg $INDEX.genome.sizes ${name}.bw
```

```
##两两比较
for i in ATHB13_col_a ATHB18_col_a ATHB6_col_a ATHB7_col_a
do
for j in ATHB13_col_a ATHB18_col_a ATHB6_col_a ATHB7_col_a
do
 if [[ "$i" != "$j" ]]; then
 bedtools intersect -a /scratch/xsli/AtDAPseq/dap_data_v4/peaks/Homeobox_tnt/$i/chr1-5/chr1-5_GEM_events.narrowPeak -b /scratch/xsli/AtDAPseq/dap_data_v4/peaks/Homeobox_tnt/$j/chr1-5/chr1-5_GEM_events.narrowPeak > /scratch/xsli/AtDAPseq/dap_data_v4/peaks/Homeobox_tnt/${i}_${j}
 sort -n /scratch/xsli/AtDAPseq/dap_data_v4/peaks/Homeobox_tnt/${i}_${j} | uniq > /scratch/xsli/AtDAPseq/dap_data_v4/peaks/Homeobox_tnt/${i}_${j}_sorted
 rm /scratch/xsli/AtDAPseq/dap_data_v4/peaks/Homeobox_tnt/${i}_${j}
 fi
done
done
```

去除矩阵中NA的行

```
 aa<- sam_all_deg[complete.cases(sam_all_deg),]##complete.cases函数
 aa<- as.data.frame(aa)
 参考https://blog.csdn.net/qq_36481674/article/details/111478445
```

去除重复的gene_id

```
bb<- aa[!duplicated(aa$aa),]
bb<- as.data.frame(bb)
参考https://www.jianshu.com/p/994f2d1a3002
```

注意：先去NA再去重，否则会造成后续分析的问题

找出差异表达基因对应的标准化后的counts

```
R中的操作-筛选连接
#semi_join(x,y)
library(dplyr)
ear_all_deg<- semi_join(ear_all,bb,by='gene_id') #根据bb和ear_all中所共有的以gene_id为列名的列在ear_all中删选bb中所对应的行
参考https://blog.csdn.net/weixin_42437924/article/details/108701702

linux中的操作-join函数
awk -F ',' '{print$5}' SAM_DEG.csv > test ##取出SAM_DEG.csv文件第五列（以','为分隔符）保存在test
sort -k 1b,1 test > test2 ##在使用join之前，文件最好使用sort进行排序处理。
join -t ',' -1 1 -2 1 test2 SAM_All_normalized.csv > test3 ##提取出对应行的数据
```

