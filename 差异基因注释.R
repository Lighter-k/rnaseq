##差异基因注释、KEGG、GO
library('biomaRt')
library("curl")
listMarts(host="plants.ensembl.org") #显示能连接的数据库
ensembl_plants=useMart("plants_mart",host="plants.ensembl.org") #选定数据库
listDatasets(ensembl_plants) #显示对应数据库的数据集
mart <- useDataset("zmays_eg_gene", useMart("plants_mart", host="plants.ensembl.org")) #指定玉米的数据集
filelist<- list.files()
for (i in filelist) {
  counts_data<- read.csv(i,header=T,stringsAsFactors=F)
  rownames(counts_data)<- counts_data[,1]
  my_ensembl_gene_id<- row.names(counts_data) #提取出要注释的基因id
  maize_symbols<- getBM(attributes = c('ensembl_gene_id','external_gene_name',"description"),filters = 'ensembl_gene_id',values = my_ensembl_gene_id,mart = mart) #对提取出的基因id进行注释
  ##head(maize_symbols) #查看基因注释结果
  name<- paste(sub(".csv","",i),"_annotation")##设置文件名字
  write.csv(maize_symbols,paste0(name,".csv"))
}


##GO
library(clusterProfiler)
go_term2gene <- read.table("..\\GO\\go_term2gene.txt",header = T,quote = "",sep = "\t")
go_term2name <- read.table("..\\GO\\go_term2name.txt",header = T,quote = "",sep = "\t")
##for (i in filelist) {
counts_data<- read.csv("C:/Users/在干嘛/Desktop/数据分析/pheatmap/所有2个重复样本/时间序列分析/sam_vs_control.csv",header=T,stringsAsFactors=F)
gene<- counts_data$X #取出所有的gene_id
go <- enricher(gene,TERM2GENE=go_term2gene,TERM2NAME=go_term2name)
barplot(go,font.size = 5)
##}
##KEGG
test <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene_id"),values= gene,mart= mart)
hak <- test$entrezgene_id
hak<- na.omit(hak)
kegg = enrichKEGG(gene=hak, organism ="zma", pvalueCutoff=1,qvalueCutoff=1)
dotplot(kegg,font.size=,label_format = 6)##label_format指定文字换行的长度

#提取go和kegg的数据
go_data<- as.data.frame(go)
kegg_data<- as.data.frame(kegg)

或许可以参考https://www.jianshu.com/p/eee2cc315f77 提取出不同cluster 与目标通路相关的富集点，做成一副图。