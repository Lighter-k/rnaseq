# RNA-seq分析流程

### 一、下载参考基因组文件和注释文件

```
mkdir -p /Data/maize_ref #创建存放参考基因组文件和注释文件的文件夹maize_ref
nohup wget -p /Data/mazie_ref 参考基因组下载地址 & #使用"nohup &"可以将命令在后台运行
nohup wget -p /Data/maize_ref 注释文件下载地址 & 
```

一般上传到网上的测序数据为sra格式（如NCBI数据库），sra格式文件需要转换成fq格式，而正常测序的原始数据都是fq格式的文件。

### fastq文件

一般我们从公司得到的原始测序文件都是fastq.gz文件，gz是一种压缩格式，如果需要对文件进行解压可使用以下命令：

```
gzip -d 19R576_combined_R1.fastq.gz 19R576_combined_R2.fastq.gz 19R577_combined_R1.fastq.gz 19R577_combined_R2.fastq.gz #解压gz文件，解压好的文件即正常的fstq文件。
```

illumina采用双端测序（paired-end），一个样本得到的是seq_R1.fastq.gz和seq_R2.fastq.gz两个文件，每个文件存放一端测序文件。在illumina的测序的cDNA短链被修饰为以下形式（图源见水印）：

![img](https://upload-images.jianshu.io/upload_images/3194654-34575c8f297fc99e.jpg)

两端的序列分别是：

- terminal sequence（两端序列，用以保护碱基）
- adapter（接头序列）
- index（索引序列，用以区分测序数据来自哪个样本）
- primer binding site（引物结合位点）

illumina公司测序所得文件经过处理以fastq文件协议存储为*.fastq格式文件。在fastq文件中每4行存储一个read。

- 第一行，以@开头接readID和其他信息，
- 第二行，read测序信息
- 第三行，必须以+开头，如果后面有内容必须与第一行@后面的内容相同
- 第四行，每个碱基的质量得分，用ASCII码表示。

![img](https://upload-images.jianshu.io/upload_images/3194654-3a962a7125288d55.jpg)

测序的质控环节还是使用gz格式的文件作为输入文件。

```
zcat filename.fastq.gz|head #查看fq.gz文件的头五行
```



### 二、数据质控

测序数据常常包含着各种接头序列以及低质量的reads，在进行后续分析之前需要对测序结果进行筛选。

Trimmomatic软件——数据质控 #软件安装使用conda一键安装

### Trimmomatic过滤

Trimmomatic工具是用于illumina二代测序数据的reads处理，主要对接头（adapter）序列和低质量序列进行过滤。根据单端和双端测序两种模式，Trimmomatic也有两种质控方法。catadapt也可用于去除接头

- **1.SE模式**

  SE模式下只有一个输入文件和一个质控后的输出文件，运行命令如下：

  ```
  Java –jar < trimmomatic的安装路径> SE –threads <线程数> <input> <output> <step1> <step2> …<step1><step2>… 表示每一步的质控参数
  ```

- **2.PE模式**

  PE模式下，有两个输入文件（正向测序reads和反向测序reads））和四个质控后的输出文件（双端序列都保留的paired序列文件和只保留一端序列的unpaired序列文件），运行命令如下：

  ```
  java -jar $trimmomatic PE -threads 12 -phred33 $R1.fq.gz $R2.fq.gz $R1.paired.fq.gz $R1.unpaired.fq.gz $R2.paired.fq.gz $R2.unpaired.fq.gz ILLUMINACLIP:$adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  参数设置说明：
  同一个命令下的不同参数用“:”来界定
  $表示软件或文件所主的路径
  R1.fq.gz、R2.fq.gz是两个输入文件；R1.paired.fq.gz 、R1.unpaired.fq.gz、 R2.paired.fq.gz 、R2.unpaired.fq.gz是对应的四个输出文件
  Phred33 设置碱基的质量格式，默认的是-phred64。
  ```

  ```
  ILLUMINACLIP:$adapter.fa:2:30:10 adapter.fa为接头文件，2表示最大mismatch数，30表示palindrome模式下碱基的匹配阈值，10表示simple模式下碱基的匹配阈值。
  LEADING: 3 表示切除reads 5’端碱基质量低于3的碱基。
  TRAILING:3 表示切除3’ 端碱基质量低于3的碱基。
  SLIDINGWINDOW:4:15 表示以4个碱基为窗口进行滑动，切除窗口内碱基平均质量小于15的。
  MINLEN:36 丢弃以上步骤处理后，序列长度小于36的reads。
  ```

  运行代码如下：

  ```
  trimmomatic PE -threads 12 -phred33 /Data/dek2014/170722A_dekWT1_R1.fq.gz /Data/dek2014/170722A_dekWT1_R2.fq.gz /Data/dek2014/170722A_dekWT1_R1.paired.fq.gz /Data/dek2014/170722A_dekWT1_R1.unpaired.fq.gz /Data/dek2014/170722A_dekWT1_R2.paired.fq.gz /Data/dek2014/170722A_dekWT1_R2.unpaired.fq.gz ILLUMINACLIP:$adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  可以使用“&”让程序在后台运行，也可使用“ctrl+z、bg”组合让程序在后台运行
  ```

  开始执行去接头

  ![image-20201206174959239](C:\Users\在干嘛\AppData\Roaming\Typora\typora-user-images\image-20201206174959239.png)

  运行结束

  ![image-20201206175258700](C:\Users\在干嘛\AppData\Roaming\Typora\typora-user-images\image-20201206175258700.png)

最终的输出结果有两个paired文件和两个unpaired文件，接下来的序列比对用文件只需要使用到两个paired文件。

### Fastqc质检

fastqc是一款基于Java的软件，它可以快速地对测序数据进行质量评估。

使用命令fastqc -o <output dir> <seqfile1,seqfile2..>

```
fastqc -f fastq -q -o /Data/dek2014/clean_data/fastqc_result /Data/dek2014/clean_data/样品名.paired.fq.gz #基础使用，-f指定文件类型默认为是fastq，这里指定为fastq，-o指定输出文件位置
```

参数定义：

```
-o --outdir fastqc 根据报告文件的储存路径，生成的报告的文件名是根据输入来设定的，一般是最后指定这个参数
-t --threads 选择程序运行的线程数
-q --quite 安静运行模式，一般不选这个选项的时候，程序会实时报告运行的状况
-f --format 指定输入的文件类型，支持bam|sam|fastqc
```

```
ls *gz | while read id; do fastqc -t 4 $id; done #脚本运行自动使用fastqc进行质检，直到读完所有以gz为后缀的文件
fastqc -o . *.fastq.gz #将所有数据进行质控，得到zip的压缩文件和html文件，-o后面的空格表示直接输入到当前文件，.后面也有空格
```

fastqc会生成一个zip文件和一个html文件，html文件为质检报告。当多个样本时，会产生多个html报告，可以使用multiqc整合多个qc的报告。支持fastqc、trimmomatic、bowtie、STAR等多种软件的结果整合。

##### multiqc整合

```
multiqc 目标文件路径/*fastqc.zip --pdf #创建一个pdf报告
multiqc . #扫描当前文件夹
multiqc --ignore *_R2* #--ignore参数，忽略某些文件
multiqc --file-list_my_file_list.txt #--file参数，使用文本指定要分析的文件路径
```

multiqc分析结果默认命名为“multiqc_report.html"，-n参数改变结果文件的名字，-o改变输出文件的位置，-f参数输出结果时会自动覆盖同名文件。

### 比对到参考基因组

#### tophat2比对

##### 构建bowtie2索引

bowtie是一个高效的端序列拼接至模板基因组的工具，适合将小序列比对到大基因组上。最长能读取1024个碱基。模板序列最短不能小于1024个碱基。**在使用bowtie前需要使用bowtie-build来构建比对模板（index）。**通常没有参考基因组的常用bowtie，有参考基因组的常用hisat2。

```
nohup bowtie2-build /Data/maize_ref/Zea_mays.B73_RefGen_v4.dna.toplevel.fa maize_B73_index >bowtie2-build.log 2>&1 & #以fa文件作为输入，maize_B73_index为索引，构建参考基因组的bowtie2索引。也有直接用gz文件比对
```

![image-20201209135853292](C:\Users\在干嘛\AppData\Roaming\Typora\typora-user-images\image-20201209135853292.png)

##### tophat2比对

```
nohup tophat2 -p 6 -G /gtf文件路径 -o /bowtie建立的参考基因组索引的前缀（公共名） /双端测序中的一端 /双端测序中的一端 >wt1.log(程序运行日志) 2>&1 &
```



#### hista2比对

##### 构建hista2索引

```
hisat2-build -p 4 /Data/maize_ref/Zea_mays.B73_RefGen_v4.dna.toplevel.fa Zea_mays #4线程构建索引，以Zea_mays为前缀，构建的index应该是8个
```

![image-20201210183004643](C:\Users\在干嘛\AppData\Roaming\Typora\typora-user-images\image-20201210183004643.png)

##### hista2比对

```
/root/miniconda3/bin/hisat2 -p 6 -x /Data/maize_ref/hisat2_index/Zea_mays -1 /Data/dek2014/clean_data/170722A_dekWT1_R1.paired.fq -2 /Data/dek2014/clean_data/170722A_dekWT1_R2.paired.fq -S /Data/dek2014/hisat2_blast/wt1_paired.hista2_1.sam > program_2.log 2>&1 &
未配置环境变量，以绝对路径打开软件
vi ~/.bashrc
PATH=/home/u883604/bio_soft/hisat2/:$PATH # hisat2所在路径
source ~/.bashrc
-p 线程数，本代码指定为6
-x 参考基因组索引文件目录和前缀（genome）
-1 双端测序中的一段文件
-2 双端测序中的另一端文件
-S 输出的sam文件，本命令指定为hista2_blast文件下的dek2014_paired.hista2_1.sam
-q 输入文件为fastq格式，fastq格式为默认参数。
-qseq 输入文件为qseq格式
-f 输入文件为fasta格式
```

