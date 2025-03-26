CCMSB
==========
CCMSB is an algorithm on Integrating GEO transcriptome big data，single-cell transcriptome data, methylome data and gene-variation data in random-combining to select differentiated genes in R in the terms specified by users. The method is based on randomly selecting fold change or p value
judged by the AI-derived or manually derived prior knowledge on each dataset on all omics in filtering differentiated genes in R. 
The differentiated genes identified by Dusim are in high probability of important biological meaning. 



CCMSB是Ｒ平台下的基于不同来源的转录组大数据, 转录组单细胞大数据，甲基化数据和基因突变数据等的多组学整合和提供用户指定的生物学概念相关差异基因的工具，适用于ncbi geo等多组学大数据的整合及其用户指定领域的差异基因筛选工作，筛选到的基因可直接用于实验验证准确率较高。

## Authors

[刘文粟] (Liu Wensu)

[周晓北] (Zhou Xiaobei)

## <a name="install"></a> Prerequisites and installation

### <a name="prerequisites"></a> Prerequisites
- **R (>=3.5.3)** 

### <a name="installlation"></a> Installlation
Temporally, use:
```r
source("~/DUSIM/R/DUSIM.R")
```

## Usage
sample_data完整处理案例展示如下：
### <a name="process"></a> Process data
The input data should be formated as following: 

    gene	pval_1
    NECAP2	0.358843944
    CLIC4	0.030455456
    DBT	0.292270083
    C1orf21	0.000382696
    COP1	0.020345337
    PRUNE1	0.325788738
    LIN9	0.017277068

The first column is gene symbol (with column name's "gene"), the other columns should contain the result of selected group comparison (pvalue (named as pval\_x) or 
fold change (fc\_x)).

读取GEO数据并形式上合并

```r
library(data.table)
dir <- "~/DUSIM/data/sample_data"
setwd(dir)
f <- list.files(pattern="^GSE",full.names=TRUE)
f1 <- lapply(f,fread)
nms <- gsub(".csv","",basename(f))
f1 <- mapply(function(x,y) {nm <- names(x)
                            idx <- !nm=="gene"
                            names(x)[idx] <- paste0(y,"_",nm[idx])
                            x},x=f1,y=nms,SIMPLIFY=FALSE)

f1 <- Reduce(function(...) merge.data.table(...,all=TRUE),f1)
#f1 <- f1[grep("[0-9][0-9]-[A-Z][a-z][a-z]|^__",gene,invert=TRUE)]
f1[is.na(f1)] <- 0
fwrite(f1,"combine.csv")
```
The result should be like:

    gene	GSE55030_p_val	GSE71797_fc_1	GSE71797_fc_2	GSE82223_fc_1	GSE86532_fc_1	GSE92574_fc_1	GSE110903_fc_1
    ISG15	0.029289252	0.475583449	3.077468801	0	0.126152927	0	0.484946133
    AGRN	0.020108192	9.767047772	0.460773615	0	0	0	1.677513801
    RP11-465B22.3	0	0	0	0	0	0	0.125557896
    RP11-465B22.8	0	0	0	0	0	0	0.395281318
    B3GALT6	0.391504072	2.840062329	0.312959305	0	0	0.334790739	1.010640313

读入文献中有报道的标志性基因和已知通路
```r
marker_gene<-c("AR","CD46","STAT5A","STAT5B","NOTCH1","FOXA1",
   "KDM6A","ASXL1","KMT2A","KMT2D","WDR5","ASH2L","EZH2",
   "KDM1A","NR3C1","PIK3CB","AKT1","MTOR","KDM5D","C17ORF49")
power_gene<-sqrt(c(16,8,5,5,8,6,6,6,10,10,10,10,8,8,12,8,8,8,8,8)) 
target_GO<-c("GO:0000380","GO:0016310","GO:0030521","GO:0008083","GO:0006914")
library("org.Hs.eg.db")
go <- as.data.table(org.Hs.egGO)
geneinfo <- as.data.table(org.Hs.egSYMBOL)
go <- merge(go,geneinfo,all.x=TRUE)
gene_go <- go[go_id%in%target_GO,unique(symbol)]


数据二次整合
```r
##"hit.csv" is the co-occurrence file from AI-derived or manually-derived file on downloaded abstracts 
library(data.table)
x <- fread("combine.csv")
gene <- x[,gene]
x[,gene:=NULL]
index <- rep("pval",ncol(x))
index[grep("fc",names(x))] <- "fc"
x[,c(1,5)] <- - x[,c(1,5)]
library(data.table)
hit <- fread("hit.csv")
hit <- merge(hit,mg,all.x=TRUE,sort=FALSE)
hit[is.na(hit)] <- 0
hit[,go:=0]
hit[gene%in%gene_go,go:=1]
```

筛选基因
```r
niter <- 100000
filter_p <- -seq(0.00032,0.32,length=1000)
filter_fc <- seq(1.003,4,length=1000)
filter <- lapply(index,function(x) {
   if(x=="pval")
      y <- sample(filter_p,niter,replace=TRUE)
   else
      y <- sample(filter_fc,niter,replace=TRUE)
      y})
filter <- do.call("cbind",filter)


for(i in seq(niter))
{
   l <- 0
   for(j in seq(dim(filter)[2]))
   {
      jump <- sample(0:dim(filter)[2],1)
      if(jump==0&l<=2)
      {
         filter[i,j] <- -1000
         l <- l+1
      }
      else if(l>2)
         break
   }
}

idx <- list()
for(i in seq(niter))
{
   fi <- filter[i,]
   y=x[,mapply(function(x,y) x>y,x=.SD,y=fi)]
   idx[[i]] <- which(rowSums(y)==dim(y)[2])
   cat("i:", i,"\n")
}
names(idx) <- seq(niter)
keep <- which(sapply(idx,length)>0)
idx <- idx[keep]
idx <- idx[!duplicated(lapply(idx,sort))]
score <- sapply(idx,function(x) 
   {y <-hit[x,]
   y[,log(sum(hits))/log(2)]+y[,sum(power)]+y[,sum(go)]})


id <- names(which(score==max(score)))
keep <- idx[[id]]
h1 <- hit[keep,]
h1 <- h1[go==1,]
h1 <- h1[order(hits,decreasing=TRUE),]
```

差异基因中选取其中的重要通路基因如下:

```r
> h1
       gene hits power go
 1:   NCOA3   36     0  1
 2:     DAP   19     0  1
 3:    KAT5   19     0  1
 4:   PIAS1   18     0  1
 5:   HBEGF   17     0  1
 6:    XBP1   13     0  1
 7:   KDM3A   12     0  1
 8:   PDGFA    7     0  1
 9:     VCP    7     0  1
10:   ATG4B    6     0  1
11: SH3GLB1    4     0  1
12:  TOLLIP    3     0  1
13:   RAB1A    2     0  1
14:   RAB1B    2     0  1
15:    MANF    1     0  1
16:   CDK13    1     0  1
17:    CLTC    1     0  1
18:    VMP1    1     0  1
19:   WIPI1    1     0  1
20:    NADK    0     0  1
21:   PANK4    0     0  1
22:   DRAM2    0     0  1
23:   SGMS2    0     0  1
24:   SGMS1    0     0  1
25:   VPS51    0     0  1
26:   RAB8A    0     0  1
```