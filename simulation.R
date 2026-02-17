
#https://cran.r-project.org/web/packages/metaRNASeq/index.html
library(metaRNASeq)
data(param)
data(dispFuncs)
idx <- which(param[["DE"]])[1:500]
param[["DE"]] <- FALSE
param[["DE"]][idx] <- TRUE
labels <- param[["DE"]] 
param[[2]] <- param[[1]]
set.seed(100)
param[[2]][param[["DE"]]] <- param[[1]][param[["DE"]]]*sample(c(0.5,0.7,rep(1,5),1.5,2),length(which(param[["DE"]])),replace=TRUE) 
classes <- list(rep(1:2,each=4))
matsim1 <- sim.function(param = param, dispFuncs = dispFuncs,classes=classes)

data(param)
data(dispFuncs)
param[["DE"]] <- labels
param[[2]] <- param[[1]]
set.seed(200)
param[[2]][param[["DE"]]] <- param[[1]][param[["DE"]]]*sample(c(0.5,0.7,rep(1,10),1.5,2),length(which(param[["DE"]])),replace=TRUE) 
classes <- list(rep(1:2,each=3))
matsim2 <- sim.function(param = param, dispFuncs = dispFuncs,classes=classes)


data(param)
data(dispFuncs)
param[["DE"]] <- labels
param[[2]] <- param[[1]]
set.seed(300)
param[[2]][param[["DE"]]] <- param[[1]][param[["DE"]]]*sample(c(0.5,0.9,rep(1,10),1.1,2),length(which(param[["DE"]])),replace=TRUE) 
classes <- list(rep(1:2,each=2))
matsim3 <- sim.function(param = param, dispFuncs = dispFuncs,classes=classes)

data(param)
data(dispFuncs)
param[["DE"]] <- labels
param[[2]] <- param[[1]]
set.seed(400)
param[[2]][param[["DE"]]] <- param[[1]][param[["DE"]]]*sample(c(0.5,0.9,rep(1,20),1.1,2),length(which(param[["DE"]])),replace=TRUE) 
classes <- list(rep(1:2,each=2))
matsim4 <- sim.function(param = param, dispFuncs = dispFuncs,classes=classes)



data(param)
data(dispFuncs)
param[["DE"]] <- labels
param[[2]] <- param[[1]]
set.seed(500)
idx <- sample(which(!labels),500)
param[["DE"]][idx] <- TRUE
param[[2]][param[["DE"]]] <- param[[1]][param[["DE"]]]*sample(c(0.5,0.9,rep(1,10),1.1,2),length(which(param[["DE"]])),replace=TRUE) 
classes <- list(rep(1:2,each=2))
matsim5 <- sim.function(param = param, dispFuncs = dispFuncs,classes=classes)


data(param)
data(dispFuncs)
param[["DE"]] <- labels
param[[2]] <- param[[1]]
set.seed(600)
idx <- sample(which(!labels),300)
param[["DE"]][idx] <- TRUE
param[[2]][param[["DE"]]] <- param[[1]][param[["DE"]]]*sample(c(0.5,0.9,rep(1,10),1.1,2),length(which(param[["DE"]])),replace=TRUE) 
classes <- list(rep(1:2,each=4))
matsim6 <- sim.function(param = param, dispFuncs = dispFuncs,classes=classes)


data(param)
data(dispFuncs)
param[["DE"]] <- labels
param[[2]] <- param[[1]]
set.seed(700)
idx <- sample(which(!labels),300)
param[["DE"]][idx] <- TRUE
param[[2]][param[["DE"]]] <- param[[1]][param[["DE"]]]*sample(c(0.5,0.9,rep(1,10),1.1,2),length(which(param[["DE"]])),replace=TRUE) 
classes <- list(rep(1:2,each=3))
matsim7 <- sim.function(param = param, dispFuncs = dispFuncs,classes=classes)



data(param)
data(dispFuncs)
param[["DE"]] <- labels
param[[2]] <- param[[1]]
set.seed(800)
param[[2]][param[["DE"]]] <- param[[1]][param[["DE"]]]*sample(c(0.5,0.9,rep(1,20),1.1,2),length(which(param[["DE"]])),replace=TRUE) 
classes <- list(rep(1:2,each=5))
matsim8 <- sim.function(param = param, dispFuncs = dispFuncs,classes=classes)


# data(param)
# data(dispFuncs)
# param[["DE"]] <- labels
# param[[2]] <- param[[1]]
# set.seed(500)
# param[[2]][param[["DE"]]] <- param[[1]][param[["DE"]]]*sample(c(0.5,0.9,rep(1,10),1.1,2),length(which(param[["DE"]])),replace=TRUE) 
# classes <- list(rep(1:2,each=2))
# matsim5 <- sim.function(param = param, dispFuncs = dispFuncs,classes=classes)



# data(param)
# data(dispFuncs)
# classes <- list(rep(1:2,each=3))
# param[[2]][param[["DE"]]] <- param[[1]][param[["DE"]]] 
# param[[3]] <- rep(FALSE,length(param[[3]])) 
# matsim5 <- sim.function(param = param, dispFuncs = dispFuncs,classes=classes)



sim.conds <- colnames(matsim1)
rownames(matsim1) <- paste("tag", 1:dim(matsim1)[1],sep="")
dim(matsim1)

sim.conds <- colnames(matsim2)
rownames(matsim2) <- paste("tag", 1:dim(matsim2)[1],sep="")
dim(matsim2)


sim.conds <- colnames(matsim3)
rownames(matsim3) <- paste("tag", 1:dim(matsim3)[1],sep="")
dim(matsim3)


sim.conds <- colnames(matsim4)
rownames(matsim4) <- paste("tag", 1:dim(matsim4)[1],sep="")
dim(matsim4)


sim.conds <- colnames(matsim5)
rownames(matsim5) <- paste("tag", 1:dim(matsim5)[1],sep="")
dim(matsim5)

sim.conds <- colnames(matsim6)
rownames(matsim6) <- paste("tag", 1:dim(matsim6)[1],sep="")
dim(matsim6)

sim.conds <- colnames(matsim7)
rownames(matsim7) <- paste("tag", 1:dim(matsim7)[1],sep="")
dim(matsim7)

sim.conds <- colnames(matsim8)
rownames(matsim8) <- paste("tag", 1:dim(matsim8)[1],sep="")
dim(matsim8)


simstudy1 <- extractfromsim(matsim1,"study1")
simstudy2 <- extractfromsim(matsim2,"study1")
simstudy3 <- extractfromsim(matsim3,"study1")
simstudy4 <- extractfromsim(matsim4,"study1")
simstudy5 <- extractfromsim(matsim5,"study1")
simstudy6 <- extractfromsim(matsim6,"study1")
simstudy7 <- extractfromsim(matsim7,"study1")
simstudy8 <- extractfromsim(matsim8,"study1")

dds1 <- DESeq2::DESeqDataSetFromMatrix(countData = simstudy1$study,
   colData = simstudy1$pheno,design = ~ condition)
res1 <- DESeq2::results(DESeq2::DESeq(dds1))

dds2 <- DESeq2::DESeqDataSetFromMatrix(countData = simstudy2$study,
   colData = simstudy2$pheno,design = ~ condition)
res2 <- DESeq2::results(DESeq2::DESeq(dds2))

dds3 <- DESeq2::DESeqDataSetFromMatrix(countData = simstudy3$study,
   colData = simstudy3$pheno,design = ~ condition)
res3 <- DESeq2::results(DESeq2::DESeq(dds3))

dds4 <- DESeq2::DESeqDataSetFromMatrix(countData = simstudy4$study,
   colData = simstudy4$pheno,design = ~ condition)
res4 <- DESeq2::results(DESeq2::DESeq(dds4))

dds5 <- DESeq2::DESeqDataSetFromMatrix(countData = simstudy5$study,
   colData = simstudy5$pheno,design = ~ condition)
res5 <- DESeq2::results(DESeq2::DESeq(dds5))

dds6 <- DESeq2::DESeqDataSetFromMatrix(countData = simstudy6$study,
   colData = simstudy6$pheno,design = ~ condition)
res6 <- DESeq2::results(DESeq2::DESeq(dds6))

dds7 <- DESeq2::DESeqDataSetFromMatrix(countData = simstudy7$study,
   colData = simstudy7$pheno,design = ~ condition)
res7 <- DESeq2::results(DESeq2::DESeq(dds7))

dds8 <- DESeq2::DESeqDataSetFromMatrix(countData = simstudy8$study,
   colData = simstudy8$pheno,design = ~ condition)
res8 <- DESeq2::results(DESeq2::DESeq(dds8))

rawpval <- list(pval1=res1[["pvalue"]],pval2=res2[["pvalue"]],pval3=res3[["pvalue"]],pval4=res4[["pvalue"]],pval5=res5[["pvalue"]],pval6=res6[["pvalue"]],pval7=res7[["pvalue"]],pval8=res8[["pvalue"]])
adjpval <- list(adjpval1=res1[["padj"]],adjpval2=res2[["padj"]],adjpval3=res3[["padj"]],adjpval4=res4[["padj"]],adjpval5=res5[["padj"]],adjpval6=res6[["padj"]],adjpval7=res7[["padj"]],adjpval8=res8[["padj"]])
FC <- list(FC1=res1[["log2FoldChange"]],FC2=res2[["log2FoldChange"]],FC3=res3[["log2FoldChange"]],FC4=res4[["log2FoldChange"]],FC5=res5[["log2FoldChange"]],FC6=res6[["log2FoldChange"]],FC7=res7[["log2FoldChange"]],FC8=res8[["log2FoldChange"]])

fishcomb <- fishercomb(rawpval, BHth = 0.05)
invnormcomb <- invnorm(rawpval,nrep=c(4,3,2,5,2,4,3,5),BHth = 0.05)


p1 <- fishcomb$adjpval
p2 <- invnormcomb$adjpval
p1[is.na(p1)] <- 1
p2[is.na(p2)] <- 1

library(ROCR)

pred1 <- prediction(1-p1,labels)
perf1 <- performance(pred1, "tpr", "fpr")

pred2 <- prediction(1-p2,labels)
perf2 <- performance(pred2, "tpr", "fpr")

plot(perf1,ylim=c(0,1))
plot(perf2,add=TRUE,col=2)



.fdXfun <- function(pval, padj, labels, thresholdX,transformation = "1-x")
{
    tf <- function(x) eval(parse(text=transformation))
    pval <- pval+(1e-20)
    score <- tf(pval)
    id <- which(labels == 1)	
    x <- 1:length(id)
    o <- order(score,decreasing=TRUE)
    y <- !o[x] %in% id
    y <- cumsum(y)
    if(!is.null(padj) & !is.null(thresholdX))
    {
        thresholdX <- thresholdX+(1e-20)
        thresholdX <- tf(thresholdX)
        padj <- padj+(1e-20)
        scoreX <- tf(padj)      
        yX <- sum(scoreX > thresholdX & labels==0)
        xX <- approx(y=x, x=y, xout=yX)$y     	
    	
    	
    }
 
    else
    {   
         yX <- xX <- NULL 
    }
    list(number=x, fd=y, numberX=xX, fdX=yX)  
}





library(data.table)
x1 <- do.call("cbind",rawpval)
#x2 <- do.call("cbind",FC)
#x <- cbind(x1,x2)
x <- x1
gene <- rownames(matsim1)
index <- rep("pval",ncol(x))
index[grep("FC",colnames(x))] <- "fc"
x[is.na(x)] <- 1

x[,index=="pval"]<- -x[,index=="pval"]
x[,!index=="pval"]<- abs(x[,!index=="pval"])
x <- data.table(x)

hit <- data.table(gene=gene,hits=sample(1:10,length(gene),replace=TRUE))
set.seed(100)
idx <- which(labels)[1:200]
hit[idx,hits:=sample(c(100:10000),length(idx),replace=TRUE)]
hit[,go:=0]
hit[,power:=0]
hit[,labels:=as.numeric(labels)]



niter <- 500000
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



library(foreach)
library(doParallel)
registerDoParallel(cores=10)

#idx <- list()
idx <- foreach(i=seq(niter)) %dopar%
{
   fi <- filter[i,]
   y=x[,mapply(function(x,y) x>y,x=.SD,y=fi)]
   which(rowSums(y)==dim(y)[2])

}
stopImplicitCluster()
names(idx) <- seq(niter)
keep <- which(sapply(idx,length)>0)
idx <- idx[keep]
idx <- idx[!duplicated(lapply(idx,sort))]


score <- sapply(idx,function(x) 
   {y <-hit[x,]
   y[,log(sum(hits))/log(2)]+y[,sum(power)]+y[,sum(go)]})


id <- names(sort(score,decreasing=TRUE))[1]
keep <- idx[[id]]
h1 <- hit[keep,]
setorder(h1,-hits)

id <- names(sort(score,decreasing=TRUE))[1:10]
keep <- idx[id]
h <- lapply(keep, function(z) hit[z,])
h <- rbindlist(h)
h <- h[,list(.N),by=gene]
h <- merge(h,hit,all.x=TRUE)
setorder(h,-N)
fwrite(h,file="sim_result.csv",row.names=FALSE)


fd1 <- .fdXfun(pval=p1,labels=as.numeric(labels),padj=NULL,thresholdX=NULL)
fd2 <- .fdXfun(pval=p2,labels=as.numeric(labels),padj=NULL,thresholdX=NULL)
fd3 <- .fdXfun(pval=h$N,labels=as.numeric(h$labels),padj=NULL,thresholdX=NULL,transformation = "x")
plot(x=fd1$number,y=fd1$fd,type="l",lwd=2,xlim=c(0,500),ylim=c(0,100),xlab="Top ranked features",ylab="Number of false discoveries")
points(x=fd2$number,y=fd2$fd,type="l",lwd=2,col=2)
points(x=fd3$number,y=fd3$fd,type="l",lwd=2,col=4)
legend("topleft",pch=13,lwd=4,col=c(1,2,4), c(" Fisher","Inverse Normal","CCMSB"))


pdf("false_discovery.pdf",w=10,h=10)
plot(x=fd1$number,y=fd1$fd,type="l",lwd=2,xlim=c(0,500),ylim=c(0,100),xlab="Top ranked features",ylab="Number of false discoveries")
points(x=fd2$number,y=fd2$fd,type="l",lwd=2,col=2)
points(x=fd3$number,y=fd3$fd,type="l",lwd=2,col=4)
legend("topleft",pch=13,lwd=4,col=c(1,2,4), c("Fisher","Inverse Normal","CCMSB"))
dev.off()









## 1 simulation





library(data.table)
x1 <- do.call("cbind",rawpval)
#x2 <- do.call("cbind",FC)
#x <- cbind(x1,x2)
x <- x1
gene <- rownames(matsim1)
index <- rep("pval",ncol(x))
index[grep("FC",colnames(x))] <- "fc"
x[is.na(x)] <- 1

x[,index=="pval"]<- -x[,index=="pval"]
x[,!index=="pval"]<- abs(x[,!index=="pval"])
x <- data.table(x)

hit <- data.table(gene=gene,hits=sample(1:10,length(gene),replace=TRUE))
set.seed(100)
idx <- which(labels)[1:200]
hit[idx,hits:=sample(c(100:10000),length(idx),replace=TRUE)]
hit[,go:=0]
hit[,power:=0]
hit[,labels:=as.numeric(labels)]



niter <- 500000
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



library(foreach)
library(doParallel)
registerDoParallel(cores=10)

#idx <- list()
idx <- foreach(i=seq(niter)) %dopar%
{
   fi <- filter[i,]
   y=x[,mapply(function(x,y) x>y,x=.SD,y=fi)]
   which(rowSums(y)==dim(y)[2])

}
stopImplicitCluster()
names(idx) <- seq(niter)
keep <- which(sapply(idx,length)>0)
idx <- idx[keep]
idx <- idx[!duplicated(lapply(idx,sort))]


score <- sapply(idx,function(x) 
   {y <-hit[x,]
   y[,log(sum(hits))/log(2)]+y[,sum(power)]+y[,sum(go)]})


id <- names(sort(score,decreasing=TRUE))[1]
keep <- idx[[id]]
h1 <- hit[keep,]
setorder(h1,-hits)

id <- names(sort(score,decreasing=TRUE))[1:10]
keep <- idx[id]
h <- lapply(keep, function(z) hit[z,])
h <- rbindlist(h)
h <- h[,list(.N),by=gene]
h <- merge(h,hit,all.x=TRUE)
setorder(h,-N)
fwrite(h,file="sim_result_1_1.csv",row.names=FALSE)




niter <- 500000
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



library(foreach)
library(doParallel)
registerDoParallel(cores=10)

#idx <- list()
idx <- foreach(i=seq(niter)) %dopar%
{
   fi <- filter[i,]
   y=x[,mapply(function(x,y) x>y,x=.SD,y=fi)]
   which(rowSums(y)==dim(y)[2])

}
stopImplicitCluster()
names(idx) <- seq(niter)
keep <- which(sapply(idx,length)>0)
idx <- idx[keep]
idx <- idx[!duplicated(lapply(idx,sort))]


score <- sapply(idx,function(x) 
   {y <-hit[x,]
   y[,log(sum(hits))/log(2)]+y[,sum(power)]+y[,sum(go)]})


id <- names(sort(score,decreasing=TRUE))[1]
keep <- idx[[id]]
h1 <- hit[keep,]
setorder(h1,-hits)

id <- names(sort(score,decreasing=TRUE))[1:10]
keep <- idx[id]
h <- lapply(keep, function(z) hit[z,])
h <- rbindlist(h)
h <- h[,list(.N),by=gene]
h <- merge(h,hit,all.x=TRUE)
setorder(h,-N)
fwrite(h,file="sim_result_1_2.csv",row.names=FALSE)







##2 simulation



library(data.table)
x1 <- do.call("cbind",rawpval)
#x2 <- do.call("cbind",FC)
#x <- cbind(x1,x2)
x <- x1
gene <- rownames(matsim1)
index <- rep("pval",ncol(x))
index[grep("FC",colnames(x))] <- "fc"
x[is.na(x)] <- 1

x[,index=="pval"]<- -x[,index=="pval"]
x[,!index=="pval"]<- abs(x[,!index=="pval"])
x <- data.table(x)

hit <- data.table(gene=gene,hits=sample(1:10,length(gene),replace=TRUE))
set.seed(100)
idx <- which(labels)[201:400]
hit[idx,hits:=sample(c(100:10000),length(idx),replace=TRUE)]
hit[,go:=0]
hit[,power:=0]
hit[,labels:=as.numeric(labels)]



niter <- 500000
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



library(foreach)
library(doParallel)
registerDoParallel(cores=10)

#idx <- list()
idx <- foreach(i=seq(niter)) %dopar%
{
   fi <- filter[i,]
   y=x[,mapply(function(x,y) x>y,x=.SD,y=fi)]
   which(rowSums(y)==dim(y)[2])

}
stopImplicitCluster()
names(idx) <- seq(niter)
keep <- which(sapply(idx,length)>0)
idx <- idx[keep]
idx <- idx[!duplicated(lapply(idx,sort))]


score <- sapply(idx,function(x) 
   {y <-hit[x,]
   y[,log(sum(hits))/log(2)]+y[,sum(power)]+y[,sum(go)]})


id <- names(sort(score,decreasing=TRUE))[1]
keep <- idx[[id]]
h1 <- hit[keep,]
setorder(h1,-hits)

id <- names(sort(score,decreasing=TRUE))[1:10]
keep <- idx[id]
h <- lapply(keep, function(z) hit[z,])
h <- rbindlist(h)
h <- h[,list(.N),by=gene]
h <- merge(h,hit,all.x=TRUE)
setorder(h,-N)
fwrite(h,file="sim_result_seed200_1.csv",row.names=FALSE)




niter <- 500000
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



library(foreach)
library(doParallel)
registerDoParallel(cores=10)

#idx <- list()
idx <- foreach(i=seq(niter)) %dopar%
{
   fi <- filter[i,]
   y=x[,mapply(function(x,y) x>y,x=.SD,y=fi)]
   which(rowSums(y)==dim(y)[2])

}
stopImplicitCluster()
names(idx) <- seq(niter)
keep <- which(sapply(idx,length)>0)
idx <- idx[keep]
idx <- idx[!duplicated(lapply(idx,sort))]


score <- sapply(idx,function(x) 
   {y <-hit[x,]
   y[,log(sum(hits))/log(2)]+y[,sum(power)]+y[,sum(go)]})


id <- names(sort(score,decreasing=TRUE))[1]
keep <- idx[[id]]
h1 <- hit[keep,]
setorder(h1,-hits)

id <- names(sort(score,decreasing=TRUE))[1:10]
keep <- idx[id]
h <- lapply(keep, function(z) hit[z,])
h <- rbindlist(h)
h <- h[,list(.N),by=gene]
h <- merge(h,hit,all.x=TRUE)
setorder(h,-N)
fwrite(h,file="sim_result_seed200_2.csv",row.names=FALSE)




library(data.table)
f1 <- fread("sim_result_1_1.csv")
f2 <- fread("sim_result_1_2.csv")
f3 <- fread("sim_result_2_1.csv")
f4 <- fread("sim_result_2_2.csv")

g1 <- f1[N==10&hits<=10,gene]
g2 <- f2[N==10&hits<=10,gene]
g3 <- f3[N==10&hits<=10,gene]
g4 <- f4[N==10&hits<=10,gene]

# Create example lists
x <- list(g1=g1,g2=g2,g3=g3,g4=g4)

library(VennDiagram)
venn.diagram(
  x = x,
  filename = "venn_diagram.png",
  output = TRUE,
  height = 2000,
  width = 2000,
  resolution = 300,
  fill = c("blue","red","purple","green"),
  alpha = 0.5
)