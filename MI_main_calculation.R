rm(list=ls())

### User define parameters
wd<-""
setwd(wd)
prefix<-"mi"

### Load function and data
source(paste0(wd,"MI_functions.R"))
ct<-read.csv("ct.csv",row.names = 1, check.names = F, stringsAsFactors = F) # gene (row) x cell (col) matrix
day<-t(read.csv("day.csv",row.names = 1)) # day (row) x cell (col) matrix
day<-day[,colnames(ct),drop=F]

### Prepare data for mutual information (MI) calculation
# For each gene, permutate expression against day for n=1000 times
# Format and output both permutated and un-permutated data to wd
mi_xy.out<-mi_formatxy(ct=ct,wd=wd,prefix=prefix,n_permute = 1000,time=day)

### Calculate MI using java script with the following command line (prefix=mi):
# java -Xmx10000M -jar interCellMI.jar -e x_mi_miTest.dat -d y_mi_miTest.dat -o mi_miOut.txt

### import MI output 
miOut<-mi_load(paste0(prefix,"_miOut.txt"),mi_xy.out)

### for each gene, calculate significance level of MI results
# and fdr correct
ecdfOut<-ecdf.cal.sig(miOut)
fdr.correction<-fdr.adjustment(miOut,0.05)

### select genes with significance level that rejects null
sorted.by.MI<-names(sort(miOut[1,],decreasing = T))
selected<-sorted.by.MI[sorted.by.MI%in%names(sort(fdr.correction$p.val))[1:fdr.correction$FDR.correction$Rejections]]

### visualize 
ylim<-max(miOut[1,names(ecdfOut$p.sig_genes.sortedMI)])
boxplot(miOut[-1,names(ecdfOut$p.sig_genes.sortedMI)],outline=T ,main=paste0(prefix,": MI (exp, day)"),xlab="Gene",ylab="MI",ylim=c(0,ylim))
points(miOut[1,names(ecdfOut$p.sig_genes.sortedMI)],col="red",cex=0.4,type = "b")
legend("topright",col=c("red","black"),legend = c("MI (exp, day)","MI (exp, permutated.day)"),pch=15,bty="n")


