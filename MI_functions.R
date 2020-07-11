### prepare data for MI calculation
mi_formatxy<-function(ct,wd,prefix,n_permute,time){
  print("permuting day...")
  day.pmt<-permute_day(pca.out_rownames=colnames(ct), day_master=time, n=n_permute)
  sum(names(day.pmt$ori)==colnames(day.pmt$pmt))
  
  print("generating data x: permutated data and y: original data...")
  x<-structure(rbind(day.pmt$ori,day.pmt$pmt), dimnames=list(c("X",paste0("NULL",1:n_permute)),names(day.pmt$ori)))
  y<-as.matrix(ct)
  print(dim(x))
  print(dim(y))
  
  print("output x, y for MI calculation...")
  setwd(wd)
  print(paste0("x_",prefix,"_miTest.dat"))
  print(paste0("y_",prefix,"_miTest.dat"))
  format.for.mi(x, paste0("x_",prefix,"_miTest.dat"))
  format.for.mi(y, paste0("y_",prefix,"_miTest.dat"))
  out<-list(x,y);names(out)<-c("x","y")
  return(out)
}

# permute data, such that cell labels are permutated agst day for n times
permute_day<-function(pca.out_rownames, day_master, n){
  require(gtools)
  day.mi<-day_master[,pca.out_rownames]
  day.tmp<-day.mi;names(day.tmp)<-NULL
  day.pmt<-replicate(n,permute(day.tmp))
  print(dim(day.pmt))
  rownames(day.pmt)<-names(day.mi)
  structure(list(day.mi,t(day.pmt)),names=c("ori","pmt"))
}

# format data for mi
format.for.mi<-function(dt, filename){
  tmp <- cbind(genes=rownames(dt), round(dt, 5))
  tmp <- rbind(colnames(tmp), tmp)
  cat(unlist(apply(tmp, 1, paste, collapse="\t"), use.names=FALSE), sep="\n", file=paste(filename, sep=""))
} 

### load MI results
mi_load<-function(miOut.path,mi_xy.out){
  print(paste0("loading miOut:",miOut.path))
  
  miOut<-as.matrix(read.csv(miOut.path, row.names=1, sep="",check.names = F))
  print(dim(miOut))
  print(length(rownames(mi_xy.out$x)))
  print(length(rownames(mi_xy.out$y)))
  miOut<-miOut[rownames(mi_xy.out$x), rownames(mi_xy.out$y)]
  print(dim(miOut))
  
  return(miOut)
}

### ecdf
ecdf.cal.sig<-function(miOut){
  print("calculating ecdf")
  ecdf_out<-apply(miOut[-1,],2,ecdf)
  p_gene<-sapply(1:ncol(miOut),function(a) ecdf_out[[a]](miOut[1,a]));names(p_gene)<-names(ecdf_out)
  
  print("selecting temporal genes")
  p.sig_genes<-p_gene[p_gene==1]
  p.sig_genes.sortedMI<-sort(miOut[1,names(p.sig_genes)],decreasing = T)
  
  genes_sorted.by.MI<-names(sort(miOut[1,],decreasing=T))
  
  ecdf.out<-list(ecdf_out,p_gene,p.sig_genes,p.sig_genes.sortedMI,genes_sorted.by.MI)
  names(ecdf.out)<-c("ecdf_out","p.val_gene","p.sig_genes","p.sig_genes.sortedMI","genes_sorted.by.MI")
  return(ecdf.out)
}

### correct p value based on FDR correction
fdr.adjustment<-function(miOut,fdr){
  library(sgof)
  print(fdr)
  ecdf.out<-ecdf(c(miOut[-1,]))
  p.val<-1-ecdf.out(miOut[1,]);names(p.val)<-colnames(miOut)
  fdr.correction_out<-BH(p.val,fdr)
  out<-list(fdr.correction_out,sort(p.val));names(out)<-c("FDR.correction","p.val")
  out
}