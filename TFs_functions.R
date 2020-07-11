#### functions for characterizing dynamics of transcriptional factors (TF) and 
#### genes encoding secretory proteins (sec genes) across the menstrual cycle
### 1. smooth
wrapper_smooth.norm<-function(ct,spar){
  o_smthed<-wrapper_smooth.spline(ct=ct, spar=spar)
  print(class(o_smthed[["y"]]))
  o_smthed_norm<-norm_minmax(o_smthed[["y"]])
  o<-list(o_smthed[["raw"]],o_smthed[["y"]],o_smthed_norm)
  names(o)<-c("raw","y","y_minmax")
  o
}
wrapper_smooth.spline<-function(ct, spar){
  raw<-apply(ct, MARGIN = 1, function(x) smooth.spline(c(1:ncol(ct)), x, spar=spar))
  y<-lapply(raw, function(x) x$y) 
  y<-as.data.frame(matrix(unlist(y), nrow = nrow(ct), byrow = TRUE));rownames(y)<-rownames(ct)
  o<-list(raw,y);names(o)<-c("raw","y")
  o
}
#[0,1] normalize
norm_minmax<-function(y){
  o<-t(apply(y, MARGIN=1, function(x) {
    x<-x-min(x)
    x<-x/max(x)
  }))
  o<-as.data.frame(o)
}

### 2. locate peak
wrapper_peak<-function(y_norm_minmax){
  library(IDPmisc)
  o<-lapply(c(1:nrow(y_norm_minmax)), function(i) {
    peak<-peaks(x=as.numeric(y_norm_minmax[i,]) 
    )
  })
  names(o)<-rownames(y_norm_minmax)
  o
}

## 3. determine the phase of peak
peak_assign.ph<-function(peaks_time, boundaries){
  # locate
  peak_ph<-lapply(peaks_time, function(a){
    xs<-a$x
    phs<-sapply(xs,function(b){
      phase.of.peak(b, boundaries)
    })
  });names(peak_ph)<-names(peaks_time)
  
  o<-list(peaks_time, peak_ph)
  names(o)<-c("peaks_time","phase")
  o
}

## determine the phase where the peak occurs 
## boundary: (n1,n2,n3,n4) where,
## (0) - PH1 - (n1) - PH2 - (n2) - PH3 - (n3) - PH4 - (n4)
phase.of.peak<-function(x.peak,boundaries){
  signs<-x.peak-boundaries
  #print(signs)
  if(signs[1]<=0){
    1
  }else{
    if(signs[2]<=0){
      2
    }else{
      if (signs[3]<=0){
        3
      }else{
        4
      }
    }
  }
}

### 4. categorize TFs by phase of peak and # of peak
get_peak.cat.unique<-function(peak_ph){
  peak.unique<-unique(peak_ph)
  peak.unique<-order_peak.cat(peak.unique)
}

order_peak.cat<-function(peak.unique){
  npeak<-unlist(lapply(peak.unique,function(x) length(x)))
  p<-lapply(sort(unique(npeak)),function(x){
    o<-peak.unique[npeak==x]
    sum<-sapply(o,function(i) sum(i))
    o<-o[order(sum)]
  })
  p<-unlist(p,recursive = F)
}

assign_peak.cat<-function(peak_ph,peak.unique){
  peak.cat<-unlist(lapply(peak_ph, function(x){
    #print(x)
    grp<-sapply(peak.unique,function(y){
      if(length(y)==length(x)){
        all(y==x)
      }else{
        FALSE
      }
    })
    grp<-which(grp==T)
  }))
}

### 5. order peak 
## first by pt(peak): pseuodotime of peak 
## then by pt0(peak): an estimate of inflection point
wrapper_order_pt.t0<-function(TFs, peak_summary, y_minmax, boundary_l=1, boundary_r=ncol(y_minmax), ph4=T){
  o.pt<-order.peak_pt(TFs, peak_summary,ph4=ph4)
  o.t0<-order.peak_t0(TFs, y_minmax, boundary_l, boundary_r)
  
  order.peak<-cbind(o.pt,o.t0[names(o.pt)])
  order.peak<-order.peak[order(order.peak[,1],order.peak[,2]),]
}

# get pt(peak)
order.peak_pt<-function(TFs, peak_summary, ph4=T){
  pt<-unlist(lapply(TFs,function(x){
    if(ph4==T){
      nr<-nrow(peak_summary[[x]])
    }else{
      nr<-which.max(peak_summary[[x]][,"y"])
    }
    o<-peak_summary[[x]][nr,1]
  }));
  names(pt)<-TFs
  pt<-sort(pt)
}

# get pt0(peak): an estimate of inflection point of a peak
order.peak_t0<-function(TFs, y_minmax, boundary_l,boundary_r){
  t0<-unlist(lapply(TFs, function(x){
    ct_po<-y_minmax[x,c(boundary_l:boundary_r)]
    o<-ct_po[which.min(abs(0-ct_po))]
    o<-as.numeric(substring(names(o),2))
  }));names(t0)<-TFs
  t0
}
