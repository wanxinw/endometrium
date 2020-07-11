# For analyzing dynamics of transcriptional regulators
# and genes encoding secretory genes
rm(list=ls())

### User define parameters
wd<-""
setwd(wd)

### Load functions and data
source("TFs_functions.R")
ct_by.pseudot<-read.csv("ct_by.pseudot.csv",row.names = 1, check.names = F, stringsAsFactors = F) # normalized gene (row) x cell (col) matrix with cells ordered by pseudot
cell.stage<-read.csv("cell.stage.csv",row.names = 1, check.names = F, stringsAsFactors = F) # Phase 1-4 assignment of each cell as in Fig. 3a.

## 1. smooth and cleanup
TF_smth.norm<-wrapper_smooth.norm(ct=ct_by.pseudot,spar=1)

## 2. locate peak
peak_t.w<-wrapper_peak(TF_smth.norm$y_minmax)

## 3. determine the phase where the peak occured
peak_ph<-peak_assign.ph(peak_t.w, boundaries=cumsum(table(cell.stage)))

## 4. categorize TFs by phase of peak and # of peak
peak.cat.unique<-get_peak.cat.unique(peak_ph$phase)

## 5. order TFs of interest
TF_order<-wrapper_order_pt.t0(TF_of.interest, peak_summary=peak_ph$peaks_time, y_minmax=TF_smth.norm$y_minmax, boundary_r=ncol(ct_by.pseudot))
