install.packages("statmod");library(statmod)
install.packages("fastICA");library(fastICA)

### User define parameters
wd<-""
setwd(wd)

### Load function and data
source("overdispersion_functions.R")
ct<-read.csv("ct.csv",row.names = 1, check.names = F, stringsAsFactors = F) # normalized gene (row) x cell (col) matrix

### calculate overdispersion for a normalized gene (row) by cell (col) matrix ct
genes_ordered<-sel.by.cv(ct)

### select the top 1000 overdispersed genes
od.genes<-rownames(ct)[genes_ordered][1:1000]
