source("./data-raw/anno.R")

library("RFpurify")
RFpurify_ABSOLUTE<-RFpurify::RFpurify_ABSOLUTE
library(ChAMP)
data(hm450.manifest.hg19)
data("AllGenes")
data("CancerGenes")
control<-readRDS("/data/SCNA/control.rds")
control@intensity<-control@intensity[,1:25]
usethis::use_data(AllGenes,CancerGenes,anno_450K,anno_epic,anno_overlap,RFpurify_ABSOLUTE,hm450.manifest.hg19, overwrite = TRUE,internal = TRUE)
