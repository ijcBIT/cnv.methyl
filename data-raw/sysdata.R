source("./data-raw/anno.R")

library("RFpurify")
RFpurify_ABSOLUTE<-RFpurify::RFpurify_ABSOLUTE
library(ChAMP)
data(hm450.manifest.hg19)
data("AllGenes")
data("CancerGenes")
control<-readRDS("/data/SCNA/control.rds")
control@intensity<-control@intensity[,1:25]
# install.packages("Polychrome")
library(Polychrome)

# build-in color palette
cols = glasbey.colors(32)
P36 = createPalette(36,  c("#ff00ff", "#f0000f", "#f0f0f0"))
cols<-c("#02055a",P36[-1],"#aa6400","#005900")
#swatch(cols)
usethis::use_data(AllGenes,CancerGenes,cols,anno_450K,anno_epic,anno_overlap,RFpurify_ABSOLUTE,hm450.manifest.hg19, overwrite = TRUE,internal = TRUE)
