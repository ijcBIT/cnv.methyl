# Deprecated as of 26/04/20222. See data-raw/anno.R

anno_epic = data.table::setDT(melt(as.data.frame(
  minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
  , id.vars = "Name"))
anno_450k = data.table::setDT(melt(as.data.frame(
  minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19))
  , id.vars = "Name"))
setkey(anno_epic,Name,variable)
overlap=merge(anno_epic,anno_450k,all=FALSE)[,value:=ifelse(value.x == value.y,value.x,paste0(value.x,";",value.y))]
anno_merged<-dcast.data.table(overlap,
                              Name~variable,value.var = "value")
anno_overlap<-as(anno_merged,"DataFrame")
sysdata_filenames <- load("R/sysdata.rda")
#save(list = c(sysdata_filenames, "anno_overlap"), file = "R/sysdata.rda")
#usethis::use_data(anno_overlap, overwrite = TRUE)
