## Download Whole Blood samples.

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GEOquery")
#install.packages("data.table")
library(data.table)
setDTthreads(0L)
library(GEOquery)
## Whole Blood:
##https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73103

n=20
gsms<-1886364:1886803
cdir="./data-raw/WholeBlood_Controls"
if(!dir.exists(cdir))dir.create(cdir)
cdirf=paste(cdir,"Females",sep = "/")
if(!dir.exists(cdirf))dir.create(cdirf)
cdirm=paste(cdir,"Males",sep = "/")
if(!dir.exists(cdirm))dir.create(cdirm)

#Fetch all 48 Females Age >= 18
cl<- parallel::makePSOCKcluster(get_ncores(ncores), outfile="")
parallel::clusterEvalQ(cl,{
  library("GEOquery")

})
doParallel::registerDoParallel(cl)

res<-foreach::foreach(k=gsms),#isplitIndices(1400,chunks=ncores),
                      .combine='c',
                      .multicombine = F,
                      .inorder=F,
                      .packages = "data.table",
                      .export = c("AllGenes","CancerGenes"),
                      .errorhandling = "pass"
)%dopar%{

  countf=length(list.dirs(cdirf))-1
  countm=length(list.dirs(cdirm))-1

  gsm <- GEOquery::getGEO(sprintf("GSM%i",k));
  Sex <- gsub(" ","",strsplit(GEOquery::Meta(gsm)[["characteristics_ch1"]][c(1,3)][1],split=":")[[1]][2]);
  Age <- gsub(" ","",strsplit(GEOquery::Meta(gsm)[["characteristics_ch1"]][c(1,3)][2],split=":")[[1]][2]);
  if(Sex=="Female" & Age >=18 & countf < n){GEOquery::getGEOSuppFiles(sprintf("GSM%i",k),fetch_files = T,baseDir = cdirf)
  if(countf>=n &countm>=n ){}
  }
  }


}
parallel::stopCluster(cl)

countf=length(list.dirs(cdirf))-1
for(k in 1886364:1886803){
  gsm <- GEOquery::getGEO(sprintf("GSM%i",k));
  Sex <- gsub(" ","",strsplit(GEOquery::Meta(gsm)[["characteristics_ch1"]][c(1,3)][1],split=":")[[1]][2]);
  Age <- gsub(" ","",strsplit(GEOquery::Meta(gsm)[["characteristics_ch1"]][c(1,3)][2],split=":")[[1]][2]);
  if(Sex=="Female" & Age >=18){GEOquery::getGEOSuppFiles(sprintf("GSM%i",k),fetch_files = T,baseDir = cdirf);countf=countf+1;
  print(count)
  if(count==n)break;
  }
}

#Fetch  48 Males Age >= 18
dir_females <- gsub("./","",list.dirs(cdirf)[-1])
GSM_all <- paste0("GSM",c(1886364:1886588))
GSM_males <- GSM_all[!GSM_all%in%dir_females]


count=length(list.dirs(cdirm))-1
for(GSM in GSM_males){
  count=count+1;
  print(count);
  if(count==n)break;
  gsm <- getGEO(GSM);
  Sex <- gsub(" ","",strsplit(Meta(gsm)[["characteristics_ch1"]][c(1,3)][1],split=":")[[1]][2]);
  Age <- gsub(" ","",strsplit(Meta(gsm)[["characteristics_ch1"]][c(1,3)][2],split=":")[[1]][2]);
  if(Sex=="Male" & Age >=18)GEOquery::getGEOSuppFiles(GSM,fetch_files = T,baseDir = cdirm)
}

# Dump all  Male/ and Female/ dirs .idat files in WholeBlood_Controls.
# Load and combine all those and perform preprocessing:
#cgs<-names(anno_450K@probes)
#control<-readRDS("/data/SCNA/control.rds")
idx<-rownames(control@intensity) %in% cgs
control@intensity<-control@intensity[idx,]
usethis::use_data(control, overwrite = TRUE)
