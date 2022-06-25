bp_stopifnot = getFromNamespace("stopifnot", "backports")
#`%dopar%` <- foreach::`%dopar%`
check_input_files<-function(Sample_Name,folder,std_file="standard_format"){

  file_list<-list.files(folder,full.names = T)
  path<-sapply(Sample_Name,function(ID){
    match=FALSE
    match<-sapply(file_list, function(file) ifelse(data.table::like(basename(file),ID),TRUE,FALSE))
    if(sum(match) == 0) warning("No input file found for sample: ",ID)
    if(sum(match) > 1){
      #std_file <- paste0(conumee.folder,"/",.folder,"/",ID,"_Segments.txt")
      match_files <- file_list[match]
      std_file<-eval(std_file)
      if(std_file %in% match_files){
        ok <- match_files[which(std_file %in% match_files)]}else{
          ok <- file_list[which.max(match)]}
      overwrite<-setdiff(match_files,ok)
      warning("multiple files match ",ID, ". \n ", paste(overwrite,collapse = ", ")," ignored.")
      return(ok)
    }# otherwise match == 1
    match_files <- file_list[match]
    return(match_files)
  })
  return(path)
}
Kc<-function(Kc=c("curated","balanced")){
  Kclist<-list()
  Kclist[["curated"]]<-list(Amp10=3.666485,Amp=1.495666,Gains=1.096756,HetLoss=-1.887569,HomDel=-5.124848)
  Kclist[["balanced"]]<-list(Amp10=4.23526646153085,	Amp=1.93402122903165,	Gains=1.02725359043027,	HetLoss=-2.03950685716771,	HomDel=-4.11934749601642)
  match.arg(Kc)

  switch(Kc,
         "curated" = Kclist[["curated"]],
         "balanced" = Kclist[["balanced"]]
  )
  K<-Kclist[[Kc]]
  KN<-list()
  KN[as.character(0:1)]=seq(K[["HomDel"]],K[["HetLoss"]],length.out=2)
  KN[as.character(1:3)]=seq(K[["HetLoss"]],K[["Gains"]],length.out=3)

  KN[as.character(3:5)]=seq(K[["Gains"]],K[["Amp"]],length.out=3)
  KN[as.character(5:10)]=seq(K[["Amp"]],K[["Amp10"]],length.out=6)
  step=((KN[["10"]]-KN[["3"]])/7)
  KN[as.character(10+5)]=KN[[10]]+(step*5)
  KN[as.character(10+10)]=KN[[10]]+(step*10)
  K$Diploid<-KN[["2"]]
  return(list(K,KN))
}

get_Tresholds<-function(ss,ID,log2rfile_folder,Kc_method="balanced"){
  ss<-data.table::setDT(ss)
  data.table::setkey(ss,"Sample_Name")
  Kcs<-Kc(Kc_method)
  K<-Kcs[[1]]
  KN<-Kcs[[2]]
  # Load log2r file and calculate baseline & Var.
  std_log2rfile <- rlang::expr(paste0(folder=log2rfile_folder,ID,'_log2r.txt'))
  log2file <- check_input_files(Sample_Name = ss$Sample_Name,folder = log2rfile_folder, std_file = std_log2rfile)
  suppressWarnings(Log2<-data.table::fread(log2file,col.names = c("probeid","log2r")))
  baseline <- mean(Log2$log2r, na.rm=T)
  purity<-names(ss)[names(ss) %ilike% "impute" & names(ss) %ilike% "purity" ]
  p<-ss[,..purity]
  Var <- p*sd(Log2$log2r,na.rm=T)
  # feature value threshold
  Tresholds<-sapply(K, function(x) baseline + unname(x) * Var)
  Tresholds<-unlist(unlist(Tresholds))
  names(Tresholds)<-names(K)
  # Numeric value threshold
  TresholdsN<-sapply(KN, function(x) baseline + unname(x) * Var)
  TresholdsN<-unlist(unlist(TresholdsN))
  names(TresholdsN)<-names(KN)

  return(list(K,KN))
}

get_int<-function(seg,ref_genes="all",cn_genes,output_vars=c("UCSC_RefGene_Name"),arraytype="450K",anno=NULL){

  if(ref_genes=="all"){
      interest_geneset<- get_anno(anno = anno,x = arraytype)@probes
      interest_genes <- unique(do.call("c",sapply(interest_geneset$Name,function(g)unique(strsplit(g, split = ";")[[1]]))))
  }else if(ref_genes == "curated"){
      #data("CancerGenes")
      interest_geneset <- CancerGenes
      interest_genes<-interest_geneset$name
  }else{
      bp_stopifnot("annotation file must be either a GRanges object or a character string with value all/curated." = class(ref_genes)=="GRanges")
      interest_geneset<-ref_genes
      interest_genes<-interest_geneset$name
  }
  if(is.null(cn_genes))cn_genes<-interest_genes
  genes <- intersect(cn_genes,interest_genes)
  #print(genes)
  if(length(genes)<1){
    df<-data.frame(ID=NULL,Int=NULL,X=NULL,Var=NULL)
    warning("No genes of interest present in cnstate: ")
    return(df)
  }
  old_names<-c("chrom","loc.start","loc.end","seg.mean")
  new_names<-c("chr","start","end","log2r")
  data.table::setnames(seg, old_names,new_names,skip_absent = TRUE)
  data.table::setcolorder(seg,new_names)
  seggr <- GenomicRanges::makeGRangesFromDataFrame(seg,keep.extra.columns=TRUE)
  grset <- GenomicRanges::makeGRangesFromDataFrame(interest_geneset,keep.extra.columns=TRUE)
  fol <- suppressWarnings(SummarizedExperiment::findOverlaps(seggr,grset))
  ag <- sapply(unique(S4Vectors::queryHits(fol)),function(x){
     idx <- S4Vectors::subjectHits(fol)[S4Vectors::queryHits(fol) == x]
     Name=paste(intersect(unique(grset[idx,]$Name),genes),collapse = ";")
     #gsub("\\s", "", Name)
   })
  seggr.matched<-seggr[unique(S4Vectors::queryHits(fol))]
  seggr.matched$gene.name<-ag
  message("CN finish")
  return(seggr.matched)
}

.default.27k.annotation  <- "ilmn12.hg19"
.default.450k.annotation <- "ilmn12.hg19"
.default.epic.annotation <- "ilm10b4.hg19"
.default.allergy.annotation <- "ilm10.hg19"
.metharray.types <- c("IlluminaHumanMethylation450k",
                      "IlluminaHumanMethylationEPIC",
                      "IlluminaHumanMethylation27k",
                      "IlluminaHumanMethylationAllergy",
                      "HorvathMammalMethylChip40")

guessArrayTypes <- function(nProbes) {
  if (nProbes >= 622000 && nProbes <= 623000) {
    arrayAnnotation <- c(
      array = "IlluminaHumanMethylation450k",
      annotation = .default.450k.annotation)
  } else if (nProbes >= 1050000 && nProbes <= 1053000) {
    # NOTE: "Current EPIC scan type"
    arrayAnnotation <- c(
      array = "IlluminaHumanMethylationEPIC",
      annotation = .default.epic.annotation)
  } else if (nProbes >= 1032000 && nProbes <= 1033000) {
    # NOTE: "Old EPIC scan type"
    arrayAnnotation <- c(
      array = "IlluminaHumanMethylationEPIC",
      annotation = .default.epic.annotation)
  } else if (nProbes >= 54000 && nProbes <= 56000) {
    arrayAnnotation <- c(
      array = "IlluminaHumanMethylation27k",
      annotation = .default.27k.annotation)
  } else if (nProbes >= 41000 & nProbes <= 41100) {
    arrayAnnotation <- c(
      array = "HorvathMammalMethylChip40",
      annotation = "test.unknown")
  } else if (nProbes >= 43650 & nProbes <= 43680) {
    arrayAnnotation <- c(
      array = "IlluminaHumanMethylationAllergy",
      annotation = .default.allergy.annotation)
  } else {
    arrayAnnotation <- c(array = "Unknown", annotation = "Unknown")
    warning("Unable to detect Array type")
  }
  arrayAnnotation
}
get_anno <- function(x = c("450K", "EPIC", "overlap"),anno=NULL) {
  if(!is.null(anno)){

    out<-tryCatch(
      {
        bp_stopifnot("annotation file must be a CNV.anno object generated by conumee." = class(anno) == "CNV.anno")
        return(anno)
      },error=function(cond){
        message(cond)
        message("\n overlaping probes between 450K and epic will be used.")

        #load("data/anno_overlap.rda")
        return(anno_overlap)
      }
    )
  }else{
    out<-tryCatch(
      {  match.arg(x)
        switch(x,
               "EPIC" = anno_epic,
               "450K" = anno_450K,
               "overlap" = anno_overlap
        )
      },error=function(cond){
        message(cond)
        message("\n overlaping probes between 450K and epic will be used.")
        #load("data/anno_overlap.rda")
        return(anno_overlap)
      }
    )
    return(out)
  }
}

get_ncores<-function(cores=NULL){
  if(is.null(cores)){
    cores<-tryCatch({
      cores<-parallel::detectCores()-2
    if (!is.na(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))))cores=as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
    return(cores)
  },error=function(e)return(1)
  )
  }
  cores<-tryCatch(
    {
      if(!cores>=1)return(1)else return(cores)
    },error=function(e){
      return(1)
    }
  )
  return(cores)
}

