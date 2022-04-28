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

get_anno<-function(anno_file=NULL,arraytype){
  if(!is.null(anno_file )){anno <- anno_file
  # }else if (arraytype == "450K"){anno<-anno_450K
  # }else if (arraytype == "EPIC"){anno<-anno_epic
  # }else anno <- anno_overlap
  }else{
  anno<-switch(arraytype,
               "EPIC" = anno_epic,
               "450K" = anno_450K,
               "overlap" = anno_overlap,
               anno_overlap
  )}
}

get_int<-function(seg,ref_genes="all",cn_genes,output_vars=c("UCSC_RefGene_Name"),arraytype="450K",anno=NULL){

  #grset <- fdata(arraytype)


  if(ref_genes=="all"){
    if(!is.null(anno)){
      interest_geneset<-tryCatch(
        {anno@probes
        },error=function(cond){
          message(cond)
          message("\n overlaping probes between 450K and epic will be used.")
          #data("anno_overlap")
          return(anno_overlap)
          #load("data/anno_overlap.rda")
        }
      )
      interest_genes <- unique(do.call("c",sapply(interest_geneset$Name,function(g)unique(strsplit(g, split = ";")[[1]]))))

    }else{

    interest_geneset <- fdata(arraytype)
    interest_genes <- unique(do.call("c",sapply(interest_geneset$Name,function(g)unique(strsplit(g, split = ";")[[1]]))))
    }

  }else if(ref_genes == "paper"){
      #data("CancerGenes")
      interest_geneset <- CancerGenes
      interest_genes<-interest_geneset$name
    }
  else{
      interest_geneset<-utils::read.table(ref_genes)
      interest_genes<-interest_geneset$name
  }
  if(is.null(cn_genes))cn_genes<-interest_genes
  genes <- intersect(cn_genes,interest_genes)
  if(length(genes)<1){
    df<-data.frame(ID=NULL,Int=NULL,X=NULL,Var=NULL)
    warning("No genes of interest present in cnstate: ")
    return(df)
  }
  old_names<-c("chrom","loc.start","loc.end","seg.mean")
  new_names<-c("chr","start","end","log2r")
  data.table::setnames(seg, old_names,new_names)
  data.table::setcolorder(seg,new_names)
  #cols<-colSums(is.na(seg))<1
  #seg<-seg[,.SD,.SDcols=cols]


  seggr <- regioneR::toGRanges(seg)
  grset <- regioneR::toGRanges(interest_geneset)
  fol <- suppressWarnings(SummarizedExperiment::findOverlaps(seggr,grset))
  ag <- sapply(unique(S4Vectors::queryHits(fol)),function(x){
     idx <- S4Vectors::subjectHits(fol)[S4Vectors::queryHits(fol) == x]
     Name=paste(intersect(unique(grset[idx,]$Name),genes),collapse = ";")
     #gsub("\\s", "", Name)
   })
  seggr.matched<-seggr[unique(S4Vectors::queryHits(fol))]
  seggr.matched$gene.mame<-ag

  #   seg[,var]
  #   for (var in output_vars){
  #     var <- grset[idx,var]
  #
  #   }
  #   genes<-sapply(RefGene, function(g) strsplit(g, split = ";")[[1]])
  #   genes<-unique(unlist(genes))
  #   return(genes)
  # })
  #mol<- suppressWarnings(SummarizedExperiment::findOverlaps(seggr,grset))
  # seggr.matched <- seggr[S4Vectors::queryHits(fol)];
  # S4Vectors::mcols(seggr.matched) <- cbind.data.frame(
  #   S4Vectors::mcols(seggr.matched),
  #   S4Vectors::mcols(grset[S4Vectors::subjectHits(fol)]));
  # int <- seggr.matched
  message("int finish")
  return(seggr.matched)
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
  return(list(K,KN))
}

fdata <- function(x = c("450K", "EPIC", "overlap"),anno=NULL) {
  if(!is.null(anno)){
    out<-tryCatch(
      {anno@probes
      },error=function(cond){
        message(cond)
        message("\n overlaping probes between 450K and epic will be used.")
        load("data/anno_overlap.rda")
        return(anno_overlap@probes)
      }
    )
  }
  out<-tryCatch(
    {  match.arg(x)

      switch(x,
             "EPIC" = anno_epic@probes,
             "450K" = anno_450K@probes,
             "overlap" = anno_overlap@probes
      )
    },error=function(cond){
      message(cond)
      message("\n overlaping probes between 450K and epic will be used.")
      load("data/anno_overlap.rda")
      return(anno_overlap@probes)
    }
  )
  return(out)

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
