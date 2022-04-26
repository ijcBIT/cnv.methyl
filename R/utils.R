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

get_int<-function(seg,ref_genes="all",cn_genes,output_vars=c("UCSC_RefGene_Name"),arraytype="450K"){

  #grset <- fdata(arraytype)

  if(ref_genes=="all"){

    interest_geneset <- fdata(arraytype)
    interest_genes <- unique(do.call("c",sapply(grset$Name,function(g)unique(strsplit(g, split = ";")[[1]]))))
    }
  else if(ref_genes == "paper"){
      data("CancerGenes")
      interest_geneset <- CancerGenes
      interest_genes<-interest_geneset$name
    }
  else{
      interest_geneset<-utils::read.table(ref_genes)
      interest_genes<-interest_geneset$name
      }
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
  }
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
  return(int)
}


fdata <- function(x = c("450K", "EPIC", "overlap")) {
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
      return(anno_overlap@probes)
    }
  )
  return(out)

}


