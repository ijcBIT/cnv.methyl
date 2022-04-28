#' Set of functions to obtain Kc from reference copy number state.
#' @param ID ID inside Sample sheet and in the Sample_Name.
#' @param ss sample sheet
#' @author Izar de Villasante
#' @inheritParams run_cnv.methyl
#' @return GRanges object of the genes of interest "ref_genes" and a cna column
#' with the copy number alterations for each of them.
#' @export
#' @examples
#' #data("TrainingSet_Sample_sheet")
#' #ss<-TrainingSet_Sample_sheet[60:200,]
#' ID="TCGA-19-A6J4-01A-11D-A33U-05"
#' ss[Sample_Name==ID,]->ss
#' cna<-Kc_get(ss=ss,ID="TCGA-19-A6J4-01A-11D-A33U-05")
#' cna




Kc_get<-function(
  ID,ss,ref_genes="all",conumee.folder="analysis/CONUMEE",seg.folder = "Segments",
  log2r.folder = "log2r",arraytype="450K",anno_file=NULL,cn_genes=NULL,Kc_method="balanced"
  ){
  Kcs<-Kc(Kc_method)
  K<-Kcs[[1]]
  KN<-Kcs[[2]]

  anno<-get_anno(anno_file,arraytype)
  #message(samp$Sample_Name)
  # Load complete log2r of whole array:
  ss<-data.table::setDT(ss)
  data.table::setkey(ss,"Sample_Name")
  segfile_folder <- paste0(conumee.folder,"/",seg.folder)
  std_segfile <- rlang::expr(paste0(folder,ID,'_Segments.txt'))#rlang::expr(paste0(folder,"/","Segments_",ID,'.txt'))
  segfile <- check_input_files(Sample_Name = ss$Sample_Name,folder = segfile_folder, std_file = std_segfile)

  log2rfile_folder <- paste0(conumee.folder,"/",log2r.folder)
  std_log2rfile <- rlang::expr(paste0(folder,ID,'_log2r.txt'))
  log2file <- check_input_files(Sample_Name = ss$Sample_Name,folder = log2rfile_folder, std_file = std_log2rfile)

  suppressWarnings(Log2<-data.table::fread(log2file,col.names = c("probeid","log2r")))
  baseline <- mean(Log2$log2r, na.rm=T)
  purity<-names(ss)[names(ss) %ilike% "impute" & names(ss) %ilike% "purity" ]
  #p<-with(ss,get(purity))
  p<-ss[,..purity]
  Var <- p*sd(Log2$log2r,na.rm=T)

  #K$Diploid<-unname(unlist(K[names(which.max(K[K<0]))])+0.00001)
  K$Diploid<-KN[["2"]]
  Tresholds<-sapply(K, function(x) baseline + unname(x) * Var)
  Tresholds<-unlist(unlist(Tresholds))
  names(Tresholds)<-names(K)

  TresholdsN<-sapply(KN, function(x) baseline + unname(x) * Var)
  TresholdsN<-unlist(unlist(TresholdsN))
  names(TresholdsN)<-names(KN)

  seg<-data.table::fread(segfile)
  cnlist<-apply(seg,1, function(x){

    rest <- as.numeric(x["seg.mean"]) - Tresholds
    restN <- as.numeric(x["seg.mean"]) - TresholdsN
    cna <- names(which.min(rest[rest>0]))
    cnv <- names(which.min(restN[restN>0]))
    if(is.null(cna)){ # In the case all Thresholds are positive and seg.mean is below
      cna<-names(which.max(rest))
      cnv<-names(which.max(restN))
    }
    return(c(cna,cnv))
  })
  seg$cna<-cnlist[1,]
  seg$cnv<-cnlist[2,]
  message("segment and log2 files loaded.")

  int<-get_int(seg,ref_genes=ref_genes,cn_genes=cn_genes,anno=anno)
  int
}
