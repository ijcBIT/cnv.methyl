#' Set of functions to obtain Kc from reference copy number state.
#' @param ID ID inside Sample sheet and in the Sample_Name.
#' @param ss sample sheet
#' @param segfile_folder Path to folder with segmentation files.
#' @param log2rfile_folder Path to folder with intensities files.
#' @author Izar de Villasante
#' @inheritParams run_cnv.methyl
#' @return GRanges object of the genes of interest "ref_genes" and a cna column
#' with the copy number alterations for each of them.
#' @export
#' @examples
#'
#' data("TrainingSet_Sample_sheet")
#' ss<-TrainingSet_Sample_sheet[1:56,]
#' ID="TCGA-19-A6J4-01A-11D-A33U-05"
#' ss[Sample_Name==ID,]->ss
#' cna<-Kc_get(ss=ss,ID="TCGA-19-A6J4-01A-11D-A33U-05")
#' cna



Kc_get<-function(
  ID,ss,ref_genes="all",segfile_folder = "analysis/CONUMEE/Segments",
  log2rfile_folder = "analysis/CONUMEE/log2r",arraytype=NULL,anno_file=NULL,cn_genes=NULL,Kc_method="balanced"
  ){
  # Prepare data:
  anno<-get_anno(anno = anno_file, x = arraytype)
  ss<-data.table::setDT(ss)
  data.table::setkey(ss,"Sample_Name")
  # Calculate Thrersholds:
  Tre_list<-get_Tresholds(ss = ss,ID = ID, log2rfile_folder = log2rfile_folder, Kc_method=Kc_method)
  Tresholds<-unlist(Tre_list[[1]])
  TresholdsN<-unlist(Tre_list[[2]])
  # Load segment file:
  std_segfile <- rlang::expr(paste0(folder = segfile_folder,ID,'_Segments.txt'))#rlang::expr(paste0(folder,"/","Segments_",ID,'.txt'))
  segfile <- check_input_files(Sample_Name = ID,folder = segfile_folder, std_file = std_segfile)
  seg<-data.table::fread(segfile)
  #Calculate cn & cna
  cnlist<-apply(seg,1, function(x){
    rest <- as.numeric(x["seg.mean"]) - Tresholds
    restN <- as.numeric(x["seg.mean"]) - TresholdsN
    cna <- names(which.min(rest[rest>0]))
    cn <- names(which.min(restN[restN>0]))
    if(is.null(cna)){ # In the case all Thresholds are positive and seg.mean is below
      cna<-names(which.max(rest))
      cn<-names(which.max(restN))
    }
    return(c(cna,cn))
  })
  seg$cna<-cnlist[1,]
  seg$cn<-cnlist[2,]
  message("segment and log2 files loaded.")
  # Calculate GRanges with matching
  int<-get_int(seg,ref_genes=ref_genes,cn_genes=cn_genes,anno=anno)
  return(int)
}
