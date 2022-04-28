#' Generate segment & log2r intensity files calculated by CONUMEE.
#' @param ref_genes Granges object with subset of genes to use for cnv analysis.
#' @param targets A sample sheet with the required fields for minfi to load files.
#' @param Kc_method Choose which Kc you want to use. Either "curated" or "balanced.
#' The first comes from a list of 94 genes manually selected by experts and an
#' uneven proportion of cancers. The second one comes from all genes and balanced
#' proportion of cancer types. default= balanced
#' @param purity Whether or not you want purity imputed by Rfpurify.
#' Default=TRUE.
#' @param out Parent working directory where you want to have your results.

#' @param folder this directory will be created as the combination of out
#'  and subf. if you save the csv file inside this folder it will be used by
#'   minfi::read.450k.sheet.
#' @param RGset Whether you want normalised "RGChannelSet" or
#' "RGChannelSetExtended" to be saved or not. path=out
#' @param arraytype Methylation array type
#' @param Sample_Name Samples to be analysed. If is NULL all columns in input
#' file will be used. Accepts numbered index and names. Default=NULL
#' @param intensities dataset with intensities. either dataframe or path to file.
#' For big datasets '.fst' file is recomended.
#' @param anno_file anno file CNV.anno object saved as .rds. default
#' CNV_Germline_GSITIC2_BROAD_SNP6.merged.151117.hg19.CNV.txt
#' @param ctrl_file CNV data object with intensities from the controls group.
#' must be compressed as '.rds' file format, default uses 96 WB samples.
#' Default='WB'
#' @param conumee.folder Parent working directory where you want to have your results.
#' @param seg.folder  Subfolder inside results where segment files are saved.
#' The segments are generated with CONUMEE CNV.segment(CNV.detail(CNV.bin(fit)))
#' default = "Segments"
#' @param log2r.folder Subfolder inside results where log2r values are saved.
#' log2r values are the log2 ratio of the intensities(_GRN + _RED channel) between the
#' query and the reference set (control) as returned by CNV.fit function from CONUMEE.
#' default = "log2r"
#' @param probeid name of column with probe ids.
#' @param ncores number of cores to use.
#' @param cn_genes Used to subset genes from the reference set. Needed in order
#' to generate new Kc from a reference set, where cn_genes are previously known
#' genes to be altered in a given cn_state. Otherwise use ref_genes.


#' @return data.frame with metadata from segmentation, genomic ranges, genes
#' and scna for each of the genes
#' @export
#'
#' @examples
#' #data("TrainingSet_Sample_sheet")
#' #ss<-TrainingSet_Sample_sheet[60:200,]
#' t0<-Sys.time()
#' cnv<-cnv.methyl(targets=ss)
#' t1<-Sys.time()
#' cnv
#' print(paste0("elapsed time: ",t1-t0))
run_cnv.methyl<-function(
  targets,ref_genes="all",cn_genes=NULL, Kc_method, out="analysis/intermediate/",
  RGset=T, purity= NULL, arraytype="450K",folder=NULL, anno_file=NULL, ctrl_file='WB',
  Sample_Name=NULL,ncores=NULL, seg.folder = "Segments", log2r.folder = "log2r",
  conumee.folder="analysis/CONUMEE/", probeid="probeid"){

  anno<-get_anno(anno_file,arraytype=arraytype)
  intensity<-pre_process(targets = targets, purity=purity, RGset = RGset,
                         out=out,folder=folder,arraytype=arraytype)
  ss<-data.table::fread(paste0(out,"Sample_Sheet.txt"))
  run_conumee(intensities = intensity,anno_file=anno, ctrl_file=ctrl_file,
              Sample_Name=Sample_Name,seg.folder = seg.folder,
              log2r.folder = log2r.folder,arraytype=arraytype,
              conumee.folder=conumee.folder, probeid=probeid)
  # Cluster:

  cl<- parallel::makePSOCKcluster(get_ncores(ncores), outfile="")
  doParallel::registerDoParallel(cl)
  res<-foreach::foreach(i=1:length(ss$Sample_Name),#isplitIndices(1400,chunks=ncores),
                        .combine='c',
                        .multicombine = F,
                        .inorder=F,
                        .packages = "data.table",
                        .export = c("AllGenes","CancerGenes"),
                        .errorhandling = "pass"
  )%dopar%{
    load("R/sysdata.rda")

    source("R/utils.R")
    source("R/Kc_get.R")
    #source()
     ss<-ss[i,]
     ID<-ss$Sample_Name
     cna <- Kc_get(ss=ss,ID=ID,ref_genes = ref_genes,arraytype=arraytype,anno=anno)
    # cna
     return(cna)
  }
  parallel::stopCluster(cl)
  return(res)
}
