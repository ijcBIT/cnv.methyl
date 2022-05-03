#' Pre-process methylation array data from (raw) .idat files to segment & log2r intensity files calculated by CONUMEE.
#'
#' @param targets A sample sheet with the required fields for ChAMP to load files.
#' @param purity Whether or not you want purity imputed by Rfpurify.
#' Default=TRUE.
#' @param query Set to TRUE when you are only interested in imputing purity.
#' Default=TRUE
#' @param out Parent working directory where you want to have your results.
#' @param subf  Subfolder inside results folder where IDAT files will be copied.
#' @param folder this directory will be created as the combination of out
#'  and subf. if you save the csv file inside this folder it will be used by
#'   minfi::read.450k.sheet.
#' @param RGset Whether you want normalised "RGChannelSet" or
#' "RGChannelSetExtended" to be saved or not. path=out
#' @param arraytype Methylation array type
#' @inheritParams run_cnv.methyl
#' @return Log2r intensities ready to be analysed with conumee.
#' @export
#'
#' @examples
#'
#' data("TrainingSet_Sample_sheet")
#' ss<-TrainingSet_Sample_sheet[1:56,]
#' pre_process(ss)

pre_process<-function(targets,purity=NULL,query=T,RGset=T,out="analysis/intermediate/",subf="IDATS/",folder=NULL,ncores=NULL,arraytype="450K"){
  #data('hm450.manifest.hg19')
  if (!is.character(folder)) folder<-paste0(out,subf)
  dir.create(folder,recursive=TRUE)
  targets$Basename<-paste0(folder, basename(targets$filenames))
  myLoad<-pre_process.myLoad(targets,folder=folder, arraytype=arraytype,nocres=ncores)
  saveRDS(myLoad,paste0(out,"myLoad.rds"))
  targets<-SummarizedExperiment::colData(myLoad)
  if(is.null(purity)){
    purity<-purify(myLoad=myLoad)
    message("Rfpurifies")
    myLoad@colData$Purity_Impute_RFPurify<-purity
    #targets$`Purity_Impute_RFPurify(Absolute)` <- purity
    utils::write.table(myLoad@colData,paste(out,"Sample_Sheet.txt",sep="/"), col.names = T, row.names = F, quote = F, sep="\t")
    message(paste("targets is saved ",out))
  }
  if(query==T){
    query <- queryfy(myLoad)
    if (RGset==T){
      saveRDS(query,paste0(out,"/intensities.rds"),compress = FALSE)
    }
    message("Pre-processing Completed successfully!")
    return(query)
  }

}

#' DEPRECATED!
#' prepare idats folder and targets csv in parallel using foreach
#' Now part of read.metharray.exp.par
#'
#' @title generate champ folder with all idats
#' @param targets data.frame representing valid targets file as is
#' @param ncores cores
#' @param folder workdir where to copy files
#' used within the minfi package
#' @return folder with idats
#' @author izar de Villasante
#' @export
#' @importFrom fs file_copy
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
idats_folder <- function(targets,folder,ncores=NULL){
  ncores<-get_ncores(ncores)
  cl<- parallel::makePSOCKcluster(ncores)
  doParallel::registerDoParallel(cl)
  res<-foreach::foreach(it=itertools::isplitIndices(nrow(targets), chunks=ncores),#isplitIndices(1400,chunks=ncores),
                        .combine='c',
                        .multicombine = F,
                        .inorder=F,
                        .errorhandling = "pass"
  )%dopar%{
    requireNamespace("fs")
    subdf<-as.data.frame(targets[it,])
    fs::file_copy(paste0(subdf$filenames,"_Grn.idat"),new_path=folder,overwrite = T)
    fs::file_copy(paste0(subdf$filenames,"_Red.idat"),new_path=folder,overwrite=T)

  }
  parallel::stopCluster(cl)
}

#' construct RGChannelSet in parallel using foreach
#'
#'
#' @title construct RGChannelSet in parallel
#'
#' @param verbose logical default TRUE
#' @param ... optional arguments to read.metharray.exp
#' @inheritParams run_cnv.methyl
#' @return RGset
#' @author izar de Villasante
#' @export
#' @import minfi
#' @importFrom utils str
#' @importFrom BiocGenerics cbind combine
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%




read.metharray.exp.par <- function(targets,folder, verbose = TRUE,arraytype="450K",ncores=NULL, ...) {
  #Make cluster:
  ncores<-get_ncores(ncores)
  message("Reading multiple idat-files in parallel. Using ",ncores," cores.")
  cl<- parallel::makePSOCKcluster(ncores)
  doParallel::registerDoParallel(cl)
  res<-foreach::foreach(it=itertools::isplitIndices(nrow(targets), chunks=ncores),#isplitIndices(1400,chunks=ncores),
               .combine='cbind',
               .multicombine = F,
               .inorder=F,
               .export= ...,
               .errorhandling = "pass"
  )%dopar%{
    requireNamespace("minfi")
    requireNamespace("fs")
    subdf<-as.data.frame(targets[it,])
    fs::file_copy(paste0(subdf$filenames,"_Grn.idat"),new_path=folder,overwrite = T)
    fs::file_copy(paste0(subdf$filenames,"_Red.idat"),new_path=folder,overwrite=T)
    rgSet<-minfi::read.metharray.exp(targets = subdf, base=folder, ...)
    if(arraytype=="EPIC") rgSet@annotation <- c(array="IlluminaHumanMethylationEPIC",annotation="ilm10b4.hg19")
    return(rgSet)
  }
  parallel::stopCluster(cl)
  return(res)
}

#' construct RGChannelSet in parallel using foreach
#'
#'
#' @title construct RGChannelSet in parallel
#' @rdname pre_process
#' @param ... optional arguments to read.metharray.exp
#' @inheritParams run_cnv.methyl
#' @return RGset
#' @author izar de Villasante
#' @export
pre_process.myLoad <-function(targets,folder,arraytype="450K",ncores=NULL, ...) {
  message("working directory: ",folder)
  file.remove(list.files(folder, full.names = TRUE))
  #idats_folder(targets,folder = folder)
  utils::write.csv(targets,paste0(folder,"/Sample_sheet.csv"))
  myLoad <- read.metharray.exp.par(targets=targets,folder = folder,arraytype=arraytype,ncores=ncores)
  message("Minfi does load data")
  return(myLoad)
}

#' @export
#' @rdname pre_process
#' @param myLoad a RGset as returned by minfi::read.metharray.exp. will use this as basedir to
#' load the idats. default is results/IDATS/.
#' @param knn number of neighbors for knn algorithm

purify <- function(myLoad,knn=5){
  requireNamespace("randomForest")
  sink(tempfile())
  on.exit(sink())

  betas<- minfi::getBeta(myLoad)
  betas<-tryCatch(
    {betas<-impute::impute.knn(data = betas, k=knn)$data},
    error={ betas<-betas[stats::complete.cases(betas),]},
    finally ={
      betas <- betas[
        match(rownames(RFpurify_ABSOLUTE$importance), rownames(betas))
        ,
        , drop = FALSE
        ]
      }
    )
  absolute<-stats::predict(RFpurify_ABSOLUTE, t(betas))
  return(absolute)
}

#' @export
#' @rdname pre_process
#' @param out Results folder to store the normalized and clean RGset.
#' Default = NULL

queryfy<-function(myLoad){
  requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  mSetSqn <- tryCatch( minfi::preprocessQuantile(myLoad),error=function(e) {message("failed")
    return(NULL)})
  ##Percentage of faulty probes: more than 10% we throw samples away.
  detP <- minfi::detectionP(myLoad)
  bad <- colnames(detP)[colSums(detP >=0.01)/nrow(detP) > 0.1]

  ## Ensure probes are in the same order in the mSetSqn and detP objects
  detP <- detP[match(Biobase::featureNames(mSetSqn),rownames(detP)),]
  ## Remove rows with at elast one 'NA' entry
  keep <- rowSums(detP < 0.01) == ncol(mSetSqn)
  mSetSqn <- mSetSqn[keep,]
  ## Remove probes with SNPs at CpG site
  mSetSqn <- dropLociWithSnps(mSetSqn)
  mSetSqn <-  maxprobes::dropXreactiveLoci(mSetSqn)
  return(mSetSqn)
  # mSetSqn <- CNV.load(mSetSqn)
  # return(mSetSqn@intensity)
}


