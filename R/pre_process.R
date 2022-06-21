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
#' @param files Files
#' @param RGset Whether you want normalised "RGChannelSet" or
#' "RGChannelSetExtended" to be saved or not. path=out
#' @param arraytype Methylation array type
#' @param idats_folder Folder containing idats. It will try to guess which idat
#' belongs to which row in the targets data.frame.
#' @param copy Change to TRUE if you want to copy files to idats folder.
#' Useful when working in HPC and idats are in ISILON.
#' @inheritParams run_cnv.methyl
#' @return Log2r intensities ready to be analysed with conumee.
#' @export
#'
#' @examples
#'
#' data("TrainingSet_Sample_sheet")
#' ss<-TrainingSet_Sample_sheet[1:56,]
#' pre_process(ss)

pre_process<-function(targets,purity=NULL,query=T,RGset=T,out="./analysis/intermediate/",idats_folder=NULL,subf="IDATS/",folder=NULL,ncores=NULL,arraytype="450K",copy=FALSE){
  #data('hm450.manifest.hg19')
  if (!is.character(folder)) folder<-paste0(out,subf)
  dir.create(folder,recursive=TRUE)
  if (!is.null(targets)) {
    if (!"Basename" %in% names(targets)) {
      warning("Need 'Basename' amongst the column names of 'targets', will construct from idats.folder")
      if(is.null(idats_folder)){
        stop("Either provide Basename or a valid idats.foler in order to continue" )
      }else{
        base<-idats_folder
        Grn.files <- list.files(path = base, pattern = "_Grn.idat$",
                                recursive = T, ignore.case = TRUE, full.names = TRUE)
        Red.files <- list.files(path = base, pattern = "_Red.idat$",
                                recursive = T, ignore.case = TRUE, full.names = TRUE)
        if (length(Grn.files) == 0 || length(Red.files) == 0) {
          stop("No IDAT files were found")
        }
        commonFiles <- intersect(x = sub("_Grn.idat$", "", Grn.files),
                                 y = sub("_Red.idat$", "", Red.files))
        if (length(commonFiles) == 0) {
          stop("No IDAT files with both Red and Green channel were found")
        }
        commonFiles.Grn <- paste(commonFiles, "_Grn.idat", sep = "")
        commonFiles.Red <- paste(commonFiles, "_Red.idat", sep = "")
        if (!setequal(commonFiles.Grn, Grn.files)) {
          warning(sprintf("the following files only exists for the green channel: %s",
                          paste(setdiff(Grn.files, commonFiles.Grn), collapse = ", ")))
        }
        if (!setequal(commonFiles.Red, Red.files)) {
          warning(sprintf("the following files only exists for the red channel: %s",
                          paste(setdiff(Red.files, commonFiles.Red), collapse = ", ")))
        }

        # Try to figure out which file belongs to which sample
        best<-NULL
        nbest<- round(dim(targets)[1]*0.9)
        counts<-sapply(names(targets),function(column){
          sum(sapply(basename(commonFiles),function(s)  targets[,..column] %like% s))
        })
        if(max(counts)>nbest){
          best<- which.max(counts)
        }else{
          warning(paste0("only ",counts, "out of ", nrow(targets), "samples in sample sheet matched" ))
          best<- which.max(counts)
        }
        targets[sapply(basename(commonFiles),function(sample){which(sapply(targets[,..best],function(x) x %like% sample))}),Basename:=commonFiles]

      }
    }
    files<-targets$Basename

    }else{
    stop("Must provide a 'targets' data.frame in order to continue")
  }


  myLoad<-pre_process.myLoad(targets,folder=folder, files=files, arraytype=arraytype,copy=copy, ncores=ncores)
  saveRDS(myLoad,paste0(out,"myLoad.rds"))
  #targets<-SummarizedExperiment::colData(myLoad)
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

#' construct RGChannelSet in parallel using foreach
#'
#'
#' @title construct RGChannelSet in parallel
#'

#' @inheritParams run_cnv.methyl
#' @inheritParams pre_process
#' @return RGset
#' @author izar de Villasante
#' @export
#' @import minfi
#' @importFrom utils str
#' @importFrom BiocGenerics cbind combine
#' @importFrom foreach foreach
#' @importFrom minfi read.metharray
#' @inheritParams minfi::read.metharray






read.metharray.exp.par <- function(targets,folder,files,copy=FALSE, verbose = TRUE, arraytype="450K",ncores=NULL, extended = FALSE, force =FALSE){
  #Make cluster:
  ncores<-get_ncores(ncores)
  message("Reading multiple idat-files in parallel. Using ",ncores," cores.")
  cl<- parallel::makePSOCKcluster(ncores)
  doParallel::registerDoParallel(cl)
  if (copy ==TRUE){
    files<-paste0(folder,basename(targets$Basename))
  }
  requireNamespace("S4Vectors")
  res<-foreach::foreach(it=itertools::isplitIndices(nrow(targets), chunks=ncores),#isplitIndices(1400,chunks=ncores),
                        .combine='cbind',
                        .multicombine = F,
                        .inorder=F,
                        .errorhandling = "pass"
  )%dopar%{
    requireNamespace(c("minfi","fs"))
    subdf<-as.data.frame(targets[it,])
    if (copy ==TRUE){
      fs::file_copy(paste0(subdf$filenames,"_Grn.idat"),new_path=folder,overwrite = T)
      fs::file_copy(paste0(subdf$filenames,"_Red.idat"),new_path=folder,overwrite=T)

    }

    rgSet<-minfi::read.metharray(basenames = files[it], extended = extended, verbose = verbose, force =force)
    if(arraytype=="EPIC") rgSet@annotation <- c(array="IlluminaHumanMethylationEPIC",annotation="ilm10b4.hg19")
    return(rgSet)
  }
  parallel::stopCluster(cl)
  pD <- data.frame(targets)
  pD$filenames <- files
  rownames(pD) <- colnames(res)
  res@colData <- methods::as(pD, "DataFrame")
  #rownames(res@colData)<-colnames(res)
  return(res)
}

#' construct RGChannelSet in parallel using foreach
#'
#'
#' @title construct RGChannelSet in parallel
#' @rdname pre_process
# #' @param ... optional arguments to read.metharray.exp
#' @inheritParams run_cnv.methyl
#' @return RGset
#' @author izar de Villasante
#' @export
pre_process.myLoad <-function(targets,folder,arraytype="450K",ncores=NULL, files=NULL, copy=FALSE) {
  message("working directory: ",folder)
  file.remove(list.files(folder, full.names = TRUE))
  #idats_folder(targets,folder = folder)
  utils::write.csv(targets,paste0(folder,"/Sample_sheet.csv"))
  myLoad <- read.metharray.exp.par(targets=targets,files=files,folder = folder,arraytype=arraytype,ncores=ncores,copy=copy)
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
  betas <- betas[
    match(rownames(RFpurify_ABSOLUTE$importance), rownames(betas))
    ,
    , drop = FALSE
  ]
  betas<-tryCatch(
    {betas<-impute::impute.knn(data = betas, k=knn)$data},
    error={ betas[is.na(betas)]<-mean(betas,na.rm=T)}

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

