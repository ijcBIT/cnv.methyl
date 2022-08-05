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
#' @param sampGroups  Groups for qc report & plots color category.
#' @inheritParams run_cnv.methyl
#' @return Log2r intensities ready to be analysed with conumee.
#' @export
#'
#' @examples
#'
#' data("TrainingSet_Sample_sheet")
#' ss<-TrainingSet_Sample_sheet[1:56,]
#' pre_process(ss)

pre_process<-function(targets,purity=NULL,query=T,RGset=T,out="./analysis/intermediate/",
                      idats_folder=NULL,subf="IDATS/",folder=NULL,ncores=NULL,arraytype=NULL,
                      copy=FALSE,frac=0.1,pval=0.01,remove_sex=TRUE,sampGroups=NULL){
  #data('hm450.manifest.hg19')
  #targets <- targets[,colSums(is.na(targets))<nrow(targets)]
  if (!is.character(folder)) folder<-paste0(out,subf)
  dir.create(folder,recursive=TRUE)
  if (!is.null(targets)) {
    if (!"Sample_Name" %in% names(targets)) {
      stop("targets must contain Sample_Names column with unique identifier of each sample")
    }else{
      if(any(duplicated(targets$Sample_Name)))stop("Sample_Name can't be repeated")
    }

    if (!"Basename" %in% names(targets)) {
      warning("Need 'Basename' amongst the column names of 'targets', will construct from idats_folder")
      if(is.null(idats_folder)){
        stop("Either provide Basename or a valid idats_folder in order to continue" )
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

  # Read idats and generate rgSet with metadata:
  myLoad<-pre_process.myLoad(targets,folder=folder, files=files, arraytype=arraytype,copy=copy, ncores=ncores)
  saveRDS(myLoad,paste0(out,"myLoad.rds"))
  # Qc:
  qc_folder<-paste0(out,"/QC/")
  dir.create(normalizePath(qc_folder),showWarnings = F)
  if(is.null(sampGroups)){
    g<-  intersect(c("Sample_Group","Group","Type","Condition","Category"), colnames(targets))
    g<-g[1]
  }else{
    g<-sampGroups
  }
  qc(rgSet = myLoad,sampGroups = g,sampNames = "Sample_Name")

  if(is.null(purity)){
    purity<-purify(myLoad=myLoad)
    message("Rfpurifies")
    myLoad@colData$Purity_Impute_RFPurify<-purity
    #targets$`Purity_Impute_RFPurify(Absolute)` <- purity
    utils::write.table(myLoad@colData,paste(out,"Sample_Sheet.txt",sep="/"), col.names = T, row.names = F, quote = F, sep="\t")
    message(paste("targets is saved ",out))
  }
  if(query==T){
    query <- queryfy(targets = targets,myLoad,arraytype=arraytype,frac=frac,pval=pval,remove_sex=remove_sex,qc_folder=qc_folder)
    if (RGset==T){
      saveRDS(query,paste0(out,"/intensities.rds"),compress = FALSE)
    }
    message("Pre-processing Completed successfully!")
    return(query)
  }

}

#' construct RGChannelSet in parallel using foreach
#'
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


read.metharray.exp.par <- function(targets,folder,files=NULL,copy=FALSE, verbose = TRUE, arraytype = NULL, ncores=NULL, extended = FALSE, force = TRUE){
  #Make cluster:
  if(is.null(files)) files<- targets$Basename
  if (copy ==TRUE){
    files<-paste0(folder,basename(targets$Basename))
  }


  ncores<-get_ncores(ncores)

  cl<- parallel::makePSOCKcluster(ncores)#,outfile="")
  parallel::clusterEvalQ(cl,{
    requireNamespace(c("minfi","S4Vectors"))
  })
  doParallel::registerDoParallel(cl)

  message("Reading multiple idat-files in parallel. Using ",ncores," cores.")
  res<-foreach::foreach(it=itertools::isplitIndices(nrow(targets), chunks=ncores),#isplitIndices(1400,chunks=ncores),
                        .combine='cbind',
                        .multicombine = F,.export = c("guessArrayTypes",".default.epic.annotation"),
                        .inorder=F,
                        .errorhandling = "pass"
  )%dopar%{

    subdf<-as.data.frame(targets[it,])
    # handle idats path
    if (copy ==TRUE){
      requireNamespace("fs")
      fs::file_copy(paste0(subdf$Basename,"_Grn.idat"),new_path=folder,overwrite = T)
      fs::file_copy(paste0(subdf$Basename,"_Red.idat"),new_path=folder,overwrite=T)
    }

    # read idats:
    rgSet<-minfi::read.metharray(basenames = files[it], extended = extended, verbose = verbose, force =force)

    # arraytype:
    if (is.null(arraytype)){
      rgSet@annotation<-guessArrayTypes(nrow(rgSet))
    }else{
      if (arraytype=="EPIC") {rgSet@annotation <- c(array="IlluminaHumanMethylationEPIC",annotation="ilm10b4.hg19")}
      else if (arraytype=="450K"){rgSet@annotation <- c(array="IlluminaHumanMethylation450k",annotation="ilmn12.hg19")
      }else{rgSet@annotation<-guessArrayTypes(nrow(rgSet))}
    }
    return(rgSet)
  }

  parallel::stopCluster(cl)
  cn<-colnames(res)
  class(res)
  data.table::setkey(targets,"barcode")
  pD <- data.frame(targets[cn,])
  pD$filenames <- files
  #rownames(pD) <- colnames(res)
  res@colData <- methods::as(pD, "DataFrame")
  rownames(res@colData)<-cn
  colnames(res)<-cn

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
  #file.remove(list.files(folder, full.names = TRUE))
  #idats_folder(targets,folder = folder)
  utils::write.csv(targets,paste0(folder,"/Sample_sheet.csv"))
  myLoad <- read.metharray.exp.par(targets=targets,files=files,folder = folder,arraytype=arraytype,ncores=ncores,copy=copy)
  message("Minfi does load data")

  return(myLoad)
}

#' Generate Qc reports
#'
#'
#' @title construct RGChannelSet in parallel
#' @rdname pre_process
#' @param qc_folder Path to location where qc_report and plot will be saved.
#' @inheritParams minfi::qcReport
#' @return Qc plots
#' @author izar de Villasante
#' @export
#'
qc <- function(rgSet,sampGroups=NULL, sampNames= "Sample_Name",qc_folder="analysis/intermediate/QC/"){
  if (!(sampNames %in% names(rgSet@colData))) sampNames <- colnames(rgSet)
  n <- ncol(rgSet)
  o <- rev(order(sampNames))
  rgSet <- rgSet[, o]
  sampNames <- sampNames[o]
  if (is.null(sampGroups))
    sampGroups <- rep(1, n)
  sampGroups <- sampGroups[o]
  dir.create(qc_folder,recursive=T,showWarnings = F)
  minfi::qcReport(rgSet = rgSet,
                  pdf = paste0(qc_folder,"Report.pdf"),
                  sampGroups = rgSet@colData[[sampGroups]],
                  sampNames = rgSet@colData[[sampNames]])
  if(length(unique(rgSet@colData[[sampGroups]]))>8){
    grDevices::pdf(file = paste0(qc_folder,"density_plot.pdf"),
                   width = 9, # The width of the plot in inches
                   height = 11) # The height of the plot in inches
    minfi::densityPlot(
      rgSet, sampGroups = rgSet@colData[[sampGroups]],main = "Beta",
      pal =  c(
        "#191919", "#F0A0FF", "#0075DC", "#993F00", "#005C31", "#5EF1F2", "#FF0010",
        "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
        "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
        "#740AFF", "#990000", "#FFFF80", "#FFE100", "#FF5005")
    )
    grDevices::dev.off()()

  }
  mSet <- minfi::preprocessRaw(rgSet)
  qc   <- minfi::getQC(mSet)
  grDevices::pdf(file = paste0(qc_folder,"mean_qc.pdf"),   # The directory you want to save the file in
                 width = 7, # The width of the plot in inches
                 height = 7) # The height of the plot in inches
  minfi::plotQC(qc)

  grDevices::dev.off()()
}
#' @export
#' @rdname pre_process
#' @param myLoad a RGset as returned by minfi::read.metharray.exp. will use this as basedir to
#' load the idats. default is results/IDATS/.
#' @param knn number of neighbors for knn algorithm

purify <- function(myLoad,knn=5){
  requireNamespace("randomForest")
  #RFpurify_ABSOLUTE<-cnv.methyl:::RFpurify_ABSOLUTE
  sink(tempfile())
  on.exit(sink())

  betas<- minfi::getBeta(myLoad)
  betas <- betas[
    match(rownames(RFpurify_ABSOLUTE$importance), rownames(betas))
    ,
    , drop = FALSE
  ]
  rownames(betas)<-rownames(RFpurify_ABSOLUTE$importance)
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
#' @param frac Fraction of faulty probes (P-value > pval ) allowed. default 0.1
#' @param pval P-value cutoff. default 0.01
#' @param remove_sex Boolean. Should sex chromosomes be removed? default = TRUE
#' @param rgSet Input RGChannelSet object with raw idats and metadata.
#' @inheritParams pre_process

#' @importFrom graphics barplot abline legend par
#' @importFrom grDevices  dev.off
#' @importFrom stats sd
#' @return Normalized & filtered RGset.

queryfy<-function(targets, rgSet,sampGroups=NULL,sampNames="Sample_Name",frac=0.1,pval=0.01,remove_sex=TRUE,arraytype=NULL,qc_folder= "analysis/intermediate/QC"){
  # requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  # requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  if (!(sampNames %in% names(rgSet@colData))) sampNames <- colnames(rgSet)
  n <- ncol(rgSet)
  #o <- rev(order(sampNames))
  #rgSet <- rgSet[, o]
  #sampNames <- sampNames[o]
  if (is.null(sampGroups)) g <- rep(1, n) else g <- targets[[sampGroups]]
  if (is.null(g)) g<-1:n
  g <- factor(g)

  pal =  c(
    "#191919", "#F0A0FF", "#0075DC", "#993F00", "#005C31", "#5EF1F2", "#FF0010",
    "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
    "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
    "#740AFF", "#990000", "#FFFF80", "#FFE100", "#FF5005")


  # Check quality of the combined probe signal (testing against negative controls)
  # 1. Remove low quality samples. Defaults more than 10% we throw samples away.
  detP <- minfi::detectionP(rgSet, type = "m+u")
  grDevices::pdf(file = paste0(qc_folder,"mean_detection_pvalues.pdf"),width = 7,height = 7)
  ylabels<-colnames(detP)
  par(mar=c(max(4.1,max(nchar(ylabels))/2.2) ,4.1 , 4.1, 2.1))

  barplot(colMeans(detP), col=pal[g], las=2,
          cex.names=0.8, ylim=c(0,max(0.002,max(colMeans(detP))*2)), main ="Mean detection p-values")
  graphics::abline(h=0.05,col="red")
  graphics::legend("topleft", legend=levels(g), fill=pal[1:length(levels(g))],
         bg="white")
  grDevices::dev.off()()

  bad <- colnames(detP)[colSums(detP >=pval)/nrow(detP) > frac]
  if(length(bad)>0){
    warning("The following samples will be discarded since they fail to pass the p-value filter ( ",
            frac*100,"% of the probes with p-val >", pval, "): \n ", paste(bad,collapse = ", " ))
    rgSet <- rgSet[,setdiff(colnames(detP),bad)]
  }else{
    cat("All samples passed detection P-value filter")
  }

  # 2. Preprocessing:
  mSets<-list()
  mSets[["noob"]] <- minfi::preprocessNoob(rgSet)
  # mSets[["pq"]] <- minfi::preprocessQuantile(rgSet)
  # mSetSqn_funn <- minfi::preprocessFunnorm(rgSet)
  # mSetSqn_noob_pq <- minfi::preprocessQuantile(mSetSqn_noob)
  # mSetSqn_funn_pq <- minfi::preprocessQuantile(mSetSqn_funn)

  #mSets<-list(mSetSqn_noob,mSetSqn_pq,mSetSqn_funn,mSetSqn_noob_pq,mSetSqn_funn_pq)
  mSets_out<-list()
  for (i in names(mSets)){
    mSetSqn<-mSets[[i]]
    # 3. Removing low-quality probes (with p-value below pval)
    mSetSqn <- mSetSqn[rowSums(detP < pval) == ncol(mSetSqn),]

    # 4. Removing probes with known SNPs at CpG site
    mSetSqn <-  minfi::mapToGenome(mSetSqn)
    mSetSqn <- minfi::dropLociWithSnps(mSetSqn)

    # 5. Removing cross reactive probes
    mSetSqn <-  maxprobes::dropXreactiveLoci(mSetSqn)

    # 6. Sex. prediction & removal
    mSetSqn$Sex_pred <- minfi::getSex(mSetSqn, cutoff = -2)$predictedSex
    if(remove_sex){
      if(!is.null(arraytype)){anno<-get_anno(arraytype)
      }else{
        anno<-minfi::getAnnotation(mSetSqn)
        anno<-anno[!(anno$chr %in% c("chrX","chrY")),]
      }
      mSets_out[[i]]<-mSetSqn[rownames(anno),]
    }

  }
  return(mSets_out)

  # mSetSqn <- CNV.load(mSetSqn)
  # return(mSetSqn@intensity)
}
