
#' For each contrast extract a set of DMPs and add gene annotation and methylation values
#'
#'
#' @title Extract DMPs, annotation and methylation difference for each contrast
#'
#' @return data.table
#' @author Izar de Villasante
#' @author Angelika Merkel
#' @export
#' @import minfi
#' @import data.table
#' @import limma
#' @param beta_normalized normalized betavalues, as produce by minfi::getBeta(grSet_noob)),
#'  where colnames(beta_normalized) == metadata$sample_Name
#' @param ContrastsDM list of contrasts as returned by limma::makeContrasts()
#' which will pass to limma topTable as input
#' @param mDiff absolute mean methylation difference between groups to filter by
#' @param ann annotation dataset from manifest with metadata such as gene info,
#' CGI, RefGene, etc. see topTable genelist arg.
#' @param writeOut save result as .csv default = TRUE.
#'
#' @inheritParams limma::topTable
#' @examples
#'
#' betas<-readRDS("data/beta_noob.rds")
#' fit<-readRDS("data/fit2.rds")
#' ann<-readRDS("data/ann.rds")
#' DMPann <- DMPextr(fit = fit,                       # linear contrast model
#'                   ContrastsDM = ContrastsDM,          # contrasts
#'                   p.value = 0.01,                      # filter significantly different probes
#'                   beta_normalized = beta_noob,        # extract mean group betas
#'                   mDiff = 0.5,                        # select mean methylation differences
#'                   ann = ann,                          # annotate positions (CGI, RefGene, etc)
#'                   writeOut = FALSE                    # write output to file


DMPextr <- function(
  fit, ContrastsDM, p.value, beta_normalized, mDiff, ann, writeOut = TRUE
  ){
  ann<-data.table::as.data.table(ann,keep.rownames = "rn")
  data.table::setkey(ann,"rn")
  # ann      = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  annSub   = ann[rownames(beta_normalized)]
  DMP.list = list()

  for (i in 1:length(ContrastsDM)){
    # Extract DMPs with significant differences (+ annotate)
    DMP_1 <- limma::topTable(fit,
                      num = Inf,
                      coef = i ,
                      genelist = annSub,
                      p.value = p.value  # = p.adj
    )

    if(nrow(DMP_1)==0){
      warning(paste("No DMP found for contrast:", ContrastsDM[i], sep=" "))
      next
    }
    else{
      DMP_1$Type <- "Hyper"
      DMP_1$Type[which(DMP_1$logFC > 0)] <- "Hypo"  # invert direction (see above0)
      #
      # DMP_1 <- DMP_1[ , c("chr", "pos", "strand", "Name",
      #                     "Type", "P.Value" , "adj.P.Val" ,
      #                     "Islands_Name", "Relation_to_Island",
      #                     "UCSC_RefGene_Name", "UCSC_RefGene_Accession","UCSC_RefGene_Group",
      #                     "Phantom4_Enhancers","Phantom5_Enhancers","X450k_Enhancer",
      #                     "Regulatory_Feature_Name", "Regulatory_Feature_Group",
      #                     "GencodeBasicV12_NAME","GencodeBasicV12_Accession", "GencodeBasicV12_Group",
      #                     "GencodeCompV12_NAME", "GencodeCompV12_Accession", "GencodeCompV12_Group"
      # )]
      DMP_1$Contrast  <- ContrastsDM[i]

      # Extract the methylation values for the respective DMP and samples
      # Get contrast i:
      c1 <-fit$contrasts[,ContrastsDM[i]]
      c1 <-c1[c1!=0]
      # Variable names in contrast i:
      vars1 <-names(c1)
      # design matrix for contrasts:
      design <-fit$design[,vars1]
      design <-t(design) * c1

      # betas
      cg         <- which(rownames(beta_normalized) %in% rownames(DMP_1))
      betas <-beta_normalized[cg,]
      # diff  _mean:
      DMP_1$diff_meanMeth<-rowSums(apply(design,1,function(x) rowSums(betas %*% diag(x))))

      # filter for absolute methylation difference
      DMP_1 <- DMP_1[which(abs(DMP_1$diff_meanMeth)>= mDiff), ]

      # save DMPs
      DMP.list[[i]] <- DMP_1

      # write output file
      if(writeOut == TRUE){

        cat(paste("writing analysis/DMP_", ContrastsDM[i], ".csv\n",sep =""))
        data.table::fwrite(DMP_1, file = paste("analysis/DMP_", ContrastsDM[i], ".csv", sep =""))
      } else {
        next
      }
    }
  }
  return(do.call(rbind, DMP.list))
}
