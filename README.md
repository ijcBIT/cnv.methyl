cnv-methyl manual: dynamic Somatic Copy Nunmber Alterations for
methylation array data.
================
Izar de Villasante
03 May 2022

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Introduction

The `cnv.methyl`package implements an automatic and dynamic thresholding
technique for **c**opy **n**umver **v**ariation analysis using Illumina
450k or EPIC DNA **methyl**ation arrays. It enhances the whole pipline
of methylation array processing from raw data to cnv calls. The main
reasons to use it are:

1.  It is fast. Runs in parallel and it is prepared to be run on HPC
    environments seamlessly.
2.  It is dynamic. It provides cnv calls and copy number values for your
    genes of interest based on each array metrics (purity, mean, sd).
3.  It is flexible. Since it accepts both **450k** and **EPIC**
    methylation arrays, different controls and genome annotations and
    although it relies on `conumee` package to calculate segmentations
    and generate log2r intensities, it can also accept input from any
    other tools.

## Context

Although the primary purpose of methylation arrays is the detection of
genome-wide DNA methylation levels \[@bibikovahigh2011\], it can
additionally be used to extract useful information about copy-number
alterations, e.g. in clinical cancer samples. The approach was initially
described in Sturm et al., 2012 \[@sturmhotspot2012\]. Some tools have
been developed for this purpose, such as `conumee`, `ChAMP` and
`CopyNumber450k` \[@conumee;@champ;@cnv450k\].

Nevertheless, all this tools require a certain level of human
interaction and interpretation in order to obtain meaningful information
from the data. Some of these tasks are automatically resolved by the
package such as providing the right annotation for the genes of interest
or setting a threshold for each cna.

Setting the threshold, may be the most decisive and challenging step. So
far, the main approach is to visually inspect the log2ratios plot in
order to get some insight of the genomic alterations. Nevertheless this
threshold may vary between arrays depending on purity of the samples and
noise. SNPs based arrays tools, such as `ASCAT` or `PURPLE`, are aware
of this problems and correct it in order to provide a copy number value.
Nevertheless, this level of precision had not been yet accomplished by
the available tools for methylation arrays. Until now.

This CNV analysis pipline can be broken into 3 main blocks:
pre-processing, segmentation & cnv calling:

# Pre-processing:

This pipeline uses an enhanced parallel version of `read.metharray.exp`
function from `minfi` package @minfi in order to load the idats and a
precalculated matrix from `RFpurify` @RFpurify for imputing the
purities.

The minimum requirements to run this pipline are the sample sheet and
its corresponding idat files. Both of them can be found within the
package as example data. Let’s have a look:

``` r
sample_sheet<-data.table::fread(system.file("extdata", "Sample_sheet_example.csv",package="cnv.methyl"))
str(sample_sheet)
#> Classes 'data.table' and 'data.frame':   20 obs. of  7 variables:
#>  $ Sample_Name : chr  "TCGA-63-A5MN-01A-22D-A27L-05" "TCGA-63-A5MP-01A-11D-A26N-05" "TCGA-63-A5MR-01A-31D-A27L-05" "TCGA-63-A5MS-01A-11D-A26N-05" ...
#>  $ filenames   : chr  "inst/extdata/sample_IDATS/9283265144_R06C01" "inst/extdata/sample_IDATS/9305216074_R05C01" "inst/extdata/sample_IDATS/9283265144_R04C02" "inst/extdata/sample_IDATS/9305216121_R02C01" ...
#>  $ Sample_Plate: logi  NA NA NA NA NA NA ...
#>  $ Sample_Group: chr  "Cancer" "Cancer" "Cancer" "Cancer" ...
#>  $ Pool_ID     : logi  NA NA NA NA NA NA ...
#>  $ Project     : chr  "TCGA-LUSC" "TCGA-LUSC" "TCGA-LUSC" "TCGA-LUSC" ...
#>  $ Sample_Well : logi  NA NA NA NA NA NA ...
#>  - attr(*, ".internal.selfref")=<externalptr>
```

This format has to be respected in order to make everything work fine.
`Sample_Name` contains the sample ids and `filenames` contain their
current path. The other parameters are required by minfi in order to
read the arrays. Basename contains the working directory where minfi
will looks for files and it depends on the folder parameter.

``` r
library(cnv.methyl)
folder="analysis/intermediate/IDATS/"
sample_sheet$Basename<-paste0(folder, basename(sample_sheet$filenames))
myLoad <- pre_process.myLoad(
  targets=sample_sheet,folder=folder,arraytype="450K"
  )
#> working directory: analysis/intermediate/IDATS/
#> Reading multiple idat-files in parallel. Using 14 cores.
#> Minfi does load data
```

``` r
myLoad
#> Loading required package: minfi
#> Loading required package: BiocGenerics
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#>     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#>     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#>     union, unique, unsplit, which.max, which.min
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: GenomeInfoDb
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
#> Loading required package: Biostrings
#> Loading required package: XVector
#> 
#> Attaching package: 'Biostrings'
#> The following object is masked from 'package:base':
#> 
#>     strsplit
#> Loading required package: bumphunter
#> Loading required package: foreach
#> Loading required package: iterators
#> Loading required package: parallel
#> Loading required package: locfit
#> locfit 1.5-9.4    2020-03-24
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
#> class: RGChannelSet 
#> dim: 622399 20 
#> metadata(0):
#> assays(2): Green Red
#> rownames(622399): 10600313 10600322 ... 74810490 74810492
#> rowData names(0):
#> colnames(20): 9283265144_R06C01 9305216074_R05C01 ... 6004791006_R06C02
#>   6004791020_R02C02
#> colData names(8): Sample_Name filenames ... Sample_Well Basename
#> Annotation
#>   array: IlluminaHumanMethylation450k
#>   annotation: ilmn12.hg19
```

Once data is loaded tumor purity in the sample is calculated with
`purify()` :

``` r
library(cnv.methyl)
purity<-purify(myLoad=myLoad)
myLoad@colData$Purity_Impute_RFPurify<-purity
```

This function uses the pre-computed matrix of absolute purities
`RFpurify_ABSOLUTE` provided in RFpurify @RFpurify and also performs
imputation of missing values with `impute.knn`.

Then data is normalized with `queryfy`:

``` r
query <- queryfy(myLoad)
#> Loading required namespace: IlluminaHumanMethylationEPICanno.ilm10b4.hg19
#> [preprocessQuantile] Mapping to genome.
#> [preprocessQuantile] Fixing outliers.
#> [preprocessQuantile] Quantile normalizing.
```

This function perforrms the following steps:

1.  Normalization by `minfi::reprocessQuantile`

2.  Filtering:

    -   `minfi::detectionP()` Remove arrays: Minimum detection threshold
        with p-value \< 0.01 per probe and a maximum of 10% failed
        probes per sample.

    -   `minfi::dropLociWithSnps()` Removes snps: probes that were
        located +/- 10 bases away from known SNPs were also filtered
        out.

    -   `maxprobes::dropXreactiveLoci()` Remove cross reactive probes.

You can use `pre_process()` function in order to run all these steps at
once:

``` r
library(cnv.methyl)
sample_sheet<-data.table::fread(system.file("extdata", "Sample_sheet_example.csv",package="cnv.methyl"))
out="analysis/intermediate/"
pre_process(targets = sample_sheet, out = out, RGset = F )
#> Warning in dir.create(folder, recursive = TRUE): 'analysis/intermediate/IDATS'
#> already exists
#> working directory: analysis/intermediate/IDATS/
#> Reading multiple idat-files in parallel. Using 14 cores.
#> Minfi does load data
#> Loading required namespace: randomForest
#> Rfpurifies
#> targets is saved  analysis/intermediate/
#> Loading required namespace: IlluminaHumanMethylationEPICanno.ilm10b4.hg19
#> [preprocessQuantile] Mapping to genome.
#> [preprocessQuantile] Fixing outliers.
#> [preprocessQuantile] Quantile normalizing.
#> Pre-processing Completed successfully!
#> class: GenomicRatioSet 
#> dim: 422646 20 
#> metadata(0):
#> assays(2): M CN
#> rownames(422646): cg16619049 cg23100540 ... cg07211220 cg02233183
#> rowData names(0):
#> colnames(20): 9283265144_R06C01 9305216074_R05C01 ... 6004791006_R06C02
#>   6004791020_R02C02
#> colData names(12): Sample_Name filenames ... yMed predictedSex
#> Annotation
#>   array: IlluminaHumanMethylation450k
#>   annotation: ilmn12.hg19
#> Preprocessing
#>   Method: Raw (no normalization or bg correction)
#>   minfi version: 1.38.0
#>   Manifest version: 0.4.0
```

``` elixir
[//]: (It is important to know the folder structure. 
out is for working directory, 
subf is the subfolder inside out where minfi reads the idats. 
Folder overwrites this parameters.  )
```

As you can see it is pretty easy to substitute any part of the analysis
to your convenience. The most interesting part here is that all packages
and data needed for the analysis are bundled together and the speedup of
`minfi` loading idats in `parallel` .

# Segmentation

Segmentation is performed using a two-step approach as described in
`conumee`. First, the combined intensity values of both ‘methylated’ and
‘unmethylated’ channel of each CpG are normalized using a set of normal
controls (i.e. with a flat genome not showing any copy-number
alterations).

By default a set of 96 Whole Blood Samples are used unless the user
specifies something else:

``` r
 cnv.methyl:::control
#> Loading required package: conumee
#> Loading required package: IlluminaHumanMethylationEPICanno.ilm10b2.hg19
#> 
#> Attaching package: 'IlluminaHumanMethylationEPICanno.ilm10b2.hg19'
#> The following objects are masked from 'package:IlluminaHumanMethylation450kanno.ilmn12.hg19':
#> 
#>     Islands.UCSC, Locations, Manifest, Other, SNPs.132CommonSingle,
#>     SNPs.135CommonSingle, SNPs.137CommonSingle, SNPs.138CommonSingle,
#>     SNPs.141CommonSingle, SNPs.142CommonSingle, SNPs.144CommonSingle,
#>     SNPs.146CommonSingle, SNPs.147CommonSingle, SNPs.Illumina
#> Loading required package: IlluminaHumanMethylationEPICmanifest
#> CNV data object
#>    created   : 
#>   @intensity : available (25 samples, 416166 probes)
```

Also annotation for **EPIC** , **450K** , or an overlap of both arrays
is also built-in:

``` r
cnv.methyl:::anno_epic
#> CNV annotation object
#>    created  : Tue Apr 26 19:32:10 2022
#>   @genome   : 22 chromosomes
#>   @gap      : 313 regions
#>   @probes   : 844316 probes
#>   @exclude  : 10042 regions (overlapping 115308 probes)
#>   @detail   : 94 regions (overlapping 2565 probes)
#>   @bins     : 21937 bins (min/avg/max size: 50/91.1/5000kb, probes: 15/32.4/437)
cnv.methyl:::anno_450K
#> CNV annotation object
#>    created  : Tue Apr 26 19:31:54 2022
#>   @genome   : 22 chromosomes
#>   @gap      : 313 regions
#>   @probes   : 470870 probes
#>   @exclude  : 10042 regions (overlapping 73549 probes)
#>   @detail   : 94 regions (overlapping 1496 probes)
#>   @bins     : 12888 bins (min/avg/max size: 50/132.4/5000kb, probes: 15/29.5/478)
cnv.methyl:::anno_overlap
#> CNV annotation object
#>    created  : Tue Apr 26 19:32:22 2022
#>   @genome   : 22 chromosomes
#>   @gap      : 313 regions
#>   @probes   : 439635 probes
#>   @exclude  : 10042 regions (overlapping 68583 probes)
#>   @detail   : 94 regions (overlapping 1403 probes)
#>   @bins     : 12272 bins (min/avg/max size: 50/137.9/5000kb, probes: 15/28.8/432)
```

The right annotation will be chosen in each case according to the
arraytype. If not specified it defaults to `anno_overlap`.

If you are using the built-in controls and epic arrays as input you
should change arraytype to overlap. The intensities can be given in the
following formats:

a genomics Ratio set or path to file. Accepted formats are: .fst, .rds,
.txt or other text formats readable by fread.

The output of `pre_process` function above with intensities from 20
samples is used as example dataset:

``` r
# 
# intensity<-readRDS(system.file("extdata", "intensities.RDS",package="cnv.methyl"))
# run_conumee(intensities=intensity,arraytype="450K")
# run_conumee(intensities = intensity,anno_file=anno, ctrl_file=ctrl_file,
#               Sample_Name=Sample_Name,seg.folder = seg.folder,
#               log2r.folder = log2r.folder,arraytype=arraytype,
#               conumee.folder=conumee.folder, probeid=probeid)
```

This step is required for correcting for probe and sample bias (e.g.
caused by GC-content, type I/II differences or technical variability).
Secondly, neighboring probes are combined in a hybrid approach,
resulting in bins of a minimum size and a minimum number of probes. This
step is required to reduce remaining technical variability and enable
meaningful segmentation results.

Plotting is not yet implemented by the pipeline. If you are interested
in this functionality, please make a request.

# cnv calling:

In order to remove biological and technical noise of each sample, the
relationship between log2ratio signal and noise \[tumor purity & signal
sd\] has been calculated for each copy number alteration, such that

![{Amp10,Amp,Gains,Diploid,HetLoss,HomDel}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7BAmp10%2CAmp%2CGains%2CDiploid%2CHetLoss%2CHomDel%7D "{Amp10,Amp,Gains,Diploid,HetLoss,HomDel}")

A precalculated constant Kcn is used in order to adjust each array’s
purity and

<!-- # Load data -->
<!-- The recommended input format are `Mset` objects generated from raw IDAT -->
<!-- files using the `minfi` package. Depending on the analysis workflow, -->
<!-- these objects might be already available. A good example 450k data-set -->
<!-- for testing the `conumee` package can be downloaded from TCGA (971.6 -->
<!-- MB): -->
<!-- [~~https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/cgcc/jhu-usc.edu/humanmethylation450/methylation/jhu-usc.edu_BRCA.HumanMethylation450.Level_1.8.8.0.tar.gz~~](https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/cgcc/jhu-usc.edu/humanmethylation450/methylation/jhu-usc.edu_BRCA.HumanMethylation450.Level_1.8.8.0.tar.gz) -->
<!-- UPDATE: The TCGA Data Portal is no longer operational. TCGA DNA -->
<!-- methylation data can now be downloaded from the NCI GDC Legacy Archive: -->
<!-- <https://gdc-portal.nci.nih.gov/legacy-archive> -->
<!-- This data-set comprises of 42 primary breast cancer samples, 2 breast -->
<!-- cancer cell lines and 16 control samples. Make sure to unpack the -->
<!-- archive. For the purpose of this vignette, only two illustrative -->
<!-- examples are downloaded using the `read.450k.url` method. -->
<!-- ```{r, echo = TRUE, collapse = TRUE, results = 'markup', warning = FALSE, message = FALSE} -->
<!-- library("minfi") -->
<!-- library("conumee") -->
<!-- #RGsetTCGA <- read.450k.exp(base = "~/conumee_analysis/jhu-usc.edu_BRCA.HumanMethylation450.Level_1.8.8.0")  # adjust path -->
<!-- RGsetTCGA <- read.450k.url()  # use default parameters for vignette examples -->
<!-- MsetTCGA <- preprocessIllumina(RGsetTCGA) -->
<!-- MsetTCGA -->
<!-- ``` -->
<!-- Alternatively, you can use the `minfiData` example data-set, comprising -->
<!-- of 3 cancer samples and 3 normal controls. -->
<!-- ```{r, echo = TRUE, collapse = TRUE, results = 'markup', warning = FALSE, message = FALSE} -->
<!-- library("minfiData") -->
<!-- data(MsetEx) -->
<!-- MsetEx -->
<!-- ``` -->
<!-- The `CopyNumber450k` package provides a large data-set of 52 control -->
<!-- samples which can be used for normalization. -->
<!-- ```{r, echo = TRUE, collapse = TRUE, results = 'markup', warning = FALSE, message = FALSE} -->
<!-- library("CopyNumber450kData") -->
<!-- data(RGcontrolSetEx) -->
<!-- MsetControls <- preprocessIllumina(RGcontrolSetEx) -->
<!-- ``` -->
<!-- If raw IDAT files are unavailable, data can also be supplied by -->
<!-- importing text-based input files, e.g. as generated by GenomeStudio or -->
<!-- downloaded from the GEO repository. More details are given later. -->
<!-- # Create annotation object -->
<!-- To begin with the CNV analysis, an annotation object, generated by the -->
<!-- `CNV.create_anno` function, is required. This object holds information -->
<!-- which only needs to be generated once, irrespective of the number of -->
<!-- samples that are processed. -->
<!-- Arguments `bin_minprobes` and `bin_minsize` define the minimum number of -->
<!-- probes per bin and the minimum bin size (default values that were -->
<!-- optimized for 450k data are 15 and 50000, respectively). Argument -->
<!-- `exclude_regions` defines regions to be excluded (e.g. polymorphic -->
<!-- regions, an example is given in `data(exclude_regions)`). Please see -->
<!-- `?CNV.create_anno` for more details. -->
<!-- Argument `detail_regions` defines regions to be examined in detail (e.g. -->
<!-- dedicated detail plots or text output, see below). For example, detail -->
<!-- regions can contain known oncogenes or tumor suppressor genes. These -->
<!-- regions should either be supplied as a path to a BED file or as a -->
<!-- GRanges object (an example is given in `data(detail_regions)`). The -->
<!-- start and stop coordinates indicate the regions in which probes are -->
<!-- analyzed in detail. The plotting range of detail plots are defined in a -->
<!-- second set of start and stop coordinates in the GRanges object values -->
<!-- (or the 7^th^ and 8^th^ column in the BED file). -->
<!-- ```{r, echo = TRUE, collapse = TRUE, results = 'markup', warning = TRUE, message = TRUE} -->
<!-- data(exclude_regions) -->
<!-- data(detail_regions) -->
<!-- head(detail_regions, n = 2) -->
<!-- anno <- CNV.create_anno(array_type = "450k", exclude_regions = exclude_regions, detail_regions = detail_regions) -->
<!-- anno -->
<!-- ``` -->
<!-- # Combine intensity values -->
<!-- Intensity values of the 'methylated' and 'unmethylated' channel are -->
<!-- combined using the `CNV.load` function. Input can be either an `Mset` -->
<!-- object (recommended, see above), or a `data.frame` or `matrix` object -->
<!-- containing intensities generated by GenomeStudio or downloaded from GEO -->
<!-- (imported using e.g. `read.table`, should work without reformatting for -->
<!-- most tables, please refer to `?CNV.load` for more details). -->
<!-- ```{r, echo = TRUE, collapse = TRUE, results = 'markup', warning = FALSE, message = FALSE} -->
<!-- data(tcgaBRCA.sentrix2name)  # TCGA sample IDs are supplied with the conumee package -->
<!-- sampleNames(MsetTCGA) <- tcgaBRCA.sentrix2name[sampleNames(MsetTCGA)] -->
<!-- tcga.data <- CNV.load(MsetTCGA) -->
<!-- tcga.controls <- grep("-11A-", names(tcga.data)) -->
<!-- names(tcga.data) -->
<!-- tcga.data -->
<!-- minfi.data <- CNV.load(MsetEx) -->
<!-- minfi.controls <- pData(MsetEx)$status == "normal" -->
<!-- controls.data <- CNV.load(MsetControls) -->
<!-- ``` -->
<!-- # Perform CNV analysis -->
<!-- The main CNV analysis is separated into four parts: -->
<!-- -   First, `CNV.fit` is used to normalize a single query sample to a set -->
<!--     of control samples by multiple linear regression. For best results -->
<!--     control samples of matched normal tissues that are profiled within -->
<!--     the same experiment are used (which are likely to have the same -->
<!--     technical bias). Essentially this regression analysis yields the -->
<!--     linear combination of control samples that most closely fits the -->
<!--     intensities of the query sample. Subsequently, the log2-ratio of -->
<!--     probe intensities of the query sample versus the combination of -->
<!--     control samples are calculated and used for further analysis. More -->
<!--     details are given in the publication. -->
<!-- -   Secondly, `CNV.bin` is used to combine probes within predefined -->
<!--     genomic bins. Bins are previously generated using `CNV.create_anno`. -->
<!--     Intensity values are shifted to minimize the median absolute -->
<!--     deviation from all bins to zero to determine the copy-number neutral -->
<!--     state.   -->
<!-- -   Thirdly, `CNV.detail` is used to analyze detail regions in detail. -->
<!--     This step is optional, but required if detail regions should be -->
<!--     outputted in plots and text files. Detail regions are defined using -->
<!--     `CNV.create_anno`. -->
<!-- -   Finally, `CNV.segment` is used to segment the genome into regions of -->
<!--     the same copy-number state. This function is a wrapper of the `CNA`, -->
<!--     `segment`, `segments.summary` and `segments.p` functions of the -->
<!--     `DNAcopy` package. Default parameters were optimized for 450k data, -->
<!--     but can be changed. See `?CNV.segment` for more details. -->
<!-- ```{r, echo = TRUE, collapse = TRUE, results = 'markup', warning = FALSE, message = FALSE} -->
<!-- x <- CNV.fit(minfi.data["GroupB_1"], minfi.data[minfi.controls], anno) -->
<!-- y <- CNV.fit(tcga.data["TCGA-AR-A1AU-01A-11D-A12R-05"], controls.data, anno)  # use TCGA control samples for better results -->
<!-- z <- CNV.fit(tcga.data["TCGA-AR-A1AY-01A-21D-A12R-05"], controls.data, anno) -->
<!-- x <- CNV.bin(x) -->
<!-- x <- CNV.detail(x) -->
<!-- x <- CNV.segment(x) -->
<!-- y <- CNV.segment(CNV.detail(CNV.bin(y))) -->
<!-- z <- CNV.segment(CNV.detail(CNV.bin(z))) -->
<!-- x -->
<!-- ``` -->
<!-- # Output plots and text files -->
<!-- The `conumee` package supports two types of plots: -->
<!-- -   The `CNV.genomeplot` method produces plots of the complete genome or -->
<!--     of one or multiple chromosomes. Intensity values of each bin are -->
<!--     plotted in colored dots. Segments are shown as blue lines. See -->
<!--     `?CNV.genomeplot` for more details. -->
<!-- -   The `CNV.detailplot` methods produces plots of individual detail -->
<!--     regions, as defined in `CNV.create_anno`. Intensity values of -->
<!--     individual probes are plotted in colored crosses. Bins are shown as -->
<!--     blue lines. `CNV.detailplot_wrap` is a wrapper function that -->
<!--     produces a single plot of all detail regions. -->
<!-- Text output is generated using the `CNV.write` method. Parameter `what` -->
<!-- specifies if "probes", "bins", "detail" or "segments" should be -->
<!-- returned. If parameter `file` is specified, the output is written into a -->
<!-- file, otherwise a `data.frame` is returned. See `?CNV.write` for more -->
<!-- details. -->
<!-- ```{r, echo = TRUE, fig.show = 'hold', collapse = TRUE, results = 'markup', warning = FALSE, message = FALSE, fig.width = 14, fig.height = 7, out.width = "800px", fig.retina = 1} -->
<!-- #pdf("~/conumee_analysis/CNVplots.pdf", height = 9, width = 18) -->
<!-- CNV.genomeplot(x) -->
<!-- CNV.genomeplot(y) -->
<!-- CNV.genomeplot(y, chr = "chr6") -->
<!-- CNV.genomeplot(z) -->
<!-- CNV.genomeplot(z, chr = "chr10") -->
<!-- CNV.detailplot(z, name = "PTEN") -->
<!-- CNV.detailplot_wrap(z) -->
<!-- #dev.off() -->
<!-- head(CNV.write(y, what = "segments"), n = 5) -->
<!-- #CNV.write(y, what = "segments", file = "~/conumee_analysis/TCGA-AR-A1AU.CNVsegments.seg")  # adjust path -->
<!-- #CNV.write(y, what = "bins", file = "~/conumee_analysis/TCGA-AR-A1AU.CNVbins.igv") -->
<!-- #CNV.write(y, what = "detail", file = "~/conumee_analysis/TCGA-AR-A1AU.CNVdetail.txt") -->
<!-- #CNV.write(y, what = "probes", file = "~/conumee_analysis/TCGA-AR-A1AU.CNVprobes.igv") -->
<!-- ``` -->
<!-- # Contact & citation -->
<!-- For bug-reports, comments and feature requests please write to -->
<!-- [conumee\@hovestadt.bio](mailto:conumee@hovestadt.bio). -->
<!-- When using `conumee` in your work, please cite as: -->
<!-- ```{r, echo = TRUE, collapse = TRUE, results = 'markup', warning = FALSE, message = FALSE} -->
<!-- citation("cnv.methyl") -->
<!-- ``` -->
<!-- # Session info -->
<!-- ```{r, echo = TRUE, collapse = TRUE, results = 'markup', warning = FALSE, message = FALSE} -->
<!-- sessionInfo() -->
<!-- ``` -->
<!-- # References -->
<!-- # Manual: -->
<!-- <!-- badges: start -->

–>

<!-- <!-- badges: end -->

–>

<!-- The goal of cnv.methyl is to ... -->
<!-- ## Installation -->
<!-- You can install the development version of cnv.methyl from -->
<!-- [GitHub](https://github.com/) with: -->
<!-- ``` r -->
<!-- # install.packages("devtools") -->
<!-- devtools::install_github("https://github.com/ijcBIT/cnv.methyl.git") -->
<!-- ``` -->
<!-- ## Example -->
<!-- This is a basic example which shows you how to get somatic copy number -->
<!-- alterations using the pre-calculated constant used in the paper -->
<!-- [Blecua,P et al.](https://academic.oup.com/bib). First load your targets -->
<!-- sample sheet: -->
<!-- ```{r samplesheetasd} -->
<!-- sample_sheet<-utils::read.csv(system.file("extdata", "Sample_sheet_example.csv",package="cnv.methyl")) -->
<!-- str(sample_sheet) -->
<!-- head(sample_sheet) -->
<!-- ``` -->
<!-- The main function is `run_cnv.methyl` that will perform the whole -->
<!-- pipline from idats to GRanges object with gene and copy number -->
<!-- alteration annotations: -->
<!-- ```{r main function,cache=TRUE} -->
<!-- library(cnv.methyl) -->
<!-- out="analysis/example/" -->
<!-- gr_list<-run_cnv.methyl(targets = sample_sheet,arraytype = "450K",Kc = "curated", out=out) -->
<!-- gr_list -->
<!-- ``` -->
<!-- If you wish to combine 450K and epic arrays you can use the variable -->
<!-- `arraytype = overlap` Here is an example using the pre-processed data -->
<!-- from the previous code chunk: **Important: don't forget to load the -->
<!-- right sample sheet with the purities. This file will be in -->
<!-- *"path-to-out-folder/Sample_Sheet.txt"*** -->
<!-- ```{r overlapping arraytypes} -->
<!-- library(cnv.methyl) -->
<!-- out="analysis/example/" -->
<!-- ss<-data.table::fread(paste0(out,"Sample_Sheet.txt")) -->
<!-- ID <- ss$Sample_Name[1] -->
<!-- ss <- ss[Sample_Name==ID,] -->
<!-- cna<-Kc_get(ss=ss,ID=ID,arraytype = "overlap") -->
<!-- cna -->
<!-- ## basic example code -->
<!-- ``` -->
