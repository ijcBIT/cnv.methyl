cnv-methyl manual: dynamic Somatic Copy Nunmber Alterations for
methylation array data.
================
Izar de Villasante
28 April 2022

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Manual:

<!-- badges: start -->
<!-- badges: end -->

The goal of cnv.methyl is to …

## Installation

You can install the development version of cnv.methyl from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("izarvillasante/cnv.methyl")
```

## Example

This is a basic example which shows you how to get somatic copy number
alterations using the pre-calculated constant used in the paper
[Blecua,P et al.](https://academic.oup.com/bib). First load your targets
sample sheet:

``` r
sample_sheet<-utils::read.csv(system.file("extdata", "Sample_sheet_example.csv",package="cnv.methyl"))
str(sample_sheet)
#> 'data.frame':    20 obs. of  9 variables:
#>  $ X           : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ Sample_Name : chr  "TCGA-63-A5MN-01A-22D-A27L-05" "TCGA-63-A5MP-01A-11D-A26N-05" "TCGA-63-A5MR-01A-31D-A27L-05" "TCGA-63-A5MS-01A-11D-A26N-05" ...
#>  $ filenames   : chr  "inst/extdata/sample_IDATS/9283265144_R06C01" "inst/extdata/sample_IDATS/9305216074_R05C01" "inst/extdata/sample_IDATS/9283265144_R04C02" "inst/extdata/sample_IDATS/9305216121_R02C01" ...
#>  $ Sample_Plate: logi  NA NA NA NA NA NA ...
#>  $ Sample_Group: chr  "Cancer" "Cancer" "Cancer" "Cancer" ...
#>  $ Pool_ID     : logi  NA NA NA NA NA NA ...
#>  $ Project     : chr  "TCGA-LUSC" "TCGA-LUSC" "TCGA-LUSC" "TCGA-LUSC" ...
#>  $ Sample_Well : logi  NA NA NA NA NA NA ...
#>  $ Basename    : chr  "data/sample_IDATS/TCGA-LUSC/TCGA-63-A5MN-01A-22D-A27L-05/9283265144_R06C01" "data/sample_IDATS/TCGA-LUSC/TCGA-63-A5MP-01A-11D-A26N-05/9305216074_R05C01" "data/sample_IDATS/TCGA-LUSC/TCGA-63-A5MR-01A-31D-A27L-05/9283265144_R04C02" "data/sample_IDATS/TCGA-LUSC/TCGA-63-A5MS-01A-11D-A26N-05/9305216121_R02C01" ...
head(sample_sheet)
#>   X                  Sample_Name                                   filenames
#> 1 1 TCGA-63-A5MN-01A-22D-A27L-05 inst/extdata/sample_IDATS/9283265144_R06C01
#> 2 2 TCGA-63-A5MP-01A-11D-A26N-05 inst/extdata/sample_IDATS/9305216074_R05C01
#> 3 3 TCGA-63-A5MR-01A-31D-A27L-05 inst/extdata/sample_IDATS/9283265144_R04C02
#> 4 4 TCGA-63-A5MS-01A-11D-A26N-05 inst/extdata/sample_IDATS/9305216121_R02C01
#> 5 5 TCGA-63-A5MT-01A-21D-A26N-05 inst/extdata/sample_IDATS/9305216126_R03C02
#> 6 6 TCGA-63-A5MU-01A-11D-A26N-05 inst/extdata/sample_IDATS/9305216121_R06C01
#>   Sample_Plate Sample_Group Pool_ID   Project Sample_Well
#> 1           NA       Cancer      NA TCGA-LUSC          NA
#> 2           NA       Cancer      NA TCGA-LUSC          NA
#> 3           NA       Cancer      NA TCGA-LUSC          NA
#> 4           NA       Cancer      NA TCGA-LUSC          NA
#> 5           NA       Cancer      NA TCGA-LUSC          NA
#> 6           NA       Cancer      NA TCGA-LUSC          NA
#>                                                                     Basename
#> 1 data/sample_IDATS/TCGA-LUSC/TCGA-63-A5MN-01A-22D-A27L-05/9283265144_R06C01
#> 2 data/sample_IDATS/TCGA-LUSC/TCGA-63-A5MP-01A-11D-A26N-05/9305216074_R05C01
#> 3 data/sample_IDATS/TCGA-LUSC/TCGA-63-A5MR-01A-31D-A27L-05/9283265144_R04C02
#> 4 data/sample_IDATS/TCGA-LUSC/TCGA-63-A5MS-01A-11D-A26N-05/9305216121_R02C01
#> 5 data/sample_IDATS/TCGA-LUSC/TCGA-63-A5MT-01A-21D-A26N-05/9305216126_R03C02
#> 6 data/sample_IDATS/TCGA-LUSC/TCGA-63-A5MU-01A-11D-A26N-05/9305216121_R06C01
```

The main function is `run_cnv.methyl` that will perform the whole
pipline from idats to GRanges object with gene and copy number
alteration annotations:

``` r
library(cnv.methyl)
out="analysis/example/"
gr_list<-run_cnv.methyl(targets = sample_sheet,arraytype = "450K",Kc = "curated", out=out)
#> Warning in dir.create(folder, recursive = TRUE): 'analysis/example/IDATS'
#> already exists
#> working directory: analysis/example/IDATS/
#> Reading multiple idat-files in parallel
#> Warning in read.metharray.exp.par(targets = targets, folder = folder, arraytype
#> = arraytype, : 14
#> Minfi does load data
#> Loading required namespace: randomForest
#> Rfpurifies
#> targets is saved  analysis/example/
#> [preprocessQuantile] Mapping to genome.
#> [preprocessQuantile] Fixing outliers.
#> [preprocessQuantile] Quantile normalizing.
#> Pre-processing Completed successfully!
#> Loading required namespace: conumee
#> anno
#> controls
#> dataset
#> Error in as.vector(x) : no method for coercing this S4 class to a vector
#> intensities
#> error calling combine function:
#> <simpleError in as(from, to_class, strict = FALSE): no method or default for coercing "simpleError" to "GRanges">
gr_list
#> GRanges object with 483 ranges and 9 metadata columns:
#>      seqnames              ranges strand |     log2r                     ID
#>         <Rle>           <IRanges>  <Rle> | <numeric>            <character>
#>    1     chr1   3232366-247756853      * |    -0.009 TCGA-63-A5MN-01A-22D..
#>    2    chr10    431113-135214495      * |     0.009 TCGA-63-A5MN-01A-22D..
#>    3    chr11    129255-134118566      * |    -0.002 TCGA-63-A5MN-01A-22D..
#>    4    chr12   1032854-133175503      * |    -0.047 TCGA-63-A5MN-01A-22D..
#>    5    chr13  19594416-114949919      * |    -0.100 TCGA-63-A5MN-01A-22D..
#>   ..      ...                 ...    ... .       ...                    ...
#>   33     chr8  97103770-102225000      * |     0.353 TCGA-97-8547-01A-11D..
#>   34     chr8 102319019-124375000      * |     0.132 TCGA-97-8547-01A-11D..
#>   35     chr8 124530104-131458854      * |     0.360 TCGA-97-8547-01A-11D..
#>   36     chr8 132144454-145193933      * |     0.132 TCGA-97-8547-01A-11D..
#>   37     chr9    855123-140723607      * |    -0.155 TCGA-97-8547-01A-11D..
#>       num.mark     bstat        pval seg.median         cna         cnv
#>      <integer> <numeric>   <numeric>  <numeric> <character> <character>
#>    1      1342        NA          NA      0.006     Diploid           2
#>    2       708        NA          NA      0.023     Diploid           2
#>    3       782        NA          NA      0.014     Diploid           2
#>    4       683        NA          NA     -0.028     Diploid           2
#>    5       341        NA          NA     -0.079     Diploid           2
#>   ..       ...       ...         ...        ...         ...         ...
#>   33        29   6.85972 3.76134e-10      0.346         Amp           7
#>   34        67   6.88306 3.16307e-10      0.125       Gains           3
#>   35        27   9.30070 9.55175e-19      0.353         Amp           7
#>   36        70        NA          NA      0.128       Gains           3
#>   37       310        NA          NA     -0.151     HetLoss           1
#>                   gene.mame
#>                 <character>
#>    1 ARHGEF16;MEGF6;TPRG1..
#>    2 DIP2C;LARP4B;GTPBP4;..
#>    3 LOC100133161;SCGB1C1..
#>    4 RAD52;LOC100292680;F..
#>    5 DKFZp686A1627;TUBA3C..
#>   ..                    ...
#>   33 GDF6;UQCRB;MTERFD1;P..
#>   34 NACAP1;GRHL2;NCALD;R..
#>   35 FBXO32;KLHL38;FAM91A..
#>   36 EFR3A;OC90;HHLA1;KCN..
#>   37 DMRT1;DMRT3;FLJ35024..
#>   -------
#>   seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

If you wish to combine 450K and epic arrays you can use the variable
`arraytype = overlap` Here is an example using the pre-processed data
from the previous code chunk: **Important: don’t forget to load the
right sample sheet with the purities. This file will be in
*“path-to-out-folder/Sample_Sheet.txt”***

``` r
library(cnv.methyl)
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
out="analysis/example/"
ss<-data.table::fread(paste0(out,"Sample_Sheet.txt"))
ID <- ss$Sample_Name[1]
ss <- ss[Sample_Name==ID,]
cna<-Kc_get(ss=ss,ID=ID,arraytype = "overlap")
#> segment and log2 files loaded.
#> int finish
cna
#> GRanges object with 24 ranges and 9 metadata columns:
#>      seqnames             ranges strand |     log2r                     ID
#>         <Rle>          <IRanges>  <Rle> | <numeric>            <character>
#>    1     chr1  3232366-247756853      * |    -0.009 TCGA-63-A5MN-01A-22D..
#>    2    chr10   431113-135214495      * |     0.009 TCGA-63-A5MN-01A-22D..
#>    3    chr11   129255-134118566      * |    -0.002 TCGA-63-A5MN-01A-22D..
#>    4    chr12  1032854-133175503      * |    -0.047 TCGA-63-A5MN-01A-22D..
#>    5    chr13 19594416-114949919      * |    -0.100 TCGA-63-A5MN-01A-22D..
#>   ..      ...                ...    ... .       ...                    ...
#>   20     chr7    416242-80526460      * |    -0.016 TCGA-63-A5MN-01A-22D..
#>   21     chr7  81108077-83692514      * |    -0.490 TCGA-63-A5MN-01A-22D..
#>   22     chr7 84772779-158367792      * |    -0.008 TCGA-63-A5MN-01A-22D..
#>   23     chr8   623652-145193933      * |    -0.041 TCGA-63-A5MN-01A-22D..
#>   24     chr9   855123-140723607      * |     0.000 TCGA-63-A5MN-01A-22D..
#>       num.mark     bstat        pval seg.median         cna         cnv
#>      <integer> <numeric>   <numeric>  <numeric> <character> <character>
#>    1      1342        NA          NA      0.006     Diploid           2
#>    2       708        NA          NA      0.023     Diploid           2
#>    3       782        NA          NA      0.014     Diploid           2
#>    4       683        NA          NA     -0.028     Diploid           2
#>    5       341        NA          NA     -0.079     Diploid           2
#>   ..       ...       ...         ...        ...         ...         ...
#>   20       370   7.31596 2.89656e-11     -0.003     Diploid           2
#>   21         5   7.29464 3.46445e-11     -0.486      HomDel           0
#>   22       394        NA          NA      0.005     Diploid           2
#>   23       564        NA          NA     -0.029     Diploid           2
#>   24       310        NA          NA      0.010     Diploid           2
#>                   gene.mame
#>                 <character>
#>    1 PRDM16;ARHGEF16;MEGF..
#>    2 DIP2C;LARP4B;GTPBP4;..
#>    3 LOC100133161;ODF3;BE..
#>    4 RAD52;ERC1;LOC100292..
#>    5 DKFZp686A1627;TUBA3C..
#>   ..                    ...
#>   20 PDGFA;PRKAR1B;HEATR2..
#>   21 HGF;CACNA2D1;PCLO;SE..
#>   22 GRM3;KIAA1324L;DMTF1..
#>   23 ERICH1;DLGAP2;CLN8;M..
#>   24 DMRT1;DMRT3;DMRT2;SM..
#>   -------
#>   seqinfo: 22 sequences from an unspecified genome; no seqlengths

## basic example code
```
