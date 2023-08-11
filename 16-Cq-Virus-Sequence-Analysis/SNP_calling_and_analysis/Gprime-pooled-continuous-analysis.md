# Gprime-pooled-continuous-analysis

Load packages

``` r
#install.packages("devtools")
#devtools::install_github("bmansfeld/QTLseqr")
# install.packages("vcfR")
#install.packages("SNPfiltR")
library(SNPfiltR)
```

    This is SNPfiltR v.1.0.1

    Detailed usage information is available at: devonderaad.github.io/SNPfiltR/ 

    If you use SNPfiltR in your published work, please cite the following papers: 

    DeRaad, D.A. (2022), SNPfiltR: an R package for interactive and reproducible SNP filtering. Molecular Ecology Resources, 22, 2443-2453. http://doi.org/10.1111/1755-0998.13618 

    Knaus, Brian J., and Niklaus J. Grunwald. 2017. VCFR: a package to manipulate and visualize variant call format data in R. Molecular Ecology Resources, 17.1:44-53. http://doi.org/10.1111/1755-0998.12549

``` r
library(QTLseqr)
library(ggplot2)
library(vcfR)
```


       *****       ***   vcfR   ***       *****
       This is vcfR 1.14.0 
         browseVignettes('vcfR') # Documentation
         citation('vcfR') # Citation
       *****       *****      *****       *****

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

#### Load in dataset

``` r
# read in vcf
DiNV_vcf <- read.vcfR("~/Desktop/KU/sequences/16Cq-DiNV-Test/freebayes/continuous-DiNV.vcf")
```

    Scanning file to determine attributes.
    File attributes:
      meta lines: 59
      header_line: 60
      variant count: 13116
      column count: 11

    Meta line 59 read in.
    All meta lines processed.
    gt matrix initialized.
    Character matrix gt created.
      Character matrix gt rows: 13116
      Character matrix gt cols: 11
      skip: 0
      nrows: 13116
      row_num: 0

    Processed variant 1000
    Processed variant 2000
    Processed variant 3000
    Processed variant 4000
    Processed variant 5000
    Processed variant 6000
    Processed variant 7000
    Processed variant 8000
    Processed variant 9000
    Processed variant 10000
    Processed variant 11000
    Processed variant 12000
    Processed variant 13000
    Processed variant: 13116
    All variants processed

``` r
# filter vcf to remove sites (SNPs) with more than two alleles present (from SNPfiltR package)
DiNV_vcf <- filter_biallelic(DiNV_vcf)
```

    2668 SNPs, 0.203% of all input SNPs, contained more than 2 alleles, and were removed from the VCF

![](Gprime-pooled-continuous-analysis_files/figure-commonmark/unnamed-chunk-2-1.png)

``` r
# extract the number of reads supporting the reference allele for each sample at each SNP 
# DiNV_vcf has a column containing RO or reference allele information
# the extract.gt function will subset from that column
# the same column contains a lot of other information separated by : and , which is why a specific function is needed 
ref <- as.data.frame(extract.gt(DiNV_vcf, element = "RO"))
# extract the number of reads supporting the alternate allele for each sample at each SNP 
# AO is alternate allele
# because I removed all SNPs that are more than biallelic, there should only be 1 alternate allele
alt <- as.data.frame(extract.gt(DiNV_vcf, element = "AO"))
DiNV_SNPs<-as.data.frame(DiNV_vcf@fix[,c(1,2,4,5)])
# add the read counts supporting the reference and alternate alleles for each of the two sequenced strains as their own columns
DiNV_SNPs$inn.ref<-ref$INN
DiNV_SNPs$inn.alt<-alt$INN
DiNV_SNPs$vir.ref<-ref$VIR
DiNV_SNPs$vir.alt<-alt$VIR

# there might be some missing data in this file, which would come out as rows where either inn or vir have an NA for alleles 
# remove any rows that have an NA 
DiNV_SNPs <- na.omit(DiNV_SNPs)
```

#### Prepare SNP dataset for Gprime format. Based off of QTLSeqr program.

``` r
# what I need is columns that say:
# chrom, position, ref base, alt base, and allele counts per sample
# but the allele counts have to be named AD_ALT.sample or AD_REF.sample 
# so I need AD_ALT.inn, AD_REF.inn, AD_ALT.vir, and AD_REF.vir

colnames(DiNV_SNPs)[5] ="AD_REF.inn"
colnames(DiNV_SNPs)[6] ="AD_ALT.inn"
colnames(DiNV_SNPs)[7] ="AD_REF.vir"
colnames(DiNV_SNPs)[8] ="AD_ALT.vir"

head(DiNV_SNPs)
```

            CHROM POS REF ALT AD_REF.inn AD_ALT.inn AD_REF.vir AD_ALT.vir
    1 NC_040699.1  43   A   C          8          0         47          1
    2 NC_040699.1  47   T   G          9          0         55          1
    3 NC_040699.1  66   T   G         10          0         80          1
    4 NC_040699.1  70   T   A         11          0         79          2
    5 NC_040699.1  79   A   C         12          0         90          1
    6 NC_040699.1  84   T   A         13          0         93          1

``` r
# looks good so far

# save this as a csv because the program doesn't seem to want it in R 
write.csv(DiNV_SNPs, file = "~/Desktop/KU/sequences/16Cq-DiNV-Test/freebayes/DiNV_SNPs_For_Gprime.csv", row.names = FALSE)
```

#### More preparing the table and loading in

``` r
# vignette says "We define the sample name for each of the bulks."
# from what I can tell, a bulk is a pooled sample or a bulked sample
# for whatever reason this program calls them a high bulk and a low bulk 
# but I'm not sure what those terms mean 

# naming the samples exactly what it says after AD_REF. 
# also want to 
HighBulk <- "inn"
LowBulk <- "vir"
# assuming this is the way they want to represent chromosomes
# they usually want you to specify which chromosomes to use
# because there is just one it doesn't matter probably 
Chroms <- "NC_040699.1"

# see if importing the table works 
Gprime_df <- importFromTable(file = "~/Desktop/KU/sequences/16Cq-DiNV-Test/freebayes/DiNV_SNPs_For_Gprime.csv",
                      highBulk = HighBulk,
                      lowBulk = LowBulk,
                      chromList = Chroms)
```

    Removing the following chromosomes: 

    Renaming the following columns: AD_REF.inn, AD_ALT.inn

    Renaming the following columns: AD_REF.vir, AD_ALT.vir

``` r
# for whatever reason it renamed my columns, so I will just have to remember that high is innubila cells and low is virilis cells
```

The import from table function has added in some columns. They mean:

- GQ.HIGH - The genotype quality score, (how confident we are in the
  genotyping)
- SNPindex.HIGH - The calculated SNP-index for the high bulk. SNP index
  is alternate allele depth over total read depth
- Same as above for the low bulk
- REF_FRQ - reference allele depth for both samples divided by total
  read depth for both samples
- deltaSNP - The change in SNP index, high bulk index minus low bulk
  index

#### Looking at SNPs to get an idea for filtering

“QTLseqr offers some options for filtering that may help reduce noise
and improve results. Filtering is mainly based on read depth for each
SNP, such that we can try to eliminate SNPs with low confidence, due to
low coverage, and SNPs that may be in repetitive regions and thus have
inflated read depth. You can also filter based on the absolute
difference in read depth between the bulks.”

``` r
#One way to assess filtering thresholds and check the quality of our data is by plotting histograms of the read depths. We can get an idea of where to draw our thresholds

# total read depth histogram:

ggplot(data = Gprime_df) + 
    geom_histogram(aes(x = DP.HIGH + DP.LOW), bins = 50) + 
    xlim(0,500)
```

    Warning: Removed 2 rows containing missing values (`geom_bar()`).

![](Gprime-pooled-continuous-analysis_files/figure-commonmark/unnamed-chunk-5-1.png)

``` r
# total reference allele frequency:
ggplot(data = Gprime_df) +
    geom_histogram(aes(x = REF_FRQ), bins = 50)
```

![](Gprime-pooled-continuous-analysis_files/figure-commonmark/unnamed-chunk-5-2.png)

``` r
# this one looks different than their example, their ref_frq was centered around 0.5... but maybe their data is really different. This looks like the vast majority of the SNPs are reference allele? Not sure exactly what that means...


# We can plot our per-bulk SNP-index to check if our data is good.We expect to find two small peaks on each end and most of the SNPs should be approximiately normally distributed arround 0.5 in an F2 population. Here is the HIGH bulk for example:

# for innubila sample
ggplot(data = Gprime_df) +
    geom_histogram(aes(x = SNPindex.HIGH, bins = 50))
```

    Warning in geom_histogram(aes(x = SNPindex.HIGH, bins = 50)): Ignoring unknown
    aesthetics: bins

    `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Gprime-pooled-continuous-analysis_files/figure-commonmark/unnamed-chunk-5-3.png)

``` r
# for vir sample
ggplot(data = Gprime_df) +
    geom_histogram(aes(x = SNPindex.LOW), bins = 50)
```

![](Gprime-pooled-continuous-analysis_files/figure-commonmark/unnamed-chunk-5-4.png)

``` r
# this is not an F2 population (not even possible with our virus) so this is maybe why these don't look as expected 
```

#### Filtering SNPs

“Now that we have an idea about our read depth distribution we can
filter out low confidence SNPS. In general we recommend filtering
extremely low and high coverage SNPs, either in both bulks
(`minTotalDepth/maxTotalDepth`) and/or in each bulk separately
(`minSampleDepth`). We have the option to filter based on reference
allele frequency (`refAlleleFreq`), this removes SNPs that for some
reason are over- or under-represented in *BOTH* bulks. We can also
filter SNPs that have large discrepancies in read depth between the
bulks(i.e. one bulk has a depth of 500 and the other has 5). Such
discrepancies can throw off the G statistic. We can also use the GATK GQ
score (Genotype Quality) to filter out low confidence SNPs. If the
`verbose` parameter is set to `TRUE` (default) the function will report
the numbers of SNPs filtered in each step.”

Additionally I added a minor allele count filter

``` r
# First filter for minor allele count 
# first determine which allele is minor
# by comparing the sum of the reference and sum of the alternates
# whichever is smaller is minor 
# and then if the minor is less than 2 reads
# get rid of that row 

# empty vector for positions of SNPs to remove
POS_to_remove <- c()

for ( i in 1:nrow(Gprime_df)) {
  if (Gprime_df$AD_ALT.LOW[i] + Gprime_df$AD_ALT.HIGH[i] <= Gprime_df$AD_REF.HIGH[i] + Gprime_df$AD_REF.LOW[i] & 
      Gprime_df$AD_ALT.LOW[i] + Gprime_df$AD_ALT.HIGH[i] <= 2) {
    POS_to_remove[i] <- 1
  }
  else if (Gprime_df$AD_REF.LOW[i] + Gprime_df$AD_REF.HIGH[i] <= Gprime_df$AD_ALT.HIGH[i] + Gprime_df$AD_ALT.LOW[i] & 
      Gprime_df$AD_REF.LOW[i] + Gprime_df$AD_REF.HIGH[i] <= 2) {
    POS_to_remove[i] <- 1
  }

  else {POS_to_remove[i] <- 0 }
}

# how many get removed 
table(POS_to_remove)
```

    POS_to_remove
       0    1 
     690 9747 

``` r
# 9747 to remove 

# keep all rows where minor allele depth is more than 2 (aka was coded 0)
trimmed_Gprime_df <- Gprime_df[POS_to_remove == 0,]



# use the filtering snps 
# minimum depth should probably be 10? This is total so both samples, basically same as min sample depth
# max depth of 200? not sure 
# don't have GATK score so not filtering by that 
# min coverage per sample of 5 reads, seems pretty low but there are inn samples that I'm sure have few reads 
# for refAlleleFreq this is keeping all SNPs with a REF_FRQ between 0.05 and 0.95
 

Gprime_df_filt <-
    filterSNPs(
        SNPset = trimmed_Gprime_df,
        refAlleleFreq = 0.05,
        minTotalDepth = 10,
        maxTotalDepth = 200, 
        minSampleDepth = 5,
        verbose = TRUE
    )
```

    Filtering by reference allele frequency: 0.05 <= REF_FRQ <= 0.95

    ...Filtered 68 SNPs

    Filtering by total sample read depth: Total DP >= 10

    ...Filtered 2 SNPs

    Filtering by total sample read depth: Total DP <= 200

    ...Filtered 2 SNPs

    Filtering by per sample read depth: DP >= 5

    ...Filtered 15 SNPs

    Original SNP number: 690, Filtered: 87, Remaining: 603

#### Running G prime analysis

``` r
# use Hampel method for outlier filter 
Gprime_df_filt_prime <- runGprimeAnalysis(Gprime_df_filt,
    windowSize = 2000,
    outlierFilter = "Hampel",
    filterThreshold = 0.1)
```

    Counting SNPs in each window...

    Calculating tricube smoothed delta SNP index...

    Calculating G and G' statistics...

    Warning: There were 41 warnings in `dplyr::mutate()`.
    The first warning was:
    ℹ In argument: `Gprime = tricubeStat(POS = POS, Stat = G, windowSize =
      windowSize, ...)`.
    ℹ In group 1: `CHROM = NC_040699.1`.
    Caused by warning in `lfproc()`:
    ! procv: no points with non-zero weight
    ℹ Run `dplyr::last_dplyr_warnings()` to see the 40 remaining warnings.

    Using Hampel's rule to filter outlier regions

    Estimating the mode of a trimmed G prime set using the 'modeest' package...

    Calculating p-values...

    Warning: There was 1 warning in `dplyr::mutate()`.
    ℹ In argument: `pvalue = getPvals(...)`.
    Caused by warning:
    ! encountered a tie, and the difference between minimal and 
                       maximal value is > length('x') * 'tie.limit'
    the distribution could be multimodal

``` r
# there are some errors with this, but I will continue anyways 


#Due to the fact that p-values are estimated from the null distribution of G', an important check is to see if the null distribution of G' values is close to log normally distributed. For this purpose we use the `plotGprimeDist` function, which plots the G' histograms of both raw and filtered G' sets alongside the log-normal null distribution (which is reported in the legend). We can also use this to test which filtering method (Hampel or DeltaSNP) estimates a more accurate null distribution. If you use the `"deltaSNP"` method plotting G' distributions with different filter thresholds might also help reveal a better G' null distribution. 

# Hampel method
plotGprimeDist(SNPset = Gprime_df_filt_prime, outlierFilter = "Hampel")
```

    Warning: encountered a tie, and the difference between minimal and 
                       maximal value is > length('x') * 'tie.limit'
    the distribution could be multimodal

    Warning: Removed 2 rows containing missing values (`geom_bar()`).
    Removed 2 rows containing missing values (`geom_bar()`).

![](Gprime-pooled-continuous-analysis_files/figure-commonmark/unnamed-chunk-7-1.png)

``` r
# deltaSNP method 
plotGprimeDist(SNPset = Gprime_df_filt_prime, outlierFilter = "deltaSNP", filterThreshold = 0.1)
```

    Warning: Removed 2 rows containing missing values (`geom_bar()`).
    Removed 2 rows containing missing values (`geom_bar()`).

![](Gprime-pooled-continuous-analysis_files/figure-commonmark/unnamed-chunk-7-2.png)

``` r
# The hampel method looks most like a log normal distribution

# plot number of SNPs in a window
p1 <- plotQTLStats(SNPset = Gprime_df_filt_prime, var = "nSNPs")
p1
```

![](Gprime-pooled-continuous-analysis_files/figure-commonmark/unnamed-chunk-7-3.png)

``` r
# plot the G' statistic
# try to plot the threshold can pass the FDR (q) of 0.01.
p3 <- plotQTLStats(SNPset = Gprime_df_filt_prime, var = "Gprime", plotThreshold = TRUE, q = 0.01)

#plot
p3
```

![](Gprime-pooled-continuous-analysis_files/figure-commonmark/unnamed-chunk-7-4.png)

``` r
# visualize with neg log ten p-value instead
QTLplots <- plotQTLStats(
    SNPset = Gprime_df_filt_prime, 
    var = "negLog10Pval", 
    plotThreshold = TRUE, 
    q = 0.01, 
    )
QTLplots
```

![](Gprime-pooled-continuous-analysis_files/figure-commonmark/unnamed-chunk-7-5.png)

From this, I feel confident about the SNPs put into the program, so I
think these 2 maybe 3 “QTLs” are accurate

#### Extracting peak data

``` r
# generate a list of datasets with each QTL
QTL <- getSigRegions(SNPset = Gprime_df_filt_prime, alpha = 0.01)
```

    Adding missing grouping variables: `qtl`

``` r
# there are 3 QTLs
# make into dfs
QTL1 <- QTL[[1]]
QTL2 <- QTL[[2]]
QTL3 <- QTL[[3]]

# combine the df
QTLS <- rbind(QTL1, QTL2, QTL3)
# 16 significant SNPs 
QTLS
```

       qtl       CHROM    POS REF ALT AD_REF.LOW AD_ALT.LOW DP.LOW SNPindex.LOW
    1    1 NC_040699.1 122041   C   A         49          6     55    0.1090909
    2    1 NC_040699.1 122377   C   T          0         48     48    1.0000000
    3    1 NC_040699.1 122463   C   T          0         48     48    1.0000000
    4    1 NC_040699.1 123223   G   A          1        105    106    0.9905660
    5    1 NC_040699.1 123319   G   T         51         32     83    0.3855422
    6    1 NC_040699.1 123601   G   T          1         73     74    0.9864865
    7    2 NC_040699.1 138144   T   C          0         56     56    1.0000000
    8    2 NC_040699.1 138524   A   T         33          0     33    0.0000000
    9    2 NC_040699.1 138891   T   C         24         27     51    0.5294118
    10   2 NC_040699.1 139038   C   T         43          0     43    0.0000000
    11   2 NC_040699.1 139171   T   G          0         42     42    1.0000000
    12   2 NC_040699.1 139398   T   C          0         41     41    1.0000000
    13   2 NC_040699.1 139665   A   G         57          0     57    0.0000000
    14   2 NC_040699.1 139848   A   C          0         44     44    1.0000000
    15   2 NC_040699.1 139876   A   G          0         46     46    1.0000000
    16   3 NC_040699.1 144408   C   A          0         45     45    1.0000000
       AD_REF.HIGH AD_ALT.HIGH DP.HIGH SNPindex.HIGH   REF_FRQ   deltaSNP nSNPs
    1           23           0      23    0.00000000 0.9230769 -0.1090909     4
    2            8           7      15    0.46666667 0.1269841 -0.5333333     5
    3           18           6      24    0.25000000 0.2500000 -0.7500000     5
    4           17           2      19    0.10526316 0.1440000 -0.8853029     5
    5           24           0      24    0.00000000 0.7009346 -0.3855422     5
    6           17           2      19    0.10526316 0.1935484 -0.8812233     3
    7           12           3      15    0.20000000 0.1690141 -0.8000000    12
    8            4           8      12    0.66666667 0.8222222  0.6666667    10
    9           18           0      18    0.00000000 0.6086957 -0.5294118    11
    10          11           5      16    0.31250000 0.9152542  0.3125000     9
    11          10           3      13    0.23076923 0.1818182 -0.7692308     9
    12          18           1      19    0.05263158 0.3000000 -0.9473684    12
    13           8           9      17    0.52941176 0.8783784  0.5294118    11
    14          20           1      21    0.04761905 0.3076923 -0.9523810    13
    15          22           0      22    0.00000000 0.3235294 -1.0000000    13
    16          14           0      14    0.00000000 0.2372881 -1.0000000     4
       tricubeDeltaSNP        G   Gprime       pvalue negLog10Pval      qvalue
    1       -0.5959731      NaN 77.32953 3.965540e-05     4.401698 0.003125050
    2       -0.5946541      NaN 76.00459 4.511755e-05     4.345654 0.003125050
    3       -0.6003328      NaN 75.66546 4.664254e-05     4.331218 0.003125050
    4       -0.6505164 78.93543 72.31448 6.509037e-05     4.186483 0.003694139
    5       -0.6568553      NaN 71.81011 6.849154e-05     4.164363 0.003694139
    6       -0.6162488 68.00551 70.32853 7.964148e-05     4.098861 0.003694139
    7       -0.6268032      NaN 66.55629 1.179157e-04     3.928428 0.004443948
    8       -0.5837769      NaN 79.14389 3.330203e-05     4.477529 0.003125050
    9       -0.5422226      NaN 81.57515 2.645196e-05     4.577542 0.003125050
    10      -0.5255783      NaN 82.54898 2.414895e-05     4.617102 0.003125050
    11      -0.5355115      NaN 83.43007 2.225087e-05     4.652653 0.003125050
    12      -0.5672328      NaN 84.93387 1.937259e-05     4.712812 0.003125050
    13      -0.6045437      NaN 83.61283 2.187765e-05     4.659999 0.003125050
    14      -0.6301164      NaN 70.38455 7.918589e-05     4.101352 0.003694139
    15      -0.6340291      NaN 68.36056 9.758525e-05     4.010616 0.004203136
    16      -0.4960274      NaN 66.59259 1.174643e-04     3.930094 0.004443948

#### Plot scatter plot of negative log p-values

``` r
# another way of visualizing the data
ggplot(Gprime_df_filt_prime, aes(x=POS, y=negLog10Pval))+
  geom_point()+
  theme_classic() + xlab("Position") + ylab("Negative log10 p-value")
```

![](Gprime-pooled-continuous-analysis_files/figure-commonmark/unnamed-chunk-9-1.png)

#### Plot G prime statistic as scatter plot

``` r
# another way of visualizing the data
# tried adding the FDR line, I'm not sure the value but it seemed like it was at 67 in the Gprime plot from the function? 
# I just really hate the grey background in that plot
# but the smooth lines are nice 
# IDK best way to visualize this 

ggplot(Gprime_df_filt_prime, aes(x=POS, y=Gprime))+
  geom_point()+
  theme_classic() + xlab("Position") + ylab("G'") +
   geom_hline(yintercept = 67, color="deeppink")
```

![](Gprime-pooled-continuous-analysis_files/figure-commonmark/unnamed-chunk-10-1.png)

save filtered dataset from From this analysis, want to use it for
Fisher’s exact test, as well as save the list of SNPs found in the QTLs

``` r
#write.csv(Gprime_df_filt, "~/Desktop/KU/sequences/16Cq-DiNV-Test/freebayes/Gprime_df_filt.csv")
#write.csv(QTLS, "~/Desktop/KU/sequences/16Cq-DiNV-Test/freebayes/QTL_SNPs.csv")
```
