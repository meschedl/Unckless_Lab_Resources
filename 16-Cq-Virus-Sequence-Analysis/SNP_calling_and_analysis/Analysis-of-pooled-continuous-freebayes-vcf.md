Analysis of pooled-continuous freebayes vcf file
================
2023-07-25

#### Load in libraries needed for analysis

``` r
# load library
# install.packages("vcfR")
library(vcfR)
```

    ## 
    ##    *****       ***   vcfR   ***       *****
    ##    This is vcfR 1.14.0 
    ##      browseVignettes('vcfR') # Documentation
    ##      citation('vcfR') # Citation
    ##    *****       *****      *****       *****

``` r
# install.packages("SNPfiltR")
library(SNPfiltR)
```

    ## This is SNPfiltR v.1.0.1
    ## 
    ## Detailed usage information is available at: devonderaad.github.io/SNPfiltR/ 
    ## 
    ## If you use SNPfiltR in your published work, please cite the following papers: 
    ## 
    ## DeRaad, D.A. (2022), SNPfiltR: an R package for interactive and reproducible SNP filtering. Molecular Ecology Resources, 22, 2443-2453. http://doi.org/10.1111/1755-0998.13618 
    ## 
    ## Knaus, Brian J., and Niklaus J. Grunwald. 2017. VCFR: a package to manipulate and visualize variant call format data in R. Molecular Ecology Resources, 17.1:44-53. http://doi.org/10.1111/1755-0998.12549

``` r
library(ggplot2)
library(ggVennDiagram)
library(units)
```

    ## udunits database from /Library/Frameworks/R.framework/Versions/4.2/Resources/library/units/share/udunits/udunits2.xml

#### Read in vcf, filter it, and generate a dataframe that includes the position, allele counts, and bases for each SNP separated by reference and alternate allele. Only considering biallelic SNPs here.

``` r
# read in vcf
DiNV_vcf <- read.vcfR("~/Desktop/KU/sequences/16Cq-DiNV-Test/freebayes/continuous-DiNV.vcf")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 59
    ##   header_line: 60
    ##   variant count: 13116
    ##   column count: 11
    ## Meta line 59 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 13116
    ##   Character matrix gt cols: 11
    ##   skip: 0
    ##   nrows: 13116
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant 3000Processed variant 4000Processed variant 5000Processed variant 6000Processed variant 7000Processed variant 8000Processed variant 9000Processed variant 10000Processed variant 11000Processed variant 12000Processed variant 13000Processed variant: 13116
    ## All variants processed

``` r
# filter vcf to remove sites (SNPs) with more than two alleles present (from SNPfiltR package)
DiNV_vcf <- filter_biallelic(DiNV_vcf)
```

    ## 2668 SNPs, 0.203% of all input SNPs, contained more than 2 alleles, and were removed from the VCF

![](Analysis-of-pooled-continuous-freebayes-vcf_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

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

# ref and alt are just SNPs and the number of reads
# But I want to also add other information with these, like what are the actual bases for the alleles, and have a column for the position 
# that information is in the fix portion of the vcf file, which is hard to look at 
# check what it looks like by naming it and looking at it 
x<-DiNV_vcf@fix
# after looking at x, I can see that the columns I want are 1: chromosome, 2: position, 4: reference allele, and 5: alternate allele
# subset out those columns 
DiNV_SNPs<-as.data.frame(DiNV_vcf@fix[,c(1,2,4,5)])
# add the read counts supporting the reference and alternate alleles for each of the two sequenced strains as their own columns
DiNV_SNPs$inn.ref<-ref$INN
DiNV_SNPs$inn.alt<-alt$INN
DiNV_SNPs$vir.ref<-ref$VIR
DiNV_SNPs$vir.alt<-alt$VIR

# there might be some missing data in this file, which would come out as rows where either inn or vir have an NA for alleles 
# remove any rows that have an NA 
DiNV_SNPs <- na.omit(DiNV_SNPs)

# additionally, the position needs to be read as a number by R
# turn it into numeric
head(DiNV_SNPs) # check numbers before
```

    ##         CHROM POS REF ALT inn.ref inn.alt vir.ref vir.alt
    ## 1 NC_040699.1  43   A   C       8       0      47       1
    ## 2 NC_040699.1  47   T   G       9       0      55       1
    ## 3 NC_040699.1  66   T   G      10       0      80       1
    ## 4 NC_040699.1  70   T   A      11       0      79       2
    ## 5 NC_040699.1  79   A   C      12       0      90       1
    ## 6 NC_040699.1  84   T   A      13       0      93       1

``` r
DiNV_SNPs$POS <- as.numeric(DiNV_SNPs$POS)
head(DiNV_SNPs) # check to make sure numbers stay the same, sometimes as.numeric can change numbers
```

    ##         CHROM POS REF ALT inn.ref inn.alt vir.ref vir.alt
    ## 1 NC_040699.1  43   A   C       8       0      47       1
    ## 2 NC_040699.1  47   T   G       9       0      55       1
    ## 3 NC_040699.1  66   T   G      10       0      80       1
    ## 4 NC_040699.1  70   T   A      11       0      79       2
    ## 5 NC_040699.1  79   A   C      12       0      90       1
    ## 6 NC_040699.1  84   T   A      13       0      93       1

#### Read in Filtered dataframe from G prime analysis \[LINK\] to see how it was filtered

``` r
filtered_SNPs <- read.csv("~/Desktop/KU/sequences/16Cq-DiNV-Test/freebayes/Gprime_df_filt.csv")

# there are a lot of columns in this that I don't need (from the G prime statistics, etc) and some columns should be labeled better
# keep all rows
# want columns 
filtered_SNPs <- filtered_SNPs[,c(2:7,10:11)]

# need to rename the columns 
colnames(filtered_SNPs)[5] ="vir.ref"
colnames(filtered_SNPs)[6] ="vir.alt"
colnames(filtered_SNPs)[7] ="inn.ref"
colnames(filtered_SNPs)[8] ="inn.alt"
```

#### From the SNPs, run a Fisher’s Exact Test on each one and determine the p-value. “Fisher’s exact test is a statistical test used to determine if there are nonrandom associations between two categorical variables.” Basically, I want to see if the number of alleles is specifically associated with one sample type or the other.

<https://mathworld.wolfram.com/FishersExactTest.html>

``` r
# this will need to be run in a loop because there are 10,437 SNPs and each one needs an individual dataframe 
# and a test run on them 
# starting by testing the process with just one SNP
# start by making the first column be the reference allele: reference=
# give reference a list of the numbers I want it to be: the first allele count entry in DiNV_SNPs for innubila reference and virilis reference 
# then the second column is the alternative allele: alt=
# and do the same thing, the first entry in DiNV_SNPs for innubila alternate and virilis alternate 
# then I give the dataframe row names, which are in the same order as the allele counts are specified 
test_df <- data.frame(reference=as.numeric(c(filtered_SNPs$inn.ref[1],filtered_SNPs$vir.ref[1])), 
                alt=as.numeric(c(filtered_SNPs$inn.alt[1],filtered_SNPs$vir.alt[1])),
                row.names = c("inn","vir"))

# use this test_df to perform a Fisher's Exact Test
test1<-fisher.test(test_df)
# look at the p-value
test1$p.value
```

    ## [1] 0.005952381

``` r
# check the test_df
head(test_df)
```

    ##     reference alt
    ## inn         3   8
    ## vir         0  45

A p-value of 0.006 kind of makes sense here, becuase inn sample has ref
reads, which vir does not. So they are different. But this isn’t very
different, and might not show up as significant when doing multiple
testing.

#### Scale up making the test dataframe for each SNP, and run the Fisher’s Exact Test, and extract p-values and add them to the DiNV_SNPs dataframe

``` r
# open an empty vector to hold each p-value
# this is so R has somewhere to put the p-values
pval<-c()
# for loop to iterate over each row in the dataframe named 'DiNV_SNPs' one after the other and record p-value for each row (ie, SNP)
# usesi in 1 through nrow because it will use the exact number of rows in the SNP dataframe
for (i in 1:nrow(filtered_SNPs)){
  # make 1 data frame for the given SNP with numbers of alleles 
  df1<-data.frame(reference=as.numeric(c(filtered_SNPs$inn.ref[i],filtered_SNPs$vir.ref[i])),
                alt=as.numeric(c(filtered_SNPs$inn.alt[i],filtered_SNPs$vir.alt[i])),
                row.names = c("inn","vir"))
  #perform a fisher's exact test for above dataframe
  test2<-fisher.test(df1)
  # store p value for this test in our vector 'pval'
  # the [i] tells R to store the value in the i-th position
  pval[i]<-test2$p.value
}

# check out my p-values

hist(pval, breaks = 50)
```

![](Analysis-of-pooled-continuous-freebayes-vcf_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# add in the pval vector to the DiNV_SNPs and name it 
filtered_SNPs$fisher_pval<-pval

# check  new DiNV_SNPs to make sure it looks right
head(filtered_SNPs)
```

    ##         CHROM  POS REF ALT vir.ref vir.alt inn.ref inn.alt  fisher_pval
    ## 1 NC_040699.1  499   G   T       0      45       3       8 5.952381e-03
    ## 2 NC_040699.1  810   T   A      20       0       9       4 1.747312e-02
    ## 3 NC_040699.1 1033   T   C       0      56      12       9 7.998252e-09
    ## 4 NC_040699.1 1069   C   A       0      49      13       3 3.410252e-11
    ## 5 NC_040699.1 1147   C   T       0      68       8       2 1.919232e-09
    ## 6 NC_040699.1 1266   T   C       0      13       7       1 6.879945e-05

Interestingly, many p-values are low here.

#### Make a scatterplot of the p-values on the y-axis and the position in the DiNV genome on the x-axis. This is basically a manhattan plot. The p-value is transformed into the negative log of the p-value. And a multiple test correction (bonferroni) line is applied to the plot to show SNPs that are significant above that line.

``` r
# first want to transform the p-value in to the negative log 10 so that really small p-values get plotted as larger 
# create a new column in the dataframe 
filtered_SNPs$neg_log_pval <- (log10(filtered_SNPs$fisher_pval))*-1

# get multiple testing corrected p value (0.05/ number of tests)
MTC <- 0.05/nrow(filtered_SNPs)

# plot scatter plot and add in the line for the MTC (must take negative log10 of that value)
ggplot(filtered_SNPs, aes(x=POS, y=neg_log_pval))+
  geom_point()+
  geom_hline(yintercept = -log10(MTC), color="blueviolet")+
  theme_classic() + xlab("Position") + ylab("Negative log10 p-value")
```

![](Analysis-of-pooled-continuous-freebayes-vcf_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# which SNPs have a p-value less than the multiple testing threshold
sig.snps<-filtered_SNPs[filtered_SNPs$fisher_pval < MTC,]
```

Most p-values here are are significant, but I only started with 603
SNPs. Out of those 404 are significant with a fisher’s exact test and
multiple testing. I might need to do a different multiple testing
threshold

#### Looking for the SNPs from Hill and Unckless (2020) (<https://elifesciences.org/articles/58931>) in this data

``` r
# list of positions of important SNPs (12)
snps <- c(14249,41210,42389,59194,59275,59276,66615,78978,78991,126118,132593,140117)

# which SNPs are present in the DiNV_SNPs (not just significant ones)
table(snps %in% filtered_SNPs$POS)
```

    ## 
    ## FALSE  TRUE 
    ##    11     1

``` r
# 6 are SNPs, and 6 were not found as SNPs

# isolate the details of the found SNPs
# separate out the 6 that are present in DiNV_SNPs
found.snps <- snps[snps %in% filtered_SNPs$POS]
# show the p-values and allele frequency differences between the samples for these SNPs of interest
filtered_SNPs[filtered_SNPs$POS %in% found.snps,]
```

    ##           CHROM   POS REF ALT vir.ref vir.alt inn.ref inn.alt fisher_pval
    ## 252 NC_040699.1 66615   G   A       0      11       3       3  0.02941176
    ##     neg_log_pval
    ## 252     1.531479

There is only one retained SNP in the filtered dataframe, and it is not
considered significantly different between the two samples with the
Fisher’s exact test. If those SNPs were in the original dataset, they
got filtered. Probably by minor allele count because that removed most
SNPs in the dataset.

#### Another thing to look at is the raw allele frequency differences between the two samples

``` r
# allele frequency is the reads for the allele divided by all alleles 
# ex. alternate allele frequency is the number of alternate alleles divided by the number of alt alleles plus reference alleles
# calculate frequency of alternate allele in innubila
freq.in <- as.numeric(filtered_SNPs$inn.alt)/(as.numeric(filtered_SNPs$inn.alt)+as.numeric(filtered_SNPs$inn.ref))
# calculate frequency of alternate allele in virilis
freq.vir <- as.numeric(filtered_SNPs$vir.alt)/(as.numeric(filtered_SNPs$vir.alt)+as.numeric(filtered_SNPs$vir.ref))
# calculate allele frequency difference between the two samples by subtracting innubila frequency from virilis 
filtered_SNPs$af.dif<-abs(freq.in-freq.vir)
# how many SNPs are fixed for the alternate allele in innubila?
table(freq.in == 1)
```

    ## 
    ## FALSE  TRUE 
    ##   590    13

``` r
# how many SNPs are fixed for the alternate allele in virilis?
table(freq.vir == 1)
```

    ## 
    ## FALSE  TRUE 
    ##   222   381

``` r
# how many SNPs are fixed different between the samples?
table(filtered_SNPs$af.dif == 1)
```

    ## 
    ## FALSE  TRUE 
    ##   506    97

``` r
# plot a histogram of the allele freq divergence from the ref genome for innubila with mean value highlighted by vertical red line
hist(freq.in, breaks=100, xlab = "allele frequency divergence from reference genome")
abline(v=mean(freq.in), col="red")
```

![](Analysis-of-pooled-continuous-freebayes-vcf_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# plot a histogram of the allele freq divergence from the ref genome for virilis with mean value highlighted by vertical red line
hist(freq.vir, breaks=100, xlab = "allele frequency divergence from reference genome")
abline(v=mean(freq.vir), col="red")
```

![](Analysis-of-pooled-continuous-freebayes-vcf_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

Interesting here that many more SNPs are fixed alternate in virilis
cells than in innubila cells. I’m not sure what to make of this. 13 SNPs
are fixed alternate in innubila cell virus, while 381 SNPs are fixed
alternate in virilis cell virus. 97 SNPs are fixed different between the
two virus samples.

#### Plot manhattan plot with fixed differences between the samples highlighted in blue

``` r
#isolate a dataframe containing only those fixed SNPs
fixed <- filtered_SNPs[filtered_SNPs$af.dif== 1,]

#plot overlaid onto same plot as above, just add in another geom_point
ggplot(filtered_SNPs, aes(x=POS, y=neg_log_pval))+
  geom_point()+
  geom_point(data=fixed,
             aes(x=POS,y=neg_log_pval), color = "deepskyblue", size=3)+
  geom_hline(yintercept = -log10(MTC), color="deeppink")+
  theme_classic() + xlab("Position") + ylab("Negative log10 p-value")
```

![](Analysis-of-pooled-continuous-freebayes-vcf_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Most fixed SNPs are above the significance line, which makes sense. But
there does not seem to be a pattern to them. And some of the highest
p-values are not fixed. (But they could be close)

#### Look at the raw read counts for highly significant SNPs and plot them

``` r
# make an adjusted pval column for my significant SNPs
sig.snps$adjusted_pval <- sig.snps$fisher_pval*603

# order the sig.snps dataframe 
sig.snps <- sig.snps[order(sig.snps$adjusted_pval),]

# separate out the lowesst 25 adjusted pvalues
top_25_sig_snps <- sig.snps[1:25,]

# separate out the table above into 4 tables, 1 for each allele, and the combine them together
# note for this all the separate dfs are going to have inn_ref as the column name, but not contents
# it will get changed below

# make empty df with column names but no rows
inn_ref <- top_25_sig_snps[0,c(2:4,7, 9:11)]
# make in.ref dataframe
  for(i in 1:nrow(top_25_sig_snps)){
  #isolate the row (SNP) you want to work on
  inn_ref[i,]<-top_25_sig_snps[i,c(2:4,7,9:11)]
  }
inn_ref$cell_type <- rep("innubila", times = nrow(inn_ref))
inn_ref$allele <- rep("reference", times = nrow(inn_ref))

inn_alt <- top_25_sig_snps[0,c(2:4,7, 9:11)]
# make in.ref dataframe
  for(i in 1:nrow(top_25_sig_snps)){
  #isolate the row (SNP) you want to work on
  inn_alt[i,]<-top_25_sig_snps[i,c(2:4, 8:11)]
  }
# make cell_type and allele type columns 
inn_alt$cell_type <- rep("innubila", times = nrow(inn_alt))
inn_alt$allele <- rep("alternate", times = nrow(inn_alt))

vir_ref <- top_25_sig_snps[0,c(2:4,7, 9:11)]
# make in.ref dataframe
  for(i in 1:nrow(top_25_sig_snps)){
  #isolate the row (SNP) you want to work on
  vir_ref[i,]<-top_25_sig_snps[i,c(2:5,9:11)]
  }
vir_ref$cell_type <- rep("Dv-1", times = nrow(vir_ref))
vir_ref$allele <- rep("reference", times = nrow(vir_ref))

vir_alt <- top_25_sig_snps[0,c(2:4,7, 9:11)]
# make in.ref dataframe
  for(i in 1:nrow(top_25_sig_snps)){
  #isolate the row (SNP) you want to work on
  vir_alt[i,]<-top_25_sig_snps[i,c(2:4,6, 9:11)]
  }
vir_alt$cell_type <- rep("Dv-1", times = nrow(vir_alt))
vir_alt$allele <- rep("alternate", times = nrow(vir_alt))

# combine the dfs
# and rename the inn_ref column (which actually has every count) as count
top_25_snps_tidy <- rbind(inn_alt, inn_ref, vir_alt, vir_ref)
colnames(top_25_snps_tidy)[4] <- "count"

top_25_snps_tidy$count <- as.numeric(top_25_snps_tidy$count)

ggplot(data=top_25_snps_tidy, aes(x=allele, y=count, fill = cell_type)) +
  geom_bar(stat="identity", position=position_dodge()) + scale_y_continuous(limits = c(0, 150))+ facet_wrap(~POS, nrow = 5) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + geom_text(aes(label=count), position=position_dodge(width=0.8), vjust=-0.35)
```

![](Analysis-of-pooled-continuous-freebayes-vcf_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Most of the most significant SNPs are where the dv-1 virus has mostly
fixed alternate reads, and the innubila virus has mostly fixed reference
reads.

#### Compare the significant SNPs found here with the QTL SNPs from the G prime analysis

``` r
# make vectors of the positions of the significant SNPs
Gprime_SNPS <- read.csv("~/Desktop/KU/sequences/16Cq-DiNV-Test/freebayes/QTL_SNPs.csv")

Gprime_SNP_POS <- Gprime_SNPS$POS
Fisher_SNP_POS <- sig.snps$POS

# combine those vectors into a list
c <- list(Gprime_SNP_POS, Fisher_SNP_POS)
# use function to make a venn diagram from that list
ggVennDiagram(c, category.names = c("Gprime QTL SNPs", "Fisher Exact Test SNPs"), label = "count") + scale_fill_gradient(low = "lavenderblush1", high = "lightcyan1") +
 theme(legend.position = "none")
```

![](Analysis-of-pooled-continuous-freebayes-vcf_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

#### What are the overlapping SNPs?

``` r
# what are the shared SNPs?
shared_SNPS <- Gprime_SNPS[Gprime_SNPS$POS %in% Fisher_SNP_POS,]

#write.csv(shared_SNPS, "~/Desktop/KU/sequences/16Cq-DiNV-Test/freebayes/shared_sig_snps.csv")
```

Go through each of these SNPs, find their locations and if they change
the amino acid codon.
