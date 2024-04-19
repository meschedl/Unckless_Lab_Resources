# 20231025-male-female-rep3-16Cq-nanoject

Load in packages needed for the analysis

``` r
library("survival")
library("survminer")
```

    Loading required package: ggplot2

    Loading required package: ggpubr


    Attaching package: 'survminer'

    The following object is masked from 'package:survival':

        myeloma

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
library(tidyr)
library(stringr)
library(AICcmodavg)
```

### Loop to convert the example data.frame ‘df’ into properly formatted data.frame ‘results’

``` r
#write a function to transform a data.frame that has the column format 'vial | treatment | D0 | D1 | D2...', with one row for each vial
#into a long version in tidy format that can be input to make a survivorship curve
convert_df<-function(df){
  #open empty data.frame to store results
  results<-data.frame(vial=character(),treatment=character(),dead=numeric(),status=numeric())
  #This loop will index out each row (one row per vial) one at a time, transform it into long format (one row per fly), and add the information to the empty data.frame called results
  for(i in 1:nrow(df)){
  #isolate the row (vial) you want to work on
  temp<-df[i,]
    #iteratively isolate each day for this vial (day 0 must be column 3, day 1 column 4, etc.). Loop stops the column before the last day
    for(j in 3:(ncol(temp)-1)){
      #assign the number of flies that died in the vial on that day (starting with day 1) to the variable 'z'
      z<-temp[1,j]-temp[1,j+1]
        #if >0 flies died add this information to the results dataframe
        if(z>0){
          #iterate over 1 through total number of dead flies
          for(k in 1:z){
            #add a new row to the 'results' data.frame for the given dead fly, specifying vial #, treatment, day died, and
            #record the current vial #
            vial<-temp[,1]
            #record the genotype of the current vial
            treatment<-temp[,2]
            #record the death date of the flies that died on this day (assumes that your input DF starts with day 0 in column 3)
            dd<-j-2
            #append this information into a new row in the 'results' data.frame, and add a '1' in the 4th column to indicate mortality
            results[nrow(results)+1,]<- c(vial,treatment,dd,1)
          } #close for loop
        } #close if loop
    } #close for loop
  
  #now assign the number of flies remaining in the vial on the last day (value in the last column of the row) to the variable 'z'
  z<-temp[1,j+1]
    #if there are any flies alive in the vial on the last day
    if(z>0){
      #iterate over 1:(number of flies alive on the last day)
      for(l in 1:z){
        #record the current vial #
        vial<-temp[,1]
        #record the genotype of the current vial
        treatment<-temp[,2]
        #record the last day we recorded this fly alive (assumes that your input DF starts with day 0 in column 3)
        dd<-j-2
        #append this information into a new row in the 'results' data.frame, and add a '0' in the 4th column to indicate that the fly made it to the end of the experiment
        results[nrow(results)+1,]<- c(vial,treatment,dd,0)
      } #close for loop
    } #close if loop
  } #close original for loop
results$dead<-as.numeric(results$dead)  #reiterate that this column must be class numeric
results$status<-as.numeric(results$status)  #reiterate that this column must be class numeric
results$vial <- as.factor(results$vial) # make sure vial is considered a factor
# gives you only the results dataframe as output from function 
return(results) 
} #close function
```

Read in raw data

``` r
#read the file from csv
df<-read.csv("/Users/maggieschedl/Desktop/Github/Unckless_Lab_Resources/Infection_survival_analyses/20231025-rep-3-male-female-nanoject-16Cq/20231025-sheet.csv")

# separate out columns needed
df<-df[,c(1,3,14:32)]
```

Convert dataframe

``` r
df.convert<-convert_df(df)
```

Plot survivial curve

``` r
# change to not have confidence intervals in this one so you can see them 
df_fit<- survfit(Surv(dead, status) ~ treatment, data=df.convert)
ggsurvplot(df_fit,
          pval = TRUE, conf.int = FALSE,
          #risk.table = TRUE, # Add risk table
          #risk.table.col = "strata", # Change risk table color by groups
          #linetype = "strata", # Change line type by groups
          #surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw(), # Change ggplot2 theme
          palette = c("orchid", "aquamarine", "blueviolet", "darkslategray3")) + ylab("Survival Proporation") + xlab("Days post injection")
```

![](rep3-and-combo-analysis_files/figure-commonmark/unnamed-chunk-5-1.png)

Combined replicates 1, 2, and 3

``` r
#read the file from csv
df2<-read.csv("/Users/maggieschedl/Desktop/Github/Unckless_Lab_Resources/Infection_survival_analyses/20231025-rep-3-male-female-nanoject-16Cq/male-female-rep-1-2-3-counts-combo.csv")

# convert dataframe

df2.convert <- convert_df(df2)

# add in block information 
# first 69 rows are block A, second 69 rows are block B, and last 68 rows are block C
df2.convert$Block <- rep(c("A","B", "C"), c(69, 68,68))

# add in sex information and DiNV information as separate columns by splitting the treatment column 
# split the columns 
df2.convert_S <- str_split_fixed(df2.convert$treatment, "-", 2)
# change column names
colnames(df2.convert_S) <- c("sex", "injection")

# add columns to df 
df2.convert_full <- cbind(df2.convert,df2.convert_S )
```

Plot all replicates as one

``` r
df2.convert_full <- df2.convert_full %>% 
  mutate(treatment = factor(treatment, levels = c("male-CCM", "female-CCM", "male-DiNV", "female-DiNV")))

df_fit_combo_1 <- survfit(Surv(dead, status) ~ treatment, data=df2.convert_full)
ggsurvplot(df_fit_combo_1, size = 5,
          pval = FALSE, conf.int = FALSE,
          legend = "bottom",
          font.tickslab = c(14),
          font.x = c(16),
          font.y = c(16),
          font.t = c(16),
          ggtheme = theme_light(),
          title = expression(paste("Comparing Male and Female",italic(" D. innubila "), " injected with 16Cq DiNV")),
          legend.title="Treatment",
          legend.labs=c("male CCM", "female CCM", "male 16Cq DiNV", "female 16Cq DiNV"),
          font.legend = c(14),
          palette = c("#ccf9ff", "#62CFF4" ,"#2C67F2",  "#191970")) + ylab("Survival Proporation") + xlab("Days post injection")
```

![](rep3-and-combo-analysis_files/figure-commonmark/unnamed-chunk-7-1.png)

Find median survival time by treatment

``` r
surv_median(df_fit_combo_1, combine = FALSE)
```

    Warning: `select_()` was deprecated in dplyr 0.7.0.
    ℹ Please use `select()` instead.
    ℹ The deprecated feature was likely used in the survminer package.
      Please report the issue at <https://github.com/kassambara/survminer/issues>.

                     strata median lower upper
    1    treatment=male-CCM     NA    NA    NA
    2  treatment=female-CCM     NA    NA    NA
    3   treatment=male-DiNV      9     9     9
    4 treatment=female-DiNV      8     8     9

Start looking at models

Model just looking at significance of block and treatment

``` r
# model including block 
df_fit_combo_2<- coxph(Surv(dead, status) ~ treatment + Block, data=df2.convert_full)
summary(df_fit_combo_2)
```

    Call:
    coxph(formula = Surv(dead, status) ~ treatment + Block, data = df2.convert_full)

      n= 205, number of events= 150 

                             coef exp(coef) se(coef)      z Pr(>|z|)    
    treatmentfemale-CCM    1.9585    7.0888   1.0801  1.813   0.0698 .  
    treatmentmale-DiNV     4.5685   96.3993   1.0162  4.496 6.94e-06 ***
    treatmentfemale-DiNV   4.8149  123.3306   1.0123  4.756 1.97e-06 ***
    BlockB                -0.4061    0.6662   0.2073 -1.959   0.0501 .  
    BlockC                 0.3178    1.3741   0.2068  1.537   0.1244    
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

                         exp(coef) exp(-coef) lower .95 upper .95
    treatmentfemale-CCM     7.0888   0.141069    0.8534    58.882
    treatmentmale-DiNV     96.3993   0.010374   13.1543   706.448
    treatmentfemale-DiNV  123.3306   0.008108   16.9590   896.894
    BlockB                  0.6662   1.500951    0.4438     1.000
    BlockC                  1.3741   0.727771    0.9162     2.061

    Concordance= 0.774  (se = 0.025 )
    Likelihood ratio test= 170  on 5 df,   p=<2e-16
    Wald test            = 68.33  on 5 df,   p=2e-13
    Score (logrank) test = 131.3  on 5 df,   p=<2e-16

Model just looking at significance of block and treatment without sex

``` r
# model including block 
df_fit_combo_3<- coxph(Surv(dead, status) ~ injection + Block, data=df2.convert_full)
summary(df_fit_combo_3)
```

    Call:
    coxph(formula = Surv(dead, status) ~ injection + Block, data = df2.convert_full)

      n= 205, number of events= 150 

                     coef exp(coef) se(coef)      z Pr(>|z|)    
    injectionDiNV  3.3847   29.5102   0.4051  8.354   <2e-16 ***
    BlockB        -0.4446    0.6411   0.2061 -2.157    0.031 *  
    BlockC         0.2448    1.2773   0.1996  1.226    0.220    
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

                  exp(coef) exp(-coef) lower .95 upper .95
    injectionDiNV   29.5102    0.03389   13.3387   65.2874
    BlockB           0.6411    1.55984    0.4280    0.9603
    BlockC           1.2773    0.78288    0.8638    1.8889

    Concordance= 0.757  (se = 0.026 )
    Likelihood ratio test= 163.2  on 3 df,   p=<2e-16
    Wald test            = 75.94  on 3 df,   p=2e-16
    Score (logrank) test = 127.6  on 3 df,   p=<2e-16

Model looking at significance of block and treatment with sex

``` r
# model including block 
df_fit_combo_4<- coxph(Surv(dead, status) ~ injection + Block + sex, data=df2.convert_full)
summary(df_fit_combo_4)
```

    Call:
    coxph(formula = Surv(dead, status) ~ injection + Block + sex, 
        data = df2.convert_full)

      n= 205, number of events= 150 

                     coef exp(coef) se(coef)      z Pr(>|z|)    
    injectionDiNV  3.3495   28.4885   0.4043  8.286   <2e-16 ***
    BlockB        -0.3949    0.6737   0.2070 -1.908   0.0564 .  
    BlockC         0.3410    1.4063   0.2063  1.653   0.0984 .  
    sexmale       -0.3177    0.7278   0.1739 -1.827   0.0677 .  
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

                  exp(coef) exp(-coef) lower .95 upper .95
    injectionDiNV   28.4885     0.0351   12.8993    62.918
    BlockB           0.6737     1.4842    0.4491     1.011
    BlockC           1.4063     0.7111    0.9385     2.107
    sexmale          0.7278     1.3740    0.5176     1.023

    Concordance= 0.77  (se = 0.026 )
    Likelihood ratio test= 166.6  on 4 df,   p=<2e-16
    Wald test            = 79.42  on 4 df,   p=2e-16
    Score (logrank) test = 130.7  on 4 df,   p=<2e-16

Model looking at significance of block and injection by sex interaction

``` r
# model including block 
df_fit_combo_5<- coxph(Surv(dead, status) ~ Block + sex*injection, data=df2.convert_full)
summary(df_fit_combo_5)
```

    Call:
    coxph(formula = Surv(dead, status) ~ Block + sex * injection, 
        data = df2.convert_full)

      n= 205, number of events= 150 

                             coef exp(coef) se(coef)      z Pr(>|z|)    
    BlockB                -0.4061    0.6662   0.2073 -1.959   0.0501 .  
    BlockC                 0.3178    1.3741   0.2068  1.537   0.1244    
    sexmale               -1.9585    0.1411   1.0801 -1.813   0.0698 .  
    injectionDiNV          2.8564   17.3981   0.4380  6.522 6.94e-11 ***
    sexmale:injectionDiNV  1.7121    5.5408   1.0948  1.564   0.1178    
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

                          exp(coef) exp(-coef) lower .95 upper .95
    BlockB                   0.6662    1.50095   0.44379     1.000
    BlockC                   1.3741    0.72777   0.91620     2.061
    sexmale                  0.1411    7.08875   0.01698     1.172
    injectionDiNV           17.3981    0.05748   7.37405    41.048
    sexmale:injectionDiNV    5.5408    0.18048   0.64816    47.366

    Concordance= 0.774  (se = 0.025 )
    Likelihood ratio test= 170  on 5 df,   p=<2e-16
    Wald test            = 68.33  on 5 df,   p=2e-13
    Score (logrank) test = 131.3  on 5 df,   p=<2e-16

Compare models, which is best?

``` r
models <- list(df_fit_combo_3, df_fit_combo_4, df_fit_combo_5)

model.names <- c( 'block and injection', 'block injection and sex', 'block and sex injection interaction')

aictab(cand.set = models, modnames = model.names)
```


    Model selection based on AICc:

                                        K    AICc Delta_AICc AICcWt Cum.Wt      LL
    block and sex injection interaction 5 1283.22       0.00   0.56   0.56 -636.46
    block injection and sex             4 1284.53       1.30   0.29   0.85 -638.16
    block and injection                 3 1285.84       2.62   0.15   1.00 -639.86

``` r
# best model is df_fit_combo_5 but not by much 
summary(df_fit_combo_5)
```

    Call:
    coxph(formula = Surv(dead, status) ~ Block + sex * injection, 
        data = df2.convert_full)

      n= 205, number of events= 150 

                             coef exp(coef) se(coef)      z Pr(>|z|)    
    BlockB                -0.4061    0.6662   0.2073 -1.959   0.0501 .  
    BlockC                 0.3178    1.3741   0.2068  1.537   0.1244    
    sexmale               -1.9585    0.1411   1.0801 -1.813   0.0698 .  
    injectionDiNV          2.8564   17.3981   0.4380  6.522 6.94e-11 ***
    sexmale:injectionDiNV  1.7121    5.5408   1.0948  1.564   0.1178    
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

                          exp(coef) exp(-coef) lower .95 upper .95
    BlockB                   0.6662    1.50095   0.44379     1.000
    BlockC                   1.3741    0.72777   0.91620     2.061
    sexmale                  0.1411    7.08875   0.01698     1.172
    injectionDiNV           17.3981    0.05748   7.37405    41.048
    sexmale:injectionDiNV    5.5408    0.18048   0.64816    47.366

    Concordance= 0.774  (se = 0.025 )
    Likelihood ratio test= 170  on 5 df,   p=<2e-16
    Wald test            = 68.33  on 5 df,   p=2e-13
    Score (logrank) test = 131.3  on 5 df,   p=<2e-16
