20221022-Food-vial-innubila
================
2022-11-09

``` r
library("survival")
library("survminer")
```

    ## Loading required package: ggplot2

    ## Loading required package: ggpubr

    ## 
    ## Attaching package: 'survminer'

    ## The following object is masked from 'package:survival':
    ## 
    ##     myeloma

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

### Read in the real raw data and make subsets

``` r
# read the file from csv
innubila_food_2<-read.csv("~/Desktop/Github/Unckless_Lab_Resources/Infection_survival_analyses/20221022/20221022_innubila_food_vial_test_2.csv")

# Subset dataframe to be only the columns needed
innubila_food_2_sub<-innubila_food_2[,c(1,4,12:26)]
```

### Convert the dataframe to long and tidy format using function defined above

``` r
innubila_food_2_sub_convert<-convert_df(innubila_food_2_sub)
```

### Make full plot

``` r
innubila_food_2_sub_convert_fit<- survfit(Surv(dead, status) ~ treatment, data=innubila_food_2_sub_convert)
ggsurvplot(innubila_food_2_sub_convert_fit,
          pval = TRUE, conf.int = TRUE,
          #risk.table = TRUE, # Add risk table
          #risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          # surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw()) # Change ggplot2 theme
```

![](2022-10-22-innubila-food-test-2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

### Separate out based on poke treatment

``` r
# Subset dataframe to be only the rows for not poked flies
innubila_food_2_nopoke<-innubila_food_2_sub[c(1,4:7,11),]

# Subset dataframe to be only the rows for poked flies
innubila_food_2_CCM<-innubila_food_2_sub[c(2:3,8:10,12),]
```

### Convert the dataframe to long and tidy format using function defined above

``` r
# no poke
innubila_food_2_nopoke_convert<-convert_df(innubila_food_2_nopoke)

# CCM
innubila_food_2_CCM_convert<-convert_df(innubila_food_2_CCM)
```

### Make full plot for No Poke

``` r
innubila_food_2_nopoke_convert_fit<- survfit(Surv(dead, status) ~ treatment, data=innubila_food_2_nopoke_convert)
ggsurvplot(innubila_food_2_nopoke_convert_fit,
          pval = TRUE, conf.int = TRUE,
          #risk.table = TRUE, # Add risk table
          #risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          # surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw()) # Change ggplot2 theme
```

![](2022-10-22-innubila-food-test-2_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Make full plot for CCM

``` r
innubila_food_2_CCM_convert_fit<- survfit(Surv(dead, status) ~ treatment, data=innubila_food_2_CCM_convert)
ggsurvplot(innubila_food_2_CCM_convert_fit,
          pval = TRUE, conf.int = TRUE,
          #risk.table = TRUE, # Add risk table
          #risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          # surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw()) # Change ggplot2 theme
```

![](2022-10-22-innubila-food-test-2_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
