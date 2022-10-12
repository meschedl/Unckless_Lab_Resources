2022-09-24-food-vial-test
================
2022-10-12

Data from this experiment is 5-7 day old D. innubila poked with sterile
cell culture medium and placed on different food vials. Vials were
either molasses food, mushroom instant food with cotton roll, sugar
agar, or sugar agar with 1/4 mushroom broth.

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
#read the file from csv
food_data<-read.csv("~/Desktop/Github/Unckless_Lab_Resources/Infection_survival_analyses/20220924/20220924_food_vial_test_pokes.csv")

#make subsets of the data frame based on males and females
food_data.male<-food_data[food_data$sex == "M",]
food_data.female<-food_data[food_data$sex == "F",]

# remove unpoked flies
food_data.male <- food_data.male[food_data.male$poked == "yes",]
food_data.female <- food_data.female[food_data.female$poked == "yes",]

#remove extraneous columns (3&4) from each to get them in proper format
food_data.male<-food_data.male[,c(1,5,14:28)]
food_data.female<-food_data.female[,c(1,5,14:28)]
```

### Convert each of these dataframes to long and tidy format using function defined above

``` r
food_data_male_convert<-convert_df(food_data.male)
food_data_female_convert<-convert_df(food_data.female)
```

### make survivorship curves

### Male

``` r
male_fit<- survfit(Surv(dead, status) ~ treatment, data=food_data_male_convert)
ggsurvplot(male_fit,
          pval = TRUE, conf.int = TRUE,
          #risk.table = TRUE, # Add risk table
          #risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw()) # Change ggplot2 theme
```

![](20220924-food-vial-innubila-test_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
\### make survivorship curves \### Female

``` r
female_fit<- survfit(Surv(dead, status) ~ treatment, data=food_data_female_convert)
ggsurvplot(female_fit,
          pval = TRUE, conf.int = TRUE,
          #risk.table = TRUE, # Add risk table
          #risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw()) # Change ggplot2 theme
```

![](20220924-food-vial-innubila-test_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
