R-Markdown-and-Survival-Curve-Template
================
2022-11-14

R markdown is a type of R file that is like a text markdown file, except
that you can include R code in it, and it will “knit” and include any
outputs or plots in the document.

Text written outside of a “chunk” of code will be just text, and follows
traditional markdown syntax. Text inside of a chunk needs to be
commented so it is not seen as code.

1.  To make a markdown file, use the dropdown menu for making a new
    file, and select R Markdown. Give your file a descriptive title.
2.  This will open up a file with a header, and some basic information
    about the file format. This is the most important piece: “For more
    details on using R Markdown see <http://rmarkdown.rstudio.com>”
3.  However you don’t want much of that text in your actual document, so
    delete everything in the file except for the title, output, date
    section at the top
4.  To make this file easily readable on Github, change “html_document”
    to “github_document” in the output section at the top
5.  The first thing to do in an R script is to load in packages needed
    for the script. To do this you first have to create a chunk. You can
    press the C+ button on the top right section of the R Markdown
    document and select R

``` r
# this has made a chunk with the {r} specification 
# this means the code inside this chunk will be read as R code 
# there are other code options, but those aren't needed right now 
# then you can put in the code you normally use for loading in a package in this space 
# notice how all the text in this chunk has to be commented 
# if these packages are not installed, install them first 
# install.packages("survival")
# install.packages("survminer")

# load in survival and surminer packages 
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

6.  You can now run the commands inside the chunk by pressing the green
    arrow button on the right side of the chunk
7.  This is basically all you need to know about R Markdown format: you
    can write descriptions outside of the code chunks, code inside the
    chunks, and you can break up your code into sections to be run
    separately by using multiple chunks. The rest of this document will
    be about how to take data collected from an infection experiment and
    make a survival curve plot

## The function for converting the dataframe

The format of the data that we collect in our infection studies is not
the correct format for making a survival curve/Kaplan Meier plot. For
the survival curve, the dataframe needs to have a row for every
individual fly, and for each fly what treatment it got, what day it died
(or the last day of the experiment), and the “status.” Because this
format treats a fly that survived to the end of the experiment (ex. 7
day experiment) as dead on day 7, then there is a “status” variable.
That is “1” if the animal actually died, and 0 if it is still alive on
the “day dead.”

This is not at all like the format of the data we have now, but we can
use this function below to go through our data and transform it into the
data we need. Below is code to make a function to do the data format
conversion. You can go line by line through the comments, but it is
pretty complicated. In short, the code takes each row (vial) and goes
through each day column to figure out how many flies died on each day.
If a fly died on a day, it’ll put that information into a new dataframe,
along with the status. The code will go through each day for that vial,
until it gets to the last day and all the alive flies left will be put
into the new dataframe as dead on the last day but the status is 0. The
the code repeats for the next row (vial) and so on.

That is a lot of code, and it would be nice not to write it out every
time you want to do the conversion in the script. So it has been made
into a function, where all you have to do is put
convert_df(dataframe-to-be-converted) and it will do all of the code
encompassed in that one function.

### Convert df function

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
            #record the treatment of the current vial
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
        #record the treatment of the current vial
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

Nothing will be printed after doing this, but you can see in your
environment after running this that you have convert_df as a function.

Now the next step is to actually load in your data. You should download
your Google Sheets file as a csv and make sure the name of the file does
not have any spaces in it, use - or \_. You should also put the file in
a folder in your github repository for survival curve analysis, so that
the data and the R files will be together in the same place, and can be
pushed to Github.

Additionally we need to subset the dataframe to be the correct columns
needed for the analysis only. Columns like replicate, species, and fly
age are not needed, and will actually mess up the function above. The
function only works on a dataset that has the vial in the first column,
the treatment in the second column, and the day 0. day 1, day 2 etc
after that. See below as an example of the data needed:

| unique_vial_name | treatment         | starting_n_number | day1 | day2 | day3 | day4 | day5 | day6 |
|------------------|-------------------|-------------------|------|------|------|------|------|------|
| 1                | D.mel_LB          | 10                | 10   | 10   | 9    | 9    | 9    | 9    |
| 2                | D.mel_LB          | 9                 | 9    | 9    | 10   | 10   | 10   | 10   |
| 3                | D.simulans_LB     | 10                | 10   | 9    | 9    | 9    | 9    | 9    |
| 4                | D.simulans_LB     | 10                | 10   | 10   | 10   | 10   | 10   | 10   |
| 5                | D.mel_P.rett      | 9                 | 10   | 10   | 9    | 8    | 7    | 6    |
| 6                | D.mel_P.rett      | 10                | 10   | 8    | 8    | 7    | 7    | 7    |
| 7                | D.simulans_P.rett | 10                | 10   | 7    | 6    | 6    | 4    | 3    |

If you have extra columns in your data you can remove them in R, but you
have to at least have these types of columns for the function to work.

### Read in the raw data and make subsets

``` r
# Note that all code from now on will be commented to be able to knit with R Markdown. If your code doesn't work, the program will not knit it. Because this code is examples, it doesn't actually run. Delete the three hashes ### before the actual line of code and modify it with your own code

# read the file from csv
###
### file_name  <- read.csv("full/path/to/your/file.csv")
###

# Now you want to check what the dataframe looks like so you can know what columns you need to keep 
# Just look at the first 10 rows with head()
###
### head(file_name)
###

# Subset dataframe to be only the columns needed, the numbers of the columns can be different for each dataset depending on what information you record 
# to subset the columns, you do file_name[,c(column numbers)]
# inside the [ ] before the comma means rows, and after the comma means columns
# you just use the number of the column to tell R what to keep 
# if you have a set of columns that are next to each other, like columns 9 through 12, you can represent that as 9:12
###
### file_name_subset <- file_name[,c(#,#,#:#)]
###

# Again you want to check that the code worked, and that your subset is the columns you actually wanted 
###
### head(file_name_subset)
###
```

Now we have to covert the subsetted dataframe from above into the format
for the survival curve, but because we made that function this is very
easy!

### Convert dataframe to correct format for plotting

``` r
# convert the subsetted dataframe using the function above 
###
### file_name_subset_convert <- convert_df(file_name_subset)
###

# check how the converted dataframe looks 
###
### head(file_name_subset_convert)
###
```

Next is to apply the functions from the survival and survminer packages
that we loaded in at the beginning. They take that modified dataframe
and make a survival model, and then you are able to plot the model.
These packages also so statistical analysis to see if the differences in
survival are statistically different between the treatment groups. These
models are rather complex because they take into account the whole
shapes of the survival curves between the treatments. For further
reading on the actual nuts and bolts of these plots and statistics, see
<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3932959/> and
<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059453/>

### Make the survivorship model and plot the model

``` r
# first you need to make the model using the survfit function. All you have to do is add in the name of the dataframe, the names of the 
# columns of the dataframe should always be the same because of the convert function we used 
###
### file_name_subset_convert_fit <- survfit(Surv(dead, status) ~ treatment, data=file_name_subset_convert)
###

# then we use a ggplot specific for the survfit model to plot it. Again the only thing you need to add in is the name of your model
###
### ggsurvplot(file_name_subset_convert_fit,
          ### pval = FALSE, conf.int = TRUE, # you can change these to remove the confidence intervals from the lines  
          ### linetype = "strata", # Change line type by groups
          ### ggtheme = theme_bw()) + # Change ggplot2 theme
          ### ylab("Survival Proportion") + xlab("Days post infection") # add in labels for your axis 
          ### palette = c() # you can add colors if you would like for each line
###
```

If you want to know more about the model statsitics you should run a
specific cox proportional hazard model on your data. You can look at
individual p values and determine other ways to run the model if needed.

``` r
### df_fit_t<- coxph(Surv(dead, status) ~ treatment, data=df.convert)
### summary(df_fit_t)
```

Once all of your code is running and it looks good, you can knit your R
Markdown by pressing the knit icon on the top left of the code file
window, and R will turn this file into a Markdown with with all your
tables and plots printed as images. This can be pushed to Github to be
viewed in your Notebook repository if you have one.
