## Data table R script

# Setting up and reading in data

getwd() # check working directory 
# if you are in the directory (folder) where your csv files are, then they can easily be read in by R 
# without needing to put in the whole path to them 

# read in data tables
# they will be "dataframes" in R language 
egg_counts_main <- read.csv("innubila_egg_counts_main_sheet.csv") # read in main datasheet
egg_counts_individual <- read.csv("innubila_egg_counts_individual_plate_counts.csv") # read in plate count datasheet 


# View dataframe in multiple ways 

head(egg_counts_individual) # look at the first 10 rows
tail(egg_counts_main) # look at the last 10 rows
print(egg_counts_main) # print the entire dataframe


### Summing Rows
# make column in the egg_counts_individual dataframe that is sum of counts a b c and d in individual counts sheet
# this should be the same column that you already have (total) that excel generated for you
# use rowSums to get the sum of the rows, and only for the columns you are interested in summing 
# call the function, then specify the dataframe you want to use, then the columns wanted 
# in [ , 5 : 8] the first space is left blank because that pertains to rows
# then you have 5 : 8 because you want to sum columns 5-8 
egg_count_R  <-  rowSums(egg_counts_individual[ , 5 : 8])
# this creates a vector that is a list of all the sums of those columns by row
# now it needs to be combined with your dataframe 
# use the function cbind (to bind columns) and tell it to use your dataframe and the list you just made 
egg_counts_individual_R  <-  cbind(egg_counts_individual, egg_count_R)

# check to see if it worked and if the total column and the egg_counts_R columns are the same
head(egg_counts_individual_R)

###########
# Now on your own make 2 new columns in egg_counts_individual_R where one is the sum of A_count and B_count
# and the over is the sum of C_count and D_count 
# use cbind to add them to egg_counts_individual_R
# then use head() or tail() to check your work 


### Averaging Rows 
# Make a column in the egg_counts_main main dataframe that is the average of your replicate counts 
# this time you can do the mean and cbind step in one line of code
# the $ symbol is a signifier for column name
# so egg_counts_main$egg_count_mean indicates to create a column in egg_counts_main called egg_count_mean
# and the <- assigns the rowMeans function output to that column 
egg_counts_main$egg_count_mean <- rowMeans(egg_counts_main[ ,10 : 12])

# check to see if it worked
head(egg_count_main)

##########
# Now on your own create a new column in egg_counts_main that is the average of egg_count_1 and egg_count_3 
# add them into the egg_counts_main dataframe in one step 
# then use head() or tail() to check your work 


### Averaging rows by column categories
# the whole point of this process is to understand how many eggs were laid based on the difference categories of plates tested
# R can separate out means by the contents of columns!
# first you need to install and load a package that has the functions that do this
# R Studio has a lot of functions built into it (like all the ones you already used), but you don't download it with every
# possible function, plus people are making new functions and packages all the time 
# so you have to install a package to your computer from the internet, and then also load it each time you open your script
# we will be using the package dplyr
intall.packages("dplyr")
# this will ask you "Do you want to install from sources the packages which need compilation? (Yes/no/cancel)" 
# you type in yes to your console and press enter
# your console will stream out a lot of text that means the package is installing, don't worry about it being in read 
# once you have installed the package to your computer, you will not need to install it again unless there is a new 
# version that comes out that you want to re-install
# when the installation is finished, you will need to load the package into this R session 
# this is because you don't always use every R package each time you use R
library(dplyr)
# now you will be able to use all the functions associated with the package dplyr
# if you want to know more about the package, type help("dplyr") into your console and the help message will come up on
# your lower right window

# for some reason when loading in the dataset, R doesn't see the average_egg_count as numeric 
# this is likely because we used an equation to make those numbers in google sheets, so it got marked as a character
# this makes it impossible to calculate the mean because it doesn't think it's a number!
# but you can tell R to reassign the things in the average_egg_count column to be classed as numeric 
# this line of code will take the average_egg_count column and code is as numeric, and assign it back to the same column
egg_counts_main$average_egg_count <- as.numeric(as.character(egg_counts_main$average_egg_count))


# now we are going to take the egg_counts_main dataframe and look at the average of the average_egg_count based off
# the category of fly_age: either two or three 
# this code uses slightly difference syntax that is specific to the dplyr package
# you take the egg_counts_main dataframe then use %>% to "pipe" it into other functions 
# group_by tells R to group the data by the categories in the fly_age column 
# then summarise is a function that you can use many calculations for, but for this purpose we want to know the mean
# and we tell R that we want to know the mean of the average_egg_count column
# the results from this function will print out on your console
egg_counts_main %>% 
  group_by(fly_age) %>% 
  summarise(average_egg_count = mean(average_egg_count))

#######
# Now on your own use the same code to find the mean of average_egg_count by the categories of yeast and plate_type
average_egg_count_by_age <- egg_counts_main %>% 
  group_by(fly_age) %>% 
  summarise(average_egg_count = mean(average_egg_count)) %>% 
  as.data.frame()

# this is great, but it would be even better to have R give us this information in a dataframe instead of just 
# printing it to the console 
# you can easily add in a function to the code above to make it output the information in a dataframe

average_egg_count_by_age <- egg_counts_main %>% 
  group_by(fly_age) %>% 
  summarise(average_egg_count = mean(average_egg_count)) %>% 
  as.data.frame()

#######
# Do this for the categories of yeast and plate_type 
# make sure to change the name of the dataframe you are assigning the data to each time 
# then print to the console the dataframes for egg count by fly age, yeast, and plate type 
