# Making a frequency distribution graph

# In this script you will learn how to make a graph of the frequency distribution of your data 

# first you will need to install and then load in a new package
# see the data_table.R script for what that means if you have forgotten
# install the package ggplot2
install.packages("ggplot2")
# remember that you only need to install the package once to your computer/R Studio, unless you download a different version of R
# load in the package 
library(ggplot2)
# usually you load in all the packages you need for a script at the beginning of the document, before any other code 

# then make sure you have your two dataframes loaded in
egg_counts_main <- read.csv("innubila_egg_counts_main_sheet.csv") # read in main datasheet
egg_counts_individual <- read.csv("innubila_egg_counts_individual_plate_counts.csv") # read in plate count datasheet

# making a plot is a more complicated piece of code than simple functions, but has a base that is simple
# start out with a basic histogram
# the plot will show up on your right in the plots window 
# the first part of the code says you are using the ggplot package, and that you want to use the egg_counts_main dataframe
# the second argument inside the ggplot is "aes" or the aesthetics (I think of it as the components) of the function
# we are telling it to use the column egg_count_1 as the x-axis
# because this is a histogram, the y-axis is just frequency of observations at that value, so you do not tell it a column for the y-axis
# that is all you tell the ggplot function, so you close the parentheses 
# then you add in a function that tells it to make this into a histogram plot 
ggplot(egg_counts_main, aes(x=egg_count_1)) + geom_histogram()
# you will get on the right a sort of sad looking histogram, it tells you a little about the spread of the data
# but it also doesn't look very informative 
# you will also get a message in your console in red below that says:
# `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
# this is not exactly an error, because the code ran ok and made the plot, it is more like a suggestion
# it is suggesting that we change the size of the "bins" that it lumps the columns in the histogram into 
# you can tell in the plot that the x-axis doesn't have increments of 1 egg, and it doesn't have increments of 50 eggs either
# it's easier to see this by actually making those plots 
# let's make the binwidth 1 and you'll see what that means
# binwidth is a specific parameter for the geom_histogram function, so you put in the specifics for that within the geom_histogram parentheses
ggplot(egg_counts_main, aes(x=egg_count_1)) + geom_histogram(binwidth = 1)
# wow this looks even worse! Almost all your plates had slightly different numbers of eggs, only a few had the exact same number 
# so this isn't very informative at all
# what if you made the binwidth much higher, so that many egg counts would fall into that bin and be counted up into the frequency?
ggplot(egg_counts_main, aes(x=egg_count_1)) + geom_histogram(binwidth = 25)
# this visually looks a little better than a bin for every egg, but you are hiding some of the variation with bins so large here 
# what about something in-between
ggplot(egg_counts_main, aes(x=egg_count_1)) + geom_histogram(binwidth = 10)
# this looks nicer, you can see more variation in the counts, and you can see that there is a little bit of a tail to the right end
# where some plates had many eggs, but not very many 

# on your own, decide what number for binwidth looks best to you. There is not a right answer.


# that was the basics of the plot
# now we can get into real aesthetics 
# what if you want to change the color of the plot?
# confusingly, when you specify "color" that really means the outline of the bars 
# and fill is the color of in the inside of the bars
# notice that the color terms are in quotes 
ggplot(egg_counts_main, aes(x=egg_count_1)) + geom_histogram(binwidth = 10, color = "black", fill = "palegreen2")

# another thing you can do is change the background 
# this is added on additionally, because it is not a specific part of just the histogram function, you can change the background/theme in all ggplot plots
ggplot(egg_counts_main, aes(x=egg_count_1)) + geom_histogram(binwidth = 10, color = "black", fill = "palegreen2") + theme_bw()

# and finally you can change the labels of the axes to be more informative 
# ylab stands for y axis label, and 
ggplot(egg_counts_main, aes(x=egg_count_1)) + geom_histogram(binwidth = 10, color = "black", fill = "palegreen2") + 
  theme_bw() + ylab("Frequency of count") + xlab("Number of eggs counted per plate")

# you can also name your plot as an object 
egg_count_1_histogram <- ggplot(egg_counts_main, aes(x=egg_count_1)) + geom_histogram(binwidth = 10, color = "blue", fill = "palegreen2") + 
  theme_bw() + ylab("Frequency of count") + xlab("Number of eggs counted per plate")
# and then use the object name to plot it 

egg_count_1_histogram 

# and there are ways you can save your plot as an external picture 
# you can either go to the bottom right window and click on export and save image
# or you can do it with code! 
# the ggsave function only works on the last plot you make, so it's easiest to put the code one after the other
ggplot(egg_counts_main, aes(x=egg_count_1)) + geom_histogram(binwidth = 10, color = "blue", fill = "palegreen2") + 
  theme_bw() + ylab("Frequency of count") + xlab("Number of eggs counted per plate")
ggsave("egg_count_1_histogram.jpeg")

#######
# now make your own histogram
# use the average_egg_count column from the egg_counts_main dataframe (the data we actually care about)
# chose your own color scheme using this document: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# and chose your own theme using this link: https://ggplot2.tidyverse.org/reference/ggtheme.html (note that not all of these apply to your type of data)
# when you are done save the plot and upload to googledrive! 
 




























