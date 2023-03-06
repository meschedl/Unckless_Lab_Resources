## Making Box Plots

# Setting up and reading in data

# load packages needed
library(ggplot2)

# read in dataset using full path to file

egg_counts_main <- read.csv("C:Your\path\to\file\innubila_egg_counts_main_sheet.csv")
# look at it 
head(egg_counts_main)


# Make a simple box plot of the average egg count separated by the two ages of the groups of flies, either "2" or "3"
# x axis is fly_age
# y axis is the average_egg_count
# geom_boxplot() is what tells R what graph to make
# it will show up on the right lower corner in the plots tab
age <- ggplot(egg_counts_main, aes(x=fly_age, y=average_egg_count)) + 
  geom_boxplot()
# plot the figure
age


# On your own make a simple box plot of the average egg count separated out by whether the plate had yeast or not
# x axis is yeast
# y axis is the average_egg_count

# We can also separate out the dataset to plot different aspects of the data separatly 
# Separating plots out by age

# First, subset the dataframe by the ages 
# You tell R to only consider rows from 1:24 and the comma says keep all columns 
egg_counts_two <- egg_counts_main[1:24,]
# You tell R to only consider rows from 25:48 and the comma says keep all columns 
egg_counts_three <- egg_counts_main[25:48,]

# Then use these separated dataframes to make the same plots as above 
# Effect of yeast on age 2 flies 
yeast_2 <- ggplot(egg_counts_two, aes(x=yeast, y=average_egg_count)) + 
  geom_boxplot()
# plot the figure
yeast_2 

# On your own make a simple box plot of the effect of yeast or not on just the age 3 flies 

# We can also plot the differences in egg count by plate type
# Make a box plot where plate type is the x axis
# Average egg count is the y axis 
# And these can be made into separate plots by age because we already subsetted the dataframe by age 

# Average egg count by plate type for age 2 flies 
plate_2 <- ggplot(egg_counts_two, aes(x=plate_type_ab, y=average_egg_count)) + 
  geom_boxplot()
# plot the figure
plate_2

# On your own make a plot for average egg count by plate type for age 3 flies 

# What if we wanted to look at both yeast and plate type? 
# This can be done with a grouped box plot 
# For this, we will still keep the ages separate 
# But we can use the aesthetic "fill" to specify another variable to separate out the plot into
# Fill will color the boxes based on an additional variable = yeast status

# Grouped box plot for age 2 flies by yeast status and plate type
age_2_grouped <- ggplot(egg_counts_two, aes(x=plate_type_ab, y=average_egg_count, fill = yeast)) + 
  geom_boxplot()

age_2_grouped 

# On your own make a grouped box plot for average egg count by plate type will yeast as a fill variable fo age 3 flies 


# Here are other combinations of plots you should try making, for fun and to gain R skills 
# Grouped plot of plate type and age, without separating by yeast 
# Grouped plot for plate type and yeast status without separating by age 
# Can you separate out specific plate types? 
# hint: use the notation below subset your dataframe by a specific value or character within the column/row
# egg_counts_MNS <- egg_counts_two[egg_counts_two$plate_type_ab == "MNS",]
# where egg_counts_two$plate_type_ab tells R to specifically look at the column plate_type_ab
# and = "MNS", tells R to only look at where the plate type is MNS
# OR 
# egg_counts_NOT_MNS <- egg_counts_two[egg_counts_two$plate_type_ab != "MNS",]
# where the != excludes any row where the plate_type_ab is MNS


# Finally, you can do aesthetic changes to your plots 
# Making the plots prettier and clearer 


# Redo the age graph with clearer labels 
# Order the boxes by creating "levels" so that the boxes appear in order on the graph 
egg_counts_main$fly_age <- factor(egg_counts_main$fly_age , levels=c("two", "three"))

# Use ylab for the y axis label and xlab for the x axis label 
# geom_jitter() adds the individual data points to your graph 
# Theme_minimal() makes the background of the plot minimal 
# And the theme(legend.position = "none") removes the legend from the graph 
# scale_fill_manual will let you select which colors to use for making the graph
# the colors are in "" and are named so they can be specificially recognized by R 

age <- ggplot(egg_counts_main, aes(x=fly_age, y=average_egg_count, fill = fly_age)) + 
  geom_boxplot() + ylab("Average egg counts") + xlab("Age of flies") + theme_minimal() + geom_jitter() + 
  theme(legend.position = "none") + scale_fill_manual(values=c("paleturquoise1", "plum2"))
# plot the figure
age

# Redo the plate type and yeast graph with clearer labels for just age 3 
# Create a value that will be the title of the legend
legend_title <- "Yeast condition"

age_3_grouped <- ggplot(egg_counts_three, aes(x=plate_type_ab, y=average_egg_count, fill = yeast)) + 
  geom_boxplot() + ylab("Average egg counts") + xlab("Plate Type") + theme_minimal() + geom_jitter() + 
  scale_fill_manual(legend_title,values=c("sienna","burlywood1"))

age_3_grouped 


# Redo the plate type and yeast graph with clearer labels for both ages 

All_grouped <- ggplot(egg_counts_main, aes(x=plate_type_ab, y=average_egg_count, fill = yeast)) + 
  geom_boxplot() + ylab("Average egg counts") + xlab("Plate Type") + theme_minimal() + geom_jitter() + 
  scale_fill_manual(legend_title,values=c("burlywood1","sienna"))

All_grouped 


# Make any plot with the aesthetic modifications you want
# you can use this link http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf to find other colors to use 

