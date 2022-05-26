# Learning how to calculate the Standard Error of the Mean and the "Quick" 95% Confidence Interval in R 


# Read in the dataframes 
egg_counts_main <- read.csv("innubila_egg_counts_main_sheet.csv") # read in main datasheet
egg_counts_individual <- read.csv("innubila_egg_counts_individual_plate_counts.csv") # read in plate count datasheet

# Calculate standard error of the mean of average_egg_count by fly age 

# separate out counts into ages 2 and 3
# use the subset function to take only the rows in the main sheet that are flies age 2 
# Then turn the average egg count column into a list of counts 
# the == double equal sign means that you are specifying the character in the fly_age column has to be two
# and the two is in quotes because it is a character
average_egg_count_2 <- subset(egg_counts_main, fly_age == "two") 
average_egg_count_2 <- average_egg_count_2$average_egg_count
# then make your own function that is the equation for standard error 
# this is the standard deviation of your data, divided by the square root of the number of observations
std <- function(x) sd(x)/sqrt(length(x)) # note that when you run this line of code once in this script you won't need to again
# then use this function on the list you just created
std(average_egg_count_2)
# the output of this is 5.7
# this means that the precision around the mean of average egg counts for age 2
# (56.3 eggs, calculated in the previous script) is + or - 5.7 eggs 
# or 56.3 ± 5.7 eggs 

# this can then be done for flies ages 3 
average_egg_count_3 <- subset(egg_counts_main, fly_age == "three")
average_egg_count_3 <- average_egg_count_3$average_egg_count

std(average_egg_count_3)
# the precision around the mean of average egg counts for age 3 is 
# 79.3 ± 6.7 eggs 

#### 
# Now on your own calculate the standard error of the mean of average egg counts for other categories:
# for yeast or no yeast (either "Y" or "N" in yeast column)
# for each plate type (either "MAS", "MNS", "AJS", or "MS" in the plate_type column)
# remember to name lists differently so you don't overwrite things you've already made 


# Calculating the "quick" version of 95% confidence interval of average_egg_counts by fly age 
# we can follow the quick and dirty rule to get an approximate 95% confidence interval of the mean 
# this is done by adding and subtracting 2 standard error of the means from the mean of your value of interest 

# flies age 2 
# create values for twice the SE of egg counts for age 2
# and the average egg counts for age 2
age_2_interval <- 5.7 * 2
mean_age_2 <- 56.3

# then you can do the math with the values
mean_age_2 + age_2_interval
mean_age_2 - age_2_interval
# The approximate 95% confidence interval for egg counts for flies age 2 is 44.9 - 67.7 eggs 

# flies age 3
age_3_interval <- 6.7 * 2
mean_age_3 <- 79.3 

mean_age_3 + age_3_interval
mean_age_3 - age_3_interval
# the approximate 95% confidence interval for egg counts for flies age 3 is 65.9 - 92.7 eggs 

### 
# Now on your own calculate the approximate 95% confidence interval of average egg counts for other categories:
# for yeast or no yeast
# for each plate type 


