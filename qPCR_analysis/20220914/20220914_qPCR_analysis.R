
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("Desktop/Github/Unckless_Lab_Resources/qPCR_analysis/20220914/")
qPCR_data <- read.csv("20220914-sheet.csv")

# remove the control samples 

qPCR_data <- qPCR_data[which(qPCR_data$type != "control"),]


# use the variance function, and calculates the variance in Cq by the unique.name (each sample/primer has 3 Cq values to 
# calculate the variance 
qPCR_data$Cq_var <- ave(qPCR_data$Cq, qPCR_data$unique.name, FUN=var)

# use the mean function, and calculates the mean of Cq values by the unique.name (each sample/primer has 3 Cq values to 
# calculate the mean 
qPCR_data$Cq_mean <- ave(qPCR_data$Cq, qPCR_data$unique.name, FUN=mean)

# Keep all rows where the replicate is 1 (or you could do 2 or 3)
# make into new Df so we keep the original with all the Cq values
qPCR_data_1rep <- qPCR_data[which(qPCR_data$replicate == "1"),]

# histogram of all variances
ggplot(qPCR_data_1rep, aes(x=Cq_var)) + geom_histogram(bins = 50)

# there is one sample where one of the replicates is a lot more different than the others, just want to remove that replicate well

# reset row names 
rownames(qPCR_data) <- 1:nrow(qPCR_data)
# remove line 17, the one replicate that is different than the other 3 
qPCR_data_1 <- qPCR_data[c(1:16, 18:24),]

# recalculate variance and means 
# variance
qPCR_data_1$Cq_var <- ave(qPCR_data_1$Cq, qPCR_data_1$unique.name, FUN=var)
# mean
qPCR_data_1$Cq_mean <- ave(qPCR_data_1$Cq, qPCR_data_1$unique.name, FUN=mean)
# Keep all rows where the replicate is 1 (or you could do 2 or 3)
# make into new Df so we keep the original with all the Cq values
qPCR_data_rep1 <- qPCR_data_1[which(qPCR_data_1$replicate == "1"),]

# histogram of all variances
ggplot(qPCR_data_rep1, aes(x=Cq_var)) + geom_histogram(bins = 50)

# 0.5 is the highest, I will keep with that for now 

# now to generate the delta Cq 
# order df by sample?
# need to get the 115 and RPL11 from the same sample together 

qPCR_data_rep1 <- qPCR_data_rep1[order(qPCR_data_rep1$Sample),]
# this does RPL11 first and 115 second 

nrow(qPCR_data_rep1)

# Separate that dataframe, incriminating by 2, every number between 1-8 (number of rows in dataframe)
qPCR_data_rep1$Cq_mean[seq(1,8,2)] # these are the RPL11 Cq means 
qPCR_data_rep1$Cq_mean[seq(2,8,2)] # these are the 115 3 primer Cq means 

# make the delta Cq by subtracting the 115 values from the RPL11 primer values
# and this is saved as a vector in R 
delta_Cqs <- qPCR_data_rep1$Cq_mean[seq(1,8,2)] - qPCR_data_rep1$Cq_mean[seq(2,8,2)]
#vector
delta_Cqs

# Keep only rows that are 115
Cq_values1rep_Delta <- qPCR_data_rep1[which(qPCR_data_rep1$Target == "115"),]
# And then add in the delta Cqs as a new column
Cq_values1rep_Delta$delta_Cq <- delta_Cqs

# do 2^ delta Cq
Cq_values1rep_Delta$delta_Cq_2 <- 2^(delta_Cqs)





