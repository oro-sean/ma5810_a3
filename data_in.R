## set up environment and import csv file

rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")
setwd("~/Documents/GitHub/ma5810_a3")


library(caret, warn.conflicts = FALSE, quietly = TRUE) # handy ml package, data splitting, training ect ect
library(cluster, warn.conflicts = FALSE, quietly = TRUE) # for clustering
library(tidyverse, warn.conflicts = FALSE, quietly = TRUE) # handy for data prep
library(ggplot2, warn.conflicts = FALSE, quietly = TRUE) # plotting
library(ggpubr, warn.conflicts = FALSE, quietly = TRUE) # added plotting function
library(ggtext, warn.conflicts = FALSE, quietly = TRUE) # more plotting
library(DataExplorer, warn.conflicts = FALSE, quietly = TRUE) # quick exploratory vis
library(corrplot, warn.conflicts = FALSE, quietly = TRUE) # plotting corrmatrix
library(alookr, warn.conflicts = F, quietly = T) # for removing correlated variables
library(proxy, warn.conflicts = FALSE, quietly = TRUE) # for computing dissimilarity
library(factoextra, warn.conflicts = FALSE, quietly = TRUE) # visualizing clustering
library(ggdendro, warn.conflicts = FALSE, quietly = TRUE) # for some clever dendrograms

file <- "210228.1s_clean.txt"
rawData <- read.csv(file, sep = "\t", header = TRUE)
rawData <- rawData %>% mutate(Utc = round(Utc*24*60*60, digits = 0) - round(rawData$Utc[1]*24*60*60, digits = 0))
names(rawData) <- c("sec", "lat", "long", "heel", "bsp", "awa", "aws", "leeway", "course", "twd", "twa", "tws", "sog", "cog", "drift", "set", "hdg")


dataCluster<- rawData %>% select(lat, long, heel, bsp, awa, twa, sog, cog, hdg) %>% drop_na()
dataCluster <- dataCluster %>%  filter(between(lat,-3351, -3348)) %>% filter(between(long, 15115, 15118))


corMatrix <- round(cor(dataCluster, method = "pearson"), 2) # calculate correlation matrix, using pearson
corrplot.mixed(corMatrix, order = "AOE") # plot correlation matrix


dataCluster_cor <- dataCluster %>% select(heel, bsp, twa, cog)

corMatrix_cor <- round(cor(dataCluster_cor, method = "pearson"), 2) # calculate correlation matrix, using pearson
corrplot.mixed(corMatrix_cor, order = "AOE") # plot correlation matrix

modelData_df <- dataCluster_cor # move data into model data frame
modelData_matrix <- as.matrix(modelData_df) # create matrix of model data

## calculate dissimilarity matrix using various methods.
dissimilarityArray_all <- array(dim = c(16362, 16362, 3))
dissimilarityArray_all[ , ,1] <- as.matrix(dist(modelData_matrix, method = "cosine")) # dissimilarity based on cosine distance and save as the first matrix in array
dissimilarityArray_all[ , ,2] <- as.matrix(dist(modelData_matrix, method = "euclidean")) # dissimilarity based on euclidean distance and save as second matrix in array
dissimilarityArray_all[ , ,3] <- as.matrix(dist(modelData_matrix, method = "manhattan"))


names_dis <- c("cosine", "eucledian", "manhattan") # vector of matrix identifiers
methods <- c("single", "complete", "average", "ward") # vector of method types
count <- 0 # set counter
ac_value <- c() # empty vector to store Ac value for quick access later
ac_label <- c() # empty vector to stor AC labels
ac_sim <- c() # empty vector to store the matrix type for each ac score

model <- agnes(dissimilarityArray_all[ , ,1], diss = TRUE, method = "ward") # generate model

for (i in 3:3) { # loop over array full of dissimilarity matrix 
  for (m in 4:4){ # loop over methods
    count <- count + 1 # increase counter by 1
    
    variable_id <- paste("cluster",names_dis[i], methods[m], sep = "_") # generate variable id
    model <- agnes(dissimilarityArray_all[ , ,i], diss = TRUE, method = methods[m]) # generate model
    assign(variable_id, model) # assign model to variable id
    ac_value[count] <- model$ac # record ac value
    ac_label[count] <- paste(names_dis[i], methods[m], sep = " ") # record label to label ac value
    ac_sim[count] <- methods[m] # record similarity matrix used
    
  }
}
