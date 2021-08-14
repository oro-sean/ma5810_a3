## set up environment and import csv file

rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")
setwd("~/GitHub/ma5810_a3")

library(cluster, warn.conflicts = FALSE, quietly = TRUE) # for clustering
library(tidyverse, warn.conflicts = FALSE, quietly = TRUE) # handy for data prep
library(ggplot2, warn.conflicts = FALSE, quietly = TRUE) # plotting
library(ggpubr, warn.conflicts = FALSE, quietly = TRUE) # added plotting function
library(ggtext, warn.conflicts = FALSE, quietly = TRUE) # more plotting
library(alookr, warn.conflicts = F, quietly = T) # for removing correlated variables
library(proxy, warn.conflicts = FALSE, quietly = TRUE) # for computing dissimilarity
library(factoextra, warn.conflicts = FALSE, quietly = TRUE) # visualizing clustering
library(doParallel, warn.conflicts = FALSE, quietly = TRUE)
library(parallel, warn.conflicts = FALSE, quietly = TRUE)
library(foreach, warn.conflicts = FALSE, quietly = TRUE)
library(dbscan, warn.conflicts = FALSE, quietly = TRUE)
library(useful, warn.conflicts = FALSE, quietly = TRUE)
library(zoo, warn.conflicts = FALSE, quietly = TRUE)

file <- "210228.1s_clean.txt"
rawData <- read.csv(file, sep = "\t", header = TRUE)
names(rawData) <- c("sec", "lat", "long", "heel", "bsp", "awa", "aws", "leeway", 
                    "course", "twd", "twa", "tws", "sog", "cog", "drift", "set",
                    "hdg")
rawData <- rawData[4500:19000, ] %>% mutate(sec = round(sec*24*60*60, digits = 0) - round(rawData$sec[1]*24*60*60, digits = 0)) %>%  
  filter(between(lat,-3351, -3348)) %>% filter(between(long, 15115, 15118))

model <- readRDS("model_ward_euclidean_.RDS")

rawData_grouped <- rawData %>% select(sec, lat, long, heel, bsp, awa, course, twa, sog, cog, hdg) # move raw data into modelData dataframe
rawData_grouped$group <- cutree(model, k = 12) # assign each observation a group

clustersToKeep <- c(3, 4, 5, 6, 8, 9, 11, 12) # select clusters to keep

##### Functions ######
## define function to change course
CHANGE_COURSE <- function(rawData_grouped,offset)
{
  rawData_grouped %>% mutate(course = ifelse(course + offset < 360 & course + offset >= 0, course 
                                             + offset, ifelse(course + offset < 0, 360 - course 
                                                              + offset, course + offset - 360)))
}

## defien function to get model data
GET_MODEL_DATA <- function(groupedData, clustersToKeep)
{
  Data <- groupedData %>% mutate(water = pol2cart(groupedData$bsp, groupedData$course, degrees = TRUE)) %>% # vector representing velocity through water
    mutate(ground = pol2cart(groupedData$sog, groupedData$cog, degrees = TRUE)) %>% # vector representing velocity over ground
    mutate(dif_grnd_wtr = ground - water) # difference between velocity over ground and velocity through water, should be equal to tide + sensor error
  Data$difX <- Data$dif_grnd_wtr$x # x component of difference for easy access
  Data$difY <- Data$dif_grnd_wtr$y # y component of difference for easy access

  ## considering each cluster individually, calculate LOF for the observation and difference and knn for position
  Data <- foreach( i = clustersToKeep, .combine = rbind) %do%
    {
      subset <- Data %>% filter(group == i) # select each group
      
      subset_metrics <- subset %>% select(sec, heel, bsp, course, twa) # select factors for lof
      preproc <- caret::preProcess(subset_metrics, method = c("center", "scale")) # set up for center and scaling
      subset_metrics <- predict(preproc, subset_metrics)
      subset$lof_metrics <- lof(x=subset_metrics, minPts = 60) # calculate lof
      
      subset_position <- subset %>% select(sec, lat, long) # select factors for knn of position
      preproc <- caret::preProcess(subset_position, method = c("center", "scale")) # set up for center and scaling
      subset_position <- predict(preproc, subset_position)
      subset$knn_position <- kNNdist(x = subset_position, k = 30) # calculate knn using position and time
      
      subset_tide <- subset %>% select(sec, difX, difY) # select factors for knn
      preproc <- caret::preProcess(subset_tide, method = c("center", "scale")) # set up for center and scaling
      subset_tide <- predict(preproc, subset_tide)
      subset$lof_tide <- lof(x = subset_tide, minPts = 180) # calculate knn using difference vector
      
      subset
    }

  ## remove outliers
  Data <- Data %>% filter(lof_metrics <= 1.2) %>% filter(lof_tide <= 1.2) %>% filter(knn_position <= 0.5)
  #lof_tide < 0.9 & knn_position < 0.5)
  
  Data
}

## define function for calculating and extracting PC's
TIDE <- function(difference,component){
  PCA <- prcomp(difference, scale = FALSE, center = FALSE)
  x <- PCA$rotation[ ,1]*(PCA$sdev[1]/(sum(PCA$sdev)))
  x <- x[component]
  x
  
}

## Define function to calculate the average tide for each group

COMPCAL <- function(modelData)
{
  tideData <- modelData %>% select(sec, difX, difY, group)
  
  foreach(g = unique(tideData$group), .combine = rbind) %do%
    {
      sample <- tideData %>% filter(group == g)
      
      tide_x <- foreach(i = seq(1:length(sample$sec)), .combine = c) %do%
        {
          s <- i-20
          e <- i+20
          
          ifelse(s < 1 | e > length(sample$sec), NA, ifelse(
            (sample$sec[e]-sample$sec[s])>50,NA, TIDE(sample[s:e, 2:3],1) 
          ))
          
        }
      
      tide_y <- foreach(i = seq(1:length(sample$sec)), .combine = c) %do%
        {
          s <- i-20
          e <- i+20
          
          ifelse(s < 1 | e > length(sample$sec), NA, ifelse(
            (sample$sec[e]-sample$sec[s])>50,NA, TIDE(sample[s:e, 2:3],2) 
          ))
          
        }
      groupAvg <- data.frame(g, mean(tide_x, na.rm = TRUE), mean(tide_y, na.rm = TRUE))
      groupAvg
      
    }
}

TIDEFIELD <- function(modelData)
{

  foreach(g = unique(modelData$group), .combine = rbind) %do%
    {
      sample <- modelData %>% filter(group == g)
      
      tide_x <- foreach(i = seq(1:length(sample$sec)), .combine = c) %do%
        {
          s <- i-20
          e <- i+20
          
          ifelse(s < 1 | e > length(sample$sec), NA, ifelse(
            (sample$sec[e]-sample$sec[s])>50,NA, TIDE(sample[s:e, 2:3],1) 
          ))
          
        }
      
      tide_y <- foreach(i = seq(1:length(sample$sec)), .combine = c) %do%
        {
          s <- i-20
          e <- i+20
          
          ifelse(s < 1 | e > length(sample$sec), NA, ifelse(
            (sample$sec[e]-sample$sec[s])>50,NA, TIDE(sample[s:e, 2:3],2) 
          ))
          
        }
      sample$tideX <- tide_x
      sample$tideY <- tide_y
      sample
      
    }
}

comp_offsets <- c(-2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2)

totalDist <- foreach(i = comp_offsets, .combine = c) %do%
  {
    rawData_offset <- CHANGE_COURSE(rawData_grouped, i) 
    
    modelData <- GET_MODEL_DATA(rawData_offset, clustersToKeep)
    
    compcal <- COMPCAL(modelData)
    
    sum(proxy::dist(compcal, method = "cosine"))
    
  }

rawData_offset <- CHANGE_COURSE(rawData_grouped, comp_offsets[which.min(totalDist)]) 

modelData <- GET_MODEL_DATA(rawData_offset, clustersToKeep)

modelData <- TIDEFIELD(modelData)

modelData <- drop_na(modelData)

data <- ggplot(modelData[seq(1,nrow(modelData), 30),], aes(x = long, y = lat)) + coord_quickmap() +
  geom_segment(aes(xend = long + 0.1 * tideX, yend = lat + .1 * tideY), colour = "blue", alpha = .5) + facet_wrap(modelData[seq(1,nrow(modelData), 30),]$group)
data


