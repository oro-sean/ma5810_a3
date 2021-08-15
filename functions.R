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
  filter(between(lat,-3351, -3348)) %>% filter(between(long, 15115, 15118)) %>% mutate(stb = twa / abs(twa)) %>% mutate(leeCalc = heel / (bsp^2))

model <- readRDS("model_ward_euclidean_.RDS")

rawData_grouped <- rawData %>% select(sec, lat, long, heel, bsp, awa, course, twa, sog, cog, hdg, stb, leeCalc) # move raw data into modelData dataframe
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

  ## remove outliers and average BSP ect and recalculate differences
  Data <- Data %>% filter(lof_metrics <= 1.2 & lof_tide <= 1.2 & knn_position <= 0.5)
  Data$heel <- rollmean(Data$heel, k = 60, fill = "extend", align = "center")
  Data$bsp <- rollmean(Data$bsp, k = 60, fill = "extend", align = "center")
  Data$course <- rollmean(Data$course, k = 60, fill = "extend", align = "center")
  Data$sog <- rollmean(Data$sog, k = 60, fill = "extend", align = "center")
  Data$cog <- rollmean(Data$cog, k = 60, fill = "extend", align = "center")
  Data$hdg <- rollmean(Data$hdg, k = 60, fill = "extend", align = "center")
  
  Data <- groupedData %>% mutate(water = pol2cart(groupedData$bsp, groupedData$course, degrees = TRUE)) %>% # vector representing velocity through water
    mutate(ground = pol2cart(groupedData$sog, groupedData$cog, degrees = TRUE)) %>% # vector representing velocity over ground
    mutate(dif_grnd_wtr = ground - water) # difference between velocity over ground and velocity through water, should be equal to tide + sensor error
  Data$difX <- Data$dif_grnd_wtr$x # x component of difference for easy access
  Data$difY <- Data$dif_grnd_wtr$y # y component of difference for easy access
  
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
          s <- i-40
          e <- i+40
          
          ifelse(s < 1 | e > length(sample$sec), NA, ifelse(
            (sample$sec[e]-sample$sec[s])>90,NA, TIDE(sample[s:e, 2:3],1) 
          ))
          
        }
      
      tide_y <- foreach(i = seq(1:length(sample$sec)), .combine = c) %do%
        {
          s <- i-40
          e <- i+40
          
          ifelse(s < 1 | e > length(sample$sec), NA, ifelse(
            (sample$sec[e]-sample$sec[s])>90,NA, TIDE(sample[s:e, 2:3],2) 
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

###

rawData_offset$bsp <- rawData_offset$bsp - (rawData_offset$bsp * bspcor + rawData_offset$heel * heelcor)

rawData_offset$course <- rawData_offset$course + rawData_offset$leeCalc * hdgcor * rawData_offset$stb

rawData_offset <- drop_na(rawData_offset)

modelData <- GET_MODEL_DATA(rawData_offset, clustersToKeep)

tideData <- modelData %>% select(sec, difX, difY, group)

tideField <- TIDEFIELD(tideData)

modelData$tideX <- tideField$tideX
modelData$tideY <- tideField$tideY
modelData$waterX <- modelData$water$x
modelData$waterY <- modelData$water$y

regresionData <- modelData %>% select(sec, lat, long, heel, bsp, course, hdg, group, stb, difX, difY, tideX, tideY, waterX, waterY, leeCalc) %>% drop_na() %>%
  mutate(errorX = difX - tideX) %>% mutate(errorY = difY - tideY) 

data <- ggplot(regresionData[seq(1,nrow(regresionData), 90),], aes(x = long, y = lat)) + coord_quickmap() +
  geom_segment(aes(xend = long + 0.1 * difX, yend = lat + .1 * difY), colour = "green", alpha = .5, arrow = arrow(length = unit(1,"mm"))) +
  geom_segment(aes(xend = long + 0.1 * tideX, yend = lat + .1 * tideY), colour = "blue", alpha = .5, arrow = arrow(length = unit(1,"mm")))+
  geom_segment(aes(xend = long + 0.1 * errorX, yend = lat + .1 * errorY), colour = "red", alpha = .5, arrow = arrow(length = unit(1,"mm")))
data

data <- ggplot(regresionData[seq(1,nrow(regresionData), 90),], aes(x = long, y = lat)) + coord_quickmap() +
  geom_segment(aes(xend = long + 0.1 * tideX, yend = lat + .1 * tideY, colour = sec), alpha = 1, arrow = arrow(length = unit(2,"mm")))
data


ERROR_ON_BSP <- function(vec)
  {
dot <- vec$waterX * vec$errorX + vec$waterY + vec$errorY
x <- dot/(vec$bsp^2)*vec$waterX
y <- dot/(vec$bsp^2)*vec$waterY
pol <- cart2pol(x, y, degrees = TRUE)
pol$r
}

regresionData$bspError <-
  foreach(i = seq(1,nrow(regresionData),1), .combine = c) %do%
  {
    vec <- regresionData[i, ]
    error <- ERROR_ON_BSP(vec)
    error
  }


HDG_ERROR <- function(vec)
{
  hdgX <- vec$waterX + vec$errorX
  hdgY <- vec$waterY + vec$errorY
  hdgError <- cart2pol(hdgX, hdgY, degrees = TRUE)
  (hdgError$theta - vec$course) * vec$stb
}

regresionData$hdgerror <-
  foreach(i = seq(1,nrow(regresionData),1), .combine = c) %do%
  {
    vec <- regresionData[i, ]
    error <- HDG_ERROR(vec)
    error
  }


plot(regresionData$bsp, regresionData$bspError)
plot(regresionData$heel, regresionData$bspError)
plot(regresionData$leeCalc, regresionData$hdgerror)

regresionData <- regresionData %>% filter(bspError/bsp < 0.2 & abs(hdgerror) < 5 )
errorModel_BSP <- lm(bspError ~ bsp + heel, data = regresionData)
errorModel_HDG <- lm(hdgerror ~ leeCalc, data = regresionData)


coef <- summary(errorModel_BSP)
coef <- coef$coefficients
bspcor <- coef[2,1]
heelcor <- coef[3,1]

coef <- summary(errorModel_HDG)
coef <- coef$coefficients
hdgcor <- coef[2,1]


