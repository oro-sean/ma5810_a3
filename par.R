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

modelData <- rawData %>% select(sec, lat, long, heel, bsp, awa, course, twa, sog, cog, hdg) # move raw data into modelData dataframe
modelData$group <- cutree(model, k = 12) # assign each observation a group

## vis final cluster for assessment
final_clusters <- ggplot(modelData, aes(x = long, y = lat)) + coord_quickmap() +
  geom_point(aes(colour = as.factor(modelData$group)), shape = ".", alpha = 0.5) + theme(legend.position = "none") + facet_wrap(modelData$group)
final_clusters

clustersToKeep = c(3, 4, 5, 6, 8, 9, 11, 12) # select clusters to keep

## Compute vector components
modelData <- modelData %>% mutate(water = pol2cart(modelData$bsp, modelData$course, degrees = TRUE)) %>% # vector representing velocity through water
  mutate(ground = pol2cart(modelData$sog, modelData$cog, degrees = TRUE)) %>% # vector representing velocity over ground
  mutate(dif_grnd_wtr = ground - water) # difference between velocity over ground and velocity through water, should be equal to tide + sensor error
modelData$difX <- modelData$dif_grnd_wtr$x # x component of difference for easy access
modelData$difY <- modelData$dif_grnd_wtr$y # y component of difference for easy access

## Inspect difference for dependence on time
dif_df <- modelData %>% select(sec, difX, difY) # create data frome for PCA
PCA <- prcomp(dif_df, scale = TRUE, center = TRUE) # calculate PCA
check_timeDiff <- fviz_pca_biplot(PCA, geom = "point") # plot data and PCA
check_timeDiff
rm(dif_df, PCA) # remove variables not required

## considering each cluster individually, calculate LOF for the observation and difference and knn for position
modelData <- foreach( i = clustersToKeep, .combine = rbind) %do%
  {
    subset <- modelData %>% filter(group == i) # select each group
    
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

## for easy visualization of cut off's create truncated lof and knn
modelData <- modelData %>% mutate(lof_metric_thresholds = ifelse(lof_metrics <= 1, 1, 
                                                                 ifelse(lof_metrics < 1.5, DescTools::RoundTo(lof_metrics, multiple = 0.1, FUN = trunc), 
                                                                        1.5)))
modelData <- modelData %>% mutate(lof_tide_threshold = ifelse(lof_tide <= 0.8, 0.8, 
                                                              ifelse(lof_tide < 1.5, DescTools::RoundTo(lof_tide, multiple = 0.1, FUN = trunc),
                                                                     1.5)))
modelData <- modelData %>% mutate(knn_position_threshold = ifelse(knn_position <= 0.3, 0.3, 
                                                                  ifelse(knn_position < 0.6, DescTools::RoundTo(knn_position, multiple = 0.1, FUN = trunc),
                                                                         0.6)))

## visualize lof and knn 
lof_metrics_plot <- ggplot(modelData, aes(x = long, y = lat)) + coord_quickmap() + 
  geom_point(aes(colour = as.factor(lof_metric_thresholds)), stroke = .001, alpha = 0.5) + theme(legend.position = "bottom") + facet_wrap(modelData$group)
lof_metrics_plot

lof_tide_plot <- ggplot(modelData, aes(x = difX, y = difY)) + geom_point(aes(colour = as.factor(lof_tide_threshold)), stroke = .001, alpha = 0.5) + theme(legend.position = "bottom") 
lof_tide_plot

knn_position_plot <- ggplot(modelData, aes(x = long, y = lat)) + coord_quickmap() + 
  geom_point(aes(colour = as.factor(knn_position_threshold)), stroke = .001, alpha = 0.5) + theme(legend.position = "bottom") + facet_wrap(modelData$group)
knn_position_plot

## remove outliers
modelData <- modelData %>% filter(lof_metrics <= 1.2 & lof_tide < 0.9 & knn_position < 0.5)

## plot cleaned data
clean_plot <- ggplot(modelData, aes(x = long, y = lat)) + coord_quickmap() +
  geom_point(aes(colour = as.factor(modelData$group)), shape = ".", alpha = 1) + theme(legend.position = "none")
clean_plot

## define function for calculating and extracting PC's
TIDE <- function(difference,component){
  PCA <- prcomp(difference, scale = FALSE, center = FALSE)
  x <- PCA$rotation[ ,1]*(PCA$sdev[1]/(sum(PCA$sdev)))
  x <- x[component]
  x
  
}


#####  do this taking the average then cosine distance   FTW
compCal <- COMPCAL(modelData)

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

###

combin <- modelData %>% group_by(group)


modelData <- modelData %>% filter(group == 3)

data <- ggplot(modelData[seq(1,nrow(modelData), 30),], aes(x = long, y = lat)) + coord_quickmap() +
   geom_segment(aes(xend = long + 0.1 * tide_x, yend = lat + .1 * tide_y), colour = "blue", alpha = .5, arrow = arrow())
data


sample$sec[-10]




comp_cor <- c(-1.5,-1,-0.5,0,0.5,1,1.5)

dif_tide <- 
  foreach(o = comp_cor, .combine = c) %do%
  {
    modelData <- modelData %>% mutate(course = ifelse(course + o < 360 & course + o >= 0, course + o, ifelse(course + o < 0, 360 - course + o, course + o - 360)))
    
    modelData <-
      foreach(g = unique(modelData$group), .combine = rbind) %do%
      {
        sample <- modelData %>% filter(group == g)
        
        tide_x <- foreach(i = seq(1:length(sample$sec)), .combine = c) %do%
          {
            s <- i-20
            e <- i+20
            
            ifelse(s < 1 | e > length(sample$sec), NA, ifelse(
              (sample$sec[e]-sample$sec[s])>50,NA, TIDE(sample[s:e, 28:29],1) 
            ))
            
          }
        
        tide_y <- foreach(i = seq(1:length(sample$sec)), .combine = c) %do%
          {
            s <- i-20
            e <- i+20
            
            ifelse(s < 1 | e > length(sample$sec), NA, ifelse(
              (sample$sec[e]-sample$sec[s])>50,NA, TIDE(sample[s:e, 28:29],2) 
            ))
            
          }
        sample$tide_x <- tide_x
        sample$tide_y <- tide_y
        sample
        
      }

    tide <- modelData[ , 34:35]
    sum(abs(dist(tide)))
    
  }

comp_cor <- comp_cor[which.min(dif_tide)] 




modelData <- modelData %>% mutate(course = ifelse(course + comp_cor < 360 & course + comp_cor >= 0, course 
                                                  + comp_cor, ifelse(course + comp_cor < 0, 360 - course 
                                                                     + comp_cor, course + comp_cor - 360)))




foreach(i = unique(modalData$group))
modelData$knn_dist <- kNNdist(x = modelData[ , 28:29], k = 180)







diference_df <- data.frame(modelData$dif_grnd_wtr$x,modelData$dif_grnd_wtr$y)
PCA <- prcomp(diference_df, scale = FALSE)
tideVector <- PCA$rotation[ ,1]*(PCA$sdev[1]/(sum(PCA$sdev)))

modelData$waterAndTide <- modelData$water + tide
modelData$modelBsp <- modelData$waterAndTide$r
modelData$bspError <- modelData$bsp_avg - modelData$modelBsp
modelData$hdgModel <- modelData$waterAndTide$theta
modelData$hdgError <- modelData$course_avg - modelData$hdgModel

errorModel_BSP <- lm(bspError ~ bsp_avg + heel, data = modelData)

errorModel_COURSE <- lm(hdgError ~ leeCoef_avg, data = modelData)

 
data_clean <- ggplot(modelData, aes(x = long, y = lat)) + coord_quickmap() + 
  geom_point(aes(colour = bspError, stroke = .001)) + theme(legend.position = "bottom")
data_clean



