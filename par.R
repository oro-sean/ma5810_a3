## set up environment and import csv file

rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")
setwd("~/GitHub/ma5810_a3")

startTime <- Sys.time()

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
library(doParallel, warn.conflicts = FALSE, quietly = TRUE)
library(parallel, warn.conflicts = FALSE, quietly = TRUE)
library(foreach, warn.conflicts = FALSE, quietly = TRUE)
library(dbscan)
library(useful, warn.conflicts = FALSE, quietly = TRUE)
library(zoo)

file <- "210228.1s_clean.txt"
rawData <- read.csv(file, sep = "\t", header = TRUE)
names(rawData) <- c("sec", "lat", "long", "heel", "bsp", "awa", "aws", "leeway", 
                    "course", "twd", "twa", "tws", "sog", "cog", "drift", "set",
                    "hdg")
rawData <- rawData[4500:19000, ] %>% mutate(sec = round(sec*24*60*60, digits = 0) - round(rawData$sec[1]*24*60*60, digits = 0)) %>%  
  filter(between(lat,-3351, -3348)) %>% filter(between(long, 15115, 15118))


clusteringFactors<- rawData %>% select(lat, long, heel, bsp, awa, twa, sog, cog, hdg) %>% 
  drop_na()

corMatrix <- round(cor(clusteringFactors, method = "pearson"), 2) # calculate correlation matrix, using pearson
corrplot.mixed(corMatrix, order = "AOE") # plot correlation matrix


clusteringFactors_clean <- clusteringFactors %>% select(heel, bsp, twa, cog)

corMatrix_cor <- round(cor(clusteringFactors_clean, method = "pearson"), 2) # calculate correlation matrix, using pearson
corrplot.mixed(corMatrix_cor, order = "AOE") # plot correlation matrix

modelData_df <- clusteringFactors_clean # move data into model data frame
modelData_matrix <- as.matrix(modelData_df) # create matrix of model data

cl <- makeCluster(10)
registerDoParallel(cl)
startTime <- Sys.time()

dis_type <- c("cosine", "euclidean", "manhattan") # vector of matrix identifiers
methods <- c("single", "complete", "average", "ward") # vector of method types

allModells_ac <- 
  foreach(i = 1:3, .combine = c) %:%
    foreach(m = 1:4, .combine = c) %dopar%
      {
        library(cluster, warn.conflicts = FALSE, quietly = TRUE) # for clustering
        library(proxy, warn.conflicts = FALSE, quietly = TRUE) # for computing dissimilarity

        dissimilarityMatrix <- as.matrix(dist(modelData_matrix, method = dis_type[i]))
        saveRDS(dissimilarityMatrix, file = paste("dis", dis_type[i],".RDS", sep = "_" ))
        
        model <- agnes(dissimilarityMatrix, diss = TRUE, method = methods[m]) # generate model
        saveRDS(model, file = paste("model", methods[m], dis_type[i],".RDS", sep = "_" ))
        
        model$ac

}

finishTime <- Sys.time()
(finishTime - startTime)

stopCluster(cl)

dissimilarityMatrix <- readRDS("dis_cosine_.RDS")
model <- readRDS("model_ward_euclidean_.RDS")

cl <- makeCluster(10)
registerDoParallel(cl)
startTime <- Sys.time()

swc_2to24 <- 
  foreach(i = 2:24, .combine = cbind) %dopar%
    { 
      library(cluster, warn.conflicts = FALSE, quietly = TRUE) # for clustering
      sil <- silhouette(cutree(model, k = i), dissimilarityMatrix) # calculate SWC
      sil <- summary(sil) # extra summary so avg can be easily accessed
      sil$avg.width # record avg
    }

finishTime <- Sys.time()
(finishTime - startTime)

stopCluster(cl)



cl <- makeCluster(10)
registerDoParallel(cl)
startTime <- Sys.time()

hdg_stDev <- 
  foreach(i = 2:24) %dopar%
  {
    clusteringFactors$group <- cutree(model, k = i)
    stDev <- aggregate(hdg ~ group, clusteringFactors, function(x) stDev = sd(x))
    median(stDev$hdg)
  
}

finishTime <- Sys.time()
(finishTime - startTime)

stopCluster(cl)

rawData$group <- cutree(model, k = 8)

cluster_1 <- modelData %>% filter(group == 3)
cluster_1$LOF <- lof(x=cluster_1, minPts = 100)
cluster_1 <- cluster_1 %>% filter(LOF < 1.5)

cluster_2 <- modelData %>% filter(group == 5)
cluster_2$LOF <- lof(x=cluster_2, minPts = 100)
cluster_2 <- cluster_2 %>% filter(LOF < 1.5)

cluster_3 <- modelData %>% filter(group == 7)
cluster_3$LOF <- lof(x=cluster_3, minPts = 100)
cluster_3 <- cluster_3 %>% filter(LOF < 1.5)

cluster_4 <- modelData %>% filter(group == 8)
cluster_4$LOF <- lof(x=cluster_4, minPts = 100)
cluster_4 <- cluster_4 %>% filter(LOF < 1.5)

modelData <- rbind(cluster_1, cluster_2, cluster_3, cluster_4)
rawData <- modelData

data <- ggplot(rawData, aes(x = long, y = lat)) + coord_quickmap() + 
  geom_point(aes(colour = hdg, stroke = .001)) + theme(legend.position = "bottom") + facet_wrap(rawData$group)
data

modelData <- rawData %>% filter(group == 3 | group == 5 | group == 7 | group == 8)

cluster_1 <- modelData %>% filter(group == 4)
cluster_1$LOF <- lof(x=cluster_1, minPts = 100)
cluster_1 <- cluster_1 %>% filter(LOF < 1.5)

cluster_2 <- modelData %>% filter(group == 5)
cluster_2$LOF <- lof(x=cluster_2, minPts = 100)
cluster_2 <- cluster_2 %>% filter(LOF < 1.5)

cluster_3 <- modelData %>% filter(group == 6)
cluster_3$LOF <- lof(x=cluster_3, minPts = 100)
cluster_3 <- cluster_3 %>% filter(LOF < 1.5)

cluster_4 <- modelData %>% filter(group == 7)
cluster_4$LOF <- lof(x=cluster_4, minPts = 100)
cluster_4 <- cluster_4 %>% filter(LOF < 1.5)

modelData <- rbind(cluster_1, cluster_2, cluster_3, cluster_4)
rawData <- modelData

dif_tide <- 
  foreach(o = c(-1.5,-1,-0.5,0,0.5,1,1.5), .combine = c) %do%
  {
    rawData <- rawData %>% mutate(course = ifelse(course + o < 360 & course + o >= 0, course + o, ifelse(course + o < 0, 360 - course + o, course + o - 360)))
    
    tide <- foreach(i  = unique(rawData$group), .combine = rbind) %do%
      {
        GetTide <- function(sample, roll_mean)
        {
          sample <- sample %>% 
            mutate(bsp_avg = rollmean(sample$bsp, k = roll_mean, fill = "extend", align = "center"),
                   course_avg = rollmean(sample$course, k = roll_mean, fill = "extend", align = "center"),
                   sog_avg = rollmean(sample$sog, k = roll_mean, fill = "extend", align = "center"),
                   cog_avg = rollmean(sample$cog, k = roll_mean, fill = "extend", align = "center")) 
          
          sample$water <- pol2cart(sample$bsp_avg, sample$course_avg, degrees = TRUE)
          sample$ground <- pol2cart(sample$sog_avg, sample$cog_avg, degrees = TRUE)
          sample$dif_grnd_wtr <- sample$ground - sample$water
          dif_df <- data.frame(sample$dif_grnd_wtr$x,sample$dif_grnd_wtr$y)
          PCA <- prcomp(dif_df, scale = FALSE)
          PCA$rotation[ ,1]*(PCA$sdev[1]/(sum(PCA$sdev)))
        }
        sample <- rawData %>% filter(group == i)
        GetTide(sample,30)
      }
    
    sum(abs(dist(tide)))
    
  }

comp_cor <- c(-1.5,-1,-0.5,0,0.5,1,1.5)
o <- comp_cor[which.min(dif_tide)] 




modelData <- rawData %>% mutate(course = ifelse(course + o < 360 & course + o >= 0, course + o, ifelse(course + o < 0, 360 - course + o, course + o - 360)))


modelData <- modelData %>% 
  mutate(bsp_avg = rollmean(modelData$bsp, k = 60, fill = "extend", align = "center"), 
         course_avg = rollmean(modelData$course, k = 30, fill = "extend", align = "center"),
         sog_avg = rollmean(modelData$sog, k = 60, fill = "extend", align = "center"),
         cog_avg = rollmean(modelData$cog, k = 60, fill = "extend", align = "center"),
         heel_avg_abs = rollmean(abs(modelData$heel), k = 30, fill = "extend", align = "center"),
         leeCoef_avg = rollmean(modelData$heel/(modelData$bsp^2), k=30, fill = "extend", align = "center")) 

modelData$water <- pol2cart(modelData$bsp_avg, modelData$course_avg, degrees = TRUE)
modelData$ground <- pol2cart(modelData$sog_avg, modelData$cog_avg, degrees = TRUE)
modelData$dif_grnd_wtr <- modelData$ground - modelData$water

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



