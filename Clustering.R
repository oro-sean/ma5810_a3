
rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")
setwd("~/GitHub/ma5810_a3")

startTime <- Sys.time()

library(tidyverse, warn.conflicts = FALSE, quietly = TRUE) # handy for data prep
library(doParallel, warn.conflicts = FALSE, quietly = TRUE)
library(parallel, warn.conflicts = FALSE, quietly = TRUE)
library(foreach, warn.conflicts = FALSE, quietly = TRUE)



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

data <- ggplot(rawData, aes(x = long, y = lat)) + coord_quickmap() + 
  geom_point(aes(colour = hdg, stroke = .001)) + theme(legend.position = "bottom") + facet_wrap(rawData$group)
data