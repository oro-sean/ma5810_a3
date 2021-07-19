## set up environment and import csv file

rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")
setwd("~/Documents/GitHub/ma5810_a3")


library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)

set.seed(123)

file <- "210228.1s_clean.txt"

rawData <- read.csv(file, sep = "\t", header = TRUE)
rawData <- rawData %>% mutate(Utc = round(Utc*24*60*60, digits = 0) - round(rawData$Utc[1]*24*60*60, digits = 0)) %>% rename("Aest" = "Utc")

dataTrack<- rawData %>% select(GPSLat,GPSLon,Heel,Boatspeed, TW_angle) %>% drop_na()
dataTrack <- dataTrack %>%  filter(between(GPSLat,-3351, -3348)) %>% filter(between(GPSLon, 15115, 15118))

Heel <- ggplot(dataTrack, aes(x = GPSLon, y = GPSLat)) + coord_quickmap() + 
  geom_point(aes(colour = Heel)) + theme(legend.position = "bottom")
Boatspeed <- ggplot(dataTrack, aes(x = GPSLon, y = GPSLat)) + coord_quickmap() + 
  geom_point(aes(colour = Boatspeed))+ theme(legend.position = "bottom") +
  scale_colour_gradient2()
TW_angle <- ggplot(dataTrack, aes(x = GPSLon, y = GPSLat)) + coord_quickmap() + 
  geom_point(aes(colour = TW_angle))+ theme(legend.position = "bottom") +
  scale_colour_gradientn(colours = terrain.colors(10))
ggarrange(Heel,Boatspeed, TW_angle + rremove("x.text"), 
          labels = c("Heel","Boatspeed", "TW_angle"), ncol = 3, nrow = 1 )
                            
dataCluster <- rawData %>% select(GPSLat,GPSLon,Heel,TW_angle,HDG, Boatspeed) %>% drop_na() %>%  filter(between(GPSLat,-3351, -3348.75)) %>% filter(between(GPSLon, 15115.75, 15118))
dataclustered <- dataCluster %>% select(TW_angle,HDG, Heel, Boatspeed)
km.out <- kmeans(x = dataclustered, centers = 12)

km.out
dataCluster$out <- km.out$cluster


dataPlotin <- dataCluster %>% filter(out == 5 | out == 7 | out == 1 | out == 4 | out ==12)
dataPlotout <- dataCluster %>% filter(out == 3 | out == 8 | out == 9 | out == 10| out == 6 | out == 2)
data <- ggplot(dataCluster, aes(x = GPSLon, y = GPSLat)) + coord_quickmap() + 
  geom_point(aes(colour = as.factor(dataCluster$out), stroke = .001)) + theme(legend.position = "bottom")
dataIn <- ggplot(dataPlotin, aes(x = GPSLon, y = GPSLat)) + coord_quickmap() + 
  geom_point(aes(colour = as.factor(dataPlotin$out), stroke = .001)) + theme(legend.position = "bottom")

dataOut <- ggplot(dataPlotout, aes(x = GPSLon, y = GPSLat)) + coord_quickmap() + 
  geom_point(aes(colour = as.factor(dataPlotout$out), stroke = .001)) + theme(legend.position = "bottom")

ggarrange(data,dataIn,dataOut + rremove("x.text"), 
          labels = c("All", "In","Out"), ncol = 3, nrow = 1)
boxplot(km.out$withinss)

ss <- km.out$withinss
plot(ss)
