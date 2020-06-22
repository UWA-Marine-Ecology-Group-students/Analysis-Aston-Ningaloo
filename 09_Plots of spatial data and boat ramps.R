############ Import data ##############
## Load libraries --
library(raster)
library(dplyr)
library(ggplot2)
library(sp)
library(rgdal)
library(rgeos)
library(spatialEco)

## Set work directory----
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # sets working directory to where this script is saved (DON'T MOVE THE SCRIPT)

## Set sub directories----
d.dir <- paste(working.dir,"Tidy data",sep="/") 
s.dir <- paste(working.dir,"Spatial",sep="/") # spatial is where I keep spatial data files, rasters and shapefiles
p.dir <- paste(working.dir,"Plots",sep="/")

## Load bathymetry and tpi data----
setwd(s.dir)
bathy <- raster(paste(s.dir, "bathy.tif", sep='/'))
bathy <- flip(bathy, direction="y")
plot(bathy) # check it
proj4string(bathy) # check the coordinate system, resolution, etc..


tpi <- raster(paste(s.dir, "tpi.tif", sep='/'))
tpi <- flip(tpi, direction="y")
plot(tpi)
proj4string(tpi) # check the coordinate system, resolution, etc..

## Load metadata
setwd(d.dir)
metadata <- read.csv("ningaloo.checked.metadata.csv")
metadata

samples <- metadata %>% 
  select("sample", "latitude", "longitude")

## Convert lat and long of samples into spatial data and assign correct coordinate system
coordinates(samples) <- ~ longitude + latitude
points(samples)
proj4string(samples) # check coordinate system and/or projection
# If no projection give assign one,
proj4string(samples) <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
samples <- spTransform(samples, CRS("+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points(samples)

## Create plots of samples on top of the bathymetry/TPI
setwd(plots)
bathy_samples <- plot(bathy)
bathy_samples <- points(samples, pch=20, cex=0.75)

tpi_samples <- plot(tpi)
tpi_samples <- points(samples, pch=20, cex=0.75)

########### Histograms of predictors ############
hist(bathy,
     main = "Distribution of Depths",
     xlab = "Depth (meters)", ylab = "Frequency",
     col = "grey")

hist(tpi,
     main = "Distribution of TPI",
     xlab = "TPI", ylab = "Frequency",
     col = "grey")

########### Checking NA in TPI and Bathy ############
NAs_bathy <- as.data.frame(bathy, na.rm=F)
summary(NAs_bathy) #There are 23076584 NAs

NAs_tpi <- as.data.frame(tpi, na.rm=F)
summary(NAs_tpi)

## Visualise NAs
bathy.na <- bathy
bathy.na[is.na(bathy.na)] <- 100
plot(bathy.na)
points(samples, pch=20, cex=0.75)

tpi.na <- tpi
tpi.na[is.na(tpi.na)] <- 100
plot(tpi.na)
points(samples, pch=20, cex=0.75)

########## Distance to boat ramps ###########

library(argosfilter)

# Calculate the distance from boat ramp for each sample----
# Ramps
ids <- factor(c("Exmouth",
                "Bundegi",
                #                 "Onslow",
                #                 "Pilbara",
                #                 "Dampier",
                #                 "WithnellBay",
                #                 "BreadonCreek",
                "Tantabiddi",
                "Coral Bay"
                #                 "Fortescue River"
                #                 "Warroora",
                #                 "Gnaraloo"
)) 
ramps <- data.frame(
  id = rep(ids, each = 1),
  
  y = c(-21.97970240551842,
        -21.84371102054585,
  #                       -21.69503107602818,
  #                       -21.0840122708165,
  #                       -20.66562571673585,
  #                       -20.539744728225,
  #                       -21.648725,
        -21.912580,
        -23.155828)
  #                        -21.028040)
  #                        -23.485963,
  #                        -23.8732820
  ,x = c(114.1461956058358,
         114.1882002543887,
  #                        114.9237761395207,
  #                        115.9306982155436,
  #                        116.6746382564552,
  #                        116.7980700539128,
  #                        115.131243,
         113.978251,
         113.767124))
  #                        116.029232))
#                          113.772874,
#                          113.497430))
head(ramps,4)

samples.ramps <- metadata %>% 
  select("sample", "latitude", "longitude") #this has OpCode,Latitude,Longitude in it

distance.to.ramp<-samples.ramps%>%
  select(sample,latitude,longitude)%>%
  mutate(To.Exmouth=distance(lat1=ramps[1,2],lat2=.$latitude,lon1=ramps[1,3],lon2=.$longitude))%>%
  mutate(To.Bundegi=distance(lat1=ramps[2,2],lat2=.$latitude,lon1=ramps[2,3],lon2=.$longitude))%>%
  #   mutate(To.Onslow=distance(lat1=ramps[3,2],lat2=.$Latitude,lon1=ramps[3,3],lon2=.$Longitude))%>%
  #   mutate(To.Pilbara=distance(lat1=ramps[4,2],lat2=.$Latitude,lon1=ramps[4,3],lon2=.$Longitude))%>%
  #   mutate(To.Dampier=distance(lat1=ramps[5,2],lat2=.$Latitude,lon1=ramps[5,3],lon2=.$Longitude))%>%
  #   mutate(To.WithnellBay=distance(lat1=ramps[6,2],lat2=.$Latitude,lon1=ramps[6,3],lon2=.$Longitude))%>%
  #   mutate(To.BreadonCreek=distance(lat1=ramps[7,2],lat2=.$Latitude,lon1=ramps[7,3],lon2=.$Longitude))%>%
  mutate(To.Tantabiddi=distance(lat1=ramps[3,2],lat2=.$latitude,lon1=ramps[3,3],lon2=.$longitude))%>%
  mutate(To.CoralBay=distance(lat1=ramps[4,2],lat2=.$latitude,lon1=ramps[4,3],lon2=.$longitude))%>%
  #   mutate(To.Fortescue=distance(lat1=ramps[10,2],lat2=.$Latitude,lon1=ramps[10,3],lon2=.$Longitude))%>%
  #   mutate(To.Warroora=distance(lat1=ramps[11,2],lat2=.$Latitude,lon1=ramps[11,3],lon2=.$Longitude))%>%
  #   mutate(To.Gnaraloo=distance(lat1=ramps[12,2],lat2=.$Latitude,lon1=ramps[12,3],lon2=.$Longitude))%>%
  mutate(Distance.to.ramp=do.call(pmin, .[,4:7]))%>%
  select(sample,Distance.to.ramp)%>%
  distinct() #need to be distinct otherwise joins dont work
head(distance.to.ramp)
names((distance.to.ramp))

distance.to.ramp

write.csv(distance.to.ramp, "distance.to.ramp.csv")






