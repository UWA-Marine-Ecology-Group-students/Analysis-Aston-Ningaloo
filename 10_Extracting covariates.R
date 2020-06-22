library(raster)
library(rstudioapi)
library(tidyverse)
library(sf)
library(dplyr)
library(sp)
library(rgdal)
library(rgeos)
library(spatialEco)

##Set working directory----
## Set work directory----
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # sets working directory to where this script is saved (DON't MOVE THE SCRIPT)

## Set sub directories----
d.dir = paste(working.dir,"Tidy data",sep="/") 
s.dir = paste(working.dir,"Spatial",sep="/") # spatial is where I keep spatial data files, rasters and shapefiles
p.dir <- paste(working.dir,"Plots",sep="/")

lat.long <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
utm <- "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"


## Load bathymetry and tpi data----
setwd(s.dir)
bathy <- raster(paste(s.dir, "bathy.tif", sep='/'))
bathy <- flip(bathy, direction="y")
plot(bathy) # check it
proj4string(bathy) # check the coordinate system, resolution, etc..


tpi <- raster(paste(s.dir, "tpi.tif", sep='/'))
plot(tpi) # check it
proj4string(tpi) # check the coordinate system, resolution, etc..
tpi <- flip(tpi, direction="y")
plot(tpi)

## Load BRUV data
bruv <-  read.csv(paste(d.dir, "ningaloo_metadata.csv", sep ='/'))
str(bruv) # check df
names(bruv) # check the names of the columns with lat and lon
# Extract bits we need
bruv <- bruv %>%
  select("sample", "longitude", "latitude")
bruv


# Turn bruv data into spatial points --
coordinates(bruv) <- ~ longitude + latitude
points(bruv) # check it - this should plot the points on top of your previously plotted tpi or bathy
proj4string(bruv) # check coordinate system and/or projection
# If no projection give assign one,
proj4string(bruv) <- lat.long
# Convert to the same projection as the bathy/tpi
# Rasters and spatial points should have the sample coordinate system and/or projection
bruv <- spTransform(bruv, CRS("+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points(bruv)

## Create all covariates wanted - already have bathymetry and TPI

Slope <- raster::terrain(bathy, opt="slope")
Aspect <- raster::terrain(bathy, opt="aspect")
TRI <- raster::terrain(bathy, opt="tri")
Roughness <- raster::terrain(bathy, opt="roughness")
FlowDir <- raster::terrain(bathy,opt="flowdir")

## Extract perdictor data at bruv locations ----

# Stack predictors --
preds <- raster::stack(bathy,tpi,Slope,Aspect,TRI,Roughness,FlowDir) # if this fails: preds may have different extent or resolution

# Extract predictors info --
covariates <- raster::extract(preds, bruv, df = T) # this will result in dataframe with bruv points and predictor values

# Give bruv uniques IDs
bruv$ID <- covariates$ID

# Fix column names
names(covariates) <- c("ID", "bathymetry", "TPI", "Slope", "Aspect","TRI", "Roughness", "FlowDir")

# Create full dataframe
fullcovariates <- merge.data.frame(covariates, bruv, by="ID")

# Save this data frame
write.csv(fullcovariates, paste(d.dir, "covariates.csv", sep ='/')) 

# Save stack of predictors if desired
writeRaster(preds, paste(s.dir, "covariatestack.tif", sep ='/'))

###### Finding NA values #######

NA.bathy <- subset(fullcovariates,is.na(fullcovariates$bathymetry))
NA.tpi <- subset(fullcovariates,is.na(fullcovariates$TPI))
NA.slope <- subset(fullcovariates,is.na(fullcovariates$Slope))
NA.aspect <- subset(fullcovariates,is.na(fullcovariates$Aspect))
NA.tri <- subset(fullcovariates,is.na(fullcovariates$TRI))
NA.roughness <- subset(fullcovariates,is.na(fullcovariates$Roughness))
NA.flowdir <- subset(fullcovariates,is.na(fullcovariates$FlowDir))

unique(NA.bathy$sample) #10.12
unique(NA.tpi$sample) # 8.05 10.09 10.12 16.03
unique(NA.slope$sample) #10.12 16.03
unique(NA.aspect$sample) #10.12 16.03
unique(NA.tri$sample) #10.12 16.03
unique(NA.roughness$sample) #10.12 16.03
unique(NA.flowdir$sample) #10.12

## Plots of NA sites 

#Bathymetry
plot(bathy)
points(bruv, pch=20, cex=0.75, col=ifelse(bruv$sample=='10.12', "red", "black"))

#TPI
plot(tpi)
points(bruv, pch=20, cex=0.75, col=ifelse(bruv$sample %in% c("8.05", "10.09", "10.12", "16.03"), "red", "black"))

#Slope
plot(Slope)
points(bruv, pch=20, cex-0.75, col=ifelse(bruv$sample %in% c("10.12", "16.03"), "red", "black"))

#Aspect
plot(Aspect)
points(bruv, pch=20, cex-0.75, col=ifelse(bruv$sample %in% c("10.12", "16.03"), "red", "black"))

#TRI
plot(TRI)
points(bruv, pch=20, cex-0.75, col=ifelse(bruv$sample %in% c("10.12", "16.03"), "red", "black"))

#Roughness
plot(Roughness)
points(bruv, pch=20, cex-0.75, col=ifelse(bruv$sample %in% c("10.12", "16.03"), "red", "black"))

#FlowDir
plot(FlowDir)
points(bruv, pch=20, cex=0.75, col=ifelse(bruv$sample=='10.12', "red", "black"))

############# Formating habitat data to use as covariate #############

setwd(d.dir)

habitat <- read.csv("ningaloo.complete.habitat.csv")

## Want corals and sponges to stay as coral/sponge when coverage is greater than 50% otherwise we want it to
## be classified as just reef but if unconsolidated sediment is >50% we want it to be classified as unconsolidated

names(habitat)
habitat.classified <- habitat %>%
  rowwise()%>%
  mutate(reef=sum(c(broad.bryozoa,broad.crinoids,broad.hydrocoral,broad.hydroids,broad.octocoral.black,broad.sponges)))%>%
  mutate(sand=broad.unconsolidated)

habitat.classified <- habitat.classified%>%
  select("sample", "mean.relief", "sd.relief", "reef", "sand")

write.csv(habitat.classified, "habitat.classified.csv")

# Add to existing covariates
covariates<-read.csv("covariates.csv")

covariates<-full_join(covariates,habitat.classified,by="sample")

############# Adding distance to boat ramp to covariates ##############
ramps<-read.csv("distance.to.ramp.csv")

covariates<-full_join(covariates,ramps,by="sample")
  
covariates<-covariates%>%
  select(!c("X.x", "X.y"))%>%
  dplyr::rename(distance.to.ramp=Distance.to.ramp)

covariates<-covariates[,c(1,9,2,3,4,5,6,7,8,10,11,12,13,14,15,16)]

write.csv(covariates, "covariates.csv")








  



         
         
        