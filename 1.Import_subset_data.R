## Load libraries --
library(raster)
library(dplyr)
library(sp)
library(rgdal)
library(rgeos)

rm(list=ls()) # clear memory 


## Set work directory----
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # sets working directory to where this script is saved (DON't MOVE THE SCRIPT)


## Set sub directories----
d.dir = paste(working.dir,"data/tidy",sep="/") 
s.dir = paste(working.dir,"data/spatial",sep="/") # spatial is where I keep spatial data files, rasters and shapefiles


## Load bathymetry and tpi data----
bathy <- raster(paste(s.dir, "bathy.tif", sep='/'))
plot(bathy) # check it
proj4string(bathy) # check the coordinate system, resolution, etc..

tpi <- raster(paste(s.dir, "tpi.tif", sep='/'))
plot(tpi) # check it
proj4string(tpi) # check the coordinate system, resolution, etc..


## Load BRUV metadata----
bruv <-  read.csv(paste(d.dir, "bruv_data.csv", sep ='/'))
str(bruv) # check df
names(bruv) # check the names of the columns with lat and lon

# Turn bruv data into spatial points --
coordinates(bruv) <- ~ lon + lat
points(bruv) # check it - this should plot the points on to of your previously plotted tpi or bathy
proj4string(bruv) # check coordinate system and/or projection
# If no projection give assign one, for example assign the one from bathy raster --
# Rasters and spatial points should have the sample coordinate system and/or projection
proj4string(bruv) <- proj4string(bathy) 



## Subset the BRUV data ----
# this is to use only a subset of the data as a first approach to run the gams

# plot bathy and bruv points --
plot(bathy)
points(bruv)

# Draw the extent from where you want to subset points
# the extent will be saved as 'e'
e <- drawExtent() # click the top left corner of the extent (yes on the plotted map :) ) and then the bottom right
# you may also just give specific coordinates using the 'extent' function if prefered

# crop the bruv data to the area that the extent --
sub.bruv <- crop(bruv, e)
# you may do the same with the raster data
sub.bathy <- crop(bathy, e)

# plot to check
plot(sub.bathy)
points(sub.bruv)

# you may repeat this process to subset other areas of your domain



## Extract perdictor data at bruv locations ----

# Stack predictors --
preds <- raster::stack(bathy, tpi) # if this fails: preds may have different extent or resolution

# Extract predictors info --
fishpreds <- raster::extract(preds, sub.bruv, df = T) # this will result in dataframe with prub points and predictor values

# Save this data frame
write.csv(fishpreds, paste(d.dir, "fishdata.csv", sep ='/')) 

# Save stack of predictors if desired
writeRaster(preds, paste(s.dir, "predstack.tif", sep ='/'))
