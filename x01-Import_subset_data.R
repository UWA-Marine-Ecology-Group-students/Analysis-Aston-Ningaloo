## Load libraries --
library(raster)
library(dplyr)
library(sp)
library(rgdal)
library(rgeos)
library(spatialEco)

rm(list=ls()) # clear memory 

lat.long <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
utm <- "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

## Set work directory----
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # sets working directory to where this script is saved (DON't MOVE THE SCRIPT)


## Set sub directories----
d.dir = paste(working.dir,"Tidy data",sep="/") 
s.dir = paste(working.dir,"Spatial",sep="/") # spatial is where I keep spatial data files, rasters and shapefiles


## Load bathymetry and tpi data----
setwd(s.dir)
bathy <- raster(paste(s.dir, "bathy.tif", sep='/'))
bathy <- flip(bathy, direction="y")
plot(bathy) # check it
proj4string(bathy) # check the coordinate system, resolution, etc..


tpi <- raster(paste(s.dir, "tpi2.tif", sep='/'))
plot(tpi) # check it
proj4string(tpi) # check the coordinate system, resolution, etc..
tpi <- flip(tpi, direction="y")
plot(tpi)

## Load BRUV data
bruv <-  read.csv(paste(d.dir, "ningaloo.complete.maxn.csv", sep ='/'))
str(bruv) # check df
names(bruv) # check the names of the columns with lat and lon
# Extract bits we need
bruv <- bruv %>%
  select("sample", "scientific", "maxn", "longitude", "latitude")
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
sub.bruv <- crop(sub.bruv, e)
# you may do the same with the raster data
sub.bathy <- crop(sub.bathy, e)
sub.tpi <- crop(tpi, e)

# plot to check
plot(sub.bathy)
points(sub.bruv)

# you may repeat this process to subset other areas of your domain



## Extract perdictor data at bruv locations ----

# Stack predictors --

preds <- raster::stack(sub.bathy, sub.tpi) # if this fails: preds may have different extent or resolution


## Extract predictors info --
fishpreds <- raster::extract(preds, sub.bruv, df = T) # this will result in dataframe with prub points and predictor values

# Give sub maxn uniques IDs
sub.bruv$ID <- fishpreds$ID

# Fix column names
names(fishpreds) <- c("ID", "bathymetry", "TPI")

# Create full dataframe
fishdata <- merge.data.frame(fishpreds, sub.bruv, by="ID")

# Save this data frame
write.csv(fishdata, paste(d.dir, "fishdata.csv", sep ='/')) 

# Save stack of predictors if desired
writeRaster(preds, paste(s.dir, "predstack.tif", sep ='/'))

