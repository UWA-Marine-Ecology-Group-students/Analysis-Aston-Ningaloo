# This script is a continuation of script 2, it uses the gam results from script 2 --

##### Load libraries if not loaded

library(sp)
library(rgdal)
library(rgeos)
library(maptools)
library(dplyr)
library(broom)
library(data.table)
library(ggplot2)
library(raster)
library(mgcv)
library(MuMIn)
library(viridis)
library(ggplot2)
library(classInt)


# Set work directory ----

working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # sets working directory to where this script is saved (DON't MOVE THE SCRIPT)

# Set sub directories----

data.dir <- paste(working.dir,"Tidy Data",sep="/")
spatial.dir <- paste(working.dir,"Spatial",sep="/") # wherever you keep your spatial data
plots.dir <- paste(working.dir,"Plots",sep="/")
model.out <- paste(working.dir,"Model Out",sep="/")


## Load predictor data ----

preds <- stack(paste(spatial.dir, "predstack.tif", sep= '/'))
plot(preds) # to visualize them

names(preds) # check the names of your preds
# rename this if needed:
names(preds)<-c('bathymetry','TPI')


## To predict in space make preds df----
predictm<-as.data.frame(preds,xy=T,na.rm=TRUE)%>%
  glimpse()
names(predictm) <- c("longitude", "latitude", 'bathymetry', 'TPI')

predictm$site <- "a"


######## Predict using the fitted model Lutjanus sebae #########
prediction.sebae<-predict(Sebae.gam, predictm, type = 'response', se.fit=T, index=1:2, progress='text', exclude="s(site)") 

## Store prediction in a df ----
prediction.sebae_df<-as.data.frame(prediction.sebae,xy=TRUE,na.rm=TRUE)%>%
  glimpse()
# add cooridinates from preds to df
prediction.sebae_df$Longitude<-predictm$longitude
prediction.sebae_df$Latitude<-predictm$latitude
glimpse(prediction.sebae_df) # check

# save this df
write.csv(prediction.sebae_df, paste(d.dir, "Sebae.prediction.csv", sep='/'))


## Create raster for mean fit and one for se fit ----

# Mean fit --
sebae.fit<-prediction.sebae_df%>%
  dplyr::select(Longitude,Latitude,fit)%>%
  glimpse()

# SE fit --
sebae.se<-prediction.sebae_df%>%
  dplyr::select(Longitude,Latitude,se.fit)%>%
  glimpse()


## Convert into a spatialPoints dataframe ----
# Mean fit--
coordinates(sebae.fit) <- ~ Longitude + Latitude
# SE fit--
coordinates(sebae.se) <- ~ Longitude + Latitude


## coerce to SpatialPixelsDataFrame ----
# Mean fit
gridded(sebae.fit) <- TRUE
# SE fit
gridded(sebae.se) <- TRUE


## Coerce to raster ----

# Set the CRS --
sr <- "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# choose colors --
library(RColorBrewer)
my.palette <- brewer.pal(n = 9, name = "OrRd")
pal <- colorRampPalette(c("yellow","pink", "red", "dark blue"))
pal <- colorRampPalette(c("yellow","pink", "red", "dark blue"))
# plot mean fit --
sebae.fit <- raster(sebae.fit)
crs(sebae.fit)<-sr
sebae.fit
plot(sebae.fit, col=my.palette)
# plot SE fit --
sebae.se <- raster(sebae.se)
crs(sebae.se) <-sr
sebae.se
plot(sebae.se, col=my.palette)

# Save rasters--
writeRaster(sebae.fit, paste(model.out, "Sebae_fit.tif", sep ='/'))
writeRaster(sebae.se, paste(model.out, "Sebae_se.tif", sep ='/'))


## Plot with ggplot ----

### Overall plot

## Fit

sebaep<-ggplot()+
  geom_tile(data=prediction.sebae_df,aes(x=Longitude,y=Latitude,fill=fit),alpha=0.8)+
  scale_fill_viridis(option = "magma",direction = -1)+
  #geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+ # can add polygon of coastline here
  scale_color_gradient()+
  coord_equal()+
  #xlim(112,155)+ #  set limits of plot if desired
  #ylim(-45,-7)+ #  set limits of plot
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))

sebaep

######## Predict using the fitted model Pristipomoides Multidens ########
prediction.multidens<-predict(Multidens.gam, predictm, type = 'response', se.fit=T, index=1:2, progress='text', exclude="s(site)") 

## Store prediction in a df ----
prediction.multidens_df<-as.data.frame(prediction.multidens,xy=TRUE,na.rm=TRUE)%>%
  glimpse()
# add cooridinates from preds to df
prediction.multidens_df$Longitude<-predictm$longitude
prediction.multidens_df$Latitude<-predictm$latitude
glimpse(prediction.multidens_df) # check

# save this df
write.csv(prediction.multidens_df, paste(d.dir, "Multidens.prediction.csv", sep='/'))


## Create raster for mean fit and one for se fit ----

# Mean fit --
multidens.fit<-prediction.multidens_df%>%
  dplyr::select(Longitude,Latitude,fit)%>%
  glimpse()

# SE fit --
multidens.se<-prediction.multidens_df%>%
  dplyr::select(Longitude,Latitude,se.fit)%>%
  glimpse()


## Convert into a spatialPoints dataframe ----
# Mean fit--
coordinates(multidens.fit) <- ~ Longitude + Latitude
# SE fit--
coordinates(multidens.se) <- ~ Longitude + Latitude


## coerce to SpatialPixelsDataFrame ----
# Mean fit
gridded(multidens.fit) <- TRUE
# SE fit
gridded(multidens.se) <- TRUE


## Coerce to raster ----

# Set the CRS --
sr <- "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# choose colors --
library(RColorBrewer)
my.palette <- brewer.pal(n = 9, name = "OrRd")
pal <- colorRampPalette(c("yellow","pink", "red", "dark blue"))
pal <- colorRampPalette(c("yellow","pink", "red", "dark blue"))
# plot mean fit --
multidens.fit <- raster(multidens.fit)
crs(multidens.fit)<-sr
multidens.fit
plot(multidens.fit, col=my.palette)
# plot SE fit --
multidens.se <- raster(multidens.se)
crs(multidens.se) <-sr
multidens.se
plot(multidens.se, col=my.palette)

# Save rasters--
writeRaster(multidens.fit, paste(model.out, "Multidens_fit.tif", sep ='/'))
writeRaster(multidens.se, paste(model.out, "Multidens_se.tif", sep ='/'))


## Plot with ggplot ----

### Overall plot

## Fit

multidensp<-ggplot()+
  geom_tile(data=prediction.multidens_df,aes(x=Longitude,y=Latitude,fill=fit),alpha=0.8)+
  scale_fill_viridis(option = "magma",direction = -1)+
  #geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+ # can add polygon of coastline here
  scale_color_gradient()+
  coord_equal()+
  #xlim(112,155)+ #  set limits of plot if desired
  #ylim(-45,-7)+ #  set limits of plot
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))

multidensp



## Calculate population estimates in your domain (fish/squared meter) ----
## Check the resolution of your data and how many cells in your raster file
yourraster # your prediction raster - check the resolution I think is 25x25m = 625m^2 /cell
# get the number of cells also: I think 383292 ni yours

# get the number of cells in your raster that are not NA
rNA <- sum(!is.na(youraster)) # 199757
notNA <- 383292 - 199757
# calculate no of 1m^2 in the actual area modeled 
domainmeters <- notNA*625 # to standarize by 1m^2

# get the total fish abundance in the domain --
suma <- cellStats(m, sum) # 62988154 this is you estmate of total abundance for the total area of domain

# get fish per meter squared -- you can calculate these estimates to other area (different from 1 m^2) 
fpm <- suma/domainmeters # 0.5491108 

# repeat for SE 




