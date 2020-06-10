############ Import and subset data ##############
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
proj4string(tpi) # check the coordinate system, resolution, etc..
tpi <- flip(tpi, direction="y")
plot(tpi)

## Load Length data
lengths <-  read.csv(paste(d.dir, "ningaloo.complete.length.csv", sep ='/'))%>%
  mutate(scientific=paste(family,genus,species,sep=" "))
str(lengths) # check df
names(lengths) # check the names of the columns with lat and lon
# Extract bits we need
lengths <- lengths %>%
  select("sample","scientific", "length", "longitude", "latitude")
names(lengths)


# Turn length data into spatial points --
coordinates(lengths) <- ~ longitude + latitude
points(lengths) # check it - this should plot the points on top of your previously plotted tpi or bathy
proj4string(lengths) # check coordinate system and/or projection
# If no projection give assign one,
proj4string(lengths) <- lat.long
# Convert to the same projection as the bathy/tpi
# Rasters and spatial points should have the sample coordinate system and/or projection
lengths <- spTransform(lengths, CRS("+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points(lengths)

## Subset the BRUV data ----
# this is to use only a subset of the data as a first approach to run the gams

# plot bathy and bruv points --
plot(bathy)
points(lengths)

# Draw the extent from where you want to subset points
# the extent will be saved as 'e'
e <- drawExtent() # click the top left corner of the extent (yes on the plotted map :) ) and then the bottom right
# you may also just give specific coordinates using the 'extent' function if prefered

# crop the bruv data to the area that the extent --
sub.lengths <- crop(lengths, e)
# you may do the same with the raster data
sub.bathy <- crop(bathy, e)
sub.tpi <- crop(tpi, e)

# plot to check
plot(sub.bathy)
points(sub.lengths)

# you may repeat this process to subset other areas of your domain

## Extract perdictor data at bruv locations ----

# Stack predictors --

preds <- raster::stack(sub.bathy, sub.tpi) # if this fails: preds may have different extent or resolution


## Extract predictors info --
fishpreds <- raster::extract(preds, sub.lengths, df = T) # this will result in dataframe with prub points and predictor values

# Give sub maxn uniques IDs
sub.lengts$ID <- fishpreds$ID

# Fix column names
names(fishpreds) <- c("ID", "bathymetry", "TPI")

# Create full dataframe
fishlengths <- merge.data.frame(fishpreds, sub.lengths, by="ID")

# Save this data frame
write.csv(fishlengths, paste(d.dir, "fishlengthscsv", sep ='/')) 

# Save stack of predictors if desired
writeRaster(preds, paste(s.dir, "predstack.tif", sep ='/'))

############ Model Data ##############
##### Load libraries

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


# Set work directory ----
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # sets working directory to where this script is saved (DON't MOVE THE SCRIPT)

# Set sub directories----

data.dir <- paste(working.dir,"Tidy data",sep="/")
plots.dir <- paste(working.dir,"Plots",sep="/")
model.out <- paste(working.dir,"Model Out",sep="/")

# Read in data----
# data.dir<- ("C:/Users/00093391/Dropbox/UWA/Research Associate/Woodside-Exmouth/fishdata")
setwd(data.dir)
dir()

metadata<-read.csv("MEG_Labsheets_2020 - 2019-08_Ningaloo-Deep_stereo-BR.csv")%>%
  dplyr::rename(sample=Sample,site=Site)%>%
  dplyr::select(sample,site)

dat<-read.csv("fishlengths.csv")%>% # this is the df with fish and predictor data
  left_join(metadata, by = "sample")%>%
  # dplyr::rename(site=sample)%>% # rename columns in data for model - BG You need site not sample 23/04/20
  # dplyr::filter(Taxa!='Synodontidae Saurida undosquamis')%>% remove Synodontidae Saurida undosquamis
  na.omit()%>% # remove NAs
  # dplyr::filter(!sample %in% c('10.13', '10.14', '10.15'))%>% Remove preferential sites 
  glimpse() # to see data
# Convert covariates that are stored as character (factor) or integers (continous) if needed

#### Check distribution of the response for species ----

plot.new()
par(mfrow=c(1,2))

levels(dat$scientific)

## Subsample data for each species needed for the analysis
# Repeat for as many species as required --


## Species 1: Lutjanus sebae --

levels(dat$scientific)

Sebae<-dat%>%
  filter(scientific=="Lutjanidae Lutjanus sebae")%>%
  glimpse()
hist(Sebae$TPI)
plot(Sebae$bathymetry)

## Species 2: Pristipomoides multidens --

levels(dat$species)

Multidens<-dat%>%
  filter(scientific=="Lutjanidae Pristipomoides multidens")%>%
  glimpse()
hist(Multidens$length)
plot(Multidens$length)



## Convert fish data into spatial object ----
coordinates(dat)<-~longitude+latitude # check the name of the columns for lon and lat
class(dat) # should be spatial object


## Set the reference system to the widely used WGS84
proj4string(dat)<-CRS("+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
dat@proj4string


###########################


##  Fit a top model ----
# First I will do marginal predictions (no random effects) and not accounting for the offset 
# (generates rates - really low numbers)
#library(mgcv)
#library(MuMIn)

## Fit model for Lutjanus sebae ---- 

names(dat)

Sebae.gam.length <- gam(length ~ s(TPI, k = 6, bs = "cr") + s(bathymetry, k = 6, bs = "cr") + # covariate effect of tpi and bathymetry
                     + s(site, k = 3, bs ='re'), # random effect of cluster or site
                   family=tw(),data=Sebae) # family tweedy
plot(Sebae.gam.2,residuals=T,all.terms = TRUE,pages=1)
summary(Sebae.gam.2) # check results
AICc(Sebae.gam.2) # check AIC of model

## Check residuals
par(mfrow=c(2,2))
gam.check(Sebae.gam)

## Fit model for Pristipomoides multidens ----

names(dat)

Multidens.gam.length <- gam(length ~ s(TPI, k = 6, bs = "cr") + s(bathymetry, k = 6, bs = "cr") +  # covariate effect of tpi and bathymetry
                         + s(site, k = 3, bs ='re'), # random effect of cluster or site
                       family=tw(),data=Multidens) # family tweedy
plot(Multidens.gam,residuals=T,all.terms = TRUE,pages=1)
summary(Multidens.gam.length) # check results
AICc(Multidens.gam.length) # check AIC of model

## Check residuals
par(mfrow=c(2,2))
gam.check(Multidens.gam.length)

############ Predict values ##############

# Load libraries if not loaded

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
glimpse(predictm)

######## Predict using the fitted model Lutjanus sebae #########
prediction.sebae.length<-predict(Sebae.gam.length, predictm, type = 'response', se.fit=T, index=1:2, progress='text', exclude="s(site)") 

## Store prediction in a df ----
prediction.sebae.length_df<-as.data.frame(prediction.sebae.length,xy=TRUE,na.rm=TRUE)%>%
  glimpse()
# add cooridinates from preds to df
prediction.sebae.length_df$Longitude<-predictm$longitude
prediction.sebae.length_df$Latitude<-predictm$latitude
glimpse(prediction.sebae.length_df) # check

# save this df
write.csv(prediction.sebae.length_df, paste(d.dir, "Sebae.prediction.length.csv", sep='/'))


## Create raster for mean fit and one for se fit ----

# Mean fit --
sebae.fit.length<-prediction.sebae.length_df%>%
  dplyr::select(Longitude,Latitude,fit)%>%
  glimpse()

# SE fit --
sebae.se.length<-prediction.sebae.length_df%>%
  dplyr::select(Longitude,Latitude,se.fit)%>%
  glimpse()


## Convert into a spatialPoints dataframe ----
# Mean fit--
coordinates(sebae.fit.length) <- ~ Longitude + Latitude
# SE fit--
coordinates(sebae.se.length) <- ~ Longitude + Latitude


## coerce to SpatialPixelsDataFrame ----
# Mean fit
gridded(sebae.fit.2) <- TRUE
# SE fit
gridded(sebae.se.2) <- TRUE


## Coerce to raster ----

# Set the CRS --
sr <- "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# choose colors --
library(RColorBrewer)
my.palette <- brewer.pal(n = 9, name = "OrRd")
pal <- colorRampPalette(c("yellow","pink", "red", "dark blue"))
pal <- colorRampPalette(c("yellow","pink", "red", "dark blue"))
# plot mean fit --
sebae.fit.length <- raster(sebae.fit.length)
crs(sebae.fit.length)<-sr
sebae.fit.length
plot(sebae.fit.length, col=my.palette)
# plot SE fit --
sebae.se.2 <- raster(sebae.se.length)
crs(sebae.se.length) <-sr
sebae.se.length
plot(sebae.se.length, col=my.palette)

# Save rasters--
setwd(model.out)
writeRaster(sebae.fit.length,"Sebae_fit_length.tif", sep ='/')
writeRaster(sebae.se.length, "Sebae_se_length.tif", sep ='/')

## Plot with ggplot ----

### Overall plot

## Fit

sebaelengthp<-ggplot()+
  geom_tile(data=prediction.sebae.length_df,aes(x=Longitude,y=Latitude,fill=fit),alpha=0.8)+
  scale_fill_viridis(option = "magma",direction = -1)+
  #geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+ # can add polygon of coastline here
  scale_color_gradient()+
  coord_equal()+
  #xlim(112,155)+ #  set limits of plot if desired
  #ylim(-45,-7)+ #  set limits of plot
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))

sebaelengthp

######## Predict using the fitted model Pristipomoides Multidens ########
prediction.multidens.length<-predict(Multidens.gam.length, predictm, type = 'response', se.fit=T, index=1:2, progress='text', exclude="s(site)") 

## Store prediction in a df ----
prediction.multidens.length_df<-as.data.frame(prediction.multidens.length,xy=TRUE,na.rm=TRUE)%>%
  glimpse()
# add cooridinates from preds to df
prediction.multidens.length_df$Longitude<-predictm$longitude
prediction.multidens.length_df$Latitude<-predictm$latitude
glimpse(prediction.multidens.length_df) # check

# save this df
write.csv(prediction.multidens.length_df, paste(d.dir, "Multidens.prediction.length.csv", sep='/'))


## Create raster for mean fit and one for se fit ----

# Mean fit --
multidens.fit.length<-prediction.multidens.length_df%>%
  dplyr::select(Longitude,Latitude,fit)%>%
  glimpse()

# SE fit --
multidens.se.length<-prediction.multidens.length_df%>%
  dplyr::select(Longitude,Latitude,se.fit)%>%
  glimpse()


## Convert into a spatialPoints dataframe ----
# Mean fit--
coordinates(multidens.fit.length) <- ~ Longitude + Latitude
# SE fit--
coordinates(multidens.se.length) <- ~ Longitude + Latitude


## coerce to SpatialPixelsDataFrame ----
# Mean fit
gridded(multidens.fit.length) <- TRUE
# SE fit
gridded(multidens.se.length) <- TRUE


## Coerce to raster ----

# Set the CRS --
sr <- "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# choose colors --
library(RColorBrewer)
my.palette <- brewer.pal(n = 9, name = "OrRd")
pal <- colorRampPalette(c("yellow","pink", "red", "dark blue"))
pal <- colorRampPalette(c("yellow","pink", "red", "dark blue"))
# plot mean fit --
multidens.fit.length <- raster(multidens.fit.length)
crs(multidens.fit.length)<-sr
multidens.fit.length
plot(multidens.fit.length, col=my.palette)
# plot SE fit --
multidens.se.length <- raster(multidens.se.length)
crs(multidens.se.length) <-sr
multidens.se.length
plot(multidens.se.length, col=my.palette)

# Save rasters--
setwd(model.out)
writeRaster(multidens.fit.length, "Multidens_fit.length.tif", sep ='/')
writeRaster(multidens.se.length, "Multidens_se.length.tif", sep ='/')


## Plot with ggplot ----

### Overall plot

## Fit

multidenslengthp<-ggplot()+
  geom_tile(data=prediction.multidens.length_df,aes(x=Longitude,y=Latitude,fill=fit),alpha=0.8)+
  scale_fill_viridis(option = "magma",direction = -1)+
  #geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+ # can add polygon of coastline here
  scale_color_gradient()+
  coord_equal()+
  #xlim(112,155)+ #  set limits of plot if desired
  #ylim(-45,-7)+ #  set limits of plot
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))

multidenslengthp
