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

dat<-read.csv("fishdata.csv")%>% # this is the df with fish and predictor data
  left_join(metadata, by = "sample")%>%
  # dplyr::rename(site=sample)%>% # rename columns in data for model - BG You need site not sample 23/04/20
  # dplyr::filter(Taxa!='Synodontidae Saurida undosquamis')%>% remove Synodontidae Saurida undosquamis
  na.omit()%>% # remove NAs
  dplyr::filter(!sample %in% c('10.13', '10.14', '10.15'))%>% # Remove preferential sites 
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
hist(Multidens$maxn)
plot(Multidens$maxn)



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

Sebae.gam.2 <- gam(maxn ~ s(TPI, k = 6, bs = "cr") + s(bathymetry, k = 6, bs = "cr") + # covariate effect of tpi and bathymetry
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

Multidens.gam.2 <- gam(maxn ~ s(TPI, k = 6, bs = "cr") + s(bathymetry, k = 6, bs = "cr") +  # covariate effect of tpi and bathymetry
                       + s(site, k = 3, bs ='re'), # random effect of cluster or site
                     family=tw(),data=Multidens) # family tweedy
plot(Multidens.gam,residuals=T,all.terms = TRUE,pages=1)
summary(Multidens.gam.2) # check results
AICc(Multidens.gam) # check AIC of model

## Check residuals
par(mfrow=c(2,2))
gam.check(Multidens.gam)

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
prediction.sebae.2<-predict(Sebae.gam.2, predictm, type = 'response', se.fit=T, index=1:2, progress='text', exclude="s(site)") 

## Store prediction in a df ----
prediction.sebae.2_df<-as.data.frame(prediction.sebae.2,xy=TRUE,na.rm=TRUE)%>%
  glimpse()
# add cooridinates from preds to df
prediction.sebae.2_df$Longitude<-predictm$longitude
prediction.sebae.2_df$Latitude<-predictm$latitude
glimpse(prediction.sebae.2_df) # check

# save this df
write.csv(prediction.sebae.2_df, paste(d.dir, "Sebae.prediction.2.csv", sep='/'))


## Create raster for mean fit and one for se fit ----

# Mean fit --
sebae.fit.2<-prediction.sebae.2_df%>%
  dplyr::select(Longitude,Latitude,fit)%>%
  glimpse()

# SE fit --
sebae.se.2<-prediction.sebae.2_df%>%
  dplyr::select(Longitude,Latitude,se.fit)%>%
  glimpse()


## Convert into a spatialPoints dataframe ----
# Mean fit--
coordinates(sebae.fit.2) <- ~ Longitude + Latitude
# SE fit--
coordinates(sebae.se.2) <- ~ Longitude + Latitude


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
sebae.fit.2 <- raster(sebae.fit.2)
crs(sebae.fit.2)<-sr
sebae.fit.2
plot(sebae.fit.2, col=my.palette)
# plot SE fit --
sebae.se.2 <- raster(sebae.se.2)
crs(sebae.se.2) <-sr
sebae.se.2
plot(sebae.se.2, col=my.palette)

# Save rasters--
setwd(model.out)
writeRaster(sebae.fit.2,"Sebae_fit_2.tif", sep ='/')
writeRaster(sebae.se.2, "Sebae_se_2.tif", sep ='/')

## Plot with ggplot ----

### Overall plot

## Fit

sebaep.2<-ggplot()+
  geom_tile(data=prediction.sebae.2_df,aes(x=Longitude,y=Latitude,fill=fit),alpha=0.8)+
  scale_fill_viridis(option = "magma",direction = -1)+
  #geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+ # can add polygon of coastline here
  scale_color_gradient()+
  coord_equal()+
  #xlim(112,155)+ #  set limits of plot if desired
  #ylim(-45,-7)+ #  set limits of plot
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))

sebaep.2

######## Predict using the fitted model Pristipomoides Multidens ########
prediction.multidens.2<-predict(Multidens.gam.2, predictm, type = 'response', se.fit=T, index=1:2, progress='text', exclude="s(site)") 

## Store prediction in a df ----
prediction.multidens.2_df<-as.data.frame(prediction.multidens.2,xy=TRUE,na.rm=TRUE)%>%
  glimpse()
# add cooridinates from preds to df
prediction.multidens.2_df$Longitude<-predictm$longitude
prediction.multidens.2_df$Latitude<-predictm$latitude
glimpse(prediction.multidens.2_df) # check

# save this df
write.csv(prediction.multidens.2_df, paste(d.dir, "Multidens.prediction.2.csv", sep='/'))


## Create raster for mean fit and one for se fit ----

# Mean fit --
multidens.fit.2<-prediction.multidens.2_df%>%
  dplyr::select(Longitude,Latitude,fit)%>%
  glimpse()

# SE fit --
multidens.se.2<-prediction.multidens.2_df%>%
  dplyr::select(Longitude,Latitude,se.fit)%>%
  glimpse()


## Convert into a spatialPoints dataframe ----
# Mean fit--
coordinates(multidens.fit.2) <- ~ Longitude + Latitude
# SE fit--
coordinates(multidens.se.2) <- ~ Longitude + Latitude


## coerce to SpatialPixelsDataFrame ----
# Mean fit
gridded(multidens.fit.2) <- TRUE
# SE fit
gridded(multidens.se.2) <- TRUE


## Coerce to raster ----

# Set the CRS --
sr <- "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# choose colors --
library(RColorBrewer)
my.palette <- brewer.pal(n = 9, name = "OrRd")
pal <- colorRampPalette(c("yellow","pink", "red", "dark blue"))
pal <- colorRampPalette(c("yellow","pink", "red", "dark blue"))
# plot mean fit --
multidens.fit.2 <- raster(multidens.fit.2)
crs(multidens.fit.2)<-sr
multidens.fit.2
plot(multidens.fit.2, col=my.palette)
# plot SE fit --
multidens.se.2 <- raster(multidens.se.2)
crs(multidens.se.2) <-sr
multidens.se.2
plot(multidens.se.2, col=my.palette)

# Save rasters--
setwd(model.out)
writeRaster(multidens.fit.2, "Multidens_fit.2.tif", sep ='/')
writeRaster(multidens.se, "Multidens_se.2.tif", sep ='/')


## Plot with ggplot ----

### Overall plot

## Fit

multidensp.2<-ggplot()+
  geom_tile(data=prediction.multidens.2_df,aes(x=Longitude,y=Latitude,fill=fit),alpha=0.8)+
  scale_fill_viridis(option = "magma",direction = -1)+
  #geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+ # can add polygon of coastline here
  scale_color_gradient()+
  coord_equal()+
  #xlim(112,155)+ #  set limits of plot if desired
  #ylim(-45,-7)+ #  set limits of plot
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))

multidensp.2
