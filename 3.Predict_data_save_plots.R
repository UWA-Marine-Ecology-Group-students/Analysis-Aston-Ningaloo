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

work.dir <- ("C:/Users/00093391/Dropbox/UWA/data.analyses/Abundances") # path to your directory

# Set sub directories----

data.dir <- paste(work.dir,"data",sep="/")
spatial.dir <- paste(work.dir,"spatial",sep="/") # wherever you keep your spatial data
plots.dir <- paste(work.dir,"plots",sep="/")
model.out <- paste(work.dir,"ModelOut",sep="/")


## Load predictor data ----

preds <- stack(paste(spatial.dir, "predstack.tif", sep= '/'))
plot(preds) # to visualize them

names(preds) # check the names of your preds
# rename this if needed:
#names(pred)<-c('Bathy','tpi')


## To predict in space make preds df----
predictm<-as.data.frame(pred,xy=T,na.rm=TRUE)%>%
  glimpse()

## Predict using the fitted model ----
prediction.nemip<-predict(Nemip, predictm, type = 'response', se.fit=T, index=1:2, progress='text') 


## Store prediction in a df ----
prediction.nemid_df<-as.data.frame(prediction.nemip,xy=TRUE,na.rm=TRUE)%>%
  glimpse()
# add cooridinates from preds to df
prediction.nemid_df$Longitude<-predictm$x 
prediction.nemid_df$Latitude<-predictm$y
glimpse(prediction.nemid_df) # check

# save this df
write.csv(prediction.nemid_df, paste(d.dir, "Nemid.predictions.csv", sep='/'))


## Create raster for mean fit and one for se fit ----

# Mean fit --
nemip.fit<-prediction.nemid_df%>%
  dplyr::select(Longitude,Latitude,fit)%>%
  glimpse()

# SE fit --
nemip.se<-prediction.nemid_df%>%
  dplyr::select(Longitude,Latitude,se.fit)%>%
  glimpse()


## Convert into a spatialPoints dataframe ----
# Mean fit--
coordinates(nemip.fit) <- ~ Longitude + Latitude
# SE fit--
coordinates(nemip.se) <- ~ Longitude + Latitude


## coerce to SpatialPixelsDataFrame ----
# Mean fit
gridded(nemip.fit) <- TRUE
# SE fit
gridded(nemip.se) <- TRUE


## Coerce to raster ----

# Set the CRS --
sr <- "+init=epsg:4326"
# choose colors --
library(RColorBrewer)
my.palette <- brewer.pal(n = 9, name = "OrRd")
pal <- colorRampPalette(c("yellow","pink", "red", "dark blue"))
pal <- colorRampPalette(c("yellow","pink", "red", "dark blue"))
# plot mean fit --
nemip.fit <- raster(nemip.fit)
crs(nemip.fit)<-sr
nemip.fit
plot(nemip.fit)
# plot SE fit --
nemip.se <- raster(nemip.se)
crs(nemip.se) <-sr
nemip.se
plot(nemip.se)

# Save rasters--
writeRaster(nemip.fit, paste(model.out, "Nemipfit.tif", sep ='/'))
writeRaster(nemip.se, paste(model.out, "Nemipse.tif", sep ='/'))


## Plot with ggplot ----

### Overall plot

## Fit

nemidp<-ggplot()+
  geom_tile(data=prediction.nemid_df,aes(x=Longitude,y=Latitude,fill=fit),alpha=0.8)+
  scale_fill_viridis(option = "magma",direction = -1)+
  #geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+ # can add polygon of coastline here
  scale_color_gradient()+
  coord_equal()+
  #xlim(112,155)+ #  set limits of plot if desired
  #ylim(-45,-7)+ #  set limits of plot
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))

nemidp

## Repeat for other species ----
