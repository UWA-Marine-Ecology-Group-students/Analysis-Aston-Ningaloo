
rm(list=ls()) # clear memory


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

data.dir <- paste(working.dir,"data",sep="/")
plots.dir <- paste(working.dir,"plots",sep="/")
model.out <- paste(working.dir,"ModelOut",sep="/")

# Read in data----
data.dir<- ("C:/Users/00093391/Dropbox/UWA/Research Associate/Woodside-Exmouth/fishdata")
setwd(data.dir)
dir()

dat<-read.csv("fishpreds.csv")%>% # this is the df with fish and predictor data
  dplyr::rename(Taxa=Species,
                response=Carangidae.Carangoides.equula)%>% # rename columns in data for model
  dplyr::filter(Taxa!='Synodontidae Saurida undosquamis')%>% # remove Synodontidae Saurida undosquamis
  na.omit()%>% # remove NAs
  glimpse() # to see data
# Convert covariates that are stored as character (factor) or integers (continous) if needed

#### Check distribution of the response for species ----

plot.new()
par(mfrow=c(1,2))

levels(dat$Taxa)

## Subsample data for each species needed for the analysis
# Repeat for as many species as required --


## Species 1: Nemip --

levels(dat$Taxa)

Nemip<-dat%>%
  filter(Taxa=="Nemipteridae Nemipterus bathybius")%>%
  glimpse()
hist(Nemip$response)
plot(Nemip$response)

## Species 2: Carangid.2 --

levels(dat$Taxa)

carangid.2<-dat%>%
  filter(Taxa=="Carangidae Decapterus tabl")%>%
  glimpse()
hist(carangid.2$response)
plot(carangid.2$response)



## Convert fish data into spatial object ----
coordinates(dat)<-~longitude+latitude # check the name of the columns for lon and lat
class(dat) # should be spatial object


## Set the reference system to the widely used WGS84
proj4string(dat)<-CRS("+init=epsg:4326")
dat@proj4string


###########################


##  Fit a top model ----
# First I will do marginal predictions (no random effects) and not accounting for the offset 
# (generates rates - really low numbers)
#library(mgcv)
#library(MuMIn)

## Fit model for Nemip ----

names(dat)

Nemip <- gam(response ~ s(tpi, k = 6, bs = "cr") + s(bathy, k = 6, bs = "cr") + # covariate effect of tpi and bathymetry
            + s(cluster, k = 3, bs ='re'), # random effect of cluster or site
          family=tw(),data=Nemip) # family tweedy
plot(Nemip,residuals=T,all.terms = TRUE,pages=1)
summary(Nemip) # check results
AICc(Nemip) # check AIC of model

## Check residuals
par(mfrow=c(2,2))
gam.check(Nemip)

## Repeat for other species as needed ----
