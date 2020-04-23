
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

data.dir <- paste(working.dir,"Tidy data",sep="/")
plots.dir <- paste(working.dir,"Plots",sep="/")
model.out <- paste(working.dir,"Model Out",sep="/")

# Read in data----
# data.dir<- ("C:/Users/00093391/Dropbox/UWA/Research Associate/Woodside-Exmouth/fishdata")
setwd(data.dir)
dir()

dat<-read.csv("fishdata.csv")%>% # this is the df with fish and predictor data
  dplyr::rename(site=sample)%>% # rename columns in data for model
  # dplyr::filter(Taxa!='Synodontidae Saurida undosquamis')%>% remove Synodontidae Saurida undosquamis
  na.omit()%>% # remove NAs
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

Sebae.gam <- gam(maxn ~ s(TPI, k = 6, bs = "cr") + s(bathymetry, k = 6, bs = "cr") + # covariate effect of tpi and bathymetry
           + s(site, k = 3, bs ='re'), # random effect of cluster or site
           family=tw(),data=Sebae) # family tweedy
plot(Sebae.gam,residuals=T,all.terms = TRUE,pages=1)
summary(Sebae.gam) # check results
AICc(Sebae.gam) # check AIC of model

## Check residuals
par(mfrow=c(2,2))
gam.check(Sebae.gam)

## Fit model for Pristipomoides multidens ----

names(dat)

Multidens.gam <- gam(maxn ~ s(TPI, k = 6, bs = "cr") + s(bathymetry, k = 6, bs = "cr") +  # covariate effect of tpi and bathymetry
                 + s(site, k = 3, bs ='re'), # random effect of cluster or site
                 family=tw(),data=Multidens) # family tweedy
plot(Multidens.gam,residuals=T,all.terms = TRUE,pages=1)
summary(Multidens.gam) # check results
AICc(Multidens.gam) # check AIC of model

## Check residuals
par(mfrow=c(2,2))
gam.check(Multidens.gam)

