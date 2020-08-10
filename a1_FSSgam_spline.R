# librarys----
detach("package:plyr", unload=TRUE)#will error - don't worry
library(tidyr)
library(dplyr)
options(dplyr.width = Inf) #enables head() to display all coloums
library(mgcv)
library(MuMIn)
library(car)
library(doBy)
library(gplots)
library(RColorBrewer)
library(doParallel) #this can removed?
library(doSNOW)
library(gamm4)
#library(RCurl) #needed to download data from GitHub
library(FSSgam)
library(spdep)
library(spatialEco)
library(nlme)
devtools::install_github("beckyfisher/FSSgam_package") #run once


rm(list=ls())

##Set working directory----
## Set work directory----
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # sets working directory to where this script is saved (DON'T MOVE THE SCRIPT)

## Set sub directories----
d.dir <- paste(working.dir,"Tidy data",sep="/") 
s.dir <- paste(working.dir,"Spatial",sep="/") # spatial is where I keep spatial data files, rasters and shapefiles
p.dir <- paste(working.dir,"Plots",sep="/")
m.dir <- paste(working.dir,"Model Out GAM", sep="/") #suggest make a model.out.gam folder to keep things seperate

# Bring in and format the data----
name <- 'ningaloo' # for the study

# Load the dataset - from github
setwd(d.dir)
dat <-read.csv('final.data.csv')

names(dat)

dat <- dat%>%
  select(!X.1)%>%
  select(!X)

# Correct Lat/Long - had been changed to UTM in a previous script
metadata <- read.csv('ningaloo_metadata.csv')

latlongs <- metadata%>%
  select('longitude', 'latitude')


dat <- cbind(dat, latlongs)
dat <- dat[,-c(3,4)]

# Add distance to -60m bathome and then convert into use for spline
distance.60m <- read.csv("distance to 60.csv")

dat <- cbind(dat,distance.60m)

#dat <- dat%>%
#  mutate(bathome.x.depth=(distance.to.60*bathymetry))%>%
#  glimpse()


# set status as a factor
dat$status <- as.factor(dat$status)

# Remove NA values  - for one of the predictor variables
# 8.05 10.09 10.12 16.03 All have NAs

dat<-dat%>%
  filter(!sample%in%c("8.05","10.09","10.12","16.03"))


#Set bathymetry to be a positive number for loop
dat <- dat%>%
  mutate(pos.bathymetry=(bathymetry*-1))%>%
  glimpse()


# Set predictor variables---
pred.vars=c("TPI","Slope","Aspect","TRI","Roughness","FlowDir","mean.relief",
            "sd.relief","reef","distance.to.ramp","pos.bathymetry")

# Check for correlation of predictor variables- remove anything highly correlated (>0.95)---
round(cor(dat[,pred.vars], use = "complete.obs"),2)


# TRI is super correlated to both slope and roughness, Roughness also super correlated with slope,
pred.vars=c("TPI","Slope","Aspect","FlowDir","mean.relief",
            "sd.relief","reef","distance.to.ramp","pos.bathymetry")

# Plot of likely transformations
par(mfrow=c(3,2))
for (i in pred.vars) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}


# Plot of likely transformations
par(mfrow=c(3,2))
for (i in pred.vars) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist((x)^2)
  plot((x)^2)
  hist((x)^3)
  plot((x)^3)
}

# Review of individual predictors - we have to make sure they have an even distribution---
# Bathymetry, distance to ramp, sd.relief, mean.relief, flow.dir, aspect, TPI leave untransformed - should use 
# sqrt for TPI but it generates a lot of NAs as some values negative and some are positive 
# Percent reef, slope use sqrt transformation
# Roughness use log+1 transformation 

# Transform variables 
dat <- dat%>%
  mutate(sqrt.reef=sqrt(reef))%>%
  mutate(sqrt.slope=sqrt(Slope))%>%
  # mutate(sqrt.TPI=sqrt(TPI))%>%
  mutate(log.roughness=log(Roughness+1))%>%
  mutate(cube.Aspect=(Aspect)^3)%>%
  glimpse()

#Reset predictor variables 
pred.vars=c("bathymetry","sqrt.slope","cube.Aspect","log.roughness","FlowDir","distance.to.ramp")

######### Legal Model #########
dat.legal<-dat%>%filter(model=="Legal")

# Bathy, TPI and Zones have to be in the null model
# Start with null model then add things one by one and look at the effect on AIC 
# First thing to look at is the Moran's I and see if spline helps 

gamm.legal.null<-gam(target.fish ~ s(bathymetry, bs="re") + s(TPI, bs="re") + s(zone, bs="re") +
                  te(distance.to.60,bathymetry, k=3, bs="cr"), family=tw(), data=dat.legal)

summary(gamm.legal.null)
gamm.legal.null$aic
# Deviance explained 23.3%
# AIC 558.821

## Check Moran's I with spline
#Create a spatial data frame
coordinates(dat.legal) <- ~ longitude + latitude
proj4string(dat.legal) # check coordinate system and/or projection
# If no projection give assign one,
proj4string(dat.legal) <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"

neighbours <- knearneigh(dat.legal, k=1, longlat = TRUE)

knn <- knn2nb(neighbours, row.names = NULL, sym = FALSE)
nb <- nb2listw(knn, glist=NULL, style="W", zero.policy=NULL)

moran.test(gamm.legal.null$residuals, nb, randomisation=T, zero.policy=NULL,
           alternative="greater", rank = FALSE)

moran.plot(gamm.legal$residuals, nb, zero.policy=NULL, spChk=NULL, labels=NULL,
           xlab=NULL, ylab=NULL, quiet=NULL, plot=TRUE)

# p-value = 0.001074

gamm.legal.no.spline <- gamm.legal.null<-gam(target.fish ~ s(bathymetry, bs="re") + s(TPI, bs="re") + s(zone, bs="re"),
                                               family=tw(), data=dat.legal)

## Check Moran's I without spline
#Create a spatial data frame
coordinates(dat.legal) <- ~ longitude + latitude
proj4string(dat.legal) # check coordinate system and/or projection
# If no projection give assign one,
proj4string(dat.legal) <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"

neighbours <- knearneigh(dat.legal, k=1, longlat = TRUE)

knn <- knn2nb(neighbours, row.names = NULL, sym = FALSE)
nb <- nb2listw(knn, glist=NULL, style="W", zero.policy=NULL)

moran.test(gamm.legal.no.spline$residuals, nb, randomisation=T, zero.policy=NULL,
           alternative="greater", rank = FALSE)

moran.plot(gamm.legal$residuals, nb, zero.policy=NULL, spChk=NULL, labels=NULL,
           xlab=NULL, ylab=NULL, quiet=NULL, plot=TRUE)

# p-value is 5.493e-05 so the spline definitely helps with spatial autocorrelation although it
# doesn't remove it entirely 

## Adding variables
gamm.legal<-gam(target.fish ~ 
                  s(distance.to.ramp, k=3, bs="cr") + 
                  status +
                 # s(log.roughness, k=3, bs="cr") +
                 # s(sqrt.slope, k=3, bs="cr") + 
                  s(cube.Aspect, k=3, bs="cr") +
                  s(bathymetry, bs="re") + s(TPI, bs="re") + s(zone, bs="re") +
                  te(distance.to.60,bathymetry, k=3, bs="cr"), family=tw(), data=dat.legal)

summary(gamm.legal)
gamm.legal$aic

# Results are on a spreadsheet in excel 

######### Sublegal Model #########
dat.sublegal<-dat%>%filter(model=="Sublegal")

# Bathy, TPI and Zones have to be in the null model
# Start with null model then add things one by one and look at the effect on AIC 
# First thing to look at is the Moran's I and see if spline helps 

######### Sublegal ##########
dat.sublegal<-dat%>%filter(model=="Sublegal")

## Adding variables
gamm.sublegal<-gam(target.fish ~ 
                     s(distance.to.ramp, k=3, bs="cr") + 
                    # status +
                    # s(log.roughness, k=3, bs="cr") +
                     s(sqrt.slope, k=3, bs="cr") + 
                    # s(cube.Aspect, k=3, bs="cr") +
                    s(bathymetry, bs="re") + s(TPI, bs="re") + s(zone, bs="re") +
                    te(distance.to.60,bathymetry, k=3, bs="cr"), family=tw(), data=dat.sublegal)

summary(gamm.sublegal)
gamm.sublegal$aic


