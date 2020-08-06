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

# Correct Lat/Long - had been changed to UTM in a previous script
metadata <- read.csv('ningaloo_metadata.csv')

latlongs <- metadata%>%
  select('longitude', 'latitude')


dat <- cbind(dat, latlongs)
dat <- dat[,-c(4,5)]

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
legal.dat<- legal.dat%>%
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
  x<-legal.dat[ ,i]
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
#Run mot likely models and look at AIC/BIC as well as Moran's I 



