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
library(FSSgam)

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
dat <-read.csv('final.data.csv')%>%
  rename(response=target.fish)%>%
  select(!X)%>%
  select(!X.1)%>%
  glimpse()

names(dat)

# Correct Lat/Long - had been changed to UTM in a previous script
metadata <- read.csv('ningaloo_metadata.csv')

latlongs <- metadata%>%
  select('longitude', 'latitude')

dat <- cbind(dat, latlongs)
glimpse(dat)
dat <- dat[,-c(3,4)]

legal.dat <- subset(dat, model=='Legal')

# Remove NA values  - for one of the predictor variables
# (can interpolate but I don't know how to do that right now...)
# 8.05 10.09 10.12 16.03 All have NAs

legal.dat<-legal.dat%>%
  filter(!sample%in%c("8.05","10.09","10.12","16.03"))

# Transform variables 
legal.dat <- legal.dat%>%
  mutate(sqrt.reef=sqrt(reef))%>%
  mutate(sqrt.slope=sqrt(Slope))%>%
  #mutate(sqrt.TPI=sqrt(TPI))%>%
  mutate(log.roughness=log(Roughness+1))%>%
  mutate(cube.Aspect=(Aspect)^3)%>%
  glimpse()

#Set predictor variables 
pred.vars=c("bathymetry","sqrt.slope","cube.Aspect","log.roughness","FlowDir",
            "distance.to.ramp")

##### With Site as a random factor ####
setwd(m.dir)
use.dat <- legal.dat
factor.vars <- c("status") # Status as a Factor with two levels
out.all <- list()
var.imp <- list()


Model1 <- gam(response~s(bathymetry,k=3,bs='cr')+ s(site,bs="re") + s(TPI,bs="re"),
           family=tw(),  data=use.dat)

model.set <- generate.model.set(use.dat=use.dat,
                             test.fit=Model1,
                             pred.vars.cont=pred.vars,
                             pred.vars.fact=factor.vars,
                             max.predictors=5,
                             k=5,
                             null.terms="s(site,bs='re')+s(TPI,bs='re')")

out.list=fit.model.set(model.set,
                       max.models=600,
                       parallel=T)
names(out.list)

out.list$failed.models # examine the list of failed models
mod.table=out.list$mod.data.out  # look at the model selection table
mod.table=mod.table[order(mod.table$AICc),]
mod.table$cumsum.wi=cumsum(mod.table$wi.AICc)
out.i=head(mod.table,10)
#out.i=mod.table[which(mod.table$delta.AICc<=3),]
out.all=c(out.all,list(out.i))
var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) 

# Model fits and importance---
all.mod.fits=do.call("rbind",out.all)
all.var.imp=do.call("rbind",var.imp)

write.csv(all.mod.fits[,-2],file=paste(name,"site.all.mod.fits.csv",sep="_"))
write.csv(all.var.imp,file=paste(name,"site.all.var.imp.csv",sep="_"))

#### With te spline and no specified k ####
setwd(m.dir)
use.dat <- legal.dat
factor.vars <- c("status") # Status as a Factor with two levels
out.all <- list()
var.imp <- list()


Model1 <- gam(response~s(bathymetry,k=3,bs='cr') + s(TPI,bs="re")
              + te(latitude, longitude, bs="cr"),
              family=tw(),  data=use.dat)

model.set <- generate.model.set(use.dat=use.dat,
                                test.fit=Model1,
                                pred.vars.cont=pred.vars,
                                pred.vars.fact=factor.vars,
                                max.predictors=5,
                                k=5,
                                null.terms="s(TPI,bs='re')+
                                te(latitude, longitude, bs='cr')")

out.list=fit.model.set(model.set,
                       max.models=600,
                       parallel=T)
names(out.list)

out.list$failed.models # examine the list of failed models
mod.table=out.list$mod.data.out  # look at the model selection table
mod.table=mod.table[order(mod.table$AICc),]
mod.table$cumsum.wi=cumsum(mod.table$wi.AICc)
out.i=head(mod.table,10)
#out.i=mod.table[which(mod.table$delta.AICc<=3),]
out.all=c(out.all,list(out.i))
var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) 

par(mar = c(1, 1,))
# Model fits and importance---
all.mod.fits=do.call("rbind",out.all)
all.var.imp=do.call("rbind",var.imp)

write.csv(all.mod.fits[,-2],file=paste(name,"site.all.mod.fits.csv",sep="_"))
write.csv(all.var.imp,file=paste(name,"site.all.var.imp.csv",sep="_"))

model.best <- gam(response ~ s(cube.Aspect, k = 5, bs = "cr") + 
                    s(bathymetry, by = status, k = 5, bs = "cr") + status + 
                    s(TPI, bs = "re") + te(latitude, longitude, bs = "cr"), data=legal.dat)
plot.gam(model.best)
