require(rstanarm)
require(tidyverse)
require(dplyr)
require(mgcv)
require(FSSgam)
require(MuMIn)
require(doBy)


rm(list=ls())

## Set working directory----
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

# Remove NA values and sites in state no-take - for one of the predictor variables
# (can interpolate but I don't know how to do that right now...)
# 8.05 10.09 10.12 16.03 All have NAs

dat<-dat%>%
  filter(!sample%in%c("8.05","10.09","10.12","16.03", "10.14", "10.15", "10.13", "14.13", "14.12", "14.02"))

# Transform variables 
dat <- dat%>%
  mutate(sqrt.reef=sqrt(reef))%>%
  mutate(sqrt.slope=sqrt(Slope))%>%
  #mutate(sqrt.TPI=sqrt(TPI))%>%
  mutate(log.roughness=log(Roughness+1))%>%
  mutate(cube.Aspect=(Aspect)^3)%>%
  glimpse()

# Subset data 
legal.dat <- subset(dat, model=='Legal')

# Add row number 
legal.dat <- legal.dat %>% 
  mutate(id = row_number())

#Set predictor variables 
pred.vars=c("sqrt.slope","cube.Aspect","log.roughness","FlowDir",
            "distance.to.ramp", "site")

# Set legal and sublegal
unique.vars=unique(as.character(dat$model))
unique.vars.use=character()
for(i in 1:length(unique.vars)){
  temp.dat=dat[which(dat$model==unique.vars[i]),]
  if(length(which(temp.dat$response==0))/nrow(temp.dat)<0.8){
    unique.vars.use=c(unique.vars.use,unique.vars[i])}
}
unique.vars.use 
modnames <- unique.vars.use

#### FSSgam using lme4 + random site ####
setwd(m.dir)
# Remove any unused columns from the dataset 
use.dat <- legal.dat%>%
  dplyr::select(response, distance.to.ramp, cube.Aspect, sqrt.slope, log.roughness, status, FlowDir,
                bathymetry, site, id)
factor.vars <- c("status")

out.all <- list()
var.imp <- list()

Model1 <- uGamm(response~s(bathymetry,k=5, bs='cr'),
                family=poisson, random=~(1|site),
                data=use.dat,
                lme4=TRUE)


model.set <- generate.model.set(use.dat=use.dat,
                                test.fit=Model1,
                                pred.vars.cont=pred.vars,
                                pred.vars.fact=factor.vars,
                                smooth.smooth.interactions=TRUE,
                                max.predictors=5,
                                k=5,
                                null.terms= "s(bathymetry, k=5, bs='cr')")



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

write.csv(all.mod.fits[,-2],file=paste(name,"lme4.random.all.mod.fits.csv",sep="_"))
write.csv(all.var.imp,file=paste(name,"lme4.all.var.imp.csv",sep="_"))

### Predict from the model ###
fits.all=list()
mod.dat.all=list()
out.fits.all=list()
mcmc.out=list()

effect.dat <- expand.grid(cube.Aspect=seq(min(legal.dat$cube.Aspect),max(legal.dat$cube.Aspect),length.out = 20),
                        distance.to.ramp=seq(min(legal.dat$distance.to.ramp),max(legal.dat$distance.to.ramp),length.out = 20),
                        bathymetry=seq(min(legal.dat$bathymetry),max(legal.dat$bathymetry), length.out=20),
                        site=lme4.site$gam$model$site)%>%
  
  distinct()%>%
  glimpse()

use.dat <- dat%>%
  dplyr::select(response, distance.to.ramp, cube.Aspect, sqrt.slope, log.roughness, status, FlowDir,
                bathymetry, site)
factor.vars <- c("status")

for(s in 1:length(modnames)){
  
  
  Model1 <- uGamm(response~s(bathymetry,k=5,bs='cr')+ s(distance.to.ramp, k=5, bs='cr')
                  + s(cube.Aspect,k=5,bs='cr'),random=~(1|site), 
                     family=poisson(), data=use.dat, lme4=TRUE)
  
  
  mod.list <- list(full=Model1)
  
  mod.data.out <- data.frame(do.call("rbind",
                                     lapply(mod.list,FUN=function(x){unlist(extract.mod.dat(x, r2.type="r2"))}))) %>%
    mutate(delta.AICc=round(AICc-min(AICc,na.rm=T),1),
           delta.BIC=round(BIC-min(BIC,na.rm=T),1),
           wi.AICc=round(wi(AICc),2),
           wi.BIC=round(wi(BIC),2)) %>%
    dplyr::select_if(~sum(!is.na(.)) > 0) %>%
    mutate(R2=round(r2.vals,3),
           edf=round(edf,1),
           #Site=site.labels.short[[s]],
           u.R2=round(max(r2.vals)-r2.vals,3)) %>%
    dplyr::select(c("u.R2","delta.AICc","wi.AICc"))
  
  # now fit with shrinkage - to generate predicted values
  Model.2 <- uGamm(response~s(bathymetry,k=5,bs='cr')+ s(distance.to.ramp, k=5, bs='cr') 
                   + s(cube.Aspect,k=5,bs='cr'), random=~(1|site), 
                     family=poisson(), data=use.dat, lme4=TRUE)
  
  stangam.s  <- stan_gamm4(response~s(bathymetry, k = 5, bs = "cr")+ s(distance.to.ramp, k=5, bs='cr') 
                           + s(cube.Aspect,k=5,bs='cr'), random=~(1|site), adapt_delta = 0.99,
                           data=use.dat, chains=3, cores=3, iter=21000, warmup=20000,
                           family=gaussian())
  
  fits.all <- c(fits.all, list(mod.list))
  mod.dat.all <- c(mod.dat.all, list(mod.data.out))
  out.fits.all <- c(out.fits.all, list(Model.2))
  mcmc.out <- c(mcmc.out,list(stangam.s))
  
}
names(fits.all) <- modnames
names(mod.dat.all) <- modnames
names(out.fits.all) <- modnames
names(mcmc.out) <- modnames


# calcualte frequentist effects
tt=lapply(out.fits.all$legal$gam,FUN="predict.gam",newdata=effect.dat)
effect.dat <- cbind(effect.dat,do.call("cbind", tt))


## use the bayesian models to calculate a posterior distribution of the effect size.-------
effect.dat.sim <- effect.dat
bayes.effects <- list()
bayes.cv <- list()
for(s in 1:length(modnames)){
  stangam.s <- mcmc.out[[s]]
  gg <- t((posterior_predict(stangam.s, effect.dat.sim)))
  colnames(gg) <- paste("sim",1:ncol(gg),sep="_")
  gg.dat <- cbind(effect.dat.sim,gg)
  
  # bathy effects
  bathy.E.log <- summaryBy(as.formula(paste(paste(paste("sim",1:ncol(gg),sep="_"),collapse="+"),
                                            "~bathymetry",sep="")), FUN=mean, data=gg.dat, keep.names = T)
  #waves.E <- as.numeric(exp(waves.E.log[2,1:ncol(gg)])-exp(waves.E.log[1,1:ncol(gg)]))
  bathy.E <- round(apply(
    apply(exp(bathy.E.log[,paste("sim",1:ncol(gg),sep="_")]), MARGIN=2,FUN="range"),
    MARGIN=2,FUN="diff"),2)
  
  # aspect effects
  aspect.E.log <- summaryBy(as.formula(paste(paste(paste("sim",1:ncol(gg),sep="_"),collapse="+"),
                                           "~cube.Aspect",sep="")), FUN=mean, data=gg.dat, keep.names = T)
  #wind.E <- as.numeric(exp(wind.E.log[2,1:ncol(gg)])-exp(wind.E.log[1,1:ncol(gg)]))
  aspect.E <- round(apply(
    apply(exp(aspect.E.log[,paste("sim",1:ncol(gg),sep="_")]), MARGIN=2,FUN="range"),
    MARGIN=2,FUN="diff"),2)
  
  # ramp effects
  ramp.E.log <- summaryBy(as.formula(paste(paste(paste("sim",1:ncol(gg),sep="_"),collapse="+"),
                                             "~distance.to.ramp",sep="")), FUN=mean, data=gg.dat, keep.names = T)
  #wind.E <- as.numeric(exp(wind.E.log[2,1:ncol(gg)])-exp(wind.E.log[1,1:ncol(gg)]))
  ramp.E <- round(apply(
    apply(exp(ramp.E.log[,paste("sim",1:ncol(gg),sep="_")]), MARGIN=2,FUN="range"),
    MARGIN=2,FUN="diff"),2)
  
  ## calculate effect as cv
  mean.NTU <- summaryBy(as.formula(paste(paste(paste("sim",1:ncol(gg),sep="_"),collapse="+"),"~Site.Code")),FUN=mean,
                        data=gg.dat)
  
  bayes.effects.s <- list(bathy.E=bathy.E,aspect.E=aspect.E,ramp.E=ramp.E)
  bayes.cv.s <- lapply(bayes.effects.s,FUN=function(x){unlist(x/exp(mean.NTU))})
  
  bayes.effects <- c(bayes.effects,list(bayes.effects.s))
  bayes.cv <- c(bayes.cv,list(bayes.cv.s))
}
