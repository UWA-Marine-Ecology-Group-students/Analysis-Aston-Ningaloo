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

#### Model 1 ####
### Predict from the model ###

effect.dat <- expand.grid(cube.Aspect=seq(min(legal.dat$cube.Aspect),max(legal.dat$cube.Aspect),length.out = 20),
                        distance.to.ramp=seq(min(legal.dat$distance.to.ramp),max(legal.dat$distance.to.ramp),length.out = 20),
                        bathymetry=seq(min(legal.dat$bathymetry),max(legal.dat$bathymetry), length.out=20))
                        #site=lme4.site$gam$model$site)%>%

use.dat <- dat%>%
  dplyr::select(response, distance.to.ramp, cube.Aspect, sqrt.slope, log.roughness, status, FlowDir,
                bathymetry, site)
factor.vars <- c("status")
  
  
Model.1 <- uGamm(response~s(bathymetry,k=5,bs='cr')+ s(distance.to.ramp, k=5, bs='cr')
                  + s(cube.Aspect,k=5,bs='cr'),random=~(1|site), 
                     family=poisson(), data=use.dat, lme4=TRUE)
  

  
# now fit with shrinkage - to generate predicted values
Model.2 <- uGamm(response~s(bathymetry,k=5,bs='cr')+ s(distance.to.ramp, k=5, bs='cr') 
                   + s(cube.Aspect,k=5,bs='cr'), random=~(1|site), 
                     family=poisson(), data=use.dat, lme4=TRUE)
  
stangam.s  <- stan_gamm4(response~s(bathymetry, k = 5, bs = "cr")+ s(distance.to.ramp, k=5, bs='cr') 
                           + s(cube.Aspect,k=5,bs='cr'), random=~(1|site), adapt_delta = 0.99,
                           data=use.dat, chains=3, cores=3, iter=21000, warmup=20000,
                           family=poisson())


# calcualte frequentist effects
predicted <- predict.gam(Model.1$gam, newdata=effect.data, type='response')
effect.dat$predicted <- cbind(effect.dat, predicted)

# bathy effect 
bathy.P <- effect.dat%>%
  group_by(bathymetry)%>%
  summarise_at(var(predicted), list(mean))

bathy.E <- diff(range(bathy.P$predicted))

# aspect effect 
aspect.P <- effect.dat%>%
  group_by(cube.Aspect)%>%
  summarise_at(var(predicted), list(mean))

aspect.E <- diff(range(aspect.P$predicted))

# ramp effect 
ramp.P <- effect.dat%>%
  group_by(distance.to.ramp)%>%
  summarise_at(var(predicted), list(mean))

ramp.E <- diff(range(ramp.P$predicted))

## use the bayesian models to calculate a posterior distribution of the effect size.-------
effect.dat.sim <- effect.dat
gg <- t(posterior_predict(stangam.s, effect.dat.sim))
#colnames(gg) <- paste("sim",1:ncol(gg),sep="_")
gg.dat <- cbind(effect.dat.sim,gg)

gg.dat.long <- gg.dat%>%
  gather(sim, predicted.sim, 4:3000)

gg.dat.long$sim <- as.factor(gg.dat.long$sim)

# bathy effect 
bathy.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.mean")
colnames(bathy.sim) <- x

for(i in 1:3000){
  bathy.sim.1 <- gg.dat.long
  bathy.sim.2 <- filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- group_by(bathy.sim.3, bathymetry)
  bathy.sim.4 <- summarise_at(bathy.sim.3, vars(predicted.sim), list(mean))
  print(bathy.sim.5 <- diff(range(bathy.sim.4$predicted.sim)))
  bathy.sim <- rbind(bathy.sim, bathy.sim.5)
  
}

bathy.sim

write.csv(bathy.sim, "bayesian.bathy.predictions")

# aspect effect 
aspect.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("aspect.mean")
colnames(aspect.sim) <- x

for(i in 1:3000){
  aspect.sim.1 <- gg.dat.long
  aspect.sim.2 <- filter(aspect.sim.1, sim==i)
  aspect.sim.3 <- group_by(aspect.sim.3, cube.Aspect)
  aspect.sim.4 <- summarise_at(aspect.sim.3, vars(predicted.sim), list(mean))
  print(aspect.sim.5 <- diff(range(aspect.sim.4$predicted.sim)))
  aspect.sim <- rbind(aspect.sim, aspect.sim.5)
  
}

aspect.sim

write.csv(aspect.sim, "bayesian.aspect.predictions")

# ramp effect 
ramp.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("ramp.mean")
colnames(ramp.sim) <- x

for(i in 1:3000){
  ramp.sim.1 <- gg.dat.long
  ramp.sim.2 <- filter(ramp.sim.1, sim==i)
  ramp.sim.3 <- group_by(ramp.sim.3, distance.to.ramp)
  aspect.sim.4 <- summarise_at(aspect.sim.3, vars(predicted.sim), list(mean))
  print(rampsim.5 <- diff(range(ramp.sim.4$predicted.sim)))
  ramp.sim <- rbind(ramp.sim, ramp.sim.5)
  
}

ramp.sim

write.csv(ramp.sim, "bayesian.ramp.predictions")

setwd(m.dir)

#### Model 2 ####
### Predict from the model ###

use.dat <- legal.dat%>%
  dplyr::select(response, distance.to.ramp, cube.Aspect, sqrt.slope, log.roughness, status, FlowDir,
                bathymetry, site)
factor.vars <- c("status")


Model.2 <- uGamm(response~s(bathymetry,k=5,bs='cr')+ factor(status)
                 + s(cube.Aspect,k=5,bs='cr'),random=~(1|site), 
                 family=poisson(), data=use.dat, lme4=TRUE)



# now fit with shrinkage - to generate predicted values
Model.2 <- uGamm(response~s(bathymetry,k=5,bs='cr')+ s(distance.to.ramp, k=5, bs='cr') 
                 + s(cube.Aspect,k=5,bs='cr'), random=~(1|site), 
                 family=poisson(), data=use.dat, lme4=TRUE)

stangam.s.2  <- stan_gamm4(response~s(bathymetry, k = 5, bs = "cr")+ status
                         + s(cube.Aspect,k=5,bs='cr'), random=~(1|site), adapt_delta = 0.99,
                         data=use.dat, chains=3, cores=3, iter=21000, warmup=20000,
                         family=poisson())



effect.dat.2 <- expand.grid(cube.Aspect=seq(min(legal.dat$cube.Aspect),max(legal.dat$cube.Aspect),length.out = 15),
                          bathymetry=seq(min(legal.dat$bathymetry),max(legal.dat$bathymetry), length.out=15),
                          status=Model.2$gam$model$status)

# calcualte frequentist effects
predicted <- predict.gam(Model.2$gam, newdata=effect.dat, type='response')
effect.dat.2 <- cbind(effect.dat.2, predicted)

# bathy effect 
bathy.P.2 <- effect.dat.2%>%
  group_by(bathymetry)%>%
  summarise_at(vars(predicted), list(mean))

bathy.E.2 <- diff(range(bathy.P.2$predicted))
bathy.E.2

# aspect effect 
aspect.P.2 <- effect.dat.2%>%
  group_by(cube.Aspect)%>%
  summarise_at(vars(predicted), list(mean))

aspect.E.2 <- diff(range(aspect.P.2$predicted))
aspect.E.2

# ramp effect 
status.P.2 <- effect.dat.2%>%
  group_by(status)%>%
  summarise_at(vars(predicted), list(mean))

status.E.2 <- diff(range(status.P.2$predicted))
status.E.2

## use the bayesian models to calculate a posterior distribution of the effect size.-------
effect.dat.sim.2 <- effect.dat.2
gg <- t(posterior_predict(stangam.s.2, effect.dat.sim.2, re.form=NA))
#colnames(gg) <- paste("sim",1:ncol(gg),sep="_")
gg.dat <- cbind(effect.dat.sim.2,gg)

head(gg.dat)
nrow(gg.dat)

gg.dat.long <- gg.dat%>%
  gather(sim, predicted.sim, 4:3003)

gg.dat.long$sim <- as.factor(gg.dat.long$sim)
head(gg.dat.long)

str(gg.dat.long)
# bathy effect 
bathy.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.mean")
colnames(bathy.sim) <- x


for(i in 1:3000){
  bathy.sim.1 <- gg.dat.long
  bathy.sim.2 <- filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- group_by(bathy.sim.2, bathymetry)
  bathy.sim.4 <- summarise_at(bathy.sim.3, vars(predicted.sim), list(mean))
  print(bathy.sim.5 <- diff(range(bathy.sim.4$predicted.sim)))
  bathy.sim <- rbind(bathy.sim, bathy.sim.5)
  
}

bathy.sim

write.csv(bathy.sim, "bayesian.bathy.predictions.2.csv")

# aspect effect 
aspect.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("aspect.mean")
colnames(aspect.sim) <- x

for(i in 1:3000){
  aspect.sim.1 <- gg.dat.long
  aspect.sim.2 <- filter(aspect.sim.1, sim==i)
  aspect.sim.3 <- group_by(aspect.sim.2, cube.Aspect)
  aspect.sim.4 <- summarise_at(aspect.sim.3, vars(predicted.sim), list(mean))
  print(aspect.sim.5 <- diff(range(aspect.sim.4$predicted.sim)))
  aspect.sim <- rbind(aspect.sim, aspect.sim.5)
  
}

aspect.sim

write.csv(aspect.sim, "bayesian.aspect.predictions.2.csv")

# ramp effect 
status.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("status.mean")
colnames(status.sim) <- x

for(i in 1:3000){
  status.sim.1 <- gg.dat.long
  status.sim.2 <- filter(status.sim.1, sim==i)
  status.sim.3 <- group_by(status.sim.2, status)
  status.sim.4 <- summarise_at(status.sim.3, vars(predicted.sim), list(mean))
  print(status.sim.5 <- diff(range(status.sim.4$predicted.sim)))
  status.sim <- rbind(status.sim, status.sim.5)
  
}

status.sim

write.csv(ramp.sim, "bayesian.status.predictions.2.csv")

## Create density plot
full.data <- cbind(status.sim, bathy.sim, aspect.sim)
full.data <- full.data[,3:1]
full.data.df <- as.data.frame(full.data)

full.data.long <- full.data.df%>%
  gather(variable, predicted, 1:3)

effect.plot <- ggplot(full.data.long, aes(x = predicted, fill = variable)) + geom_density(alpha = 0.5)
effect.plot




