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
require(rstanarm)
library(ggplot2)

rm(list=ls())

## Set working directory----
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # sets working directory to where this script is saved (DON'T MOVE THE SCRIPT)

## Set sub directories----
d.dir <- paste(working.dir,"Tidy data",sep="/") 
s.dir <- paste(working.dir,"Spatial",sep="/") # spatial is where I keep spatial data files, rasters and shapefiles
p.dir <- paste(working.dir,"Plots",sep="/")
m.dir <- paste(working.dir,"Model Out GAM", sep="/") #suggest make a model.out.gam folder to keep things seperate

Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend.background = element_rect(fill="white"),
    legend.background = element_blank(),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=15),
    legend.title = element_blank(),
    legend.position = c(0.75, 0.8),
    text=element_text(size=15),
    strip.text.y = element_text(size = 15,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=15),
    axis.title.y=element_text(vjust=0.6, angle=90, size=15),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())

#### Species richness models ####

# Bring in and format the data
name <- 'ningaloo' # for the study

setwd(d.dir)
richness.dat <-read.csv("ningaloo.total.abundance.and.species.richness.csv")
richness.dat <- richness.dat[116:230,c(1,3)]

variables <- read.csv("final.data.csv")

metadata <- read.csv('ningaloo_metadata.csv')

latlongs <- metadata%>%
  select('longitude', 'latitude')

variables <- cbind(variables, latlongs)
variables <- variables[,-c(5,6)]

variables <- variables%>%
  filter(model=="Legal")%>%
  mutate(sqrt.reef=sqrt(reef))%>%
  mutate(sqrt.slope=sqrt(Slope))%>%
  #mutate(sqrt.TPI=sqrt(TPI))%>%
  mutate(log.roughness=log(Roughness+1))%>%
  mutate(cube.Aspect=(Aspect)^3)%>%
  glimpse()

variables <- variables%>%
  dplyr::select('sample', 'site', 'bathymetry','distance.to.ramp', 'status','log.roughness', 'cube.Aspect', 'sqrt.slope', 
                 'latitude', 'longitude')

full.data.rich <- merge(richness.dat,variables, all=TRUE)
full.data.rich[is.na(full.data.rich)] <- 0

full.data.rich <- full.data.rich%>%
  dplyr::rename(response=maxn)

# Remove NA values and sites in state no-take - for one of the predictor variables
# (can interpolate but I don't know how to do that right now...)
# 8.05 10.09 10.12 16.03 All have NAs

full.data.rich<-full.data.rich%>%
  filter(!sample%in%c("8.05","10.09","10.12","16.03", "10.14", "10.15", "10.13", "14.13", "14.12", "14.02"))

# Add row ID
full.data.rich <- full.data.rich%>%
  mutate(ID = rownames(full.data.rich))

#Set predictor variables 
pred.vars=c("sqrt.slope","cube.Aspect","log.roughness",
            "distance.to.ramp")

## Run full subsets 
use.dat <- full.data.rich%>%
  dplyr::select(response, distance.to.ramp, cube.Aspect, sqrt.slope, log.roughness, status,
                bathymetry,site)
factor.vars <- c("status")

out.all <- list()
var.imp <- list()

Model.rich <- uGamm(response~s(bathymetry,k=5, bs='cr'),
                family=poisson, random=~(1|site),
                data=use.dat,
                lme4=TRUE)


model.set <- generate.model.set(use.dat=use.dat,
                                test.fit=Model.rich,
                                pred.vars.cont=pred.vars,
                                pred.vars.fact=factor.vars,
                                factor.smooth.interactions=NA,
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

setwd(d.dir)
write.csv(all.mod.fits, "species.richness.models.csv")

## Top models
# fit with shrinkage
effect.dat.richness <- expand.grid(cube.Aspect=seq(min(full.data.rich$cube.Aspect),max(full.data.rich$cube.Aspect),length.out = 20),
                          bathymetry=seq(min(full.data.rich$bathymetry),max(full.data.rich$bathymetry), length.out=20),
                          status=unique(full.data.rich$status))

Model.richness <- uGamm(response~s(bathymetry,k=5,bs='cs')+ factor(status)
                 + s(cube.Aspect,k=5,bs='cs'),
                 random=~(1|site), 
                 family=poisson(), data=use.dat, lme4=TRUE)

gam.check(Model.richness$gam)

stangam.s.richness  <- stan_gamm4(response~s(bathymetry, k = 5, bs = "cs") + status
                           + s(cube.Aspect, k=5,bs='cs'),
                           random=~(1|site), adapt_delta = 0.99,
                           data=use.dat, chains=3, cores=3, iter=41000, warmup=40000, thin=3,
                           family=poisson)

## calcualte frequentist effects
# bathy effect 
predicted <- predict.gam(Model.richness$gam, newdata=effect.dat.richness, type='response')
effect.dat.richness <- cbind(effect.dat.richness, predicted)

bathy.P.rich <- effect.dat.richness%>%
  group_by(bathymetry)%>%
  summarise_at(vars(predicted), list(mean))

bathy.E.rich <- diff(range(bathy.P.rich$predicted))
bathy.E.rich

# Plot predicted
predicted <- predict.gam(Model.richness$gam, newdata=effect.dat.richness, type='response', se.fit=T)
effect.dat.richness <- cbind(effect.dat.richness, predicted)

predicts.legal.bathy = effect.dat.richness%>%data.frame(predicted)%>%
  group_by(bathymetry)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.legal.bathy<- ggplot() +
  ylab("Predicted Species Richness")+
  xlab('Bathymetry (m)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.data.rich,aes(x=bathymetry,y=response), colour="lightblue", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response), colour='darkblue', alpha=0.75)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response - se.fit), colour='darkblue', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response + se.fit), colour='darkblue', linetype="dashed",alpha=0.75)+
  theme_classic()+
  Theme1
  #annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.legal.bathy


# aspect by status effect 
effect.dat.richness <- expand.grid(cube.Aspect=seq(min(full.data.rich$cube.Aspect),max(full.data.rich$cube.Aspect),length.out = 20),
                                   bathymetry=seq(min(full.data.rich$bathymetry),max(full.data.rich$bathymetry), length.out=20),
                                   status=unique(full.data.rich$status))
predicted <- predict.gam(Model.richness$gam, newdata=effect.dat.richness, type='response')
effect.dat.richness <- cbind(effect.dat.richness, predicted)

aspect.P <- effect.dat.richness%>%
  group_by(cube.Aspect)%>%
  summarise_at(vars(predicted), list(mean))

aspect.E <- diff(range(aspect.P$predicted))
aspect.E


# Plot it 
predicted <- predict.gam(Model.richness$gam, newdata=effect.dat.richness, type='response', se.fit=T)
effect.dat.richness <- cbind(effect.dat.richness, predicted)

predicts.aspect = effect.dat.richness%>%data.frame(predicted)%>%
  group_by(cube.Aspect)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.aspect<- ggplot() +
  ylab("Predicted Species Richness")+
  xlab('Aspect (cubed)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.data.rich,aes(x=cube.Aspect,y=response), colour="lightgreen", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.aspect,aes(x=cube.Aspect,y=response), colour='darkgreen', alpha=0.75)+
  geom_line(data=predicts.aspect,aes(x=cube.Aspect,y=response - se.fit), colour='darkgreen', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.aspect,aes(x=cube.Aspect,y=response + se.fit), colour='darkgreen', linetype="dashed",alpha=0.75)+
  theme_classic()+
  Theme1
  #annotate("text", x = -Inf, y=Inf, label = "NT",vjust = 2, hjust = -.1,size=5)
ggmod.aspect

# status effect 
effect.dat.richness <- expand.grid(cube.Aspect=seq(min(full.data.rich$cube.Aspect),max(full.data.rich$cube.Aspect),length.out = 20),
                                   bathymetry=seq(min(full.data.rich$bathymetry),max(full.data.rich$bathymetry), length.out=20),
                                   status=unique(full.data.rich$status))
predicted <- predict.gam(Model.richness$gam, newdata=effect.dat.richness, type='response')
effect.dat.richness <- cbind(effect.dat.richness, predicted)

status.rich <- effect.dat.richness%>%
  group_by(status)%>%
  summarise_at(vars(predicted), list(mean))

status.rich <- diff(range(status.rich$predicted))
status.rich

# Plot it 
effect.dat.richness <- expand.grid(cube.Aspect=seq(min(full.data.rich$cube.Aspect),max(full.data.rich$cube.Aspect),length.out = 20),
                                   bathymetry=seq(min(full.data.rich$bathymetry),max(full.data.rich$bathymetry), length.out=20),
                                   status=unique(full.data.rich$status))
predicted <- predict.gam(Model.richness$gam, newdata=effect.dat.richness, type='response', se.fit=T)
effect.dat.richness <- cbind(effect.dat.richness, predicted)

predicts.rich.status = effect.dat.richness%>%data.frame(predicted)%>%
  group_by(status)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

status.order <- c("F", "NT")
ggmod.gam.status<- ggplot(aes(x=status,y=response,colour=status), data=full.data.rich) +
  ylab("Predicted Species Richness")+
  xlab('Status')+
  #   ggtitle(substitute(italic(name)))+
  #scale_fill_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  scale_x_discrete(limits = status.order)+
  #geom_bar(stat = "identity")+
  geom_errorbar(data=predicts.rich.status,aes(ymin = response-se.fit,ymax = response+se.fit),colour="lightcoral",width = 0.5) +
  theme_classic()+
  Theme1
ggmod.gam.status

## Plot bayesian effects 
effect.dat.rich.sim <- effect.dat.richness
gg.richness <- t(posterior_predict(stangam.s.richness, effect.dat.rich.sim, re.form=NA))
#colnames(gg) <- paste("sim",1:ncol(gg),sep="_")
gg.rich <- cbind(effect.dat.rich.sim,gg.richness)

ncol(gg.rich)

gg.rich.long <- gg.rich%>%
  gather(sim, predicted.sim, 4:1005)

gg.rich.long$sim <- as.factor(gg.rich.long$sim)

# bathy effect 
bathy.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.mean")
colnames(bathy.sim) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.rich.long
  bathy.sim.2 <- filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- group_by(bathy.sim.2, bathymetry)
  bathy.sim.4 <- summarise_at(bathy.sim.3, vars(predicted.sim), list(mean))
  print(bathy.sim.5 <- diff(range(bathy.sim.4$predicted.sim)))
  bathy.sim <- rbind(bathy.sim, bathy.sim.5)
  
}

bathy.sim

write.csv(bathy.sim, "bayesian.bathy.richness.csv")

# aspect effect by status
aspect.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("aspect.mean")
colnames(aspect.sim) <- x

for(i in 1:1002){
  aspect.sim.1 <- gg.rich.long
  aspect.sim.2 <- filter(aspect.sim.1, sim==i)
  aspect.sim.3 <- group_by(aspect.sim.2, cube.Aspect)
  aspect.sim.4 <- summarise_at(aspect.sim.3, vars(predicted.sim), list(mean))
  print(aspect.sim.5 <- diff(range(aspect.sim.4$predicted.sim)))
  aspect.sim <- rbind(aspect.sim, aspect.sim.5)
  
}

aspect.sim
write.csv(aspect.sim, "bayesian.aspect.richness.csv")

# Status effect
status.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("status.mean")
colnames(status.sim) <- x

for(i in 1:1002){
  status.sim.1 <- gg.rich.long
  status.sim.2 <- filter(status.sim.1, sim==i)
  status.sim.3 <- group_by(status.sim.2, status)
  status.sim.4 <- summarise_at(status.sim.3, vars(predicted.sim), list(mean))
  print(status.sim.5 <- diff(range(status.sim.4$predicted.sim)))
  status.sim <- rbind(status.sim, status.sim.5)
  
}

status.sim

write.csv(status.sim, "bayesian.status.richness.csv")

## Create density plot
status.sim <- as.data.frame(status.sim[-1,])
x <- 'Status'
colnames(status.sim) <- x

bathy.sim <- as.data.frame(bathy.sim[-1,])
x <- 'Bathymetry (m)'
colnames(bathy.sim) <- x

aspect.sim <- as.data.frame(aspect.sim[-1,])
x <- 'Cube Aspect (degrees)'
colnames(aspect.sim) <- x

full.data <- as.data.frame(cbind(bathy.sim, aspect.sim, status.sim))

full.data.long <- full.data%>%
  gather(variable, predicted, 1:3)

colours <- c('#619CFF', '#00BA38', '#F8766D')
effect.plot <- ggplot(full.data.long, aes(x = predicted, fill=variable, colour=variable))+ 
  geom_density(alpha = 0.5)+
  scale_color_manual(values=colours)+
  scale_fill_manual(values=colours)+
  geom_vline(xintercept = 7.496373, color = "steelblue", size=0.75)+
  geom_vline(xintercept = 1.924486, color = "tomato3", size=0.75)+
  geom_vline(xintercept = 6.250768, color = "springgreen4", size=0.75)+
  labs(y="Density", x="Effect Size")+
  xlim(-1,15)+
  theme_classic()+
  Theme1
effect.plot


## Plot the predicted results for each of the variables 
# Bathy effect
bathy.means <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.means")
colnames(bathy.means) <- x

for(i in unique(gg.rich.long$bathymetry)){
  bathy.sim.1 <- gg.rich.long
  bathy.sim.2 <- dplyr::filter(bathy.sim.1, bathymetry==i)
  print(bathy.sim.3 <- dplyr::summarise_at(bathy.sim.2, vars(predicted.sim), list(mean)))
  bathy.means <- rbind(bathy.means, bathy.sim.3$predicted.sim)
}

bathy.means <- bathy.means[-1,]
depths <- unique(gg.rich$bathymetry)
bathy.means <- as.data.frame(cbind(bathy.means,depths))

# Get means for each value of bathymetry for each simulation 
bathy.mean.sim <- data.frame(matrix(ncol=1, nrow=20))
x <- c("bathy.mean.sim")
colnames(bathy.mean.sim) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.rich.long
  bathy.sim.2 <- filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- group_by(bathy.sim.2, bathymetry)
  print(bathy.sim.4 <- summarise_at(bathy.sim.3, vars(predicted.sim), list(mean)))
  bathy.mean.sim <- cbind(bathy.mean.sim, bathy.sim.4)
}

bathy.mean.full <- bathy.mean.sim[ ,seq(3, ncol(bathy.mean.sim), 2)]
bathy.mean.full$bathymetry <- bathy.mean.sim[,2]

bathy.mean.long <- bathy.mean.full%>%
  gather(sim, mean, 1:1002)


# Plot the results bathy
predict.plot.bathy <- ggplot() +
  geom_line(data=bathy.mean.long, aes(x = bathymetry, y = mean, group=sim), colour='lightblue'
            , alpha=0.2)+
  geom_line(data=bathy.means, aes(x = depths, y = bathy.means), colour='darkblue') +
  labs(y="Predicted Species Richness", x="Bathymetry (m)")+
  theme_classic()+
  xlim(-180,-51)+
  Theme1
predict.plot.bathy

# aspect effect
aspect.means <- data.frame(matrix(ncol=1, nrow=1))
x <- c("aspect.mean")
colnames(aspect.means) <- x

for(i in unique(gg.rich.long$cube.Aspect)){
  aspect.sim.1 <- gg.rich.long
  aspect.sim.2 <- dplyr::filter(aspect.sim.1, cube.Aspect==i)
  print(aspect.sim.3 <- dplyr::summarise_at(aspect.sim.2, vars(predicted.sim), list(mean)))
  aspect.means <- rbind(aspect.means, aspect.sim.3$predicted.sim)
}

aspect.means <- aspect.means[-1,]
aspect <- unique(gg.rich$cube.Aspect)
aspect.means <- as.data.frame(cbind(aspect.means,aspect))

# Get means for each value of aspect for each simulation 
aspect.mean.sim <- data.frame(matrix(ncol=1, nrow=20))
x <- c("aspect.mean.sim")
colnames(aspect.mean.sim) <- x

for(i in 1:1002){
  aspect.sim.1 <- gg.rich.long
  aspect.sim.2 <- filter(aspect.sim.1, sim==i)
  aspect.sim.3 <- group_by(aspect.sim.2, cube.Aspect)
  print(aspect.sim.4 <- summarise_at(aspect.sim.3, vars(predicted.sim), list(mean)))
  aspect.mean.sim <- cbind(aspect.mean.sim, aspect.sim.4)
}

aspect.mean.full <- aspect.mean.sim[ ,seq(3, ncol(aspect.mean.sim), 2)]
aspect.mean.full$cube.Aspect <- aspect.mean.sim[,2]

aspect.mean.long <- aspect.mean.full%>%
  gather(sim, mean, 1:1002)

# Plot the results aspect NT
predict.plot.aspect <- ggplot() +
  geom_line(data=aspect.mean.long, aes(x = cube.Aspect, y = mean, group=sim), colour='lightgreen',
            alpha=0.2) +
  geom_line(data=aspect.means, aes(x = aspect, y = aspect.means), colour='darkgreen') +
  labs(y="Predicted Species Richness", x="Cube Aspect (degrees)")+
  theme_classic()+
  Theme1
predict.plot.aspect


# status effect 
status.means <- data.frame(matrix(ncol=1, nrow=1))
x <- c("status.mean")
colnames(status.means) <- x

for(i in unique(gg.rich.long$status)){
  status.sim.1 <- gg.rich.long
  status.sim.2 <- dplyr::filter(status.sim.1, status==i)
  print(status.sim.3 <- dplyr::summarise_at(status.sim.2, vars(predicted.sim), list(mean)))
  status.means <- rbind(status.means, status.sim.3$predicted.sim)
}

status.means <- status.means[-1,]
status <- unique(gg.rich$status)
status.means <- as.data.frame(cbind(status.means,status))

# Get means for each value of aspect for each simulation 
status.mean.sim <- data.frame(matrix(ncol=1, nrow=20))
x <- c("status.mean.sim")
colnames(status.mean.sim) <- x

for(i in 1:1002){
  status.sim.1 <- gg.rich.long
  status.sim.2 <- filter(status.sim.1, sim==i)
  status.sim.3 <- group_by(status.sim.2, status)
  print(status.sim.4 <- summarise_at(status.sim.3, vars(predicted.sim), list(mean)))
  status.mean.sim <- cbind(status.mean.sim, status.sim.4)
}

status.mean.full <- status.mean.sim[ ,seq(3, ncol(status.mean.sim), 2)]
status.mean.full$status <- status.mean.sim[,2]

status.mean.long <- status.mean.full%>%
  gather(sim, mean, 1:1002)

# Plot the results status
jitter <- position_jitter(width = 0.06, height = 0.1)
predict.plot.status <- ggplot() +
  geom_point(data=status.mean.long,  position = jitter, aes(x = status, y = mean, group=sim), colour='lightcoral',
             alpha=0.2) +
  geom_boxplot(data=status.mean.long, aes(x=status, y=mean), colour='black', fill=NA, width=0.14)+
  geom_point(data=status.means, aes(x = status, y = status.means), colour='darkred', shape=17, size=3) +
  labs(y="PredictedSpecies Richness", x="Status")+
  #geom_jitter()+
  theme_classic()+
  Theme1
predict.plot.status

#### Species richness model 2 ####
# fit with shrinkage
effect.dat.richness <- expand.grid(cube.Aspect=seq(min(full.data.rich$cube.Aspect),max(full.data.rich$cube.Aspect),length.out = 20),
                                   bathymetry=seq(min(full.data.rich$bathymetry),max(full.data.rich$bathymetry), length.out=20))

Model.richness.2 <- uGamm(response~s(bathymetry,k=5,bs='cs')+
                        + s(cube.Aspect,k=5,bs='cs'),
                        random=~(1|site), 
                        family=poisson(), data=use.dat, lme4=TRUE)

gam.check(Model.richness.2$gam)

stangam.s.richness.2  <- stan_gamm4(response~s(bathymetry, k = 5, bs = "cs")
                                  + s(cube.Aspect, k=5,bs='cs'),
                                  random=~(1|site), adapt_delta = 0.99,
                                  data=use.dat, chains=3, cores=3, iter=41000, warmup=40000, thin=3,
                                  family=poisson)

## calcualte frequentist effects
# bathy effect 
predicted <- predict.gam(Model.richness.2$gam, newdata=effect.dat.richness, type='response')
effect.dat.richness <- cbind(effect.dat.richness, predicted)

bathy.P.rich <- effect.dat.richness%>%
  group_by(bathymetry)%>%
  summarise_at(vars(predicted), list(mean))

bathy.E.rich <- diff(range(bathy.P.rich$predicted))
bathy.E.rich

# Plot predicted
predicted <- predict.gam(Model.richness.2$gam, newdata=effect.dat.richness, type='response', se.fit=T)
effect.dat.richness <- cbind(effect.dat.richness, predicted)

predicts.legal.bathy = effect.dat.richness%>%data.frame(predicted)%>%
  group_by(bathymetry)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.legal.bathy<- ggplot() +
  ylab("Predicted Species Richness")+
  xlab('Bathymetry (m)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.data.rich,aes(x=bathymetry,y=response), colour="lightblue", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response), colour='darkblue', alpha=0.75)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response - se.fit), colour='darkblue', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response + se.fit), colour='darkblue', linetype="dashed",alpha=0.75)+
  theme_classic()+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.legal.bathy


# aspect by status effect 
effect.dat.richness <- expand.grid(cube.Aspect=seq(min(full.data.rich$cube.Aspect),max(full.data.rich$cube.Aspect),length.out = 20),
                                   bathymetry=seq(min(full.data.rich$bathymetry),max(full.data.rich$bathymetry), length.out=20),
                                   status=unique(full.data.rich$status))
predicted <- predict.gam(Model.richness.2$gam, newdata=effect.dat.richness, type='response')
effect.dat.richness <- cbind(effect.dat.richness, predicted)

aspect.P <- effect.dat.richness%>%
  group_by(cube.Aspect)%>%
  summarise_at(vars(predicted), list(mean))

aspect.E <- diff(range(aspect.P$predicted))
aspect.E


# Plot it 
effect.dat.richness <- expand.grid(cube.Aspect=seq(min(full.data.rich$cube.Aspect),max(full.data.rich$cube.Aspect),length.out = 20),
                                   bathymetry=seq(min(full.data.rich$bathymetry),max(full.data.rich$bathymetry), length.out=20))
predicted <- predict.gam(Model.richness.2$gam, newdata=effect.dat.richness, type='response', se.fit=T)
effect.dat.richness <- cbind(effect.dat.richness, predicted)

predicts.aspect = effect.dat.richness%>%data.frame(predicted)%>%
  group_by(cube.Aspect)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.aspect<- ggplot() +
  ylab("Predicted Species Richness")+
  xlab('Aspect (cubed)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.data.rich,aes(x=cube.Aspect,y=response), colour="lightgreen", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.aspect,aes(x=cube.Aspect,y=response), colour='darkgreen', alpha=0.75)+
  geom_line(data=predicts.aspect,aes(x=cube.Aspect,y=response - se.fit), colour='darkgreen', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.aspect,aes(x=cube.Aspect,y=response + se.fit), colour='darkgreen', linetype="dashed",alpha=0.75)+
  theme_classic()+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "NT",vjust = 2, hjust = -.1,size=5)
ggmod.aspect

## Plot bayesian effects 
effect.dat.rich.sim <- effect.dat.richness
gg.richness <- t(posterior_predict(stangam.s.richness.2, effect.dat.rich.sim, re.form=NA))
#colnames(gg) <- paste("sim",1:ncol(gg),sep="_")
gg.rich <- cbind(effect.dat.rich.sim,gg.richness)

ncol(gg.rich)

gg.rich.long <- gg.rich%>%
  gather(sim, predicted.sim, 3:1004)

gg.rich.long$sim <- as.factor(gg.rich.long$sim)

# bathy effect 
bathy.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.mean")
colnames(bathy.sim) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.rich.long
  bathy.sim.2 <- filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- group_by(bathy.sim.2, bathymetry)
  bathy.sim.4 <- summarise_at(bathy.sim.3, vars(predicted.sim), list(mean))
  print(bathy.sim.5 <- diff(range(bathy.sim.4$predicted.sim)))
  bathy.sim <- rbind(bathy.sim, bathy.sim.5)
  
}

bathy.sim

write.csv(bathy.sim, "bayesian.bathy.richness.csv")

# aspect effect by status
aspect.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("aspect.mean")
colnames(aspect.sim) <- x

for(i in 1:1002){
  aspect.sim.1 <- gg.rich.long
  aspect.sim.2 <- filter(aspect.sim.1, sim==i)
  aspect.sim.3 <- group_by(aspect.sim.2, cube.Aspect)
  aspect.sim.4 <- summarise_at(aspect.sim.3, vars(predicted.sim), list(mean))
  print(aspect.sim.5 <- diff(range(aspect.sim.4$predicted.sim)))
  aspect.sim <- rbind(aspect.sim, aspect.sim.5)
  
}

aspect.sim
write.csv(aspect.sim, "bayesian.aspect.richness.csv")


## Create density plot
bathy.sim <- as.data.frame(bathy.sim[-1,])
x <- 'Bathymetry (m)'
colnames(bathy.sim) <- x

aspect.sim <- as.data.frame(aspect.sim[-1,])
x <- 'Cube Aspect (degrees)'
colnames(aspect.sim) <- x

full.data <- as.data.frame(cbind(bathy.sim, aspect.sim))

full.data.long <- full.data%>%
  gather(variable, predicted, 1:2)

colours2 <- c('#619CFF', '#00BA38')
effect.plot <- ggplot(full.data.long, aes(x = predicted, fill=variable, colour=variable))+ 
  geom_density(alpha = 0.5)+
  scale_color_manual(values=colours2)+
  scale_fill_manual(values=colours2)+
  geom_vline(xintercept = 7.80625, color = "steelblue", size=0.75)+
  geom_vline(xintercept = 6.112132, color = "springgreen4", size=0.75)+
  labs(y="Density", x="Effect Size")+
  xlim(-1,18)+
  theme_classic()+
  Theme1+
  theme(legend.position=c(0.79, 0.8))
effect.plot


## Plot the predicted results for each of the variables 
# Bathy effect
bathy.means <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.means")
colnames(bathy.means) <- x

for(i in unique(gg.rich.long$bathymetry)){
  bathy.sim.1 <- gg.rich.long
  bathy.sim.2 <- dplyr::filter(bathy.sim.1, bathymetry==i)
  print(bathy.sim.3 <- dplyr::summarise_at(bathy.sim.2, vars(predicted.sim), list(mean)))
  bathy.means <- rbind(bathy.means, bathy.sim.3$predicted.sim)
}

bathy.means <- bathy.means[-1,]
depths <- unique(gg.rich$bathymetry)
bathy.means <- as.data.frame(cbind(bathy.means,depths))

# Get means for each value of bathymetry for each simulation 
bathy.mean.sim <- data.frame(matrix(ncol=1, nrow=20))
x <- c("bathy.mean.sim")
colnames(bathy.mean.sim) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.rich.long
  bathy.sim.2 <- filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- group_by(bathy.sim.2, bathymetry)
  print(bathy.sim.4 <- summarise_at(bathy.sim.3, vars(predicted.sim), list(mean)))
  bathy.mean.sim <- cbind(bathy.mean.sim, bathy.sim.4)
}

bathy.mean.full <- bathy.mean.sim[ ,seq(3, ncol(bathy.mean.sim), 2)]
bathy.mean.full$bathymetry <- bathy.mean.sim[,2]

bathy.mean.long <- bathy.mean.full%>%
  gather(sim, mean, 1:1002)


# Plot the results bathy
predict.plot.bathy <- ggplot() +
  geom_line(data=bathy.mean.long, aes(x = bathymetry, y = mean, group=sim), colour='lightblue'
            , alpha=0.2)+
  geom_line(data=bathy.means, aes(x = depths, y = bathy.means), colour='darkblue') +
  labs(y="Predicted Species Richness", x="Bathymetry (m)")+
  theme_classic()+
  xlim(-180,-51)+
  Theme1
predict.plot.bathy

# aspect effect
aspect.means <- data.frame(matrix(ncol=1, nrow=1))
x <- c("aspect.mean")
colnames(aspect.means) <- x

for(i in unique(gg.rich.long$cube.Aspect)){
  aspect.sim.1 <- gg.rich.long
  aspect.sim.2 <- dplyr::filter(aspect.sim.1, cube.Aspect==i)
  print(aspect.sim.3 <- dplyr::summarise_at(aspect.sim.2, vars(predicted.sim), list(mean)))
  aspect.means <- rbind(aspect.means, aspect.sim.3$predicted.sim)
}

aspect.means <- aspect.means[-1,]
aspect <- unique(gg.rich$cube.Aspect)
aspect.means <- as.data.frame(cbind(aspect.means,aspect))

# Get means for each value of aspect for each simulation 
aspect.mean.sim <- data.frame(matrix(ncol=1, nrow=20))
x <- c("aspect.mean.sim")
colnames(aspect.mean.sim) <- x

for(i in 1:1002){
  aspect.sim.1 <- gg.rich.long
  aspect.sim.2 <- filter(aspect.sim.1, sim==i)
  aspect.sim.3 <- group_by(aspect.sim.2, cube.Aspect)
  print(aspect.sim.4 <- summarise_at(aspect.sim.3, vars(predicted.sim), list(mean)))
  aspect.mean.sim <- cbind(aspect.mean.sim, aspect.sim.4)
}

aspect.mean.full <- aspect.mean.sim[ ,seq(3, ncol(aspect.mean.sim), 2)]
aspect.mean.full$cube.Aspect <- aspect.mean.sim[,2]

aspect.mean.long <- aspect.mean.full%>%
  gather(sim, mean, 1:1002)

# Plot the results aspect 
predict.plot.aspect <- ggplot() +
  geom_line(data=aspect.mean.long, aes(x = cube.Aspect, y = mean, group=sim), colour='lightgreen',
            alpha=0.2) +
  geom_line(data=aspect.means, aes(x = aspect, y = aspect.means), colour='darkgreen') +
  labs(y="Predicted Species Richness", x="Cube Aspect (degrees)")+
  theme_classic()+
  Theme1
predict.plot.aspect

#### Species Ricness model FOV 1 ####
# Bring in and format the data
name <- 'ningaloo' # for the study

setwd(d.dir)
richness.dat <-read.csv("ningaloo.total.abundance.and.species.richness.csv")
richness.dat <- richness.dat[116:230,c(1,3)]

variables <- read.csv("final.data.csv")

metadata <- read.csv('ningaloo_metadata.csv')

latlongs <- metadata%>%
  select('longitude', 'latitude')

variables <- cbind(variables, latlongs)
variables <- variables[,-c(5,6)]

variables <- variables%>%
  filter(model=="Legal")%>%
  mutate(sqrt.reef=sqrt(reef))%>%
  mutate(sqrt.slope=sqrt(Slope))%>%
  #mutate(sqrt.TPI=sqrt(TPI))%>%
  mutate(log.roughness=log(Roughness+1))%>%
  mutate(cube.Aspect=(Aspect)^3)%>%
  glimpse()

variables <- variables%>%
  dplyr::select('sample', 'site', 'bathymetry','distance.to.ramp', 'status','log.roughness', 'cube.Aspect', 'sqrt.slope', 
                 'latitude', 'longitude', "sqrt.reef", "mean.relief", "sd.relief")

full.data.rich <- merge(richness.dat,variables, all=TRUE)
full.data.rich[is.na(full.data.rich)] <- 0

full.data.rich <- full.data.rich%>%
  dplyr::rename(response=maxn)

# Remove NA values and sites in state no-take - for one of the predictor variables
# (can interpolate but I don't know how to do that right now...)
# 8.05 10.09 10.12 16.03 All have NAs

full.data.rich<-full.data.rich%>%
  filter(!sample%in%c("8.05","10.09","10.12","16.03", "10.14", "10.15", "10.13", "14.13", "14.12", "14.02"))

# Add row ID
full.data.rich <- full.data.rich%>%
  mutate(ID = rownames(full.data.rich))

#Set predictor variables 
pred.vars=c("sqrt.slope","cube.Aspect","log.roughness",
            "distance.to.ramp", "sd.relief", "mean.relief", "sqrt.reef")

## Run full subsets 
use.dat <- full.data.rich%>%
  dplyr::select(response, distance.to.ramp, cube.Aspect, sqrt.slope, log.roughness, status,
                bathymetry, site, sd.relief, mean.relief, sqrt.reef)
factor.vars <- c("status")

out.all <- list()
var.imp <- list()

Model.rich <- uGamm(response~s(bathymetry,k=5, bs='cr'),
                    family=poisson, random=~(1|site),
                    data=use.dat,
                    lme4=TRUE)


model.set <- generate.model.set(use.dat=use.dat,
                                test.fit=Model.rich,
                                pred.vars.cont=pred.vars,
                                pred.vars.fact=factor.vars,
                                factor.smooth.interactions=NA,
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

setwd(d.dir)
write.csv(all.mod.fits, "species.richness.models.fov.csv")
write.csv(all.var.imp, "species.richness.importance.fov.csv")

## Top models
# fit with shrinkage
effect.dat.richness <- expand.grid(sd.relief=seq(min(full.data.rich$sd.relief),max(full.data.rich$sd.relief),length.out = 20),
                                   bathymetry=seq(min(full.data.rich$bathymetry),max(full.data.rich$bathymetry), length.out=20))

Model.richness.fov.1 <- uGamm(response~s(bathymetry,k=5,bs='cs')
                        + s(sd.relief,k=5,bs='cs'),
                        random=~(1|site), 
                        family=poisson(), data=use.dat, lme4=TRUE)

gam.check(Model.richness.fov.1$gam)

stangam.s.richness.fov  <- stan_gamm4(response~s(bathymetry, k = 5, bs = "cs")
                                  + s(sd.relief, k=5,bs='cs'),
                                  random=~(1|site), adapt_delta = 0.99,
                                  data=use.dat, chains=3, cores=3, iter=41000, warmup=40000, thin=3,
                                  family=poisson)

## calcualte frequentist effects
# bathy effect 
predicted <- predict.gam(Model.richness.fov.1$gam, newdata=effect.dat.richness, type='response')
effect.dat.richness <- cbind(effect.dat.richness, predicted)

bathy.P.rich <- effect.dat.richness%>%
  group_by(bathymetry)%>%
  summarise_at(vars(predicted), list(mean))

bathy.E.rich <- diff(range(bathy.P.rich$predicted))
bathy.E.rich

# Plot predicted
predicted <- predict.gam(Model.richness.fov.1$gam, newdata=effect.dat.richness, type='response', se.fit=T)
effect.dat.richness <- cbind(effect.dat.richness, predicted)

predicts.legal.bathy = effect.dat.richness%>%data.frame(predicted)%>%
  group_by(bathymetry)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.legal.bathy<- ggplot() +
  ylab("Predicted Species Richness")+
  xlab('Bathymetry (m)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.data.rich,aes(x=bathymetry,y=response), colour="lightblue", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response), colour='darkblue', alpha=0.75)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response - se.fit), colour='darkblue', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response + se.fit), colour='darkblue', linetype="dashed",alpha=0.75)+
  theme_classic()+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.legal.bathy


# relief by status effect 
effect.dat.richness <- expand.grid(sd.relief=seq(min(full.data.rich$sd.relief),max(full.data.rich$sd.relief),length.out = 20),
                                   bathymetry=seq(min(full.data.rich$bathymetry),max(full.data.rich$bathymetry), length.out=20))
predicted <- predict.gam(Model.richness.fov.1$gam, newdata=effect.dat.richness, type='response')
effect.dat.richness <- cbind(effect.dat.richness, predicted)

relief.P <- effect.dat.richness%>%
  group_by(sd.relief)%>%
  summarise_at(vars(predicted), list(mean))

relief.E <- diff(range(relief.P$predicted))
relief.E


# Plot it 
predicted <- predict.gam(Model.richness.fov.1$gam, newdata=effect.dat.richness, type='response', se.fit=T)
effect.dat.richness <- cbind(effect.dat.richness, predicted)

predicts.relief = effect.dat.richness%>%data.frame(predicted)%>%
  group_by(sd.relief)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.relief<- ggplot() +
  ylab("Predicted Species Richness")+
  xlab('relief (cubed)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.data.rich,aes(x=sd.relief,y=response), colour="lightgreen", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.relief,aes(x=sd.relief,y=response), colour='darkgreen', alpha=0.75)+
  geom_line(data=predicts.relief,aes(x=sd.relief,y=response - se.fit), colour='darkgreen', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.relief,aes(x=sd.relief,y=response + se.fit), colour='darkgreen', linetype="dashed",alpha=0.75)+
  theme_classic()+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "NT",vjust = 2, hjust = -.1,size=5)
ggmod.relief

## Plot bayesian effects 
effect.dat.rich.sim <- effect.dat.richness
gg.richness <- t(posterior_predict(stangam.s.richness.fov, effect.dat.rich.sim, re.form=NA))
#colnames(gg) <- paste("sim",1:ncol(gg),sep="_")
gg.rich <- cbind(effect.dat.rich.sim,gg.richness)

ncol(gg.rich)

gg.rich.long <- gg.rich%>%
  gather(sim, predicted.sim, 3:1004)

gg.rich.long$sim <- as.factor(gg.rich.long$sim)

# bathy effect 
bathy.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.mean")
colnames(bathy.sim) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.rich.long
  bathy.sim.2 <- filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- group_by(bathy.sim.2, bathymetry)
  bathy.sim.4 <- summarise_at(bathy.sim.3, vars(predicted.sim), list(mean))
  print(bathy.sim.5 <- diff(range(bathy.sim.4$predicted.sim)))
  bathy.sim <- rbind(bathy.sim, bathy.sim.5)
  
}

bathy.sim

write.csv(bathy.sim, "bayesian.bathy.richness.csv")

# relief effect 
relief.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("relief.mean")
colnames(relief.sim) <- x

for(i in 1:1002){
  relief.sim.1 <- gg.rich.long
  relief.sim.2 <- filter(relief.sim.1, sim==i)
  relief.sim.3 <- group_by(relief.sim.2, sd.relief)
  relief.sim.4 <- summarise_at(relief.sim.3, vars(predicted.sim), list(mean))
  print(relief.sim.5 <- diff(range(relief.sim.4$predicted.sim)))
  relief.sim <- rbind(relief.sim, relief.sim.5)
  
}

relief.sim
write.csv(aspect.sim, "bayesian.aspect.richness.csv")

## Create density plot
bathy.sim <- as.data.frame(bathy.sim[-1,])
x <- 'Bathymetry (m)'
colnames(bathy.sim) <- x

relief.sim <- as.data.frame(relief.sim[-1,])
x <- 'SD Relief'
colnames(relief.sim) <- x

full.data <- as.data.frame(cbind(bathy.sim, relief.sim))

full.data.long <- full.data%>%
  gather(variable, predicted, 1:2)

colours2 <- c('#619CFF', '#00BA38')
effect.plot <- ggplot(full.data.long, aes(x = predicted, fill=variable, colour=variable))+ 
  geom_density(alpha = 0.5)+
  scale_color_manual(values=colours2)+
  scale_fill_manual(values=colours2)+
  geom_vline(xintercept = 9.047447, color = "steelblue", size=0.75)+
  geom_vline(xintercept = 15.81526, color = "springgreen4", size=0.75)+
  labs(y="Density", x="Effect Size")+
  xlim(-1,35)+
  theme_classic()+
  Theme1
effect.plot


## Plot the predicted results for each of the variables 
# Bathy effect
bathy.means <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.means")
colnames(bathy.means) <- x

for(i in unique(gg.rich.long$bathymetry)){
  bathy.sim.1 <- gg.rich.long
  bathy.sim.2 <- dplyr::filter(bathy.sim.1, bathymetry==i)
  print(bathy.sim.3 <- dplyr::summarise_at(bathy.sim.2, vars(predicted.sim), list(mean)))
  bathy.means <- rbind(bathy.means, bathy.sim.3$predicted.sim)
}

bathy.means <- bathy.means[-1,]
depths <- unique(gg.rich$bathymetry)
bathy.means <- as.data.frame(cbind(bathy.means,depths))

# Get means for each value of bathymetry for each simulation 
bathy.mean.sim <- data.frame(matrix(ncol=1, nrow=20))
x <- c("bathy.mean.sim")
colnames(bathy.mean.sim) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.rich.long
  bathy.sim.2 <- filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- group_by(bathy.sim.2, bathymetry)
  print(bathy.sim.4 <- summarise_at(bathy.sim.3, vars(predicted.sim), list(mean)))
  bathy.mean.sim <- cbind(bathy.mean.sim, bathy.sim.4)
}

bathy.mean.full <- bathy.mean.sim[ ,seq(3, ncol(bathy.mean.sim), 2)]
bathy.mean.full$bathymetry <- bathy.mean.sim[,2]

bathy.mean.long <- bathy.mean.full%>%
  gather(sim, mean, 1:1002)


# Plot the results bathy
predict.plot.bathy <- ggplot() +
  geom_line(data=bathy.mean.long, aes(x = bathymetry, y = mean, group=sim), colour='lightblue'
            , alpha=0.2)+
  geom_line(data=bathy.means, aes(x = depths, y = bathy.means), colour='darkblue') +
  labs(y="Predicted Species Richness", x="Bathymetry (m)")+
  theme_classic()+
  xlim(-180,-51)+
  Theme1
predict.plot.bathy

# relief effect
relief.means <- data.frame(matrix(ncol=1, nrow=1))
x <- c("relief.mean")
colnames(relief.means) <- x

for(i in unique(gg.rich.long$sd.relief)){
  relief.sim.1 <- gg.rich.long
  relief.sim.2 <- dplyr::filter(relief.sim.1, sd.relief==i)
  print(relief.sim.3 <- dplyr::summarise_at(relief.sim.2, vars(predicted.sim), list(mean)))
  relief.means <- rbind(relief.means, relief.sim.3$predicted.sim)
}

relief.means <- relief.means[-1,]
relief <- unique(gg.rich$sd.relief)
relief.means <- as.data.frame(cbind(relief.means,relief))

# Get means for each value of sd.relief for each simulation 
relief.mean.sim <- data.frame(matrix(ncol=1, nrow=20))
x <- c("relief.mean.sim")
colnames(relief.mean.sim) <- x

for(i in 1:1002){
  relief.sim.1 <- gg.rich.long
  relief.sim.2 <- filter(relief.sim.1, sim==i)
  relief.sim.3 <- group_by(relief.sim.2, sd.relief)
  print(relief.sim.4 <- summarise_at(relief.sim.3, vars(predicted.sim), list(mean)))
  relief.mean.sim <- cbind(relief.mean.sim, relief.sim.4)
}

relief.mean.full <- relief.mean.sim[ ,seq(3, ncol(relief.mean.sim), 2)]
relief.mean.full$sd.relief <- relief.mean.sim[,2]

relief.mean.long <- relief.mean.full%>%
  gather(sim, mean, 1:1002)

# Plot the results relief
predict.plot.relief <- ggplot() +
  geom_line(data=relief.mean.long, aes(x = sd.relief, y = mean, group=sim), colour='lightgreen',
            alpha=0.2) +
  geom_line(data=relief.means, aes(x = relief, y = relief.means), colour='darkgreen') +
  labs(y="Predicted Species Richness", x="SD Relief")+
  theme_classic()+
  Theme1
predict.plot.relief






#### Species Richness model FOV 2 ####
## Top models
# fit with shrinkage
effect.dat.richness <- expand.grid(sd.relief=seq(min(full.data.rich$sd.relief),max(full.data.rich$sd.relief),length.out = 20),
                                   bathymetry=seq(min(full.data.rich$bathymetry),max(full.data.rich$bathymetry), length.out=20),
                                   status=unique(full.data.rich$status))

Model.richness.fov.2 <- uGamm(response~s(bathymetry,k=5,bs='cs')+ factor(status)
                        + s(sd.relief,k=5,bs='cs'),
                        random=~(1|site), 
                        family=poisson(), data=use.dat, lme4=TRUE)

gam.check(Model.richness.fov.2$gam)

stangam.s.richness.fov.2  <- stan_gamm4(response~s(bathymetry, k = 5, bs = "cs") + status
                                  + s(sd.relief, k=5,bs='cs'),
                                  random=~(1|site), adapt_delta = 0.99,
                                  data=use.dat, chains=3, cores=3, iter=41000, warmup=40000, thin=3,
                                  family=poisson)

## calcualte frequentist effects
# bathy effect 
predicted <- predict.gam(Model.richness.fov.2$gam, newdata=effect.dat.richness, type='response')
effect.dat.richness <- cbind(effect.dat.richness, predicted)

bathy.P.rich <- effect.dat.richness%>%
  group_by(bathymetry)%>%
  summarise_at(vars(predicted), list(mean))

bathy.E.rich <- diff(range(bathy.P.rich$predicted))
bathy.E.rich

# Plot predicted
predicted <- predict.gam(Model.richness.fov.2$gam, newdata=effect.dat.richness, type='response', se.fit=T)
effect.dat.richness <- cbind(effect.dat.richness, predicted)

predicts.legal.bathy = effect.dat.richness%>%data.frame(predicted)%>%
  group_by(bathymetry)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.legal.bathy<- ggplot() +
  ylab("Predicted Species Richness")+
  xlab('Bathymetry (m)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.data.rich,aes(x=bathymetry,y=response), colour="lightblue", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response), colour='darkblue', alpha=0.75)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response - se.fit), colour='darkblue', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response + se.fit), colour='darkblue', linetype="dashed",alpha=0.75)+
  theme_classic()+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.legal.bathy


# relief effect 
effect.dat.richness <- expand.grid(sd.relief=seq(min(full.data.rich$sd.relief),max(full.data.rich$sd.relief),length.out = 20),
                                   bathymetry=seq(min(full.data.rich$bathymetry),max(full.data.rich$bathymetry), length.out=20),
                                   status=unique(full.data.rich$status))
predicted <- predict.gam(Model.richness.fov.2$gam, newdata=effect.dat.richness, type='response')
effect.dat.richness <- cbind(effect.dat.richness, predicted)

relief.P <- effect.dat.richness%>%
  group_by(sd.relief)%>%
  summarise_at(vars(predicted), list(mean))

relief.E <- diff(range(relief.P$predicted))
relief.E


# Plot it 
predicted <- predict.gam(Model.richness.fov.2$gam, newdata=effect.dat.richness, type='response', se.fit=T)
effect.dat.richness <- cbind(effect.dat.richness, predicted)

predicts.relief = effect.dat.richness%>%data.frame(predicted)%>%
  group_by(sd.relief)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.relief<- ggplot() +
  ylab("Predicted Species Richness")+
  xlab('relief (cubed)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.data.rich,aes(x=sd.relief,y=response), colour="lightgreen", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.relief,aes(x=sd.relief,y=response), colour='darkgreen', alpha=0.75)+
  geom_line(data=predicts.relief,aes(x=sd.relief,y=response - se.fit), colour='darkgreen', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.relief,aes(x=sd.relief,y=response + se.fit), colour='darkgreen', linetype="dashed",alpha=0.75)+
  theme_classic()+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "NT",vjust = 2, hjust = -.1,size=5)
ggmod.relief

# status effect 
effect.dat.richness <- expand.grid(sd.relief=seq(min(full.data.rich$sd.relief),max(full.data.rich$sd.relief),length.out = 20),
                                   bathymetry=seq(min(full.data.rich$bathymetry),max(full.data.rich$bathymetry), length.out=20),
                                   status=unique(full.data.rich$status))
predicted <- predict.gam(Model.richness.fov.2$gam, newdata=effect.dat.richness, type='response')
effect.dat.richness <- cbind(effect.dat.richness, predicted)

status.rich <- effect.dat.richness%>%
  group_by(status)%>%
  summarise_at(vars(predicted), list(mean))

status.rich <- diff(range(status.rich$predicted))
status.rich

# Plot it 
effect.dat.richness <- expand.grid(sd.relief=seq(min(full.data.rich$sd.relief),max(full.data.rich$sd.relief),length.out = 20),
                                   bathymetry=seq(min(full.data.rich$bathymetry),max(full.data.rich$bathymetry), length.out=20),
                                   status=unique(full.data.rich$status))
predicted <- predict.gam(Model.richness.fov.2$gam, newdata=effect.dat.richness, type='response', se.fit=T)
effect.dat.richness <- cbind(effect.dat.richness, predicted)

predicts.rich.status = effect.dat.richness%>%data.frame(predicted)%>%
  group_by(status)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

status.order <- c("F", "NT")
ggmod.gam.status<- ggplot(aes(x=status,y=response,colour=status), data=full.data.rich) +
  ylab("Predicted Species Richness")+
  xlab('Status')+
  #   ggtitle(substitute(italic(name)))+
  #scale_fill_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  scale_x_discrete(limits = status.order)+
  #geom_bar(stat = "identity")+
  geom_errorbar(data=predicts.rich.status,aes(ymin = response-se.fit,ymax = response+se.fit),colour="lightcoral",width = 0.5) +
  theme_classic()+
  Theme1
ggmod.gam.status

## Plot bayesian effects 
effect.dat.rich.sim <- effect.dat.richness
gg.richness <- t(posterior_predict(stangam.s.richness.fov.2, effect.dat.rich.sim, re.form=NA))
#colnames(gg) <- paste("sim",1:ncol(gg),sep="_")
gg.rich <- cbind(effect.dat.rich.sim,gg.richness)

ncol(gg.rich)

gg.rich.long <- gg.rich%>%
  gather(sim, predicted.sim, 4:1005)

gg.rich.long$sim <- as.factor(gg.rich.long$sim)

# bathy effect 
bathy.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.mean")
colnames(bathy.sim) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.rich.long
  bathy.sim.2 <- filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- group_by(bathy.sim.2, bathymetry)
  bathy.sim.4 <- summarise_at(bathy.sim.3, vars(predicted.sim), list(mean))
  print(bathy.sim.5 <- diff(range(bathy.sim.4$predicted.sim)))
  bathy.sim <- rbind(bathy.sim, bathy.sim.5)
  
}

bathy.sim

write.csv(bathy.sim, "bayesian.bathy.richness.csv")

# relief effect by status
relief.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("relief.mean")
colnames(relief.sim) <- x

for(i in 1:1002){
  relief.sim.1 <- gg.rich.long
  relief.sim.2 <- filter(relief.sim.1, sim==i)
  relief.sim.3 <- group_by(relief.sim.2, sd.relief)
  relief.sim.4 <- summarise_at(relief.sim.3, vars(predicted.sim), list(mean))
  print(relief.sim.5 <- diff(range(relief.sim.4$predicted.sim)))
  relief.sim <- rbind(relief.sim, relief.sim.5)
  
}

relief.sim
write.csv(relief.sim, "bayesian.relief.richness.csv")

# Status effect
status.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("status.mean")
colnames(status.sim) <- x

for(i in 1:1002){
  status.sim.1 <- gg.rich.long
  status.sim.2 <- filter(status.sim.1, sim==i)
  status.sim.3 <- group_by(status.sim.2, status)
  status.sim.4 <- summarise_at(status.sim.3, vars(predicted.sim), list(mean))
  print(status.sim.5 <- diff(range(status.sim.4$predicted.sim)))
  status.sim <- rbind(status.sim, status.sim.5)
  
}

status.sim

write.csv(status.sim, "bayesian.status.richness.csv")

## Create density plot
status.sim <- as.data.frame(status.sim[-1,])
x <- 'Status'
colnames(status.sim) <- x

bathy.sim <- as.data.frame(bathy.sim[-1,])
x <- 'Bathymetry (m)'
colnames(bathy.sim) <- x

relief.sim <- as.data.frame(relief.sim[-1,])
x <- 'SD Relief (degrees)'
colnames(relief.sim) <- x

full.data <- as.data.frame(cbind(bathy.sim, relief.sim, status.sim))

full.data.long <- full.data%>%
  gather(variable, predicted, 1:3)

colours <- c('#619CFF', '#00BA38', '#F8766D')
effect.plot <- ggplot(full.data.long, aes(x = predicted, fill=variable, colour=variable))+ 
  geom_density(alpha = 0.5)+
  scale_color_manual(values=colours)+
  scale_fill_manual(values=colours)+
  geom_vline(xintercept = 8.96545, color = "steelblue", size=0.75)+
  geom_vline(xintercept = 0.6136397, color = "tomato3", size=0.75)+
  geom_vline(xintercept = 15.55588, color = "springgreen4", size=0.75)+
  labs(y="Density", x="Effect Size")+
  xlim(-2,33)+
  theme_classic()+
  Theme1
effect.plot


## Plot the predicted results for each of the variables 
# Bathy effect
bathy.means <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.means")
colnames(bathy.means) <- x

for(i in unique(gg.rich.long$bathymetry)){
  bathy.sim.1 <- gg.rich.long
  bathy.sim.2 <- dplyr::filter(bathy.sim.1, bathymetry==i)
  print(bathy.sim.3 <- dplyr::summarise_at(bathy.sim.2, vars(predicted.sim), list(mean)))
  bathy.means <- rbind(bathy.means, bathy.sim.3$predicted.sim)
}

bathy.means <- bathy.means[-1,]
depths <- unique(gg.rich$bathymetry)
bathy.means <- as.data.frame(cbind(bathy.means,depths))

# Get means for each value of bathymetry for each simulation 
bathy.mean.sim <- data.frame(matrix(ncol=1, nrow=20))
x <- c("bathy.mean.sim")
colnames(bathy.mean.sim) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.rich.long
  bathy.sim.2 <- filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- group_by(bathy.sim.2, bathymetry)
  print(bathy.sim.4 <- summarise_at(bathy.sim.3, vars(predicted.sim), list(mean)))
  bathy.mean.sim <- cbind(bathy.mean.sim, bathy.sim.4)
}

bathy.mean.full <- bathy.mean.sim[ ,seq(3, ncol(bathy.mean.sim), 2)]
bathy.mean.full$bathymetry <- bathy.mean.sim[,2]

bathy.mean.long <- bathy.mean.full%>%
  gather(sim, mean, 1:1002)


# Plot the results bathy
predict.plot.bathy <- ggplot() +
  geom_line(data=bathy.mean.long, aes(x = bathymetry, y = mean, group=sim), colour='lightblue'
            , alpha=0.2)+
  geom_line(data=bathy.means, aes(x = depths, y = bathy.means), colour='darkblue') +
  labs(y="Predicted Species Richness", x="Bathymetry (m)")+
  theme_classic()+
  xlim(-180,-51)+
  Theme1
predict.plot.bathy

# relief effect
relief.means <- data.frame(matrix(ncol=1, nrow=1))
x <- c("relief.mean")
colnames(relief.means) <- x

for(i in unique(gg.rich.long$sd.relief)){
  relief.sim.1 <- gg.rich.long
  relief.sim.2 <- dplyr::filter(relief.sim.1, sd.relief==i)
  print(relief.sim.3 <- dplyr::summarise_at(relief.sim.2, vars(predicted.sim), list(mean)))
  relief.means <- rbind(relief.means, relief.sim.3$predicted.sim)
}

relief.means <- relief.means[-1,]
relief <- unique(gg.rich$sd.relief)
relief.means <- as.data.frame(cbind(relief.means,relief))

# Get means for each value of relief for each simulation 
relief.mean.sim <- data.frame(matrix(ncol=1, nrow=20))
x <- c("relief.mean.sim")
colnames(relief.mean.sim) <- x

for(i in 1:1002){
  relief.sim.1 <- gg.rich.long
  relief.sim.2 <- filter(relief.sim.1, sim==i)
  relief.sim.3 <- group_by(relief.sim.2, sd.relief)
  print(relief.sim.4 <- summarise_at(relief.sim.3, vars(predicted.sim), list(mean)))
  relief.mean.sim <- cbind(relief.mean.sim, relief.sim.4)
}

relief.mean.full <- relief.mean.sim[ ,seq(3, ncol(relief.mean.sim), 2)]
relief.mean.full$sd.relief <- relief.mean.sim[,2]

relief.mean.long <- relief.mean.full%>%
  gather(sim, mean, 1:1002)

# Plot the results relief NT
predict.plot.relief <- ggplot() +
  geom_line(data=relief.mean.long, aes(x = sd.relief, y = mean, group=sim), colour='lightgreen',
            alpha=0.2) +
  geom_line(data=relief.means, aes(x = relief, y = relief.means), colour='darkgreen') +
  labs(y="Predicted Species Richness", x="SD Relief")+
  theme_classic()+
  Theme1
predict.plot.relief


# status effect 
status.means <- data.frame(matrix(ncol=1, nrow=1))
x <- c("status.mean")
colnames(status.means) <- x

for(i in unique(gg.rich.long$status)){
  status.sim.1 <- gg.rich.long
  status.sim.2 <- dplyr::filter(status.sim.1, status==i)
  print(status.sim.3 <- dplyr::summarise_at(status.sim.2, vars(predicted.sim), list(mean)))
  status.means <- rbind(status.means, status.sim.3$predicted.sim)
}

status.means <- status.means[-1,]
status <- unique(gg.rich$status)
status.means <- as.data.frame(cbind(status.means,status))

# Get means for each value of relief for each simulation 
status.mean.sim <- data.frame(matrix(ncol=1, nrow=20))
x <- c("status.mean.sim")
colnames(status.mean.sim) <- x

for(i in 1:1002){
  status.sim.1 <- gg.rich.long
  status.sim.2 <- filter(status.sim.1, sim==i)
  status.sim.3 <- group_by(status.sim.2, status)
  print(status.sim.4 <- summarise_at(status.sim.3, vars(predicted.sim), list(mean)))
  status.mean.sim <- cbind(status.mean.sim, status.sim.4)
}

status.mean.full <- status.mean.sim[ ,seq(3, ncol(status.mean.sim), 2)]
status.mean.full$status <- status.mean.sim[,2]

status.mean.long <- status.mean.full%>%
  gather(sim, mean, 1:1002)

# Plot the results status
jitter <- position_jitter(width = 0.06, height = 0.1)
predict.plot.status <- ggplot() +
  geom_point(data=status.mean.long,  position = jitter, aes(x = status, y = mean, group=sim), colour='lightcoral',
             alpha=0.2) +
  geom_boxplot(data=status.mean.long, aes(x=status, y=mean), colour='black', fill=NA, width=0.14)+
  geom_point(data=status.means, aes(x = status, y = status.means), colour='darkred', shape=17, size=3) +
  labs(y="PredictedSpecies Richness", x="Status")+
  #geom_jitter()+
  theme_classic()+
  Theme1
predict.plot.status

