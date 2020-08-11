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

############### Part 3 Plotting the most parsimonious models #################

library(gridExtra)
library(grid)
# Theme-
Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend.background = element_rect(fill="white"),
    legend.background = element_blank(),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=15),
    legend.title = element_blank(),
    legend.position = c(0.2, 0.8),
    text=element_text(size=15),
    strip.text.y = element_text(size = 15,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=15),
    axis.title.y=element_text(vjust=0.6, angle=90, size=15),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())

######### Predict and plot legal fish #########
dat.legal<-dat%>%
  dplyr::filter(model=="Legal")%>%
  dplyr::rename(response=target.fish)

gamm.legal<-gam(response ~ 
                  s(distance.to.ramp, k=3, bs="cr") + 
                  status +
                  # s(log.roughness, k=3, bs="cr") +
                  # s(sqrt.slope, k=3, bs="cr") + 
                  s(cube.Aspect, k=3, bs="cr") +
                  te(distance.to.60,bathymetry, k=3, bs="cr") +
                  s(bathymetry, bs="re") + s(TPI, bs="re") + s(zone, bs="re"),
                  family=tw(), data=dat.legal)


# predict cube.Aspect from model
mod<-gamm.legal
testdata <- expand.grid(cube.Aspect=seq(min(dat$cube.Aspect),max(dat$cube.Aspect),length.out = 5),
                        distance.to.ramp=seq(min(dat$distance.to.ramp),max(dat$distance.to.ramp),length.out=5),
                        bathymetry=seq(min(dat$bathymetry), max(dat$bathymetry), length.out=5),
                        distance.to.60=seq(min(dat$distance.to.60), max(dat$distance.to.60), length.out=5),
                        TPI=seq(min(dat$TPI), max(dat$TPI), length.out=5),
                        zone=(mod$model$zone),
                        status=(mod$model$status))%>%
                        
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.cube.Aspect <- testdata%>%data.frame(fits)%>%
  group_by(cube.Aspect)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.legal.bathy,"predict.legal.bathy.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.bathy<-read.csv("predict.legal.bathy.csv")%>%
  glimpse()

# Plot cube aspect from model
ggmod.legal.cube.Aspect<- ggplot() +
  ylab("Predicted Abundance of Legal Sized Fish")+
  xlab('Aspect (cubed)')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=dat.legal,aes(x=cube.Aspect,y=response,colour=status),  alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.legal.cube.Aspect,aes(x=cube.Aspect,y=response),alpha=0.75)+
  geom_line(data=predicts.legal.cube.Aspect,aes(x=cube.Aspect,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.cube.Aspect,aes(x=cube.Aspect,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.legal.cube.Aspect

#predict distance to ramp
mod<-gamm.legal
testdata <- expand.grid(cube.Aspect=seq(min(dat$cube.Aspect),max(dat$cube.Aspect),length.out = 5),
                        distance.to.ramp=seq(min(dat$distance.to.ramp),max(dat$distance.to.ramp),length.out=5),
                        bathymetry=seq(min(dat$bathymetry), max(dat$bathymetry), length.out=5),
                        distance.to.60=seq(min(dat$distance.to.60), max(dat$distance.to.60), length.out=5),
                        TPI=seq(min(dat$TPI), max(dat$TPI), length.out=5),
                        zone=(mod$model$zone),
                        status=(mod$model$status))%>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.ramp <- testdata%>%data.frame(fits)%>%
  group_by(distance.to.ramp)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.legal.bathy,"predict.legal.bathy.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.bathy<-read.csv("predict.legal.bathy.csv")%>%
  glimpse()


# Plot Distance to ramps
ggmod.legal.ramps<- ggplot() +
  ylab("Predicted Abundance of Legal Sized Fish")+
  xlab('Distance to Ramp')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=dat.legal,aes(x=distance.to.ramp,y=response,colour=status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.legal.ramp,aes(x=distance.to.ramp,y=response),alpha=0.75)+
  geom_line(data=predicts.legal.ramp,aes(x=distance.to.ramp,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.ramp,aes(x=distance.to.ramp,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(d)",vjust = 1, hjust = -.1,size=5)
ggmod.legal.ramps

#predict status
mod<-gamm.legal
testdata <- expand.grid(cube.Aspect=seq(min(dat$cube.Aspect),max(dat$cube.Aspect),length.out = 5),
                        distance.to.ramp=seq(min(dat$distance.to.ramp),max(dat$distance.to.ramp),length.out=5),
                        bathymetry=seq(min(dat$bathymetry), max(dat$bathymetry), length.out=5),
                        distance.to.60=seq(min(dat$distance.to.60), max(dat$distance.to.60), length.out=5),
                        TPI=seq(min(dat$TPI), max(dat$TPI), length.out=5),
                        zone=(mod$model$zone),
                        status=(mod$model$status))%>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.status <- testdata%>%data.frame(fits)%>%
  group_by(status)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.legal.bathy,"predict.legal.bathy.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.bathy<-read.csv("predict.legal.bathy.csv")%>%
  glimpse()

# Plot status
ggmod.legal.status<- ggplot(aes(x=status,y=response,fill=status,colour=status), data=predicts.legal.status) +
  ylab("Predicted Abundance of Legal Sized Fish")+
  xlab('Status')+
  #   ggtitle(substitute(italic(name)))+
  scale_fill_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  scale_colour_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  scale_x_discrete(limits = rev(levels(predicts.legal.status$status)))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  theme_classic()+
  Theme1
ggmod.legal.status


######## Predict and plot Sublegal ramp and slope ######
dat.sublegal<-dat%>%
  filter(model=="Sublegal")%>%
  dplyr::rename(response=target.fish)

gamm.sublegal<-gam(response ~ 
                     s(distance.to.ramp, k=3, bs="cr") + 
                     # status +
                     # s(log.roughness, k=3, bs="cr") +
                     s(sqrt.slope, k=3, bs="cr") + 
                     # s(cube.Aspect, k=3, bs="cr") +
                     s(bathymetry, bs="re") + s(TPI, bs="re") + s(zone, bs="re") +
                     te(distance.to.60,bathymetry, k=3, bs="cr"), family=tw(), data=dat.sublegal)


#predict slope
mod<-gamm.sublegal
testdata <- expand.grid(sqrt.slope=seq(min(dat$sqrt.slope),max(dat$sqrt.slope),length.out = 5),
                        distance.to.ramp=seq(min(dat$distance.to.ramp),max(dat$distance.to.ramp),length.out=5),
                        bathymetry=seq(min(dat$bathymetry), max(dat$bathymetry), length.out=5),
                        distance.to.60=seq(min(dat$distance.to.60), max(dat$distance.to.60), length.out=5),
                        TPI=seq(min(dat$TPI), max(dat$TPI), length.out=5),
                        zone=(mod$model$zone))%>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.sublegal.slope <- testdata%>%data.frame(fits)%>%
  group_by(sqrt.slope)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()


# Plot Bathymetry
ggmod.sublegal.slope<- ggplot() +
  ylab("Predicted Abundance of Sublegal Sized Fish")+
  xlab('Slope (sqrt)')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  #   geom_jitter(width = 0.25,height = 0)+
  geom_point(data=dat.sublegal,aes(x=sqrt.slope,y=response,colour=status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.sublegal.slope,aes(x=sqrt.slope,y=response),alpha=0.75)+
  geom_line(data=predicts.sublegal.slope,aes(x=sqrt.slope,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.sublegal.slope,aes(x=sqrt.slope,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.sublegal.slope


#Predict distance to ramp 
mod<-gamm.sublegal
testdata <- expand.grid(sqrt.slope=seq(min(dat$sqrt.slope),max(dat$sqrt.slope),length.out = 5),
                        distance.to.ramp=seq(min(dat$distance.to.ramp),max(dat$distance.to.ramp),length.out=5),
                        bathymetry=seq(min(dat$bathymetry), max(dat$bathymetry), length.out=5),
                        distance.to.60=seq(min(dat$distance.to.60), max(dat$distance.to.60), length.out=5),
                        TPI=seq(min(dat$TPI), max(dat$TPI), length.out=5),
                        zone=(mod$model$zone))%>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.sublegal.ramp = testdata%>%data.frame(fits)%>%
  group_by(distance.to.ramp)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.legal.ramp,"predict.sublegal.ramp.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.ramp<-read.csv("predict.sublegal.ramp.csv")%>%
  glimpse()

# Distance to ramps
ggmod.sublegal.ramps<- ggplot() +
  ylab("Predicted Abundance of Sublegal Sized Fish")+
  xlab('Distance to Ramps')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  #   geom_jitter(width = 0.25,height = 0)+
  geom_point(data=dat.sublegal,aes(x=distance.to.ramp,y=response,colour=status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.sublegal.ramp,aes(x=distance.to.ramp,y=response),alpha=0.5)+
  geom_line(data=predicts.sublegal.ramp,aes(x=distance.to.ramp,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.sublegal.ramp,aes(x=distance.to.ramp,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.sublegal.ramps


# combined.plot using grid() and gridExtra()------
blank <- grid.rect(gp=gpar(col="white"))

# To see what they will look like use grid.arrange() - make sure Plot window is large enough! - or will error!
grid.arrange(ggmod.sublegal.bathy,ggmod.sublegal.relief)

# Use arrangeGrob ONLY - as we can pass this to ggsave! Note use of raw ggplot's
combine.plot<-arrangeGrob(ggmod.sublegal.bathy,ggmod.sublegal.relief, ggmod.sublegal.ramps)
ggsave(combine.plot,file="Ningaloo_sublegalgamm.plot.png", width = 30, height = 30,units = "cm")




