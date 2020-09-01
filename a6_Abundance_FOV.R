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
dat <- subset(dat, model=='Legal')

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


#### FSSgam using lme4 + random site ####
setwd(m.dir)
# Remove any unused columns from the dataset 
use.dat <- dat%>%
  dplyr::select(response, distance.to.ramp, cube.Aspect, sqrt.slope, log.roughness, status,
                bathymetry, site, mean.relief, sd.relief, sqrt.reef)
                #reef, sd.relief, mean.relief)
#Set predictor variables 
pred.vars=c("sqrt.slope","cube.Aspect","log.roughness",
            "distance.to.ramp", "sd.relief", "mean.relief", "sqrt.reef")
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

write.csv(all.mod.fits[,-2],file=paste(name,"lme4.random.all.mod.fits.csv",sep="_"))
write.csv(all.var.imp,file=paste(name,"lme4.all.var.imp.fov.csv",sep="_"))

### Plot Importance Scores separates by whether it using FOV or not ####
setwd(m.dir)
dat.var.imp <-read.csv("ningaloo_lme4.all.var.imp.nofov.csv")%>% #from local copy
  rename(resp.var=X)%>%
  gather(key=predictor,value=importance,2:ncol(.))%>%
  glimpse()

dat.var.imp.fov  <- read.csv("ningaloo_lme4.all.var.imp.fov.csv")%>% #from local copy
  rename(resp.var=X)%>%
  gather(key=predictor,value=importance,2:ncol(.))%>%
  glimpse()

# Add column to say what model was used 
dat.var.imp <- dat.var.imp%>%
  mutate(model="No FOV")

dat.var.imp.fov <- dat.var.imp.fov%>%
  mutate(model="FOV")

# Stick importance scores together
all.var.imp <- rbind(dat.var.imp,dat.var.imp.fov)

# Plotting defaults----
library(ggplot2)

# Theme-
Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill="white"),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=8),
    legend.title = element_text(size=8, face="bold"),
    legend.position = "top",
    legend.direction="horizontal",
    text=element_text(size=10),
    strip.text.y = element_text(size = 10,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=10),
    axis.title.y=element_text(vjust=0.6, angle=90, size=10),
    axis.text.x=element_text(size=10,angle = 90, hjust=1,vjust=0.5),
    axis.text.y=element_text(size=10,face="italic"),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())

# colour ramps-
re <- colorRampPalette(c("lightskyblue1","royalblue3"))(200)

# Labels-
legend_title<-"Importance"

# Annotations-
dat.var.label<-all.var.imp%>%
  mutate(label=NA)%>%
  mutate(label=ifelse(predictor=="distance.to.ramp"&resp.var=="Legal","X",ifelse(predictor=="bathymetry"&resp.var=="Legal","X",ifelse(predictor=="sqrt.reef"&resp.var=="Legal","X",label))))%>%
  mutate(label=ifelse(predictor=="distance.to.ramp"&resp.var=="Sublegal","X",ifelse(predictor=="bathymetry"&resp.var=="Sublegal","X",ifelse(predictor=="sd.relief"&resp.var=="Sublegal","X",label))))%>%
  glimpse()

# Plot gg.importance.scores
gg.importance.scores <- ggplot(dat.var.label, aes(x=predictor,y=model,fill=importance))+
  geom_tile(show.legend=T) +
  scale_fill_gradientn(legend_title,colours=c("white", re), na.value = "grey98",
                       limits = c(0, max(dat.var.label$importance)))+
  scale_x_discrete(limits=c("sqrt.slope","cube.Aspect","log.roughness","distance.to.ramp", "status",
                            "mean.relief", "sd.relief","sqrt.reef"),
                   labels=c(
                     "Slope (sqrt)","Aspect (cubed)","Roughness (log)","Distance to Ramp", "Status",
                     "Mean Relief", "SD Relief","% Reef (sqrt)"
                   ))+
  scale_y_discrete(limits = c("FOV",
                              "No FOV"),
                   labels=c("FOV",
                            "No FOV"))+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  Theme1+
  geom_text(aes(label=label))
gg.importance.scores
