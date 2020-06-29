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

rm(list=ls())

##Set working directory----
## Set work directory----
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # sets working directory to where this script is saved (DON'T MOVE THE SCRIPT)

## Set sub directories----
d.dir <- paste(working.dir,"Tidy data",sep="/") 
s.dir <- paste(working.dir,"Spatial",sep="/") # spatial is where I keep spatial data files, rasters and shapefiles
p.dir <- paste(working.dir,"Plots",sep="/")
m.dir <- paste(working.dir,"Model Out", sep="/") #suggest make a model.out.gam folder to keep things seperate

# Bring in and format the data----
name <- 'ningaloo' # for the study

# Load the dataset - from github
setwd(d.dir)
dat <-read.csv('final.data.csv')%>%
  rename(response=target.fish)%>%
  glimpse()

names(dat)

#Set bathymetry to be a positive number for loop
dat<- dat%>%
  mutate(pos.bathymetry=(bathymetry*-1))%>%
  glimpse()

# Set predictor variables---
pred.vars=c("TPI","Slope","Aspect","TRI","Roughness","FlowDir","mean.relief",
            "sd.relief","reef","distance.to.ramp","pos.bathymetry")

# Check for correlation of predictor variables- remove anything highly correlated (>0.95)---
round(cor(dat[,pred.vars], use = "complete.obs"),2)


# TRI is super correlated to both slope and roughness, "Roughness",
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
pred.vars=c("bathymetry","TPI","sqrt.slope","cube.Aspect","log.roughness","FlowDir","mean.relief",
            "sd.relief","sqrt.reef","distance.to.ramp")


# Check to make sure Response vector has not got more than 80% zeros----
glimpse(dat)
str(dat)

unique.vars=unique(as.character(dat$model))
unique.vars.use=character()
for(i in 1:length(unique.vars)){
  temp.dat=dat[which(dat$model==unique.vars[i]),]
  if(length(which(temp.dat$response==0))/nrow(temp.dat)<0.5){
    unique.vars.use=c(unique.vars.use,unique.vars[i])}
}
unique.vars.use     

zero.legal<-dat%>%
  filter(model=='Legal')%>%
  filter(response=='0') #29 which is 21% - there are a lot of legal size fish compared to shallow water surveys

zero.sub<-dat%>%
  filter(model=='Sublegal')%>%
  filter(response=='0') #85 which is 64%

# Remove NA values  - for one of the predictor variables
# (can interpolate but I don't know how to do that right now...)
# 8.05 10.09 10.12 16.03 All have NAs

dat<-dat%>%
  filter(!sample%in%c("8.05","10.09","10.12","16.03"))



# Before we run a model of the data  - check the resiuduals with a 'likley' model - and investigate error distribution and potentially spatial auotcorrelation.

# Likely model 
m1=gam(response~s(bathymetry,k=3,bs='cr')+ s(distance.to.ramp,k=3,bs='cr')+s(site,bs="re"),
           family=tw(),  data=dat%>%filter(model=="Legal"))
summary(m1)
plot(m1)
dev.off()
gam.check(m1)

# We should re-vist spatial auotcorrelation - one approach is to use a 
# https://www.researchgate.net/post/How_can_I_apply_GLM_and_GAM_to_spatially_autocorrelated_data
# 1. Specfic test for spatial autocorrelation - Moran's I
# 2. or can explicitlly include - by including spatial covariance matrix?



# Run the full subset model selection----
setwd(m.dir)
resp.vars=unique.vars.use
use.dat=dat
factor.vars=c("status")# Status as a Factor with two levels
out.all=list()
var.imp=list()


# Check out the model set----

Model1=gam(response~s(bathymetry,k=3,bs='cr')+ s(site,bs="re"),
           family=tw(),  data=use.dat)

model.set=generate.model.set(use.dat=use.dat,
                             test.fit=Model1,
                             pred.vars.cont=pred.vars,
                             pred.vars.fact=factor.vars,
                             # linear.vars="distance.to.ramp",
                             k=3,
                             null.terms="s(site,bs='re')")



# Loop through the FSS function for legal and sublegal
# Loop through the FSS function for each Taxa----
for(i in 1:length(resp.vars)){
  use.dat=dat[which(dat$model==resp.vars[i]),]
  
  Model1=gam(response~s(bathymetry,k=3,bs='cr')+ s(site,bs="re"),
             family=tw(),  data=use.dat)
  
  model.set=generate.model.set(use.dat=use.dat,
                               test.fit=Model1,
                               pred.vars.cont=pred.vars,
                               pred.vars.fact=factor.vars,
                               # linear.vars="distance.to.ramp",
                               k=3,
                               null.terms="s(site,bs='re')")
  
  
  out.list=fit.model.set(model.set,
                         max.models=600,
                         parallel=T)
  names(out.list)
  
  out.list$failed.models # examine the list of failed models
  mod.table=out.list$mod.data.out  # look at the model selection table
  mod.table=mod.table[order(mod.table$AICc),]
  mod.table$cumsum.wi=cumsum(mod.table$wi.AICc)
  out.i=mod.table[which(mod.table$delta.AICc<=3),]
  out.all=c(out.all,list(out.i))
  var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) 
  
  # plot the best models
  for(m in 1:nrow(out.i)){
    best.model.name=as.character(out.i$modname[m])
    
    png(file=paste(name,m,resp.vars[i],"mod_fits.png",sep="_"))
    if(best.model.name!="null"){
      par(mfrow=c(3,1),mar=c(9,4,3,1))
      best.model=out.list$success.models[[best.model.name]]
      plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
      mtext(side=2,text=resp.vars[i],outer=F)}  
    dev.off()
  }
}

# Model fits and importance---
names(out.all)=resp.vars
names(var.imp)=resp.vars
all.mod.fits=do.call("rbind",out.all)
all.var.imp=do.call("rbind",var.imp)
write.csv(all.mod.fits[,-2],file=paste(name,"all.mod.fits.csv",sep="_"))
write.csv(all.var.imp,file=paste(name,"all.var.imp.csv",sep="_"))

# Generic importance plots-
heatmap(all.var.imp,notecex=0.4,  dendrogram ="none",
          col=colorRampPalette(c("white","yellow","red"))(10),
          trace="none",key.title = "",keysize=2,
          notecol="black",key=T,
          sepcolor = "black",margins=c(12,8), lhei=c(4,15),Rowv=FALSE,Colv=FALSE)

############## Part 2 - custom plot of importance scores ###############

# Load the importance score dataset produced above
dat.var.imp <-read.csv("ningaloo_all.var.imp.csv")%>% #from local copy
  rename(resp.var=X)%>%
  gather(key=predictor,value=importance,2:ncol(.))%>%
  glimpse()

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
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

# Labels-
legend_title<-"Importance"

# Annotations-
dat.var.label<-dat.var.imp%>%
  mutate(label=NA)%>%
  mutate(label=ifelse(predictor=="distance.to.ramp"&resp.var=="Legal","X",ifelse(predictor=="bathymetry"&resp.var=="Legal","X",ifelse(predictor=="sqrt.reef"&resp.var=="Legal","X",label))))%>%
  mutate(label=ifelse(predictor=="sd.relief"&resp.var=="Sublegal","X",ifelse(predictor=="bathymetry"&resp.var=="Sublegal","X",label)))
  glimpse()

# Plot gg.importance.scores ----
  # NEED TO ADD IN STATUS
gg.importance.scores <- ggplot(dat.var.label, aes(x=predictor,y=resp.var,fill=importance))+
  geom_tile(show.legend=T) +
  scale_fill_gradientn(legend_title,colours=c("white", re), na.value = "grey98",
                       limits = c(0, max(dat.var.label$importance)))+
  scale_x_discrete(limits=c("bathymetry","TPI","sqrt.slope","Aspect","log.roughness","FlowDir","mean.relief",
                            "sd.relief","sqrt.reef","distance.to.ramp"),
                   labels=c(
                     "bathymetry","TPI","sqrt.slope","Aspect","log.roughness","FlowDir","mean.relief",
                     "sd.relief","sqrt.reef","distance.to.ramp"
                   ))+
  scale_y_discrete(limits = c("Legal",
                              "Sublegal"),
                   labels=c("Legal",
                            "Sublegal"))+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  Theme1+
  geom_text(aes(label=label))
gg.importance.scores


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

# MODEL Legal size bathymetry, distance and sqrt.reef
dat.legal<-dat%>%filter(model=="Legal")
gamm.legal=gam(response~s(bathymetry,k=3,bs='cr')+s(sqrt.reef,k=3,bs='cr')+ s(distance.to.ramp,k=3,bs='cr')+
           s(site,bs="re"), family=tw(),data=dat.legal)

# predict bathymetry from model
mod<-gamm.legal
testdata <- expand.grid(sqrt.reef=seq(min(dat$sqrt.reef),max(dat$sqrt.reef),length.out = 20),
                        bathymetry=seq(min(dat$bathymetry),max(dat$bathymetry), length.out=20),
                        distance.to.ramp=seq(min(dat$distance.to.ramp),max(dat$distance.to.ramp),length.out=20),
                        site=(mod$model$site))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.bathy = testdata%>%data.frame(fits)%>%
  group_by(bathymetry)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.legal.bathy,"predict.legal.bathy.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.bathy<-read.csv("predict.legal.bathy.csv")%>%
  glimpse()

#predict sqrt.reef
mod<-gamm.legal
testdata <- expand.grid(sqrt.reef=seq(min(dat$sqrt.reef),max(dat$sqrt.reef),length.out = 20),
                        bathymetry=seq(min(dat$bathymetry),max(dat$bathymetry), length.out=20),
                        distance.to.ramp=seq(min(dat$distance.to.ramp),max(dat$distance.to.ramp),length.out=20),
                        site=(mod$model$site))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.sqrt.reef = testdata%>%data.frame(fits)%>%
  group_by(sqrt.reef)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.legal.sqrt.reef,"predict.legal.sqrt.reef.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.sqrt.reef<-read.csv("predict.legal.sqrt.reef.csv")%>%
  glimpse()

#predict distance to ramp
mod<-gamm.legal
testdata <- expand.grid(sqrt.reef=seq(min(dat$sqrt.reef),max(dat$sqrt.reef),length.out = 20),
                        bathymetry=seq(min(dat$bathymetry),max(dat$bathymetry), length.out=20),
                        distance.to.ramp=seq(min(dat$distance.to.ramp),max(dat$distance.to.ramp),length.out=20),
                        site=(mod$model$site))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.ramp = testdata%>%data.frame(fits)%>%
  group_by(distance.to.ramp)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.legal.ramp,"predict.legal.ramp.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.ramp<-read.csv("predict.legal.ramp.csv")%>%
  glimpse()

# MODEL Sublegal bathymetry, sd.reef
dat.sublegal<-dat%>%filter(model=="Sublegal")
gamm.sublegal=gam(response~s(bathymetry,k=3,bs='cr')+ s(sd.relief,k=3,bs='cr')+
                 s(site,bs="re"), family=tw(),data=dat.sublegal)

#predict bathymetry
mod<-gamm.sublegal
testdata <- expand.grid(sd.relief=seq(min(dat$sd.relief),max(dat$sd.relief),length.out = 20),
                        bathymetry=seq(min(dat$bathymetry),max(dat$bathymetry), length.out=20),
                        site=(mod$model$site))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.sublegal.bathy = testdata%>%data.frame(fits)%>%
  group_by(bathymetry)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.sublegal.bathy,"predict.sublegal.bathy.csv") #there is some BUG in dplyr - that this fixes
predicts.sublegal.bathy<-read.csv("predict.sublegal.bathy.csv")%>%
  glimpse()


# Predict sd.relief
mod<-gamm.sublegal
testdata <- expand.grid(sd.relief=seq(min(dat$sd.relief),max(dat$sd.relief),length.out = 20),
                        bathymetry=seq(min(dat$bathymetry),max(dat$bathymetry), length.out=20),
                        site=(mod$model$site))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predict.sublegal.sd.relief = testdata%>%data.frame(fits)%>%
  group_by(sd.relief)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predict.sublegal.sd.relief,"predict.sublegal.sd.relief.csv") #there is some BUG in dplyr - that this fixes
predict.sublegal.sd.relief<-read.csv("predict.sublegal.sd.relief.csv")%>%
  glimpse()

## PLOTS for Legal sized fish ----

# Distance to ramps
ggmod.legal.ramps<- ggplot() +
  ylab(" ")+
  xlab('Distance to Ramps')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #   geom_jitter(width = 0.25,height = 0)+
  geom_point(data=dat.legal,aes(x=distance.to.ramp,y=response,colour=status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.legal.ramp,aes(x=distance.to.ramp,y=response),alpha=0.5)+
  geom_line(data=predicts.legal.ramp,aes(x=distance.to.ramp,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.ramp,aes(x=distance.to.ramp,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.legal.ramps

# Bathymetry
ggmod.legal.bathy<- ggplot() +
  ylab(" ")+
  xlab('Bathymetry (m)')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #   geom_jitter(width = 0.25,height = 0)+
  geom_point(data=dat.legal,aes(x=bathymetry,y=response,colour=status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response),alpha=0.5)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.legal.bathy

# sqrt.reef
ggmod.legal.reef<- ggplot() +
  ylab(" ")+
  xlab('% Reef (sqrt)')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #   geom_jitter(width = 0.25,height = 0)+
  geom_point(data=dat.legal,aes(x=sqrt.reef,y=response,colour=status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.legal.sqrt.reef,aes(x=sqrt.reef,y=response),alpha=0.5)+
  geom_line(data=predicts.legal.sqrt.reef,aes(x=sqrt.reef,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.sqrt.reef,aes(x=sqrt.reef,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.legal.reef

# combined.plot using grid() and gridExtra()------
blank <- grid.rect(gp=gpar(col="white"))

# To see what they will look like use grid.arrange() - make sure Plot window is large enough! - or will error!
grid.arrange(ggmod.legal.ramps,ggmod.legal.bathy,ggmod.legal.reef,blank)

# Use arrangeGrob ONLY - as we can pass this to ggsave! Note use of raw ggplot's
combine.plot<-arrangeGrob(ggmod.legal.ramps,ggmod.legal.bathy,ggmod.legal.reef,blank)
ggsave(combine.plot,file="Ningaloo_legalgamm.plot.png", width = 30, height = 30,units = "cm")

## PLOTS for Sublegal sized fish ----

#Bathymetry
ggmod.sublegal.bathy<- ggplot() +
  ylab(" ")+
  xlab('Bathymetry (m)')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #   geom_jitter(width = 0.25,height = 0)+
  geom_point(data=dat.sublegal,aes(x=bathymetry,y=response,colour=status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.sublegal.bathy,aes(x=bathymetry,y=response),alpha=0.5)+
  geom_line(data=predicts.sublegal.bathy,aes(x=bathymetry,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.sublegal.bathy,aes(x=bathymetry,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.sublegal.bathy

#sd.relief
ggmod.sublegal.relief<- ggplot() +
  ylab(" ")+
  xlab('SD Relief')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #   geom_jitter(width = 0.25,height = 0)+
  geom_point(data=dat.sublegal,aes(x=sd.relief,y=response,colour=status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predict.sublegal.sd.relief,aes(x=sd.relief,y=response),alpha=0.5)+
  geom_line(data=predict.sublegal.sd.relief,aes(x=sd.relief,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predict.sublegal.sd.relief,aes(x=sd.relief,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.sublegal.relief

# combined.plot using grid() and gridExtra()------
blank <- grid.rect(gp=gpar(col="white"))

# To see what they will look like use grid.arrange() - make sure Plot window is large enough! - or will error!
grid.arrange(ggmod.sublegal.bathy,ggmod.sublegal.relief)

# Use arrangeGrob ONLY - as we can pass this to ggsave! Note use of raw ggplot's
combine.plot<-arrangeGrob(ggmod.sublegal.bathy,ggmod.sublegal.relief)
ggsave(combine.plot,file="Ningaloo_sublegalgamm.plot.png", width = 30, height = 30,units = "cm")









