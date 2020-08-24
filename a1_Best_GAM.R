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
library(ggplot2)
library(vcd)

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

legal.dat <- subset(dat, model=='Legal')

# Remove NA values and sites in state no-take - for one of the predictor variables
# (can interpolate but I don't know how to do that right now...)
# 8.05 10.09 10.12 16.03 All have NAs

legal.dat<-legal.dat%>%
  filter(!sample%in%c("8.05","10.09","10.12","16.03", "10.14", "10.15", "10.13", "14.13", "14.12", "14.02"))

# Transform variables 
legal.dat <- legal.dat%>%
  mutate(sqrt.reef=sqrt(reef))%>%
  mutate(sqrt.slope=sqrt(Slope))%>%
  #mutate(sqrt.TPI=sqrt(TPI))%>%
  mutate(log.roughness=log(Roughness+1))%>%
  mutate(cube.Aspect=(Aspect)^3)%>%
  glimpse()

# Add row number 
legal.dat <- legal.dat %>% 
  mutate(id = row_number())

#Set predictor variables 
pred.vars=c("sqrt.slope","cube.Aspect","log.roughness","FlowDir",
            "distance.to.ramp", "site")


##### With Site as a random factor ####
setwd(m.dir)
use.dat <- legal.dat
factor.vars <- c("status") # Status as a Factor with two levels
out.all <- list()
var.imp <- list()


Model1 <- gam(response~s(distance.to.ramp,k=5,bs='cr')+ 
                s(site,bs="re") + s(bathymetry, k=5, bs="cr"),
              family=tw(),  data=use.dat)

model.set <- generate.model.set(use.dat=use.dat,
                                test.fit=Model1,
                                pred.vars.cont=pred.vars,
                                pred.vars.fact=factor.vars,
                                #smooth.smooth.interactions=TRUE,
                                max.predictors=5,
                                k=5,
                                null.terms="s(site,bs='re') + s(bathymetry, k=5, bs='cr')")

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
dat.legal<-dat%>%filter(model=="Legal")
gam.site<-gam(response~s(bathymetry,k=3,bs='cr')+s(cube.Aspect,k=5,bs='cr') + status +
                 s(site,bs="re"), family=tw(),data=legal.dat)

# predict bathymetry from model
mod<-gam.site
testdata <- expand.grid(cube.Aspect=seq(min(legal.dat$cube.Aspect),max(legal.dat$cube.Aspect),length.out = 20),
                        bathymetry=seq(min(legal.dat$bathymetry),max(legal.dat$bathymetry), length.out=20),
                        status=(mod$model$status),
                        site=(mod$model$site))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.bathy = testdata%>%data.frame(fits)%>%
  dplyr::group_by(bathymetry)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  dplyr::ungroup()
write.csv(predicts.legal.bathy,"predict.legal.bathy.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.bathy<-read.csv("predict.legal.bathy.csv")%>%
  glimpse()

# Plot Bathymetry from model
ggmod.gam.bathy<- ggplot() +
  ylab("Predicted Abundance of Legal Sized Fish")+
  xlab('Bathymetry (m)')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=legal.dat,aes(x=bathymetry,y=response,colour=status),  alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response),alpha=0.75)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.gam.bathy

#predict cube.Aspect
mod<-gam.site
testdata <- expand.grid(cube.Aspect=seq(min(legal.dat$cube.Aspect),max(legal.dat$cube.Aspect),length.out = 20),
                        bathymetry=seq(min(legal.dat$bathymetry),max(legal.dat$bathymetry), length.out=20),
                        status=(mod$model$status),
                        site=(mod$model$site))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.aspect = testdata%>%data.frame(fits)%>%
  dplyr::group_by(cube.Aspect)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  dplyr::ungroup()
write.csv(predicts.legal.sqrt.reef,"predict.legal.sqrt.reef.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.sqrt.reef<-read.csv("predict.legal.sqrt.reef.csv")%>%
  glimpse()

# plot cube.Aspect
ggmod.gam.aspect<- ggplot() +
  ylab("Predicted Abundance of Legal Sized Fish")+
  xlab('Aspect (cubed)')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=legal.dat,aes(x=cube.Aspect,y=response,colour=status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response),alpha=0.75)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.gam.aspect

# predict status
mod<-gam.site
testdata <- expand.grid(cube.Aspect=seq(min(legal.dat$cube.Aspect),max(legal.dat$cube.Aspect),length.out = 20),
                        bathymetry=seq(min(legal.dat$bathymetry),max(legal.dat$bathymetry), length.out=20),
                        status=(mod$model$status),
                        site=(mod$model$site))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.gam.status = testdata%>%data.frame(fits)%>%
  dplyr::group_by(status)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  dplyr::ungroup()
write.csv(predicts.legal.sqrt.reef,"predict.legal.sqrt.reef.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.sqrt.reef<-read.csv("predict.legal.sqrt.reef.csv")%>%
  glimpse()

# Plot status
ggmod.gam.status<- ggplot(aes(x=status,y=response,colour=status), data=legal.dat) +
  ylab("Predicted Abundance of Legal Sized Fish")+
  xlab('Status')+
  #   ggtitle(substitute(italic(name)))+
  scale_fill_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  scale_x_discrete(limits = rev(levels(predicts.gam.status$status)))+
  #geom_bar(stat = "identity")+
  geom_errorbar(data=predicts.gam.status,aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  theme_classic()+
  Theme1
ggmod.gam.status

#### FSSgam gamm4 + random site ####
setwd(m.dir)
# Remove any unused columns from the dataset 
use.dat <- legal.dat%>%
  dplyr::select(response, distance.to.ramp, cube.Aspect, sqrt.slope, log.roughness, status, FlowDir,
                bathymetry, site)
factor.vars <- c("status")

out.all <- list()
var.imp <- list()

Model1 <- uGamm(response~s(bathymetry,k=5,bs='cr'),random = list(site = ~1),
                family=poisson(), 
                #correlation=corGaus(form = ~ latitude + longitude),
                data=use.dat,
                lme4=FALSE)


model.set <- generate.model.set(use.dat=use.dat,
                                test.fit=Model1,
                                pred.vars.cont=pred.vars,
                                pred.vars.fact=factor.vars,
                                max.predictors=5,
                                k=5,
                                null.terms= "s(bathymetry, bs='cr')")


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

write.csv(all.mod.fits[,-2],file=paste(name,"gamm.random.all.mod.fits.csv",sep="_"))
write.csv(all.var.imp,file=paste(name,"gamm.all.var.imp.csv",sep="_"))


#### Predict and plot best model ####
ugamm.site <-uGamm(response~s(bathymetry,k=5,bs='cr')+ status + s(cube.Aspect,k=5,bs='cr'),
                     random = list(site = ~1), family=poisson(), data=use.dat, lme4=FALSE)




# predict bathymetry from model
mod<-ugamm.site$gam
testdata <- expand.grid(cube.Aspect=seq(min(legal.dat$cube.Aspect),max(legal.dat$cube.Aspect),length.out = 20),
                        bathymetry=seq(min(legal.dat$bathymetry),max(legal.dat$bathymetry), length.out=20),
                        status=(mod$model$status))%>%
                       
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.bathy = testdata%>%data.frame(fits)%>%
  dplyr::group_by(bathymetry)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  dplyr::ungroup()
write.csv(predicts.legal.bathy,"predict.legal.bathy.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.bathy<-read.csv("predict.legal.bathy.csv")%>%
  glimpse()

# Plot Bathymetry from model
ggmod.gam.bathy<- ggplot() +
  ylab("Predicted Abundance of Legal Sized Fish")+
  xlab('Bathymetry (m)')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=legal.dat,aes(x=bathymetry,y=response,colour=status),  alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response),alpha=0.75)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.gam.bathy

#predict cube.Aspect
mod<-ugamm.site$gam
testdata <- expand.grid(cube.Aspect=seq(min(legal.dat$cube.Aspect),max(legal.dat$cube.Aspect),length.out = 20),
                        bathymetry=seq(min(legal.dat$bathymetry),max(legal.dat$bathymetry), length.out=20),
                        status=(mod$model$status),
                        site=(mod$model$site))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.aspect = testdata%>%data.frame(fits)%>%
  dplyr::group_by(cube.Aspect)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  dplyr::ungroup()
write.csv(predicts.legal.sqrt.reef,"predict.legal.sqrt.reef.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.sqrt.reef<-read.csv("predict.legal.sqrt.reef.csv")%>%
  glimpse()

# plot cube.Aspect
ggmod.gam.aspect<- ggplot() +
  ylab("Predicted Abundance of Legal Sized Fish")+
  xlab('Aspect (cubed)')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=legal.dat,aes(x=cube.Aspect,y=response,colour=status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response),alpha=0.75)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.gam.aspect

# predict status
mod<-ugamm.site$gam
testdata <- expand.grid(cube.Aspect=seq(min(legal.dat$cube.Aspect),max(legal.dat$cube.Aspect),length.out = 20),
                        bathymetry=seq(min(legal.dat$bathymetry),max(legal.dat$bathymetry), length.out=20),
                        status=(mod$model$status),
                        site=(mod$model$site))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.gam.status = testdata%>%data.frame(fits)%>%
  dplyr::group_by(status)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  dplyr::ungroup()
write.csv(predicts.legal.sqrt.reef,"predict.legal.sqrt.reef.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.sqrt.reef<-read.csv("predict.legal.sqrt.reef.csv")%>%
  glimpse()

# Plot status
ggmod.gam.status<- ggplot(aes(x=status,y=response,colour=status), data=legal.dat) +
  ylab("Predicted Abundance of Legal Sized Fish")+
  xlab('Status')+
  #   ggtitle(substitute(italic(name)))+
  scale_fill_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  scale_x_discrete(limits = rev(levels(predicts.gam.status$status)))+
  #geom_bar(stat = "identity")+
  geom_errorbar(data=predicts.gam.status,aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  theme_classic()+
  Theme1
ggmod.gam.status

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


#### Predict and plot best model ####
lme4.site <- uGamm(response~s(bathymetry,k=5,bs='cr')+ s(distance.to.ramp, k=5, bs='cr') + s(cube.Aspect,k=5,bs='cr'),
                  random=~(1|site), 
                  family=poisson(), data=use.dat, lme4=TRUE)

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp <- as.data.frame(overdisp_fun(lme4.site$gam))

fit <- goodfit(legal.dat$response, type='poisson')
summary(fit)
rootogram(fit)

# predict bathymetry from model
mod<-lme4.site$gam
testdata <- expand.grid(cube.Aspect=seq(min(legal.dat$cube.Aspect),max(legal.dat$cube.Aspect),length.out = 20),
                        distance.to.ramp=seq(min(legal.dat$distance.to.ramp),max(legal.dat$distance.to.ramp),length.out = 20),
                        bathymetry=seq(min(legal.dat$bathymetry),max(legal.dat$bathymetry), length.out=20))%>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.bathy = testdata%>%data.frame(fits)%>%
  dplyr::group_by(bathymetry)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  dplyr::ungroup()
write.csv(predicts.legal.bathy,"predict.legal.bathy.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.bathy<-read.csv("predict.legal.bathy.csv")%>%
  glimpse()

# Plot Bathymetry from model
ggmod.gam.bathy<- ggplot() +
  ylab("Predicted Abundance of Legal Sized Fish")+
  xlab('Bathymetry (m)')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=legal.dat,aes(x=bathymetry,y=response,colour=status),  alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response),alpha=0.75)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.gam.bathy

#predict cube.Aspect
mod<-lme4.site$gam
testdata <- expand.grid(cube.Aspect=seq(min(legal.dat$cube.Aspect),max(legal.dat$cube.Aspect),length.out = 20),
                        distance.to.ramp=seq(min(legal.dat$distance.to.ramp),max(legal.dat$distance.to.ramp),length.out = 20),
                        bathymetry=seq(min(legal.dat$bathymetry),max(legal.dat$bathymetry), length.out=20))%>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.aspect = testdata%>%data.frame(fits)%>%
  dplyr::group_by(cube.Aspect)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  dplyr::ungroup()
write.csv(predicts.legal.sqrt.reef,"predict.legal.sqrt.reef.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.sqrt.reef<-read.csv("predict.legal.sqrt.reef.csv")%>%
  glimpse()

# plot cube.Aspect
ggmod.gam.aspect<- ggplot() +
  ylab("Predicted Abundance of Legal Sized Fish")+
  xlab('Aspect (cubed)')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=legal.dat,aes(x=cube.Aspect,y=response,colour=status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response),alpha=0.75)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.gam.aspect

# predict plot distance.to.ramp
mod<-lme4.site$gam
testdata <- expand.grid(cube.Aspect=seq(min(legal.dat$cube.Aspect),max(legal.dat$cube.Aspect),length.out = 20),
                        distance.to.ramp=seq(min(legal.dat$distance.to.ramp),max(legal.dat$distance.to.ramp),length.out = 20),
                        bathymetry=seq(min(legal.dat$bathymetry),max(legal.dat$bathymetry), length.out=20))%>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.distance = testdata%>%data.frame(fits)%>%
  dplyr::group_by(distance.to.ramp)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  dplyr::ungroup()
write.csv(predicts.legal.sqrt.reef,"predict.legal.sqrt.reef.csv") #there is some BUG in dplyr - that this fixes
predicts.legal.sqrt.reef<-read.csv("predict.legal.sqrt.reef.csv")%>%
  glimpse()

# plot cube.Aspect
ggmod.gam.distance<- ggplot() +
  ylab("Predicted Abundance of Legal Sized Fish")+
  xlab('Distance to Ramp')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=legal.dat,aes(x=distance.to.ramp,y=response,colour=status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.legal.distance,aes(x=distance.to.ramp,y=response),alpha=0.75)+
  geom_line(data=predicts.legal.distance,aes(x=distance.to.ramp,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.distance,aes(x=distance.to.ramp,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.gam.distance

#### Simulating my model for verification lme4 ####

