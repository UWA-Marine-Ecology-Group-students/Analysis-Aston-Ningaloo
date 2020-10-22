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
m.dir <- paste(working.dir,"Model Out GAM", sep="/")
b.dir <- paste(working.dir, "Model Out StanGam", sep="/")

# Bring in and format the data
name <- 'ningaloo' # for the study

# Names of target species
setwd(d.dir)
target <- read.csv("target.species.csv")
target.names <- as.data.frame(unique(target$scientific))
x <- "scientific"
colnames(target.names) <- x

# All maxm
setwd(d.dir)
maxn.dat <-read.csv("ningaloo.complete.maxn.csv")
unique(maxn.dat$scientific)

#Select species from maxn that do not appear in target species list 
non.target.data <- anti_join(maxn.dat,target.names, by='scientific')

#Select only demersal species
non.target.data <- non.target.data%>%
  dplyr::filter(!genus %in% c("Aipys","Alepes","Atule","Carangoides","Caranx","Carcharhinus",
                               "Cybiosarda","Decapterus","Emydocephalus","Euthynnus","Galeocerdo",
                               "Gnathanodon","Hemipristis","Hydrophis","Loxodon","Megalaspis",
                               "Negaprion","Scombero","Selar","Seriola","Sphyraena","Sphyrna",
                               "Stegostoma","Ulua","Unk", "Pterocaesio"))

#Sum for maxn of non-target at each site 
non.target <- non.target.data%>%
  group_by(sample)%>%
  summarise_at(vars(maxn), list(sum))

#Merge with variables 
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
                'latitude', 'longitude', "sd.relief", "mean.relief", "sqrt.reef")

full.non.target <- merge(non.target,variables, all=TRUE)
full.non.target[is.na(full.non.target)] <- 0

full.non.target <- full.non.target%>%
  dplyr::rename(response=maxn)

# Remove NA values and sites in state no-take - for one of the predictor variables
# (can interpolate but I don't know how to do that right now...)
# 8.05 10.09 10.12 16.03 All have NAs

full.non.target<-full.non.target%>%
  filter(!sample%in%c("8.05","10.09","10.12","16.03", "10.14", "10.15", "10.13", "14.13", "14.12", "14.02"))

full.non.target <- full.non.target%>%
  mutate(ID = rownames(full.non.target))

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

## Run full subsets 
use.dat <- full.non.target%>%
  dplyr::select(response, distance.to.ramp, cube.Aspect, sqrt.slope, log.roughness, status,
                bathymetry,site, sd.relief, mean.relief, sqrt.reef, ID)
pred.vars=c("sqrt.slope","cube.Aspect","log.roughness",
            "distance.to.ramp","sd.relief", "mean.relief", "sqrt.reef")
factor.vars <- c("status")

out.all <- list()
var.imp <- list()

Model.nontarget.NOFOV <- uGamm(response~s(bathymetry,k=5, bs='cr'),
                    family=poisson, random=~(1|site)+(1|ID),
                    data=use.dat,
                    lme4=TRUE)


model.set <- generate.model.set(use.dat=use.dat,
                                test.fit=Model.nontarget.NOFOV,
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

setwd(m.dir)
write.csv(all.mod.fits, "nontarget.models.nofov.csv")
write.csv(all.var.imp, "nontarget.imp.nofov.csv")

#### Model 1 no FOV ####
effect.dat <- expand.grid(cube.Aspect=seq(min(full.non.target$cube.Aspect),max(full.non.target$cube.Aspect),length.out = 20),
                          bathymetry=seq(min(full.non.target$bathymetry),max(full.non.target$bathymetry), length.out=20))

Model.NFOV1 <- uGamm(response~s(bathymetry,k=5,bs='cr')
                        + s(sd.relief,k=5,bs='cr'),
                        random=~(1|site)+(1|ID), 
                        family=poisson(), data=use.dat, lme4=TRUE)

gam.check(Model.NFOV1$gam)

stangam.s.NFOV1  <- stan_gamm4(response~s(bathymetry, k = 5, bs = "cr")
                                  + s(cube.Aspect, k=5,bs='cr'),
                                  random=~(1|site)+(1|ID), adapt_delta = 0.99,
                                  data=use.dat, chains=3, cores=3, iter=41000, warmup=40000, thin=3,
                                  family=poisson)

# calculate frequentist effects - link function is log 
predicted <- predict.gam(Model.NFOV1$gam, newdata=effect.dat, type='response')
effect.dat <- cbind(effect.dat, predicted)

?predict.gam
# bathy effect 
bathy.P <- effect.dat%>%
  dplyr::group_by(bathymetry)%>%
  summarise_at(vars(predicted), list(mean))

bathy.E <- diff(range(bathy.P$predicted)) 
bathy.E

# plot bathy effect 
predicted.fit <- predict.gam(Model.NFOV1$gam, newdata=effect.dat, type='response', se.fit=T)
effect.dat.fit <- cbind(effect.dat, predicted)

predicts.legal.bathy = effect.dat.fit%>%data.frame(predicted.fit)%>%
  group_by(bathymetry)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.legal.bathy<- ggplot() +
  ylab("Predicted Abundance of Non-target Species")+
  xlab('Bathymetry (m)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.non.target,aes(x=bathymetry,y=response), colour="lightblue", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response), colour='darkblue', alpha=0.75)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response - se.fit), colour='darkblue', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.legal.bathy,aes(x=bathymetry,y=response + se.fit), colour='darkblue', linetype="dashed",alpha=0.75)+
  theme_classic()+
  xlim(-180,-52)+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.legal.bathy


# aspect effect 
aspect.P <- effect.dat%>%
  group_by(cube.Aspect)%>%
  summarise_at(vars(predicted), list(mean))

aspect.E <- diff(range(aspect.P$predicted))
aspect.E

predicts.legal.aspect = effect.dat%>%data.frame(predicted.fit)%>%
  group_by(cube.Aspect)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.legal.aspect<- ggplot() +
  ylab("Predicted Abundance of Non-target Species")+
  xlab('Cubed Aspect (degrees)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.non.target,aes(x=cube.Aspect,y=response), colour="lightgreen", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response), colour='darkgreen', alpha=0.75)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response - se.fit), colour='darkgreen', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response + se.fit), colour='darkgreen', linetype="dashed",alpha=0.75)+
  theme_classic()+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.legal.aspect

# use the bayesian models to calculate a posterior distribution of the effect size
setwd(b.dir)
effect.dat <- expand.grid(cube.Aspect=seq(min(full.non.target$cube.Aspect),max(full.non.target$cube.Aspect),length.out = 20),
                          bathymetry=seq(min(full.non.target$bathymetry),max(full.non.target$bathymetry), length.out=20))
effect.dat.sim <- effect.dat
gg <- t(posterior_predict(stangam.s.NFOV1, effect.dat.sim, re.form=NA)) # Returning in fish
#colnames(gg) <- paste("sim",1:ncol(gg),sep="_")
gg.dat <- cbind(effect.dat.sim,gg)


gg.dat.long <- gg.dat%>%
  gather(sim, predicted.sim, 3:1004)

gg.dat.long$sim <- as.factor(gg.dat.long$sim)
str(gg.dat.long)

# bathy effect 
bathy.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("Bathymetry (m)")
colnames(bathy.sim) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.dat.long
  bathy.sim.2 <- dplyr::filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- dplyr::group_by(bathy.sim.2, bathymetry)
  bathy.sim.4 <- dplyr::summarise_at(bathy.sim.3, vars(predicted.sim), list(mean))
  print(bathy.sim.5 <- diff(range(bathy.sim.4$predicted.sim)))
  bathy.sim <- rbind(bathy.sim, bathy.sim.5)
  
}

bathy.sim


write.csv(bathy.sim, "bayesian.bathy.predictions.nontarget.NFOV1.csv")

# aspect effect 
aspect.sim <- data.frame(matrix(ncol=1, nrow=1))
x <- c("Cubed Aspect (degrees)")
colnames(aspect.sim) <- x

for(i in 1:1002){
  aspect.sim.1 <- gg.dat.long
  aspect.sim.2 <- filter(aspect.sim.1, sim==i)
  aspect.sim.3 <- group_by(aspect.sim.2, cube.Aspect)
  aspect.sim.4 <- summarise_at(aspect.sim.3, vars(predicted.sim), list(mean))
  print(aspect.sim.5 <- diff(range(aspect.sim.4$predicted.sim)))
  aspect.sim <- rbind(aspect.sim, aspect.sim.5)
  
}

aspect.sim

write.csv(aspect.sim, "bayesian.aspect.predictions.nontarget.NFOV1.csv")

# Plot effect sizes 
bathy.sim <- bathy.sim[-1,]
bathy.sim <- as.data.frame(bathy.sim)
x <- c("Bathymetry (m)")
colnames(bathy.sim) <- x

aspect.sim <- aspect.sim[-1,]
aspect.sim <- as.data.frame(aspect.sim)
x <- c("Cube Aspect (Degrees)")
colnames(aspect.sim) <- x

full.data <- as.data.frame(cbind(aspect.sim, bathy.sim))

full.data.long <- full.data%>%
  gather(variable, effect.size, 1:2)

colours <- c('#619CFF', '#00BA38') #F8766D
effect.plot <- ggplot(full.data.long, aes(x = effect.size, fill = variable, colour = variable)) + geom_density(alpha = 0.5) +
  scale_fill_manual(values=colours)+
  scale_colour_manual(values=colours)+
  geom_vline(xintercept = 5.13502, color = "steelblue", size=0.75)+
  #geom_vline(xintercept = 0, color = "tomato3", size=0.75)+
  geom_vline(xintercept = 10.52997, color = "springgreen4", size=0.75)+
  xlim(-1,35)+
  theme_classic()+
  labs(y="Density", x="Effect Size")+
  Theme1
effect.plot
aspect.E

## Make prediction plots

# bathy effect 
bathy.means <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.mean")
colnames(bathy.means) <- x

for(i in unique(gg.dat.long$bathymetry)){
  bathy.sim.1 <- gg.dat.long
  bathy.sim.2 <- dplyr::filter(bathy.sim.1, bathymetry==i)
  print(bathy.sim.3 <- dplyr::summarise_at(bathy.sim.2, vars(predicted.sim), list(mean)))
  bathy.means <- rbind(bathy.means, bathy.sim.3$predicted.sim)
}

bathy.means <- bathy.means[-1,]
depths <- unique(gg.dat$bathymetry)
bathy.means <- as.data.frame(cbind(bathy.means,depths))

# Get means for each value of bathymetry for each simulation 
bathy.mean.sim <- data.frame(matrix(ncol=1, nrow=20))
x <- c("bathy.mean.sim")
colnames(bathy.mean.sim) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.dat.long
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
  labs(y="Predicted Abundance Non-Target Species", x="Bathymetry (m)")+
  theme_classic()+
  xlim(-180,-51)+
  Theme1
predict.plot.bathy

# aspect effect 
aspect.means <- data.frame(matrix(ncol=1, nrow=1))
x <- c("aspect.mean")
colnames(aspect.means) <- x

for(i in unique(gg.dat.long$cube.Aspect)){
  aspect.sim.1 <- gg.dat.long
  aspect.sim.2 <- dplyr::filter(aspect.sim.1, cube.Aspect==i)
  print(aspect.sim.3 <- dplyr::summarise_at(aspect.sim.2, vars(predicted.sim), list(mean)))
  aspect.means <- rbind(aspect.means, aspect.sim.3$predicted.sim)
}

aspect.means <- aspect.means[-1,]
aspect <- unique(gg.dat$cube.Aspect)
aspect.means <- as.data.frame(cbind(aspect.means,aspect))

# Get means for each value of aspect for each simulation 
aspect.mean.sim <- data.frame(matrix(ncol=1, nrow=20))
x <- c("aspect.mean.sim")
colnames(aspect.mean.sim) <- x

for(i in 1:1002){
  aspect.sim.1 <- gg.dat.long
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
  labs(y="Predicted Abundance of Non-target Species", x="Cube Aspect (degrees)")+
  theme_classic()+
  Theme1
predict.plot.aspect


#### Model 2 no FOV ####
effect.dat <- expand.grid(cube.Aspect=seq(min(full.non.target$cube.Aspect),max(full.non.target$cube.Aspect),length.out = 20),
                          bathymetry=seq(min(full.non.target$bathymetry),max(full.non.target$bathymetry), length.out=20),
                          status=unique(full.non.target$status))

Model.NFOV2 <- uGamm(response~s(bathymetry,k=5,bs='cs')+
                       + s(cube.Aspect,k=5,bs='cs')+factor(status),
                     random=~(1|site)+(1|ID), 
                     family=poisson(), data=use.dat, lme4=TRUE)

gam.check(Model.NFOV2$gam)

stangam.s.NFOV2  <- stan_gamm4(response~s(bathymetry, k = 5, bs = "cs") +
                                 + s(cube.Aspect, k=5,bs='cs')+status,
                               random=~(1|site)+(1|ID), adapt_delta = 0.99,
                               data=use.dat, chains=3, cores=3, iter=41000, warmup=40000, thin=3,
                               family=poisson)

# calculate frequentist effects - link function is log 
predicted <- predict.gam(Model.NFOV2$gam, newdata=effect.dat, type='response')
effect.dat <- cbind(effect.dat, predicted)


# bathy effect 
bathy.P <- effect.dat%>%
  dplyr::group_by(bathymetry)%>%
  summarise_at(vars(predicted), list(mean))

bathy.E.NOFV2 <- diff(range(bathy.P$predicted)) 
bathy.E.NOFV2

# plot bathy effect 
predicted.fit <- predict.gam(Model.NFOV2$gam, newdata=effect.dat, type='response', se.fit=T)
effect.dat.fit <- cbind(effect.dat, predicted)

predicts.NFOV2.bathy = effect.dat.fit%>%data.frame(predicted.fit)%>%
  group_by(bathymetry)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.nofov2.bathy<- ggplot() +
  ylab("Predicted Abundance of Non-target Species")+
  xlab('Bathymetry (m)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.non.target,aes(x=bathymetry,y=response), colour="lightblue", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.NFOV2.bathy,aes(x=bathymetry,y=response), colour='darkblue', alpha=0.75)+
  geom_line(data=predicts.NFOV2.bathy,aes(x=bathymetry,y=response - se.fit), colour='darkblue', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.NFOV2.bathy,aes(x=bathymetry,y=response + se.fit), colour='darkblue', linetype="dashed",alpha=0.75)+
  theme_classic()+
  xlim(-180,-52)+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.nofov2.bathy


# aspect effect 
aspect.P <- effect.dat%>%
  group_by(cube.Aspect)%>%
  summarise_at(vars(predicted), list(mean))

aspect.E.NFOV2 <- diff(range(aspect.P$predicted))
aspect.E.NFOV2

predicts.legal.aspect = effect.dat%>%data.frame(predicted.fit)%>%
  group_by(cube.Aspect)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.NFOV2.aspect<- ggplot() +
  ylab("Predicted Abundance of Non-target Species")+
  xlab('Cubed Aspect (degrees)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.non.target,aes(x=cube.Aspect,y=response), colour="lightgreen", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response), colour='darkgreen', alpha=0.75)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response - se.fit), colour='darkgreen', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.legal.aspect,aes(x=cube.Aspect,y=response + se.fit), colour='darkgreen', linetype="dashed",alpha=0.75)+
  theme_classic()+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.NFOV2.aspect

# status effect 
status.P.2 <- effect.dat%>%
  group_by(status)%>%
  summarise_at(vars(predicted), list(mean))

status.E.NFOV2 <- diff(range(status.P.2$predicted))
status.E.NFOV2

# Plot it 
predicts.legal.status = effect.dat%>%data.frame(predicted.fit)%>%
  group_by(status)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

status.order <- c("F", "NT")
ggmod.gam.status<- ggplot(aes(x=status,y=response,colour=status), data=full.non.target) +
  ylab("Predicted Abundance of Non-target Species")+
  xlab('Status')+
  #   ggtitle(substitute(italic(name)))+
  #scale_fill_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  scale_x_discrete(limits = status.order)+
  #geom_bar(stat = "identity")+
  geom_errorbar(data=predicts.legal.status,aes(ymin = response-se.fit,ymax = response+se.fit),colour="lightcoral",width = 0.5) +
  theme_classic()+
  Theme1
ggmod.gam.status


# use the bayesian models to calculate a posterior distribution of the effect size
setwd(b.dir)
effect.dat <- expand.grid(cube.Aspect=seq(min(full.non.target$cube.Aspect),max(full.non.target$cube.Aspect),length.out = 20),
                          bathymetry=seq(min(full.non.target$bathymetry),max(full.non.target$bathymetry), length.out=20),
                          status=unique(full.non.target$status))
effect.dat.sim <- effect.dat
gg <- t(posterior_predict(stangam.s.NFOV2, effect.dat.sim, re.form=NA)) # Returning in fish
#colnames(gg) <- paste("sim",1:ncol(gg),sep="_")
gg.dat <- cbind(effect.dat.sim,gg)
str(gg.dat)

gg.dat.long <- gg.dat%>%
  gather(sim, predicted.sim, 4:1005)

gg.dat.long$sim <- as.factor(gg.dat.long$sim)
str(gg.dat.long)

# bathy effect 
bathy.sim.NFOV2 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("Bathymetry (m)")
colnames(bathy.sim.NFOV2) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.dat.long
  bathy.sim.2 <- dplyr::filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- dplyr::group_by(bathy.sim.2, bathymetry)
  bathy.sim.4 <- dplyr::summarise_at(bathy.sim.3, vars(predicted.sim), list(mean))
  print(bathy.sim.5 <- diff(range(bathy.sim.4$predicted.sim)))
  bathy.sim.NFOV2 <- rbind(bathy.sim.NFOV2, bathy.sim.5)
  
}

#bathy.sim

write.csv(bathy.sim, "bayesian.bathy.predictions.nontarget.NFOV2.csv")

# aspect effect 
aspect.sim.NFOV2 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("Cubed Aspect (degrees)")
colnames(aspect.sim.NFOV2) <- x

for(i in 1:1002){
  aspect.sim.1 <- gg.dat.long
  aspect.sim.2 <- filter(aspect.sim.1, sim==i)
  aspect.sim.3 <- group_by(aspect.sim.2, cube.Aspect)
  aspect.sim.4 <- summarise_at(aspect.sim.3, vars(predicted.sim), list(mean))
  print(aspect.sim.5 <- diff(range(aspect.sim.4$predicted.sim)))
  aspect.sim.NFOV2 <- rbind(aspect.sim.NFOV2, aspect.sim.5)
  
}

#aspect.sim

write.csv(aspect.sim, "bayesian.aspect.predictions..nontarget.NFOV2.csv")

# status effect 
status.sim.NFOV2 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("Status")
colnames(status.sim.NFOV2) <- x

for(i in 1:1002){
  status.sim.1 <- gg.dat.long
  status.sim.2 <- filter(status.sim.1, sim==i)
  status.sim.3 <- group_by(status.sim.2, status)
  status.sim.4 <- summarise_at(status.sim.3, vars(predicted.sim), list(mean))
  print(status.sim.5 <- diff(range(status.sim.4$predicted.sim)))
  status.sim.NFOV2 <- rbind(status.sim.NFOV2, status.sim.5)
  
}

#aspect.sim

write.csv(aspect.sim, "bayesian.status.predictions.nontarget.NFOV2.csv")


# Plot effect sizes 
bathy.sim.NFOV2 <- bathy.sim.NFOV2[-1,]
bathy.sim.NFOV2 <- as.data.frame(bathy.sim.NFOV2)
x <- c("Bathymetry (m)")
colnames(bathy.sim.NFOV2) <- x

aspect.sim.NFOV2 <- aspect.sim.NFOV2[-1,]
aspect.sim.NFOV2 <- as.data.frame(aspect.sim.NFOV2)
x <- c("Cube Aspect (Degrees)")
colnames(aspect.sim.NFOV2) <- x

status.sim.NFOV2 <- status.sim.NFOV2[-1,]
status.sim.NFOV2 <- as.data.frame(status.sim.NFOV2)
x <- c("Status")
colnames(status.sim.NFOV2) <- x

full.data.NFOV2 <- as.data.frame(cbind(aspect.sim.NFOV2, bathy.sim.NFOV2, status.sim.NFOV2))

full.data.long.NFOV2 <- full.data.NFOV2%>%
  gather(variable, effect.size, 1:3)

colours <- c('#619CFF', '#00BA38', '#F8766D')  #9590FF
effect.plot.NFOV2 <- ggplot(full.data.long.NFOV2, aes(x = effect.size, fill = variable, colour = variable)) + geom_density(alpha = 0.5) +
  scale_fill_manual(values=colours)+
  scale_colour_manual(values=colours)+
  geom_vline(xintercept = 0.8285865, color = "steelblue", size=0.75)+
  geom_vline(xintercept = 1.543569, color = "tomato3", size=0.75)+
  geom_vline(xintercept = 9.955939, color = "springgreen4", size=0.75)+
  xlim(-10,40)+
  theme_classic()+
  labs(y="Density", x="Effect Size")+
  Theme1
effect.plot.NFOV2

## Make prediction plots

# bathy effect 
bathy.means.NFOV2 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.mean")
colnames(bathy.means.NFOV2) <- x

for(i in unique(gg.dat.long$bathymetry)){
  bathy.sim.1 <- gg.dat.long
  bathy.sim.2 <- dplyr::filter(bathy.sim.1, bathymetry==i)
  print(bathy.sim.3 <- dplyr::summarise_at(bathy.sim.2, vars(predicted.sim), list(mean)))
  bathy.means.NFOV2 <- rbind(bathy.means.NFOV2, bathy.sim.3$predicted.sim)
}

bathy.means.NFOV2 <- bathy.means.NFOV2[-1,]
depths <- unique(gg.dat$bathymetry)
bathy.means.NFOV2 <- as.data.frame(cbind(bathy.means.NFOV2,depths))

# Get means for each value of bathymetry for each simulation 
bathy.mean.sim.NFOV2 <- data.frame(matrix(ncol=1, nrow=20))
x <- c("bathy.mean.sim")
colnames(bathy.mean.sim.NFOV2) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.dat.long
  bathy.sim.2 <- filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- group_by(bathy.sim.2, bathymetry)
  print(bathy.sim.4 <- summarise_at(bathy.sim.3, vars(predicted.sim), list(mean)))
  bathy.mean.sim.NFOV2 <- cbind(bathy.mean.sim.NFOV2, bathy.sim.4)
}

bathy.mean.full.NFOV2 <- bathy.mean.sim.NFOV2[ ,seq(3, ncol(bathy.mean.sim.NFOV2), 2)]
bathy.mean.full.NFOV2$bathymetry <- bathy.mean.sim.NFOV2[,2]

bathy.mean.long.NFOV2 <- bathy.mean.full.NFOV2%>%
  gather(sim, mean, 1:1002)


# Plot the results bathy
predict.plot.bathy.NFOV2 <- ggplot() +
  geom_line(data=bathy.mean.long.NFOV2, aes(x = bathymetry, y = mean, group=sim), colour='lightblue'
            , alpha=0.2)+
  geom_line(data=bathy.means.NFOV2, aes(x = depths, y = bathy.means.NFOV2), colour='darkblue') +
  labs(y="Predicted Abundance of Non-Target Fish", x="Bathymetry (m)")+
  theme_classic()+
  xlim(-180,-51)+
  Theme1
predict.plot.bathy.NFOV2

# aspect effect 
aspect.means.NFOV2 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("aspect.mean")
colnames(aspect.means.NFOV2) <- x

for(i in unique(gg.dat.long$cube.Aspect)){
  aspect.sim.1 <- gg.dat.long
  aspect.sim.2 <- dplyr::filter(aspect.sim.1, cube.Aspect==i)
  print(aspect.sim.3 <- dplyr::summarise_at(aspect.sim.2, vars(predicted.sim), list(mean)))
  aspect.means.NFOV2 <- rbind(aspect.means.NFOV2, aspect.sim.3$predicted.sim)
}

aspect.means.NFOV2 <- aspect.means.NFOV2[-1,]
aspect <- unique(gg.dat$cube.Aspect)
aspect.means.NFOV2 <- as.data.frame(cbind(aspect.means.NFOV2,aspect))

# Get means for each value of aspect for each simulation 
aspect.mean.sim.NFOV2 <- data.frame(matrix(ncol=1, nrow=20))
x <- c("aspect.mean.sim")
colnames(aspect.mean.sim.NFOV2) <- x

for(i in 1:1002){
  aspect.sim.1 <- gg.dat.long
  aspect.sim.2 <- filter(aspect.sim.1, sim==i)
  aspect.sim.3 <- group_by(aspect.sim.2, cube.Aspect)
  print(aspect.sim.4 <- summarise_at(aspect.sim.3, vars(predicted.sim), list(mean)))
  aspect.mean.sim.NFOV2 <- cbind(aspect.mean.sim.NFOV2, aspect.sim.4)
}

aspect.mean.full.NFOV2 <- aspect.mean.sim.NFOV2[ ,seq(3, ncol(aspect.mean.sim.NFOV2), 2)]
aspect.mean.full.NFOV2$cube.Aspect <- aspect.mean.sim.NFOV2[,2]

aspect.mean.long.NFOV2 <- aspect.mean.full.NFOV2%>%
  gather(sim, mean, 1:1002)

# Plot the results aspect
predict.plot.aspect <- ggplot() +
  geom_line(data=aspect.mean.long.NFOV2, aes(x = cube.Aspect, y = mean, group=sim), colour='lightgreen',
            alpha=0.2) +
  geom_line(data=aspect.means.NFOV2, aes(x = aspect, y = aspect.means.NFOV2), colour='darkgreen') +
  labs(y="Predicted Abundance Non-Target Species", x="Cube Aspect (degrees)")+
  theme_classic()+
  Theme1
predict.plot.aspect

# status effect 
status.means.NFOV2 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("status.mean")
colnames(status.means.NFOV2) <- x

for(i in unique(gg.dat.long$status)){
  status.sim.1 <- gg.dat.long
  status.sim.2 <- dplyr::filter(status.sim.1, status==i)
  print(status.sim.3 <- dplyr::summarise_at(status.sim.2, vars(predicted.sim), list(mean)))
  status.means.NFOV2 <- rbind(status.means.NFOV2, status.sim.3$predicted.sim)
}

status.means.NFOV2 <- status.means.NFOV2[-1,]
status <- unique(gg.dat$status)
status.means.NFOV2 <- as.data.frame(cbind(status.means.NFOV2,status))

# Get means for each value of status for each simulation 
status.mean.sim.NFOV2 <- data.frame(matrix(ncol=1, nrow=20))
x <- c("status.mean.sim")
colnames(status.mean.sim.NFOV2) <- x

for(i in 1:1002){
  status.sim.1 <- gg.dat.long
  status.sim.2 <- filter(status.sim.1, sim==i)
  status.sim.3 <- group_by(status.sim.2, status)
  print(status.sim.4 <- summarise_at(status.sim.3, vars(predicted.sim), list(mean)))
  status.mean.sim.NFOV2 <- cbind(status.mean.sim.NFOV2, status.sim.4)
}

status.mean.full.NFOV2 <- status.mean.sim.NFOV2[ ,seq(3, ncol(status.mean.sim.NFOV2), 2)]
status.mean.full.NFOV2$status <- status.mean.sim.NFOV2[,2]

status.mean.long.NFOV2 <- status.mean.full.NFOV2%>%
  gather(sim, mean, 1:1002)

# Plot the results status
jitter <- position_jitter(width = 0.06, height = 0.1)
predict.plot.status <- ggplot() +
  geom_point(data=status.mean.long.NFOV2,  position = jitter, aes(x = status, y = mean, group=sim), colour='lightcoral',
             alpha=0.2) +
  geom_boxplot(data=status.mean.long.NFOV2, aes(x=status, y=mean), colour='black', fill=NA, width=0.14)+
  geom_point(data=status.means.NFOV2, aes(x = status, y = status.means.NFOV2), colour='darkred', shape=17, size=3) +
  labs(y="Predicted Abundance Legal Sized Fish", x="Status")+
  #geom_jitter()+
  theme_classic()+
  Theme1
predict.plot.status

#### Model 1 FOV ####
effect.dat <- expand.grid(bathymetry=seq(min(full.non.target$bathymetry),max(full.non.target$bathymetry), length.out=20),
                          sd.relief=seq(min(full.non.target$sd.relief),max(full.non.target$sd.relief), length.out=20))

Model.FOV1 <- uGamm(response~s(bathymetry,k=5,bs='cs')+ s(sd.relief, k=5, bs='cs'),
                     random=~(1|site)+(1|ID), 
                     family=poisson(), data=use.dat, lme4=TRUE)

gam.check(Model.FOV1$gam)

stangam.s.FOV1  <- stan_gamm4(response~s(bathymetry, k = 5, bs = "cs") + s(sd.relief, k=5, bs='cs'),
                               random=~(1|site)+(1|ID), adapt_delta = 0.99,
                               data=use.dat, chains=3, cores=3, iter=41000, warmup=40000, thin=3,
                               family=poisson)

# calculate frequentist effects - link function is log 
predicted <- predict.gam(Model.FOV1$gam, newdata=effect.dat, type='response')
effect.dat <- cbind(effect.dat, predicted)


# bathy effect 
bathy.P <- effect.dat%>%
  dplyr::group_by(bathymetry)%>%
  summarise_at(vars(predicted), list(mean))

bathy.E.FOV1 <- diff(range(bathy.P$predicted)) 
bathy.E.FOV1

# plot bathy effect 
predicted.fit <- predict.gam(Model.FOV1$gam, newdata=effect.dat, type='response', se.fit=T)
effect.dat.fit <- cbind(effect.dat, predicted.fit)

predicts.FOV1.bathy = effect.dat.fit%>%data.frame(predicted.fit)%>%
  group_by(bathymetry)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.FOV1.bathy<- ggplot() +
  ylab("Predicted Abundance of Non-fished Species")+
  xlab('Bathymetry (m)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.non.target,aes(x=bathymetry,y=response), colour="skyblue", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.FOV1.bathy,aes(x=bathymetry,y=response), colour='deepskyblue4')+
  geom_line(data=predicts.FOV1.bathy,aes(x=bathymetry,y=response - se.fit), colour='deepskyblue4', linetype="dashed")+
  geom_line(data=predicts.FOV1.bathy,aes(x=bathymetry,y=response + se.fit), colour='deepskyblue4', linetype="dashed")+
  theme_classic()+
  geom_rug(data=full.non.target, aes(x=bathymetry),colour="slategrey")+
  xlim(-180,-52)+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.FOV1.bathy


# relief effect 
relief.P <- effect.dat%>%
  group_by(sd.relief)%>%
  summarise_at(vars(predicted), list(mean))

relief.E.FOV1 <- diff(range(relief.P$predicted))
relief.E.FOV1

predicts.FOV1.relief = effect.dat.fit%>%data.frame(predicted.fit)%>%
  group_by(sd.relief)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.FOV1.relief<- ggplot() +
  ylab("Predicted Abundance of Non-fished Species")+
  xlab('Standard Deviation of Relief')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.non.target,aes(x=sd.relief,y=response), colour='goldenrod1', alpha=0.5, size=2,show.legend=F)+
  geom_line(data=predicts.FOV1.relief,aes(x=sd.relief,y=response), colour='darkgoldenrod3')+
  geom_line(data=predicts.FOV1.relief,aes(x=sd.relief,y=response - se.fit), colour='darkgoldenrod3', linetype="dashed")+
  geom_line(data=predicts.FOV1.relief,aes(x=sd.relief,y=response + se.fit), colour='darkgoldenrod3', linetype="dashed")+
  geom_rug(data=full.non.target, aes(x=sd.relief),colour="slategrey")+
  theme_classic()+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.FOV1.relief


# use the bayesian models to calculate a posterior distribution of the effect size
setwd(b.dir)
effect.dat <- expand.grid(bathymetry=seq(min(full.non.target$bathymetry),max(full.non.target$bathymetry), length.out=20),
                          sd.relief=seq(min(full.non.target$sd.relief),max(full.non.target$sd.relief), length.out=20))
effect.dat.sim <- effect.dat
gg <- t(posterior_predict(stangam.s.FOV1, effect.dat.sim, re.form=NA)) # Returning in fish
#colnames(gg) <- paste("sim",1:ncol(gg),sep="_")
gg.dat <- cbind(effect.dat.sim,gg)


gg.dat.long <- gg.dat%>%
  gather(sim, predicted.sim, 3:1004)

gg.dat.long$sim <- as.factor(gg.dat.long$sim)
str(gg.dat.long)

# bathy effect 
bathy.sim.FOV1 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("Bathymetry (m)")
colnames(bathy.sim.FOV1) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.dat.long
  bathy.sim.2 <- dplyr::filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- dplyr::group_by(bathy.sim.2, bathymetry)
  bathy.sim.4 <- dplyr::summarise_at(bathy.sim.3, vars(predicted.sim), list(mean))
  print(bathy.sim.5 <- diff(range(bathy.sim.4$predicted.sim)))
  bathy.sim.FOV1 <- rbind(bathy.sim.FOV1, bathy.sim.5)
  
}

#bathy.sim

write.csv(bathy.sim.FOV1, "bayesian.bathy.predictions.nontarget.FOV1.csv")

# relief effect 
relief.sim.FOV1 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("Standard Deviation of Relief")
colnames(relief.sim.FOV1) <- x

for(i in 1:1002){
  relief.sim.1 <- gg.dat.long
  relief.sim.2 <- filter(relief.sim.1, sim==i)
  relief.sim.3 <- group_by(relief.sim.2, sd.relief)
  relief.sim.4 <- summarise_at(relief.sim.3, vars(predicted.sim), list(mean))
  print(relief.sim.5 <- diff(range(relief.sim.4$predicted.sim)))
  relief.sim.FOV1 <- rbind(relief.sim.FOV1, relief.sim.5)
  
}

#relief.sim

write.csv(relief.sim.FOV1, "bayesian.relief.predictions.nontarget.FOV1.csv")


# Plot effect sizes 
relief.sim.FOV1 <- relief.sim.FOV1[-1,]
relief.sim.FOV1 <- as.data.frame(relief.sim.FOV1)
x <- c("Standard Deviation of Relief")
colnames(relief.sim.FOV1) <- x

bathy.sim.FOV1 <- bathy.sim.FOV1[-1,]
bathy.sim.FOV1 <- as.data.frame(bathy.sim.FOV1)
x <- c("Bathymetry (m)")
colnames(bathy.sim.FOV1) <- x

full.data.FOV1 <- as.data.frame(cbind(relief.sim.FOV1, bathy.sim.FOV1))

full.data.long.FOV1 <- full.data.FOV1%>%
  gather(variable, effect.size, 1:2)

colours <- c('skyblue', NA) #00BA3, #F8766D
effect.plot.FOV1 <- ggplot(full.data.long.FOV1, aes(x = effect.size, fill = variable, colour = variable)) + geom_density(alpha = 0.5) +
  scale_fill_manual(values=colours)+
  scale_colour_manual(values=colours)+
  geom_vline(xintercept = 6.170516, color = "deepskyblue4", size=0.75)+
  #geom_vline(xintercept = 0.8757575, color = "tomato3", size=0.75)+
  #geom_vline(xintercept = 9.650483, color = "springgreen4", size=0.75)+
  #geom_vline(xintercept= 19.93631, color= "darkgoldenrod3", size=0.75)+
  xlim(-10,110)+
  theme_classic()+
  labs(y="Density", x="Effect Size")+
  Theme1
effect.plot.FOV1

bathy.E.FOV1

## Make prediction plots

# bathy effect 
bathy.means.FOV1 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.mean")
colnames(bathy.means.FOV1) <- x

for(i in unique(gg.dat.long$bathymetry)){
  bathy.sim.1 <- gg.dat.long
  bathy.sim.2 <- dplyr::filter(bathy.sim.1, bathymetry==i)
  print(bathy.sim.3 <- dplyr::summarise_at(bathy.sim.2, vars(predicted.sim), list(mean)))
  bathy.means.FOV1 <- rbind(bathy.means.FOV1, bathy.sim.3$predicted.sim)
}

bathy.means.FOV1 <- bathy.means.FOV1[-1,]
depths <- unique(gg.dat$bathymetry)
bathy.means.FOV1 <- as.data.frame(cbind(bathy.means.FOV1,depths))

# Get means for each value of bathymetry for each simulation 
bathy.mean.sim.FOV1 <- data.frame(matrix(ncol=1, nrow=20))
x <- c("bathy.mean.sim")
colnames(bathy.mean.sim.FOV1) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.dat.long
  bathy.sim.2 <- filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- group_by(bathy.sim.2, bathymetry)
  print(bathy.sim.4 <- summarise_at(bathy.sim.3, vars(predicted.sim), list(mean)))
  bathy.mean.sim.FOV1 <- cbind(bathy.mean.sim.FOV1, bathy.sim.4)
}

bathy.mean.full.FOV1 <- bathy.mean.sim.FOV1[ ,seq(3, ncol(bathy.mean.sim.FOV1), 2)]
bathy.mean.full.FOV1$bathymetry <- bathy.mean.sim.FOV1[,2]

bathy.mean.long.FOV1 <- bathy.mean.full.FOV1%>%
  gather(sim, mean, 1:1002)


# Plot the results bathy
predict.plot.bathy.FOV1 <- ggplot() +
  geom_line(data=bathy.mean.long.FOV1, aes(x = bathymetry, y = mean, group=sim), colour='skyblue'
            , alpha=0.2)+
  geom_line(data=bathy.means.FOV1, aes(x = depths, y = bathy.means.FOV1), colour='deepskyblue4') +
  labs(y="Predicted Abundance Non-fished Species", x="Bathymetry (m)")+
  theme_classic()+
  geom_rug(data=full.non.target, aes(x=bathymetry),colour="slategrey")+
  xlim(-180,-51)+
  Theme1
predict.plot.bathy.FOV1


# relief effect 
relief.means.FOV1 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("relief.mean")
colnames(relief.means.FOV1) <- x

for(i in unique(gg.dat.long$sd.relief)){
  relief.sim.1 <- gg.dat.long
  relief.sim.2 <- dplyr::filter(relief.sim.1, sd.relief==i)
  print(relief.sim.3 <- dplyr::summarise_at(relief.sim.2, vars(predicted.sim), list(mean)))
  relief.means.FOV1 <- rbind(relief.means.FOV1, relief.sim.3$predicted.sim)
}

relief.means.FOV1 <- relief.means.FOV1[-1,]
relief <- unique(gg.dat$sd.relief)
relief.means.FOV1 <- as.data.frame(cbind(relief.means.FOV1,relief))

# Get means for each value of relief for each simulation 
relief.mean.sim.FOV1 <- data.frame(matrix(ncol=1, nrow=20))
x <- c("relief.mean.sim")
colnames(relief.mean.sim.FOV1) <- x

for(i in 1:1002){
  relief.sim.1 <- gg.dat.long
  relief.sim.2 <- filter(relief.sim.1, sim==i)
  relief.sim.3 <- group_by(relief.sim.2, sd.relief)
  print(relief.sim.4 <- summarise_at(relief.sim.3, vars(predicted.sim), list(mean)))
  relief.mean.sim.FOV1 <- cbind(relief.mean.sim.FOV1, relief.sim.4)
}

relief.mean.full.FOV1 <- relief.mean.sim.FOV1[ ,seq(3, ncol(relief.mean.sim.FOV1), 2)]
relief.mean.full.FOV1$sd.relief <- relief.mean.sim.FOV1[,2]

relief.mean.long.FOV1 <- relief.mean.full.FOV1%>%
  gather(sim, mean, 1:1002)

# Plot the results relief
predict.plot.relief.FOV1 <- ggplot() +
  geom_line(data=relief.mean.long.FOV1, aes(x = sd.relief, y = mean, group=sim), colour='goldenrod1',
            alpha=0.07) +
  geom_line(data=relief.means.FOV1, aes(x = relief, y = relief.means.FOV1), colour='darkgoldenrod3', size=0.85) +
  labs(y="Predicted Abundance Non-fished Species", x="Standard Deviation of Relief")+
  theme_classic()+
  geom_rug(data=full.non.target, aes(x=sd.relief),colour="slategrey")+
  Theme1
predict.plot.relief.FOV1


#### Model 2 FOV ####
effect.dat <- expand.grid(sqrt.reef=seq(min(full.non.target$sqrt.reef),max(full.non.target$sqrt.reef),length.out = 20),
                          bathymetry=seq(min(full.non.target$bathymetry),max(full.non.target$bathymetry), length.out=20),
                          sd.relief=seq(min(full.non.target$sd.relief),max(full.non.target$sd.relief), length.out=20))

Model.FOV2 <- uGamm(response~s(bathymetry,k=5,bs='cs')+ s(sd.relief, k=5, bs='cs')+
                      + s(sqrt.reef,k=5,bs='cs'),
                    random=~(1|site)+(1|ID), 
                    family=poisson(), data=use.dat, lme4=TRUE)

gam.check(Model.FOV2$gam)


stangam.s.FOV2  <- stan_gamm4(response~s(bathymetry, k = 5, bs = "cs")
                                + s(sqrt.reef, k=5,bs='cs'),
                              random=~(1|site)+(1|ID), adapt_delta = 0.99,
                              data=use.dat, chains=3, cores=3, iter=41000, warmup=40000, thin=3,
                              family=poisson)

# calculate frequentist effects - link function is log 
predicted <- predict.gam(Model.FOV2$gam, newdata=effect.dat, type='response')
effect.dat <- cbind(effect.dat, predicted)


# bathy effect 
bathy.P <- effect.dat%>%
  dplyr::group_by(bathymetry)%>%
  summarise_at(vars(predicted), list(mean))

bathy.E.FOV2 <- diff(range(bathy.P$predicted)) 
bathy.E.FOV2

# plot bathy effect 
predicted.fit <- predict.gam(Model.FOV2$gam, newdata=effect.dat, type='response', se.fit=T)
effect.dat.fit <- cbind(effect.dat, predicted)

predicts.FOV2.bathy = effect.dat.fit%>%data.frame(predicted.fit)%>%
  group_by(bathymetry)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.FOV2.bathy<- ggplot() +
  ylab("Predicted Abundance of Non-fished Species")+
  xlab('Bathymetry (m)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.non.target,aes(x=bathymetry,y=response), colour="skyblue", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.FOV2.bathy,aes(x=bathymetry,y=response), colour='deepskyblue4')+
  geom_line(data=predicts.FOV2.bathy,aes(x=bathymetry,y=response - se.fit), colour='deepskyblue', linetype="dashed")+
  geom_line(data=predicts.FOV2.bathy,aes(x=bathymetry,y=response + se.fit), colour='deepskyblue', linetype="dashed")+
  theme_classic()+
  geom_rug(data=full.non.target,aes(x=bathymetry), colour="slategrey")+
  xlim(-180,-52)+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.FOV2.bathy


# reef effect 
reef.P <- effect.dat%>%
  group_by(sqrt.reef)%>%
  summarise_at(vars(predicted), list(mean))

reef.E.FOV2 <- diff(range(reef.P$predicted))
reef.E.FOV2

predicts.FOV2.reef = effect.dat.fit%>%data.frame(predicted.fit)%>%
  group_by(sqrt.reef)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.FOV2.reef<- ggplot() +
  ylab("Predicted Abundance of Non-fished Species")+
  xlab('Percent Reef Cover (square root)')+
  #   ggtitle(substitute(italic(name)))+
  #scale_color_manual(labels = c("Fished", "No-take"),values=c("royalblue2", "slategrey"))+
  geom_jitter(width = 0.25,height = 0)+
  geom_point(data=full.non.target,aes(x=sqrt.reef,y=response), colour="darkslateblue", alpha=0.5, size=2,show.legend=F)+
  geom_line(data=predicts.FOV2.reef,aes(x=sqrt.reef,y=response), colour='mediumpurple4')+
  geom_line(data=predicts.FOV2.reef,aes(x=sqrt.reef,y=response - se.fit), colour='mediumpurple4', linetype="dashed")+
  geom_line(data=predicts.FOV2.reef,aes(x=sqrt.reef,y=response + se.fit), colour='mediumpurple4', linetype="dashed")+
  geom_rug(data=full.non.target,aes(x=sqrt.reef), colour="slategrey")+
  theme_classic()+
  Theme1
#annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.FOV2.reef


# use the bayesian models to calculate a posterior distribution of the effect size
setwd(b.dir)
effect.dat <- expand.grid(bathymetry=seq(min(full.non.target$bathymetry),max(full.non.target$bathymetry), length.out=20),
                          sqrt.reef=seq(min(full.non.target$sqrt.reef),max(full.non.target$sqrt.reef), length.out=20))
effect.dat.sim <- effect.dat
gg <- t(posterior_predict(stangam.s.FOV2, effect.dat.sim, re.form=NA)) # Returning in fish
#colnames(gg) <- paste("sim",1:ncol(gg),sep="_")
gg.dat <- cbind(effect.dat.sim,gg)


gg.dat.long <- gg.dat%>%
  gather(sim, predicted.sim, 3:1004)

gg.dat.long$sim <- as.factor(gg.dat.long$sim)
str(gg.dat.long)

# bathy effect 
bathy.sim.FOV2 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("Bathymetry (m)")
colnames(bathy.sim.FOV2) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.dat.long
  bathy.sim.2 <- dplyr::filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- dplyr::group_by(bathy.sim.2, bathymetry)
  bathy.sim.4 <- dplyr::summarise_at(bathy.sim.3, vars(predicted.sim), list(mean))
  print(bathy.sim.5 <- diff(range(bathy.sim.4$predicted.sim)))
  bathy.sim.FOV2 <- rbind(bathy.sim.FOV2, bathy.sim.5)
  
}

#bathy.sim

write.csv(bathy.sim.FOV2, "bayesian.bathy.predictions.nontarget.FOV2csv")


# reef effect 
reef.sim.FOV2 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("Standard Deviation of reef")
colnames(reef.sim.FOV2) <- x

for(i in 1:1002){
  reef.sim.1 <- gg.dat.long
  reef.sim.2 <- filter(reef.sim.1, sim==i)
  reef.sim.3 <- group_by(reef.sim.2, sqrt.reef)
  reef.sim.4 <- summarise_at(reef.sim.3, vars(predicted.sim), list(mean))
  print(reef.sim.5 <- diff(range(reef.sim.4$predicted.sim)))
  reef.sim.FOV2 <- rbind(reef.sim.FOV2, reef.sim.5)
  
}

#reef.sim

write.csv(reef.sim.FOV2, "bayesian.ramp.predictions.nontarget.FOV2.csv")

# Plot effect sizes 
reef.sim.FOV2 <- reef.sim.FOV2[-1,]
reef.sim.FOV2 <- as.data.frame(reef.sim.FOV2)
x <- c("Standard Deviation of reef")
colnames(reef.sim.FOV2) <- x

bathy.sim.FOV2 <- bathy.sim.FOV2[-1,]
bathy.sim.FOV2 <- as.data.frame(bathy.sim.FOV2)
x <- c("Bathymetry (m)")
colnames(bathy.sim.FOV2) <- x

full.data.FOV2 <- as.data.frame(cbind(reef.sim.FOV2, bathy.sim.FOV2))

full.data.long.FOV2 <- full.data.FOV2%>%
  gather(variable, effect.size, 1:2)

colours <- c(NA, 'darkslateblue') #
effect.plot.FOV2 <- ggplot(full.data.long.FOV2, aes(x = effect.size, fill = variable, colour = variable)) + geom_density(alpha = 0.5) +
  scale_fill_manual(values=colours)+
  scale_colour_manual(values=colours)+
  #geom_vline(xintercept = 4.124249, color = "deepskyblue4", size=0.75)+
  geom_vline(xintercept = 6.686027, color = "mediumpurple4", size=0.75)+
  #geom_vline(xintercept = 9.783041, color = "springgreen4", size=0.75)+
  xlim(-10,100)+
  theme_classic()+
  labs(y="Density", x="Effect Size")+
  Theme1
effect.plot.FOV2

reef.E.FOV2

## Make prediction plots

# bathy effect 
bathy.means.FOV2 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("bathy.mean")
colnames(bathy.means.FOV2) <- x

for(i in unique(gg.dat.long$bathymetry)){
  bathy.sim.1 <- gg.dat.long
  bathy.sim.2 <- dplyr::filter(bathy.sim.1, bathymetry==i)
  print(bathy.sim.3 <- dplyr::summarise_at(bathy.sim.2, vars(predicted.sim), list(mean)))
  bathy.means.FOV2 <- rbind(bathy.means.FOV2, bathy.sim.3$predicted.sim)
}

bathy.means.FOV2 <- bathy.means.FOV2[-1,]
depths <- unique(gg.dat$bathymetry)
bathy.means.FOV2 <- as.data.frame(cbind(bathy.means.FOV2,depths))

# Get means for each value of bathymetry for each simulation 
bathy.mean.sim.FOV2 <- data.frame(matrix(ncol=1, nrow=20))
x <- c("bathy.mean.sim")
colnames(bathy.mean.sim.FOV2) <- x

for(i in 1:1002){
  bathy.sim.1 <- gg.dat.long
  bathy.sim.2 <- filter(bathy.sim.1, sim==i)
  bathy.sim.3 <- group_by(bathy.sim.2, bathymetry)
  print(bathy.sim.4 <- summarise_at(bathy.sim.3, vars(predicted.sim), list(mean)))
  bathy.mean.sim.FOV2 <- cbind(bathy.mean.sim.FOV2, bathy.sim.4)
}

bathy.mean.full.FOV2 <- bathy.mean.sim.FOV2[ ,seq(3, ncol(bathy.mean.sim.FOV2), 2)]
bathy.mean.full.FOV2$bathymetry <- bathy.mean.sim.FOV2[,2]

bathy.mean.long.FOV2 <- bathy.mean.full.FOV2%>%
  gather(sim, mean, 1:1002)


# Plot the results bathy
predict.plot.bathy.FOV2 <- ggplot() +
  geom_line(data=bathy.mean.long.FOV2, aes(x = bathymetry, y = mean, group=sim), colour='skyblue'
            , alpha=0.2)+
  geom_line(data=bathy.means.FOV2, aes(x = depths, y = bathy.means.FOV2), colour='deepskyblue4') +
  labs(y="Predicted Abundance of Non-fished Species", x="Bathymetry (m)")+
  theme_classic()+
  geom_rug(data=full.non.target, aes(x = bathymetry), colour="slategrey")+
  xlim(-180,-51)+
  Theme1
predict.plot.bathy.FOV2

# reef effect 
reef.means.FOV2 <- data.frame(matrix(ncol=1, nrow=1))
x <- c("reef.mean")
colnames(reef.means.FOV2) <- x

for(i in unique(gg.dat.long$sqrt.reef)){
  reef.sim.1 <- gg.dat.long
  reef.sim.2 <- dplyr::filter(reef.sim.1, sqrt.reef==i)
  print(reef.sim.3 <- dplyr::summarise_at(reef.sim.2, vars(predicted.sim), list(mean)))
  reef.means.FOV2 <- rbind(reef.means.FOV2, reef.sim.3$predicted.sim)
}

reef.means.FOV2 <- reef.means.FOV2[-1,]
reef <- unique(gg.dat$sqrt.reef)
reef.means.FOV2 <- as.data.frame(cbind(reef.means.FOV2,reef))

# Get means for each value of reef for each simulation 
reef.mean.sim.FOV2 <- data.frame(matrix(ncol=1, nrow=20))
x <- c("reef.mean.sim")
colnames(reef.mean.sim.FOV2) <- x

for(i in 1:1002){
  reef.sim.1 <- gg.dat.long
  reef.sim.2 <- filter(reef.sim.1, sim==i)
  reef.sim.3 <- group_by(reef.sim.2, sqrt.reef)
  print(reef.sim.4 <- summarise_at(reef.sim.3, vars(predicted.sim), list(mean)))
  reef.mean.sim.FOV2 <- cbind(reef.mean.sim.FOV2, reef.sim.4)
}

reef.mean.full.FOV2 <- reef.mean.sim.FOV2[ ,seq(3, ncol(reef.mean.sim.FOV2), 2)]
reef.mean.full.FOV2$sqrt.reef <- reef.mean.sim.FOV2[,2]

reef.mean.long.FOV2 <- reef.mean.full.FOV2%>%
  gather(sim, mean, 1:1002)

# Plot the results reef
predict.plot.reef.FOV2 <- ggplot() +
  geom_line(data=reef.mean.long.FOV2, aes(x = sqrt.reef, y = mean, group=sim), colour='mediumpurple3',
            alpha=0.02) +
  geom_line(data=reef.means.FOV2, aes(x = reef, y = reef.means.FOV2), colour='mediumpurple4') +
  labs(y="Predicted Abundance Non-fished Species", x="Percent Cover Reef (square root)")+
  theme_classic()+
  geom_rug(data=full.non.target, aes(x=sqrt.reef), colour="slategrey")+
  Theme1
predict.plot.reef.FOV2

### Plot Importance Scores separates by whether it using FOV or not ####
setwd(m.dir)
dat.var.imp <-read.csv("nontarget.imp.nofov.csv")%>% #from local copy
  rename(resp.var=X)%>%
  gather(key=predictor,value=importance,2:ncol(.))%>%
  glimpse()

dat.var.imp.fov  <- read.csv("nontarget.imp.fov.csv")%>% #from local copy
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

# Plotting defaults
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
  mutate(label=ifelse(predictor=="distance.to.ramp"&resp.var=="Sublegal","X",ifelse(predictor=="bathymetry"&resp.var=="Sublegal","X",ifelse(predictor=="sqrt.reef"&resp.var=="Sublegal","X",label))))%>%
  glimpse()

# Plot gg.importance.scores
gg.importance.scores <- ggplot(dat.var.label, aes(x=predictor,y=model,fill=importance))+
  geom_tile(show.legend=T) +
  scale_fill_gradientn(legend_title,colours=c("white", re), na.value = "grey98",
                       limits = c(0, max(dat.var.label$importance)))+
  scale_x_discrete(limits=c("sqrt.slope","cube.Aspect","log.roughness","distance.to.ramp", "status",
                            "mean.reef", "sd.relief","sqrt.reef"),
                   labels=c(
                     "Slope (sqrt)","Aspect (cubed)","Roughness (log)","Distance to Ramp", "Status",
                     "Mean % Reef Cover", "SD reef","% Reef (sqrt)"
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

