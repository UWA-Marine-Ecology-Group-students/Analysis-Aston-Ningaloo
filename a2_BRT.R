library(dismo)
library(RCurl)
library(gbm)
remotes::install_github("JBjouffray/ggBRT", dependencies = FALSE)
library(ggBRT)

working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # sets working directory to where this script is saved (DON'T MOVE THE SCRIPT)

## Set sub directories----
d.dir <- paste(working.dir,"Tidy data",sep="/") 
s.dir <- paste(working.dir,"Spatial",sep="/") # spatial is where I keep spatial data files, rasters and shapefiles
p.dir <- paste(working.dir,"Plots",sep="/")
m.dir <- paste(working.dir,"Model Out BRT", sep="/")

setwd(d.dir)

dat <-read.csv('final.data.csv')

legal.dat <- subset(dat, model=='Legal')
sublegal.dat <- subset(dat, model=='Sublegal')

#Set co-variates. There are no distributional assumptions on these. Make those you want to be factors into factors
factor.vars=c("status","model")
cont.vars=c("bathymetry","TPI", "Slope", "Aspect", "FlowDir",
            "distance.to.ramp", "mean.relief", "sd.relief", "reef")

################## Legal Target Species Model ##################

#STEP 1: you need to find the optimal learning rate, bag fraction and tree complexity for your data

#To find optimal settings we are going to loop through a range of candidates and select the option with the 
#lowest cross validated residual deviance - i.e. the model that fits the data the best. 

# making a grid for loop
res.dev<-c() ## make a res.dev variable to store in grid
grid<-expand.grid(lr=c(0.1, 0.05, 0.01), tc=c(1, 3, 5), bf = c(0.6, 0.7, 1)) ## make grid of lr and tc that gbm.step will run through 
grid$res.dev<-NA ## add red.dev to grid
grid ## check grid

# use gbm.step to find optimal trees with a range of lr and tc combinations 
#this will take a while...
for(i in 1:nrow(grid)) {
  dp.stepW_dem<-dismo::gbm.step(data=legal.dat,
                                gbm.x=c(factor.vars, cont.vars), ## explanatory
                                gbm.y="target.fish",## response
                                lr=grid[i,"lr"], ## ref to grid 
                                tc=grid[i,"tc"], ## ref to grid 
                                family="poisson", ## distribution family  - see what else you can use?
                                bag.fraction=grid[i,"bf"])
  grid[i, "res.dev"]<-dp.stepW_dem$self.statistics$mean.resid} ### store res.dev in grid 

#STEP 2 re-fit the best model

best<-grid[grid$res.dev == min(grid$res.dev), ] #shows parameters for best fitting model

model <- dismo::gbm.step(data=legal.dat,
                         gbm.x=c(factor.vars, cont.vars), ## explanatory
                         gbm.y="target.fish",## response
                         lr=best$lr, ## ref to grid 
                         tc=best$tc, ## ref to grid 
                         family="poisson",
                         bag.fraction=best$tc)


#STEP 3 plots
summary(model) #Shows importance of different variables
dismo::gbm.plot(model) #Shows variable effects



###################### Sublegal Target Species Model #########################
#STEP 1: you need to find the optimal learning rate, bag fraction and tree complexity for your data

#To find optimal settings we are going to loop through a range of candidates and select the option with the 
#lowest cross validated residual deviance - i.e. the model that fits the data the best. 

# making a grid for loop
res.dev<-c() ## make a res.dev variable to store in grid
grid<-expand.grid(lr=c(0.01, 0.005, 0.001), tc=c(1, 3, 5), bf = c(0.6, 0.7, 1)) ## make grid of lr and tc that gbm.step will run through 
grid$res.dev<-NA ## add red.dev to grid
grid ## check grid

# use gbm.step to find optimal trees with a range of lr and tc combinations 
#this will take a while...
for(i in 1:nrow(grid)) {
  dp.stepW_dem<-dismo::gbm.step(data=sublegal.dat,
                                gbm.x=c(factor.vars, cont.vars), ## explanatory
                                gbm.y="target.fish",## response
                                lr=grid[i,"lr"], ## ref to grid 
                                tc=grid[i,"tc"], ## ref to grid 
                                bag.fraction=grid[i,"bf"], # ref to grid
                                family="poisson",
                                n.trees=50,
                                step.size = 5) 
  grid[i, "res.dev"]<-dp.stepW_dem$self.statistics$mean.resid} ### store res.dev in grid 

#STEP 2 re-fit the best model

best<-grid[grid$res.dev == min(grid$res.dev), ] #shows parameters for best fitting model

model <- dismo::gbm.step(data=sublegal.dat,
                         gbm.x=c(factor.vars, cont.vars), ## explanatory
                         gbm.y="target.fish",## response
                         lr=0.005, ## ref to grid 
                         tc=1, ## ref to grid 
                         family="poisson",
                         bag.fraction=0.6,
                         n.trees=50,
                         step.size = 5)


#STEP 3 plots
summary(model) #Shows importance of different variables
dismo::gbm.plot(model) #Shows variable effects



############# All data together #############
#STEP 1: you need to find the optimal learning rate, bag fraction and tree complexity for your data

#To find optimal settings we are going to loop through a range of candidates and select the option with the 
#lowest cross validated residual deviance - i.e. the model that fits the data the best. 

# making a grid for loop
res.dev<-c() ## make a res.dev variable to store in grid
grid<-expand.grid(lr=c(0.1, 0.05, 0.01), tc=c(1, 3, 5), bf = c(0.6, 0.7, 1)) ## make grid of lr and tc that gbm.step will run through 
grid$res.dev<-NA ## add red.dev to grid
grid ## check grid

# use gbm.step to find optimal trees with a range of lr and tc combinations 
#this will take a while...
for(i in 1:nrow(grid)) {
  dp.stepW_dem<-dismo::gbm.step(data=dat,
                                gbm.x=c(factor.vars, cont.vars), ## explanatory
                                gbm.y="target.fish",## response
                                lr=grid[i,"lr"], ## ref to grid 
                                tc=grid[i,"tc"], ## ref to grid 
                                family="poisson", ## distribution family 
                                bag.fraction=grid[i,"bf"])
                                #n.trees=50,
                                #step.size=10)
  grid[i, "res.dev"]<-dp.stepW_dem$self.statistics$mean.resid} ### store res.dev in grid 

#STEP 2 re-fit the best model

best<-grid[grid$res.dev == min(grid$res.dev), ] #shows parameters for best fitting model

model <- dismo::gbm.step(data=dat,
                         gbm.x=c(factor.vars, cont.vars), ## explanatory
                         gbm.y="target.fish",## response
                         lr=best$lr, ## ref to grid 
                         tc=best$tc, ## ref to grid 
                         family="poisson",
                         bag.fraction=best$tc)

#STEP 3 plots
summary(model) #Shows importance of different variables
dismo::gbm.plot(model) #Shows variable effects

#Bootstrapping to get CIs
brt1.prerun<- plot.gbm.4list(model)
brt1.boot <- gbm.bootstrap.functions(model, list.predictors=brt1.prerun, n.reps=1000)

bathy <- ggPD_boot(model, predictor="bathymetry", list.4.preds=brt1.prerun, 
                   booted.preds=brt1.boot$function.preds, type.ci = "ribbon",
                   col.line = "royalblue2", rug = T)

legal.sublegal <- ggPD_boot(model, predictor="model", list.4.preds=brt1.prerun, 
                            booted.preds=brt1.boot$function.preds, type.ci = "lines",
                            col.line = "royalblue2", rug = T)

sd.relief <- ggPD_boot(model, predictor="sd.relief", list.4.preds=brt1.prerun, 
                       booted.preds=brt1.boot$function.preds, type.ci = "ribbon",rug = T)

slope <- ggPD_boot(model, predictor="Slope", list.4.preds=brt1.prerun, 
                   booted.preds=brt1.boot$function.preds, type.ci = "ribbon",
                   col.line = "royalblue2", rug = T)

aspect <- ggPD_boot(model, predictor="Aspect", list.4.preds=brt1.prerun, 
                    booted.preds=brt1.boot$function.preds, type.ci = "ribbon",
                    col.line = "royalblue2", rug = T)

status <- ggPD_boot(model, predictor="status", list.4.preds=brt1.prerun,
                    booted.preds=brt.boot$function.preds, type.ci = "ribbon",
                    col.line = "royalblue2", rug = T)

#ramp <- ggPD_boot(model, predictor="distance.to.ramp", list.4.preds=brt1.prerun, 
#                  booted.preds=brt1.boot$function.preds, type.ci = "ribbon",rug = T)
