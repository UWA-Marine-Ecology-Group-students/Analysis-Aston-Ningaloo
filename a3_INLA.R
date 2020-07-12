library(raster)
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)
library(dplyr)
library(sp) 
library(fields)
library(ggplot2)
library(devtools) #From GitHub
install_github('timcdlucas/INLAutils')
library(INLAutils)
library(rgdal)


rm(list=ls())

##Set working directory----
## Set work directory----
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # sets working directory to where this script is saved (DON'T MOVE THE SCRIPT)

## Set sub directories----
d.dir <- paste(working.dir,"Tidy data",sep="/") 
s.dir <- paste(working.dir,"Spatial",sep="/") 
p.dir <- paste(working.dir,"Plots",sep="/")
m.dir <- paste(working.dir,"Model Out INLA", sep="/") 

######## Loading data and setting covariates ########

## Load data
setwd(d.dir)
data<- read.csv('final.data.csv')

data.legal <- subset(data, model=='Legal')

#remove NA sites 
data<- data%>%
  filter(!sample%in%c("8.05","10.09","10.12","16.03"))%>%
  dplyr::select(!sand)%>%
  dplyr::select(!TRI)%>%
  dplyr::select(!Roughness)

# Set your covariates/spatial data
covariates <- c("bathymetry","TPI","Slope","Aspect","FlowDir","mean.relief",
                "sd.relief","reef","distance.to.ramp")
factors <- c("status")
covariates <- na.exclude(covariates)


########## Setting up a mesh #########
# This creates a mesh of discrete sampling locations that allows the model to estimate the spatial
# autocorrelation in the data 

# changing the mesh to reflect the coastline of WA
# data <- data%>%
#   rename(easting=latitude)%>%
#   rename(northing=longitude)
# 
# Locations <- data%>%
#   dplyr::select("easting", "northing")
# 
# coordinates(Locations) <- ~ northing + easting
# 
# setwd(s.dir)
# 
# interest <- readOGR(dsn = ".", layer = "Area of Interest")
# coastline <- readOGR(dsn=".", layer = "Coastline")
# 
# plot(coastline)
# plot(interest)
# 
# interest <- spTransform(interest, CRS("+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# coastline <- spTransform(coastline, CRS("+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# pll <- Polygon(interest@polygons[[1]]@Polygons[[1]]@coords, hole = F)
# hole <- Polygon(coastline@polygons[[1]]@Polygons[[1]]@coords, hole = T)
# new_sp <- SpatialPolygons(list(Polygons(append(hole, pll),'1')))

# plot(new_sp)


# Creating the mesh
# my.mesh.land <- inla.mesh.2d(Locations, boundary=new_sp, max.edge = c(3500))

my.mesh <- inla.mesh.2d(Locations, max.edge=c(6000))

plot(my.mesh)
points(Locations, col = 2, pch = 16, cex = 0.2)

######### Setting the priors ########

spat.priors <- inla.spde2.pcmatern(
  mesh=my.mesh, alpha=2, ### alpha can usually be left at 2 unless using a 3d mesh
  constr=TRUE,### so that random effects are constrained to sum to zero
  prior.range=c(3000, 0.05), ### P(practic.range<3000)=0.05
  #This is saying that the probability the range is less than 10000 is 5%
  prior.sigma=c(5, 0.1)) ### P( sigma>5)=0.01  sigma is sd.
  #This is saying that probability the variance is greater than 1 is less than 1%

######### Making the A matrix ########
# This translates the spatial locations on the mesh to vectors in the model
A.matrix <- inla.spde.make.A(my.mesh, loc = as.matrix(data[,c("easting","northing")]))

######## Creating a stack #######
# This is a combination of the A. matrix and the mesh 

mesh.index <- inla.spde.make.index(name = "iSpat", n.spde = spat.priors$n.spde)
my.stack <- inla.stack(data=list(y=data$target.fish), A=list(A.matrix, 1),
                       tag = "preds",
                        effects=list(c(mesh.index, list( intercept=1)),#for the spatial effects and intercept
                                      list(bathymetry=data$bathymetry, TPI=data$TPI, Slope=data$Slope, 
                                           Aspect=data$Aspect, FlowDir=data$FlowDir,
                                           mean.relief=data$mean.relief, sd.relief=data$sd.relief, 
                                           reef=data$reef, distance.to.ramp=data$distance.to.ramp, 
                                           status=data$status))) #covariates

###### Run Full Model #########
f.s <- y ~ -1 + intercept + bathymetry + TPI + Slope + Aspect + FlowDir + mean.relief + sd.relief + 
  reef + distance.to.ramp + status + f(data$site, model='iid') + f(iSpat, model=spat.priors)

fm <- inla(f.s,
            family = "zeroinflatedpoisson1", #might use zeroinflated poisson as you've got a lot of zeros, link is the same
            data = inla.stack.data(my.stack),
            verbose=FALSE,
            control.predictor=list(A=inla.stack.A(my.stack), compute=TRUE, link=1),
            control.fixed = list(mean=0, prec=0.2),
            control.results = list(return.marginals.random = TRUE, return.marginals.predictor = TRUE), 
            control.compute=list(config = TRUE, dic=TRUE)
)

summary(fm)

plot(fm)

######## Plotting the results ###########
#There are no p-values with INLA you look at significance by finding the overlap of the 0.025 and the 0.975
# Quantiles with 0 - those that don't overlap may be significant 
summary <- as.data.frame(summary(fm)$fixed)

summary <- summary%>%
  rename(lower=`0.025quant`)%>%
  rename(upper=`0.975quant`)%>%
  rename(estimate=mean)%>%
  mutate(factor=c("intercept", "bathy", "TPI", "slope", "aspect", "flowdir", "mean.relief", "sd.relief",
           "reef", "ramp", "fished", "NTZ"))

effects <- ggplot(summary,
            aes(x = estimate,
                y = factor, color=factor)) +
            geom_point() +
              geom_errorbar(aes(xmin=lower, xmax=upper), width=.5)+
              geom_vline(xintercept = 0)+
  theme_bw()
effects

#So SD.reef, reef, ramp(?), flowdir, bathy look promising 
fm$summary.fixed

###### Projecting the posterior mean and SD of the spatial field########
rang <- apply(my.mesh$loc[, c(1, 2)], 2, range)
proj <- inla.mesh.projector(my.mesh,
                            xlim = rang[, 1], ylim = rang[, 2],
                            dims = c(300, 300))

mean_s <- inla.mesh.project(proj, fm$summary.random$iSpat$mean)
sd_s <- inla.mesh.project(proj, fm$summary.random$iSpat$sd)


df <- expand.grid(x = proj$x, y = proj$y)
df$mean_s <- as.vector(mean_s)
df$sd_s <- as.vector(sd_s)

library(viridis)
library(cowplot)

gmean <- ggplot(df, aes(x = x, y = y, fill = mean_s)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

gsd <- ggplot(df, aes(x = x, y = y, fill = sd_s)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

plot_grid(gmean, gsd)


####### Plotting the residuals ######
# If the model is well calibrated then the bins should be the same size. The convex model suggests that
# the model is overconfident in it's predictions and the line graph suggests we've got both over and 
# under fitting 
ggplot_inla_residuals(fm, data$target.fish, binwidth = 0.1)

####### Stepwise Model selection ########
bestmodel <- INLAstep(fam1="zeroinflatedpoisson1", data, spde = spat.priors,  in_stack=my.stack, 
                  invariant = "-1+f(iSpat, model=spat.priors)", direction="backwards", 
                  include=6:14, y='y', powerl = 1, inter=3)

autoplot(bestmodel$best_model, which = c(1, 5), CI = TRUE)


######### Predicting fish abundance using subset of data to test priors ########
# Randomly selected 25% of the sites to be used as predictive sites rather than used in the model
training <- data%>%
  filter(!sample%in%c("19.03","10.11","13.03","14.03", "4.01", "1.06", "12.06", "10.13", "19.02",
                      "4.12", "6.12", "30.02", "15.02", "30.01", "16.05", "3.06", "1.04", "16.02",
                      "9.01", "23.06", "13.08", "2.02", "14.06", "6.03", "14.11", "6.02", "23.02",
                      "14.12", "30.06", "5.08", "18.03", "18.07")) 
predictive <- data%>%
  filter(sample%in%c("19.03","10.11","13.03","14.03", "4.01", "1.06", "12.06", "10.13", "19.02",
                      "4.12", "6.12", "30.02", "15.02", "30.01", "16.05", "3.06", "1.04", "16.02",
                      "9.01", "23.06", "13.08", "2.02", "14.06", "6.03", "14.11", "6.02", "23.02",
                      "14.12", "30.06", "5.08", "18.03", "18.07")) 

# Set your covariates/spatial data
covariates <- c("bathymetry","TPI","Slope","Aspect","FlowDir","mean.relief",
                "sd.relief","reef","distance.to.ramp")
factors <- c("status")
covariates <- na.exclude(covariates)


####### Make a mesh using the training data set ########
Locations.train <- training%>%
  dplyr::select("easting", "northing")

pred.mesh <- inla.mesh.2d(Locations.train, max.edge=c(6000))

plot(pred.mesh)
points(Locations.train, col = 2, pch = 16, cex = 0.2)

######### Setting the priors for prediction ########
spat.priors.predict <- inla.spde2.pcmatern(
  mesh=pred.mesh, alpha=2, ### alpha can usually be left at 2 unless using a 3d mesh
  constr=TRUE,### so that random effects are constrained to sum to zero
  prior.range=c(3000, 0.05), ### P(practic.range<3000)=0.05
  #This is saying that the probability the range is less than 10000 is 5%
  prior.sigma=c(5, 0.1)) ### P( sigma>5)=0.01  sigma is sd.
#This is saying that probability the variance is greater than 1 is less than 1%

######### Making the A matrix for prediction ########
# Need one that has the data points and one that doesn't so we can predict over the area
Pred.matrix <- inla.spde.make.A(pred.mesh, loc = as.matrix(training[,c("easting","northing")]))
Pred.matrix.2 <- inla.spde.make.A(mesh = pred.mesh)

######## Creating the stacks #######
# Again need two stacks, one that has all the points in the training data and one that has NA for response

mesh.index <- inla.spde.make.index(name = "iSpat", n.spde = spat.priors.predict$n.spde)

stack.train <- inla.stack(data=list(y=training$target.fish), A=list(Pred.matrix, 1),
                       tag = "train",
                       effects=list(c(mesh.index, list( intercept=1)),#for the spatial effects and intercept
                                    list(bathymetry=training$bathymetry, TPI=training$TPI, Slope=training$Slope, 
                                         Aspect=training$Aspect, FlowDir=training$FlowDir,
                                         mean.relief=training$mean.relief, sd.relief=training$sd.relief, 
                                         reef=training$reef, distance.to.ramp=training$distance.to.ramp, 
                                         status=training$status))) 

stack.pred <- inla.stack(data=list(y=NA), A=list(Pred.matrix.2),
                          tag = "preds",
                          effects=list(c(mesh.index, list(intercept=1))))

StackJoin <- inla.stack(stack.train, stack.pred)

##### Run model ######

f.s <- y ~ -1 + intercept + bathymetry + TPI + Slope + Aspect + FlowDir + mean.relief + sd.relief + 
  reef + distance.to.ramp + status + f(iSpat, model=spat.priors)

fm <- inla(f.s,
           family = "zeroinflatedpoisson1", #might use zeroinflated poisson as you've got a lot of zeros, link is the same
           data = inla.stack.data(StackJoin),
           verbose=FALSE,
           control.predictor=list(A=inla.stack.A(StackJoin), compute=TRUE, link=1),
           control.fixed = list(mean=0, prec=0.2),
           control.results = list(return.marginals.random = TRUE, return.marginals.predictor = TRUE), 
           control.compute=list(config = TRUE, dic=TRUE))


###### Selecting the posterior mean and SD of the response from the model #######
index.pred <- inla.stack.index(StackJoin, "Pred")$data

post.mean.pred <- fm$summary.linear.predictor[index.pred, "mean"]
post.sd.pred <- fm$summary.linear.predictor[index.pred, "sd"]

proj.grid <- inla.mesh.projector(pred.mesh, 
                                 xlim = range(pred.grid[,1]), 
                                 ylim = range(pred.grid[,2]), 
                                 dims = c(ncol,nrow))






