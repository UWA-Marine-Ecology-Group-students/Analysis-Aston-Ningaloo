library(raster)
#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
library(INLA)
inla.setOption(mkl=FALSE)
inla.setOption(pardiso.license = "/Users/Charlotte/Documents/University/Masters/Research/pardiso.lic.rtf")
inla.pardiso.check()
library(dplyr)
library(sp) 
library(fields)
library(ggplot2)
library(devtools) #From GitHub
install_github('timcdlucas/INLAutils')
library(INLAutils)
library(rgdal)
library(RColorBrewer)
install_github('oswaldosantos/INLAOutputs')
library(INLAOutputs)

rm(list=ls())

names(inla.models()$latent$iid$name)
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

str(data)
data$status <- as.factor(data$status)

#remove NA sites 
data<- data%>%
  filter(!sample%in%c("8.05","10.09","10.12","16.03"))%>%
  dplyr::select(!sand)%>%
  dplyr::select(!TRI)%>%
  dplyr::select(!Roughness)%>%
  dplyr::select(!X)%>%
  dplyr::select(!X.1)

# Set your covariates/spatial data
covariates <- c("bathymetry","TPI","Slope","Aspect","FlowDir","mean.relief",
                "sd.relief","reef","distance.to.ramp")
factors <- c("status")
covariates <- na.exclude(covariates)

# Subet data
data.legal <- subset(data, model=='Legal')
data.sublegal <- subset(data, model=='Sublegal')

###################### Legal Model #####################

########## Setting up a mesh #########
# This creates a mesh of discrete sampling locations that allows the model to estimate the spatial
# autocorrelation in the data 

# Rename coordinates
data.legal <- data.legal%>%
  dplyr::rename(easting=latitude)%>%
  dplyr::rename(northing=longitude)

Locations <- data.legal%>%
  dplyr::select("easting", "northing")

# Turn into locations if using the modified coastline
coordinates(Locations) <- ~ northing + easting

setwd(s.dir)
 
interest <- readOGR(dsn = ".", layer = "Area of Interest")
coastline <- readOGR(dsn=".", layer = "Coastline")

plot(coastline)
plot(interest)
 
interest <- spTransform(interest, CRS("+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
coastline <- spTransform(coastline, CRS("+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

pll <- Polygon(interest@polygons[[1]]@Polygons[[1]]@coords, hole = F)
hole <- Polygon(coastline@polygons[[1]]@Polygons[[1]]@coords, hole = T)
new_sp <- SpatialPolygons(list(Polygons(append(hole, pll),'1')))

plot(new_sp)


# Creating the mesh
my.mesh.land <- inla.mesh.2d(Locations, boundary=new_sp, max.edge = c(6000)) #This is the modified one

my.mesh <- inla.mesh.2d(Locations, max.edge=c(3000))

plot(my.mesh)
points(Locations, col = 2, pch = 16, cex = 0.2)

######### Setting the priors ########

spat.priors.legal <- inla.spde2.pcmatern(
  mesh=my.mesh, alpha=2, ### alpha can usually be left at 2 unless using a 3d mesh
  constr=TRUE,### so that random effects are constrained to sum to zero
  prior.range=c(500, 0.05), ### 300 for legal, not sure for sub-legal
  #This is saying that the probability the range is less than 300 is 5%
  prior.sigma=c(0.1, 0.1)) ### 0.5,0.08 for legal, not sure for sub-legal
  #This is saying that probability the variance is greater than 0.5 is less than 8%

######### Making the A matrix ########
# This translates the spatial locations on the mesh to vectors in the model
A.matrix.legal <- inla.spde.make.A(my.mesh, loc = as.matrix(data.legal[,c("easting","northing")]))

######## Creating a stack #######
# This is a combination of the A. matrix and the mesh 

mesh.index <- inla.spde.make.index(name = "iSpat", n.spde = spat.priors.legal$n.spde)
my.stack.legal <- inla.stack(data=list(y=data.legal$target.fish), A=list(A.matrix.legal, 1),
                       tag = "preds",
                       effects=list(c(mesh.index, list( intercept=1)),#for the spatial effects and intercept
                                    list(bathymetry=data.legal$bathymetry, TPI=data.legal$TPI, 
                                         Slope=data.legal$Slope, Aspect=data.legal$Aspect, 
                                         FlowDir=data.legal$FlowDir,distance.to.ramp=
                                           data.legal$distance.to.ramp, 
                                         status=data.legal$status))) #covariates

###### Run Full Model #########
f.s <- y ~ -1 + intercept + bathymetry + TPI + Slope + Aspect + FlowDir +
  distance.to.ramp + status + f(iSpat, model=spat.priors.legal) 
#+ f(data.legal$site, model="iid")

fm <- inla(f.s,
           family = "zeroinflatedpoisson1",
           data = inla.stack.data(my.stack.legal),
           verbose=FALSE,
           control.predictor=list(A=inla.stack.A(my.stack.legal), compute=TRUE, link=1),
           control.fixed = list(mean=0, prec=0.001),
           control.results = list(return.marginals.random = TRUE, return.marginals.predictor = TRUE), 
           control.compute=list(config = TRUE, dic=TRUE, waic=TRUE)
)

summary(fm)

plot(fm)

####### Plotting the residuals and checking the model ######
# If the model is well calibrated then the bins should be the same height. The convex model suggests that
# the model is overconfident in it's predictions and the line graph suggests we've got both over and 
# under fitting 
ggplot_inla_residuals(fm, data.legal$target.fish, binwidth = 0.1)

# This plots the values predeicted by the model against the actual values 
index.pred <- inla.stack.index(my.stack.legal, "preds")$data

post.mean.pred <- fm$summary.fitted.values[index.pred, "mean"]
post.sd.pred <- fm$summary.fitted.values[index.pred, "sd"]

plot(post.mean.pred,data.legal$target.fish)


######## Plotting the results ###########
#There are no p-values with INLA you look at significance by finding the overlap of the 0.025 and the 0.975
# Quantiles with 0 - those that don't overlap may be significant 
summary <- as.data.frame(summary(fm)$fixed)

summary <- summary%>%
  dplyr::rename(lower=`0.025quant`)%>%
  dplyr::rename(upper=`0.975quant`)%>%
  dplyr::rename(estimate=mean)%>%
  dplyr::mutate(factor=c("intercept", "bathy", "TPI", "slope", "aspect", "flowdir",
                  "ramp", "fished", "NTZ"))

effects <- ggplot(summary,
                  aes(x = estimate,
                      y = factor, color=factor)) +
  geom_point() +
  geom_errorbar(aes(xmin=lower, xmax=upper), width=.5)+
  geom_vline(xintercept = 0)+
  theme_bw()
effects

#So reef, ramp(?), flowdir, bathy look promising 
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


####### Stepwise Model selection ########
bestmodel <- INLAstep(fam1="zeroinflatedpoisson1", data.legal, spde = spat.priors,  in_stack=my.stack, 
                      invariant = "-1+f(iSpat, model=spat.priors)", direction="backwards", 
                      include=6:14, y='y', powerl = 1, inter=1)

autoplot(bestmodel$best_model, which = c(1, 5), CI = TRUE)

####### Making a predictive map - Legal #########

##### Loading raster files ######
setwd(s.dir)

bathy <- raster(paste(s.dir, "bathy.tif", sep='/'))
bathy <- flip(bathy, direction="y")
plot(bathy) # check it
proj4string(bathy) # check the coordinate system, resolution, etc..

tpi <- raster(paste(s.dir, "tpi.tif", sep='/'))
tpi <- flip(tpi,direction='y')
plot(tpi) # check it
proj4string(tpi) 

aspect <- raster(paste(s.dir, "aspect.tif", sep='/'))
plot(aspect)
proj4string(aspect)

slope <- raster(paste(s.dir, "slope.tif", sep='/'))
plot(slope)
proj4string(slope)

flowdir <- raster(paste(s.dir, "flowdir.tif", sep='/'))
plot(flowdir)
proj4string(flowdir)

ramps <- raster(paste(s.dir, "d_ramp_m.tif", sep='/'))
plot(ramps)
proj4string(ramps)

zones <- raster(paste(s.dir,"Rasterised NTZ.tif", sep='/'))
plot(zones)
proj4string(zones)

# Define an empty raster to hold the coordinates of the prediction must be same extent as 
# covariate rasters
r <- raster(xmn = min(754687.5), xmx = max(851312.5),
            ymn = min(7460312), ymx = max(7619938),
            resolution = 25)

# Create data frame with everything in it
preds <- as.data.frame(xyFromCell(r, cell = 1:ncell(r)))
preds$bathymetry <- values(bathy)
preds$TPI <- values(tpi)
preds$Aspect <- values(aspect)
preds$Slope <- values(slope)
preds$FlowDir <- values(flowdir)
preds$distance.to.ramps <- values(ramps)
preds$status <- values(zones)

preds <- preds %>% 
          mutate(status = replace(status, status == '1', "NT"))
preds <- preds %>%
          mutate(status = replace(status, status== '0', "F"))
preds$status <- as.factor(preds$status)


#Remove NAs
preds <- na.omit(preds)
glimpse(preds)

#### Create mesh as usual ####
mesh.pred <- inla.mesh.2d(Locations, max.edge=c(3000))


#### Create SPDE as usual and set priors ####

pred.spat.priors <- inla.spde2.pcmatern(
  mesh=mesh.pred, alpha=2, ### alpha can usually be left at 2 unless using a 3d mesh
  constr=TRUE,### so that random effects are constrained to sum to zero
  prior.range=c(500, 0.05), ### 300 for legal, not sure for sub-legal
  #This is saying that the probability the range is less than 300 is 5%
  prior.sigma=c(0.1, 0.1)) ### 0.5,0.08 for legal, not sure for sub-legal
  #This is saying that probability the variance is greater than 0.5 is less than 8%

#### Create index ####
s.index.p <- inla.spde.make.index(name = "predict", 
                                  n.spde = pred.spat.priors$n.spde) 

#### Need two different A matrices one for prediciton and one with the existing points on it #####

# Normal A matrix
A.matrix <- inla.spde.make.A(mesh.pred, loc = as.matrix(data.legal[,c("easting","northing")]))

# Predictive matrix 
A.pred <- inla.spde.make.A(mesh=mesh.pred)

#### Make the stacks ####

# Stack with estimates 
stackEst <- inla.stack(data=list(y=data.legal$target.fish), A=list(A.matrix, 1),
                       tag = "Est",
                       effects=list(c(mesh.index, list( intercept=1)),#for the spatial effects and intercept
                                    list(bathymetry=data.legal$bathymetry, TPI=data.legal$TPI, 
                                         Slope=data.legal$Slope, Aspect=data.legal$Aspect, 
                                         FlowDir=data.legal$FlowDir, distance.to.ramp=
                                           data.legal$distance.to.ramp, 
                                         status=data.legal$status))) #covariates

# Stack for prediction
stackPred <- inla.stack(data = list(y = NA),  # NAs in the response variable  
                        A = list(A.pred),
                        effects = list(c(s.index.p, list(Intercept = 1))),
                        tag = "Pred")

# Combine these stacks into one
StackJoin <- inla.stack(stackEst, stackPred)

#### Run the model as usual but using the joint stack ####
f.s <- y ~ -1 + intercept + bathymetry + TPI + Slope + Aspect + FlowDir +
  distance.to.ramp + status + f(predict, model=pred.spat.priors) 
#+ f(data.legal$site, model="iid")

fm <- inla(f.s,
           family = "zeroinflatedpoisson1",
           data = inla.stack.data(StackJoin),
           verbose=FALSE,
           control.predictor=list(A=inla.stack.A(StackJoin), compute=TRUE, link=1),
           control.fixed = list(mean=0, prec=0.2),
           control.results = list(return.marginals.random = TRUE, return.marginals.predictor = TRUE), 
           control.compute=list(config = TRUE, dic=TRUE))

summary(fm)
plot(fm)

#### Select relevant posterior mean and SD ####
index.pred <- inla.stack.index(StackJoin, "Pred")$data

post.mean.pred <- fm$summary.linear.predictor[index.pred, "mean"]
post.sd.pred <- fm$summary.linear.predictor[index.pred, "sd"]

post.mean.pred <- fm$summary.fitted.values[index.pred, "mean"]


#### Create a grid to project on ####


Seq.X.grid <- seq(from = min(preds$x),
                  to = max(preds$x))

Seq.Y.grid <- seq(from = min(preds$y),
                   to = max(preds$y))

pred.grid <- as.matrix(expand.grid(x = Seq.X.grid,
                                   y = Seq.Y.grid)) 

str()


#### Project onto grid #####
str(preds)
nrow(preds)

proj.grid <- inla.mesh.projector(mesh.pred,
                                 xlim = range(preds[,1]), 
                                 ylim = range(preds[,2]), 
                                 dims = c(1443537,1443537))

post.mean.pred.grid <- inla.mesh.project(proj.grid, post.mean.pred)
post.sd.pred.grid <- inla.mesh.project(proj.grid, post.sd.pred)

#### Plot mean ####
my.palette.post <- rev(brewer.pal(n = 9, name = "YlGnBu"))

predmean <- t(post.mean.pred.grid)
predmean2 <- predmean[rev(1:length(predmean[,1])),]
predmean_ras <- raster(predmean2,
                       crs = CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 
                                                           +datum=OSGB36 +units=m +no_defs +ellps=airy +towgs84=446.448,
                                                           -125.157,542.060,0.1502,0.2470,0.8421,-20.4894"))

plot(predmean_ras, asp = 1, col = my.palette.post)

#### Plot SD ####
my.palette.var <- brewer.pal(n = 9, name = "BuPu")

predsd <- t(post.sd.pred.grid)
predsd2 <- predsd[rev(1:length(predsd[,1])),]
predsd_ras <- raster(predsd2,
                     crs = CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 
                                                           +datum=OSGB36 +units=m +no_defs +ellps=airy +towgs84=446.448,
                                                           -125.157,542.060,0.1502,0.2470,0.8421,-20.4894"))

par(mfrow = c(1,1), mar = c(2,2, 1,1))
plot(predsd_ras, asp = 1, col = my.palette.var)


###################### Sublegal Model #########################
########## Setting up a mesh #########
# This creates a mesh of discrete sampling locations that allows the model to estimate the spatial
# autocorrelation in the data 

# Rename coordinates
data.sublegal <- data.sublegal%>%
  rename(easting=latitude)%>%
  rename(northing=longitude)

Locations <- data.sublegal%>%
  dplyr::select("easting", "northing")

# Turn into locations if using the modified coastline
coordinates(Locations) <- ~ northing + easting

setwd(s.dir)

interest <- readOGR(dsn = ".", layer = "Area of Interest")
coastline <- readOGR(dsn=".", layer = "Coastline")

plot(coastline)
plot(interest)

interest <- spTransform(interest, CRS("+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
coastline <- spTransform(coastline, CRS("+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

pll <- Polygon(interest@polygons[[1]]@Polygons[[1]]@coords, hole = F)
hole <- Polygon(coastline@polygons[[1]]@Polygons[[1]]@coords, hole = T)
new_sp <- SpatialPolygons(list(Polygons(append(hole, pll),'1')))

plot(new_sp)


# Creating the mesh
my.mesh.land <- inla.mesh.2d(Locations, boundary=new_sp, max.edge = c(6000)) #This is the modified one

my.mesh <- inla.mesh.2d(Locations, max.edge=c(3000))

plot(my.mesh)
points(Locations, col = 2, pch = 16, cex = 0.2)

######### Setting the priors ########

spat.priors.sublegal <- inla.spde2.pcmatern(
  mesh=my.mesh, alpha=2, ### alpha can usually be left at 2 unless using a 3d mesh
  constr=TRUE,### so that random effects are constrained to sum to zero
  prior.range=c(900, 0.05), ### 300 for legal, not sure for sub-legal
  #This is saying that the probability the range is less than 300 is 5%
  prior.sigma=c(0.1, 0.1)) ### 0.5,0.08 for legal, not sure for sub-legal
#This is saying that probability the variance is greater than 0.5 is less than 8%

######### Making the A matrix ########
# This translates the spatial locations on the mesh to vectors in the model
A.matrix.sublegal <- inla.spde.make.A(my.mesh, loc = as.matrix(data.sublegal[,c("easting","northing")]))

######## Creating a stack #######
# This is a combination of the A. matrix and the mesh 

mesh.index.sublegal <- inla.spde.make.index(name = "iSpat", n.spde = spat.priors.sublegal$n.spde)
my.stack.sublegal <- inla.stack(data=list(y=data.sublegal$target.fish), A=list(A.matrix.sublegal, 1),
                       tag = "preds",
                       effects=list(c(mesh.index, list( intercept=1)),#for the spatial effects and intercept
                                    list(bathymetry=data.sublegal$bathymetry, TPI=data.sublegal$TPI, 
                                         Slope=data.sublegal$Slope, Aspect=data.sublegal$Aspect, 
                                         FlowDir=data.sublegal$FlowDir,distance.to.ramp=
                                           data.sublegal$distance.to.ramp, 
                                         status=data.sublegal$status))) #covariates

###### Run Full Model #########
f.s <- y ~ -1 + intercept + bathymetry + TPI + Slope + Aspect + FlowDir +
  distance.to.ramp + status + f(iSpat, model=spat.priors.sublegal) 
#+ f(data.legal$site, model="iid")

fm <- inla(f.s,
           family = "zeroinflatedpoisson1",
           data = inla.stack.data(my.stack.sublegal),
           verbose=FALSE,
           control.predictor=list(A=inla.stack.A(my.stack.sublegal), compute=TRUE, link=1),
           control.fixed = list(mean=0, prec=0.2), #This is where you put priors for fixed effects
           control.results = list(return.marginals.random = TRUE, return.marginals.predictor = TRUE), 
           control.compute=list(config = TRUE, dic=TRUE)
)

summary(fm)

plot(fm)


####### Plotting the residuals and checking the model ######
# If the model is well calibrated then the bins should be the same height. The convex model suggests that
# the model is overconfident in it's predictions and the line graph suggests we've got both over and 
# under fitting 
ggplot_inla_residuals(fm, data.sublegal$target.fish, binwidth = 0.1)

# This plots the values predeicted by the model against the actual values 
index.pred <- inla.stack.index(my.stack, "preds")$data

post.mean.pred <- fm$summary.fitted.values[index.pred, "mean"]
post.sd.pred <- fm$summary.fitted.values[index.pred, "sd"]

plot(post.mean.pred,data.sublegal$target.fish, ylim=c(0,13), xlim=c(0,12))


######## Plotting the results ###########
#There are no p-values with INLA you look at significance by finding the overlap of the 0.025 and the 0.975
# Quantiles with 0 - those that don't overlap may be significant 
summary <- as.data.frame(summary(fm)$fixed)

summary <- summary%>%
  rename(lower=`0.025quant`)%>%
  rename(upper=`0.975quant`)%>%
  rename(estimate=mean)%>%
  mutate(factor=c("intercept", "bathy", "TPI", "slope", "aspect", "flowdir",
                  "ramp", "fished", "NTZ"))

effects <- ggplot(summary,
                  aes(x = estimate,
                      y = factor, color=factor)) +
  geom_point() +
  geom_errorbar(aes(xmin=lower, xmax=upper), width=.5)+
  geom_vline(xintercept = 0)+
  theme_bw()
effects

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

####### Stepwise Model selection ########
bestmodel <- INLAstep(fam1="zeroinflatedpoisson1", data.legal, spde = spat.priors,  in_stack=my.stack, 
                      invariant = "-1+f(iSpat, model=spat.priors)", direction="backwards", 
                      include=6:14, y='y', powerl = 1, inter=1)

autoplot(bestmodel$best_model, which = c(1, 5), CI = TRUE)








