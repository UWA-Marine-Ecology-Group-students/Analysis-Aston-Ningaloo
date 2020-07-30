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

str(data)
data$status <- as.factor(data$status)

#remove NA sites 
data<- data%>%
  filter(!sample%in%c("8.05","10.09","10.12","16.03"))%>%
  dplyr::select(!sand)%>%
  dplyr::select(!TRI)%>%
  dplyr::select(!Roughness)%>%
  dplyr::select(!X)

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
  rename(easting=latitude)%>%
  rename(northing=longitude)

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

spat.priors <- inla.spde2.pcmatern(
  mesh=my.mesh, alpha=2, ### alpha can usually be left at 2 unless using a 3d mesh
  constr=TRUE,### so that random effects are constrained to sum to zero
  prior.range=c(300, 0.1), ### 300 for legal, not sure for sub-legal
  #This is saying that the probability the range is less than 300 is 5%
  prior.sigma=c(0.1, 0.1)) ### 0.5,0.08 for legal, not sure for sub-legal
  #This is saying that probability the variance is greater than 0.5 is less than 8%

######### Making the A matrix ########
# This translates the spatial locations on the mesh to vectors in the model
A.matrix <- inla.spde.make.A(my.mesh, loc = as.matrix(data.legal[,c("easting","northing")]))

######## Creating a stack #######
# This is a combination of the A. matrix and the mesh 

mesh.index <- inla.spde.make.index(name = "iSpat", n.spde = spat.priors$n.spde)
my.stack <- inla.stack(data=list(y=data.legal$target.fish), A=list(A.matrix, 1),
                       tag = "preds",
                       effects=list(c(mesh.index, list( intercept=1)),#for the spatial effects and intercept
                                    list(bathymetry=data.legal$bathymetry, TPI=data.legal$TPI, 
                                         Slope=data.legal$Slope, Aspect=data.legal$Aspect, 
                                         FlowDir=data.legal$FlowDir,distance.to.ramp=
                                           data.legal$distance.to.ramp, 
                                         status=data.legal$status))) #covariates

###### Run Full Model #########
f.s <- y ~ -1 + intercept + bathymetry + TPI + Slope + Aspect + FlowDir +
  distance.to.ramp + status + f(iSpat, model=spat.priors) 
#+ f(data.legal$site, model="iid")

fm <- inla(f.s,
           family = "zeroinflatedpoisson1",
           data = inla.stack.data(my.stack),
           verbose=FALSE,
           control.predictor=list(A=inla.stack.A(my.stack), compute=TRUE, link=1),
           control.fixed = list(mean=0, prec=0.2),
           control.results = list(return.marginals.random = TRUE, return.marginals.predictor = TRUE), 
           control.compute=list(config = TRUE, dic=TRUE)
)

summary(fm)

plot(fm)


####### Plotting the residuals and checking the model ######
# If the model is well calibrated then the bins should be the same height. The convex model suggests that
# the model is overconfident in it's predictions and the line graph suggests we've got both over and 
# under fitting 
ggplot_inla_residuals(fm, data.legal$target.fish, binwidth = 0.1)

# This plots the values predeicted by the model against the actual values 
index.pred <- inla.stack.index(my.stack, "preds")$data

post.mean.pred <- fm$summary.fitted.values[index.pred, "mean"]
post.sd.pred <- fm$summary.fitted.values[index.pred, "sd"]

plot(post.mean.pred,data.legal$target.fish, ylim=c(0,13), xlim=c(0,12))


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
proj4string(tpi) # check the coordinate system, resolution, etc..

aspect <- raster(paste(s.dir, "aspect.tif", spe='/'))
plot(aspect)





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

spat.priors <- inla.spde2.pcmatern(
  mesh=my.mesh, alpha=2, ### alpha can usually be left at 2 unless using a 3d mesh
  constr=TRUE,### so that random effects are constrained to sum to zero
  prior.range=c(900, 0.05), ### 300 for legal, not sure for sub-legal
  #This is saying that the probability the range is less than 300 is 5%
  prior.sigma=c(0.1, 0.1)) ### 0.5,0.08 for legal, not sure for sub-legal
#This is saying that probability the variance is greater than 0.5 is less than 8%

######### Making the A matrix ########
# This translates the spatial locations on the mesh to vectors in the model
A.matrix <- inla.spde.make.A(my.mesh, loc = as.matrix(data.sublegal[,c("easting","northing")]))

######## Creating a stack #######
# This is a combination of the A. matrix and the mesh 

mesh.index <- inla.spde.make.index(name = "iSpat", n.spde = spat.priors$n.spde)
my.stack <- inla.stack(data=list(y=data.sublegal$target.fish), A=list(A.matrix, 1),
                       tag = "preds",
                       effects=list(c(mesh.index, list( intercept=1)),#for the spatial effects and intercept
                                    list(bathymetry=data.sublegal$bathymetry, TPI=data.sublegal$TPI, 
                                         Slope=data.sublegal$Slope, Aspect=data.sublegal$Aspect, 
                                         FlowDir=data.sublegal$FlowDir,distance.to.ramp=
                                           data.sublegal$distance.to.ramp, 
                                         status=data.sublegal$status))) #covariates

###### Run Full Model #########
f.s <- y ~ -1 + intercept + bathymetry + TPI + Slope + Aspect + FlowDir +
  distance.to.ramp + status + f(iSpat, model=spat.priors) 
#+ f(data.legal$site, model="iid")

fm <- inla(f.s,
           family = "zeroinflatedpoisson1",
           data = inla.stack.data(my.stack),
           verbose=FALSE,
           control.predictor=list(A=inla.stack.A(my.stack), compute=TRUE, link=1),
           control.fixed = list(mean=0, prec=0.2),
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








