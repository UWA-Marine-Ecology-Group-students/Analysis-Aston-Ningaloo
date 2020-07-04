library(raster)
library(INLA)
library(dplyr)
library(sp) 
library(fields)
library(ggplot2)
library(INLAutils)


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
  select(!sand)%>%
  select(!TRI)%>%
  select(!Roughness)

# Set your covariates/spatial data
covariates <- c("bathymetry","TPI","Slope","Aspect","FlowDir","mean.relief",
                "sd.relief","reef","distance.to.ramp")
factors <- c("status")
covariates <- na.exclude(covariates)


########## Setting up a mesh #########
# This creates a mesh of discrete sampling locations that allows the model to estimate the spatial
# autocorrelation in the data 

data <- data%>%
  rename(easting=latitude)%>%
  rename(northing=longitude)

Locations <- data%>%
  dplyr::select("easting", "northing")

my.mesh <- inla.mesh.2d(Locations, max.edge=c(6000,11000), offset = c(7000, 7000))

plot(my.mesh)

######### Setting the priors ########

spat.priors <- inla.spde2.pcmatern(
  mesh=my.mesh, alpha=2, ### alpha can usually be left at 2 unless using a 3d mesh
  constr=TRUE,### so that random effects are constrained to sum to zero
  prior.range=c(5000, 0.05), ### P(practic.range<150)=0.05
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
                                           status=data$status)))#covariates

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
           "reef", "ramp", "statusF", "statusNT"))

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

####### Model selection ########
model <- INLAstep(fam1="zeroinflatedpoisson1", data, in_stack=my.stack, 
                  invariant = "0-1", direction="backwards", 
                  include=6:15, y='y',
                  inter=3)




