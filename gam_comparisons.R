## The first part of the code just formats the data you need for the gams correctly 
## I've also written code for the top two models for each of the spatial methods that FSSgam produced
## For the corGaus, as I can't figure out how to make work in FSSgam I've just put the two models with the 
## "lowest AIC" from the different options that I tried 


detach("package:plyr", unload=TRUE)
library(dplyr)
library(mgcv)
library(FSSgam)
library(gamm4)
library(doBy)
library(doParallel)

#### Reading in and formatting the data ####
legal.dat <- read.csv("final.data.legal.csv")%>%
  rename(response=target.fish)

# Remove NA values and sites in state no-take zones
legal.dat<-legal.dat%>%
  filter(!sample%in%c("8.05","10.09","10.12","16.03", "10.14", "10.15", "10.13", "14.13", "14.12", "14.02"))

# Transform variables - these are the transformations we worked out would be best for the gam initially
legal.dat <- legal.dat%>%
  mutate(sqrt.reef=sqrt(reef))%>%
  mutate(sqrt.slope=sqrt(Slope))%>%
  mutate(log.roughness=log(Roughness+1))%>%
  mutate(cube.aspect=(Aspect)^3)%>%
  glimpse()

#Set predictor variables 
pred.vars=c("sqrt.slope","cube.aspect","log.roughness","FlowDir",
            "distance.to.ramp")


#### Top models with no accounting for spatial structure ####
model.no.spatial.1 <- gam(response~
                          s(log.roughness, k=3, bs='cr') + status + s(bathymetry, bs="cr"),
                          family=tw(), data=legal.dat)
summary(model.no.spatial.1)

model.no.spatial.2 <- gam(response~
                          s(sqrt.slope, k=3, bs='cr') + status + s(bathymetry, bs="cr"),
                          family=tw(), data=legal.dat)
summary(model.no.spatial.2)


#### Top models with site as a random factor ####
model.site.1 <- gam(response~
                    status + s(bathymetry, bs="cr") 
                    + s(site, bs="re"),
                    family=tw(), data=legal.dat)
summary(model.site.1)

model.site.2 <- gam(response~
                    status + s(cube.aspect, k=3, bs='cr') + s(bathymetry, bs="cr") 
                    + s(site, bs="re"),
                    family=tw(), data=legal.dat)
summary(model.site.2)

#### Top models with spline k=3 ####
model.splinek3.1 <- gam(response~
                        s(distance.to.ramp,k=3,bs='cr') + s(bathymetry, bs="cr") 
                        + te(latitude, longitude, k=3, bs="cr"),
                        family=tw(),  data=legal.dat)
summary(model.splinek3.1)

model.splinek3.2 <- gam(response~
                        s(distance.to.ramp,k=3,bs='cr') + s(sqrt.slope, k=3, bs="cr") + s(bathymetry, bs="cr") 
                        + te(latitude, longitude, k=3, bs="cr"),
                        family=tw(),  data=legal.dat)
summary(model.splinek3.2)

#### Top models with spline k=5/unrestricted spline ####
model.splinek5.1 <- gam(response~
                        s(distance.to.ramp,k=3,bs='cr') + s(cube.aspect, k=3, br="cr") + s(bathymetry, bs="cr") 
                        + te(latitude, longitude, k=5, bs="cr"),
                        family=tw(),  data=legal.dat)
summary(model.splinek5.1)

model.splinek5.2 <- gam(response~
                        s(distance.to.ramp,k=3,bs='cr') + s(FlowDir, k=3, br="cr") + s(bathymetry, bs="cr") 
                        + te(latitude, longitude, k=5, bs="cr"),
                        family=tw(),  data=legal.dat)
summary(model.splinek5.2)

#### Top models using corGaus ####
model.corgaus.1 <- gamm(response~
                          s(sqrt.slope, k=3, bs='cr') + status + s(bathymetry, bs="cr"),
                        family=poisson(), correlation = corGaus(form = ~ latitude + longitude),  data=legal.dat)

AIC(model.corgaus.1)
summary(model.corgaus.1$gam)

model.corgaus.2 <- gamm(response~
                        s(log.roughness, k=3, bs='cr') + status + s(bathymetry, bs="cr"),
                        family=poisson(), correlation = corGaus(form = ~ latitude + longitude),  data=legal.dat)

AIC(model.corgaus.2)
summary(model.corgaus.2$gam)


