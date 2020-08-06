# a newdata to get the predictions
newdata <- expand.grid(bathymetry = seq(min(data.legal$bathymetry), 
                                      max(data.legal$bathymetry), 
                                      length.out = 10),
                       TPI = seq(min(data.legal$TPI),
                                 max(data.legal$TPI),
                                 length.out = 10),
                       Aspect = seq(min(data.legal$Aspect),
                                    max(data.legal$Aspect),
                                    length.out=10),
                       Slope = seq(min(data.legal$Slope),
                                   max(data.legal$Slope),
                                   length.out=10),
                       FlowDir = seq(min(data.legal$FlowDir),
                                     max(data.legal$FlowDir),
                                     length.out=10),
                       distance.to.ramp = seq(min(data.legal$distance.to.ramp),
                                           max(data.legal$distance.to.ramp),
                                           length.out=10),
                      status = unique(data.legal$status))

# the stack for these predictions
pred_stack_fixef <- inla.stack(data = list(y = NA),
                               A = list(1, 1, 1, 1, 1, 1, 1,1),
                               effects = list(Intercept = rep(1, nrow(newdata)),
                                              bathymetry=newdata$bathymetry, TPI=newdata$TPI, 
                                              Slope=newdata$Slope, Aspect=newdata$Aspect, 
                                              FlowDir=newdata$FlowDir,distance.to.ramp=
                                                newdata$distance.to.ramp, 
                                              status=newdata$status),
                               tag = "prd_fixef")
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

# create a new projection matrix for the points
Apred <- inla.spde.make.A(mesh.pred,
                          loc = as.matrix(preds[,c("x", "y")]))

# put this in a new stack
pred_stack_alleff <- inla.stack(data = list(calcium = NA),
                                A = list(Apred, 1, 1, 1),
                                effects = list(i = 1:spde$n.spde,
                                               Intercept = rep(1, nrow(newdat)),
                                               elevation = newdat$elevation,
                                               region = factor(newdat$region)),
                                tag = "prd_alleff")

#######

