# Set directories----
rm(list=ls())

# Study name ----
study <- "ningaloo" 

# Libraries required
install_github("UWAMEGFisheries/GlobalArchive") #to check for updates
library(GlobalArchive)

library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(googlesheets4)

library(ggmap)
library(rgdal)
library(raster)
library(png)

library("ggspatial")
library("rnaturalearth")
library("rnaturalearthdata")
library(cowplot)
library(scatterpie)

## Set your working directory ----
working.dir<-dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

setwd(working.dir)
dir()

# Set sub directories----
plots.dir=paste(working.dir,"plots",sep="/")
tidy.dir=paste(working.dir,"Tidy data",sep="/")
images.dir=paste(working.dir,"images",sep="/")
spatial.dir=paste(working.dir,"Spatial",sep="/")

# functions for summarising data on plots----
se <- function(x) sd(x) / sqrt(length(x))
se.min <- function(x) (mean(x)) - se(x)
se.max <- function(x) (mean(x)) + se(x)

theme.larger.text<-theme(
  strip.text.x = element_text(size = 5,angle = 0),
  strip.text.y = element_text(size = 5),
  axis.title.x=element_text(vjust=-0.0, size=10),
  axis.title.y=element_text(vjust=0.0,size=10),
  axis.text.x=element_text(size=6),
  axis.text.y=element_text(size=6),
  legend.title = element_text(family="TN",size=8),
  legend.text = element_text(family="TN",size=8))

theme.species<-theme(
  strip.text.x = element_text(size = 8,angle = 0),
  strip.text.y = element_text(size = 8),
  axis.title.x=element_text(vjust=-0.0, size=12),
  axis.title.y=element_text(vjust=0.0,size=12),
  axis.text.x=element_text(size=8),
  axis.text.y=element_text(size=8),
  legend.title = element_text(family="TN",size=11),
  legend.text = element_text(family="TN",size=11))


theme_collapse<-theme(      ## the commented values are from theme_grey
  panel.grid.major=element_line(colour = "white"), ## element_line(colour = "white")
  panel.grid.minor=element_line(colour = "white", size = 0.25),
  strip.text.x = element_text(size = 4,angle = 0),
  strip.text.y = element_text(size = 4),
  axis.title.x=element_text(vjust=-0.0, size=10),
  axis.title.y=element_text(vjust=0.0,size=10),
  axis.text.x=element_text(size=6),
  axis.text.y=element_text(size=6),
  legend.title = element_text(family="TN",size=6),
  legend.text = element_text(family="TN",size=4))

# Load fish pictures for plotting ----
setwd(images.dir)
dir()

n.b <- readPNG("Nemipterus_bathybius_nb_GIBBONS.png")
n.b <- as.raster(n.b)

n.v <- readPNG("Nemipteridae.png")
n.v <- as.raster(n.v)

d.spp <- readPNG("Decapterus_spp.png")
d.spp <- as.raster(d.spp)

s.u <- readPNG("Synodus_variegatus_nb.png")
s.u <- as.raster(s.u)

c.e <- readPNG("Carangoides_equula_nb_GIBBONS.png")
c.e <- as.raster(c.e)

d.c <- readPNG("Dentex_carpenteri_nb_GIBBONS.png")
d.c <- as.raster(d.c)

s.l <- readPNG("Sphyrna_lewini_nb_GIBBONS.png")
s.l <- as.raster(s.l)

# Read in shapefile ----
setwd(spatial.dir)
dir()

commonwealth.marine.parks <- readOGR(spatial.dir, "AustraliaNetworkMarineParks")
commonwealth.marine.parks<- fortify(commonwealth.marine.parks)

state.marine.parks <- readOGR(spatial.dir, "WA_MPA_2018")
state.marine.parks <- fortify(state.marine.parks)

# read in maxn
setwd(tidy.dir)
dir()

maxn <- read.csv("ningaloo.complete.maxn.csv")%>%
  mutate(sample=as.character(sample))

metadata <- read.csv("ningaloo.checked.metadata.csv")%>%
  mutate(sample=as.character(sample))

habitat<-read_csv("2020-01_Guardian-Ningaloo_stereoBRUVs._habitat.csv" )%>% # change this
  ga.clean.names()

# workout total abundance and species richness
maxn.ta.sr <- maxn%>%
  group_by(scientific,sample)%>%
  dplyr::summarise(maxn = sum(maxn))%>%
  spread(scientific,maxn, fill = 0)%>%
  mutate(total.abundance=rowSums(.[,2:(ncol(.))],na.rm = TRUE ))%>% #Add in Totals
  mutate(species.richness=rowSums(.[,2:17] > 0))%>%
  dplyr::select(sample,total.abundance,species.richness)%>%
  gather(.,"scientific","maxn",2:3)%>%
  left_join(metadata)

# Format maxn for species specific plots
unique(maxn$scientific)



# TOTAL ABUNDANCE ----
spatial.ta<-ggplot() +
  geom_polygon(data = commonwealth.marine.parks, aes(x = long, y = lat, group = group),color = 'black', fill = 'grey90', size = .2)+ # change colours and add wa state reserves
  geom_point(data=filter(maxn.ta.sr,scientific%in%c("total.abundance")&maxn==0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="white",alpha=0.75)+
  geom_point(data=filter(maxn.ta.sr,scientific%in%c("total.abundance")&maxn>0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="dodgerblue2",alpha=0.75)+
  xlab('Longitude')+
  ylab('Latitude')+
  labs(size = "Total\nabundance")+
  # annotate("text",x=193000, y=7606500,label="Total abundance",color="Black",hjust=0,family="TN",cex=3.5,fontface="italic")+ # change this to a different lat and lon x=lon, y=lat
  theme_bw()+
  theme_collapse+
  theme.larger.text

spatial.ta

# SPECIES RICHNESS ----
spatial.sr<-ggplot() +
  geom_polygon(data = shapefile_df, aes(x = long, y = lat, group = group),color = 'black', fill = 'grey90', size = .2)+
  geom_point(data=filter(maxn.ta.sr,scientific%in%c("species.richness")&maxn==0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="white",alpha=0.75)+
  geom_point(data=filter(maxn.ta.sr,scientific%in%c("species.richness")&maxn>0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="dodgerblue2",alpha=0.75)+
  xlab('Longitude')+
  ylab('Latitude')+
  labs(size = "Species\nrichness")+
  annotate("text",x=193000, y=7606500,label="Species richness",color="Black",hjust=0,family="TN",cex=3.5,fontface="italic")+
  
  theme_bw()+
  theme_collapse+
  theme.larger.text

spatial.sr

# Nemipterus bathybius ----
species <- c("Nemipterus bathybius")

spatial.bathybius<-ggplot() +
  geom_polygon(data = shapefile_df, aes(x = long, y = lat, group = group),color = 'black', fill = 'grey90', size = .2)+
  geom_point(data=filter(top.5.maxn,genus.species%in%c("Nemipterus bathybius")&maxn==0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="white",alpha=0.75)+
  geom_point(data=filter(top.5.maxn,genus.species%in%c("Nemipterus bathybius")&maxn>0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="dodgerblue2",alpha=0.75)+
  xlab('Longitude')+
  ylab('Latitude')+
  labs(size = "Relative \nabundance")+
  #labs(size = " ")+
  annotate("text",x=193000, y=7606500,label=species,color="Black",hjust=0,family="TN",cex=3.5,fontface="italic")+
  annotation_raster(n.b, xmin=196800, xmax=198250, ymin=7605800, ymax=7606600)+
  theme_bw()+
  theme_collapse+
  theme.species

spatial.bathybius


# Carangoides equula ----
species <- c("Carangoides equula")

spatial.equula<-ggplot() +
  geom_polygon(data = shapefile_df, aes(x = long, y = lat, group = group),color = 'black', fill = 'grey90', size = .2)+
  geom_point(data=filter(top.5.maxn,genus.species%in%c("Carangoides equula")&maxn==0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="white",alpha=0.75)+
  geom_point(data=filter(top.5.maxn,genus.species%in%c("Carangoides equula")&maxn>0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="dodgerblue2",alpha=0.75)+
  xlab('Longitude')+
  ylab('Latitude')+
  labs(size = "Relative \nabundance")+
  #labs(size = " ")+
  annotate("text",x=193000, y=7606500,label=species,color="Black",hjust=0,family="TN",cex=3.5,fontface="italic")+
  annotation_raster(c.e, xmin=196800, xmax=198250, ymin=7605800, ymax=7606600)+
  theme_bw()+
  theme_collapse+
  theme.species


spatial.equula

# Dentex carpenteri ----
species <- c("Dentex carpenteri")

spatial.carpenteri<-ggplot() +
  geom_polygon(data = shapefile_df, aes(x = long, y = lat, group = group),color = 'black', fill = 'grey90', size = .2)+
  geom_point(data=filter(top.5.maxn,genus.species%in%c("Dentex carpenteri")&maxn==0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="white",alpha=0.75)+
  geom_point(data=filter(top.5.maxn,genus.species%in%c("Dentex carpenteri")&maxn>0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="dodgerblue2",alpha=0.75)+
  xlab('Longitude')+
  ylab('Latitude')+
  labs(size = "Relative \nabundance")+
  #labs(size = " ")+
  annotate("text",x=193000, y=7606500,label=species,color="Black",hjust=0,family="TN",cex=3.5,fontface="italic")+
  annotation_raster(d.c, xmin=196800, xmax=198250, ymin=7605800, ymax=7606750)+
  theme_bw()+
  theme_collapse+
  theme.species


spatial.carpenteri

# Synodus variegatus ----
species <- c("Synodus variegatus")

spatial.variegatus<-ggplot() +
  geom_polygon(data = shapefile_df, aes(x = long, y = lat, group = group),color = 'black', fill = 'grey90', size = .2)+
  geom_point(data=filter(top.5.maxn,genus.species%in%c("Synodus variegatus")&maxn==0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="white",alpha=0.75)+
  geom_point(data=filter(top.5.maxn,genus.species%in%c("Synodus variegatus")&maxn>0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="dodgerblue2",alpha=0.75)+
  xlab('Longitude')+
  ylab('Latitude')+
  labs(size = "Relative \nabundance")+
  #labs(size = " ")+
  annotate("text",x=193000, y=7606500,label=species,color="Black",hjust=0,family="TN",cex=3.5,fontface="italic")+
  annotation_raster(s.u,xmin=196500, xmax=198250, ymin=7606000, ymax=7606800)+
  theme_bw()+
  theme_collapse+
  theme.species

spatial.variegatus

# Decapterus tabl ----
species <- c("Decapterus tabl")

spatial.tabl<-ggplot() +
  geom_polygon(data = shapefile_df, aes(x = long, y = lat, group = group),color = 'black', fill = 'grey90', size = .2)+
  geom_point(data=filter(top.5.maxn,genus.species%in%c("Decapterus tabl")&maxn==0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="white",alpha=0.75)+
  geom_point(data=filter(top.5.maxn,genus.species%in%c("Decapterus tabl")&maxn>0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="dodgerblue2",alpha=0.75)+
  xlab('Longitude')+
  ylab('Latitude')+
  labs(size = "Relative \nabundance")+
  annotate("text",x=193000, y=7606500,label=species,color="Black",hjust=0,family="TN",cex=3.5,fontface="italic")+
  annotation_raster(d.spp, xmin=196250, xmax=198300, ymin=7605500, ymax=7606700)+
  theme_bw()+
  theme_collapse+
  theme.species

spatial.tabl


# Nemipterus virgatus ----
species <- c("Nemipterus virgatus")

spatial.virgatus<-ggplot() +
  geom_polygon(data = shapefile_df, aes(x = long, y = lat, group = group),color = 'black', fill = 'grey90', size = .2)+
  geom_point(data=filter(top.5.maxn,genus.species%in%c("Nemipterus virgatus")&maxn==0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="white",alpha=0.75)+
  geom_point(data=filter(top.5.maxn,genus.species%in%c("Nemipterus virgatus")&maxn>0),aes(longitude,latitude,size=maxn),shape=21,colour="dodgerblue4",fill="dodgerblue2",alpha=0.75)+
  xlab('Longitude')+
  ylab('Latitude')+
  labs(size = "Relative \nabundance")+
  annotate("text",x=193000, y=7606500,label=species,color="Black",hjust=0,family="TN",cex=3.5,fontface="italic")+
  annotation_raster(n.v, xmin=196800, xmax=198250, ymin=7605800, ymax=7606600)+
  theme_bw()+
  theme_collapse+
  theme.species

spatial.virgatus



# SAVE - Species richness and total abundance combined ----
setwd(plots.dir)
ta.sr<-plot_grid(spatial.ta, spatial.sr, labels = c('A', 'B'), label_size = 12,ncol=1)
ggsave("total.abundance.and.species.richness.potrait.png",ta.sr,dpi=300,width=11,height=17.5,unit="cm")
# 
# ta.sr<-plot_grid(spatial.ta, spatial.sr, labels = c('A', 'B'), label_size = 12,ncol=2)
# ggsave("total.abundance.and.species.richness.landscape1.png",ta.sr,dpi=300,width = 20, height = 5.5,unit="cm") 

# Species specific ----
species.combined<-plot_grid(spatial.equula, spatial.tabl, 
                            spatial.carpenteri, spatial.bathybius, 
                            spatial.virgatus,spatial.variegatus, 
                            labels = c('A', 'B','C','D','E','F'), label_size = 12,ncol=2)
ggsave("spatial.species.png",species.combined,dpi=500, width = 21, height = 23,units = "cm")


species.combined<-plot_grid(spatial.bathybius, spatial.carpenteri, 
                            spatial.equula, spatial.tabl, 
                            spatial.variegatus, spatial.virgatus,
                            labels = c('A', 'B','C','D','E','F'), label_size = 12,ncol=2)
ggsave("spatial.species.png",species.combined,dpi=500, width = 21, height = 23,units = "cm")




# Species specific
species.combined.land<-plot_grid(spatial.equula, spatial.tabl, 
                            spatial.carpenteri, spatial.bathybius, 
                            spatial.virgatus,spatial.variegatus, 
                            labels = c('A', 'B','C','D','E','F'), label_size = 12,ncol=3)
ggsave("spatial.species.landscape.png",species.combined.land,dpi=100, width = 33, height = 16,units = "cm")

# species.combined1<-plot_grid(spatial.equula, spatial.tabl, 
#                             labels = c('A', 'B'), label_size = 12,ncol=2)
# ggsave("spatial.species1.png",species.combined1,dpi=300,width=10,height=3.75)
# 
# species.combined2<-plot_grid(spatial.carpenteri, spatial.bathybius, 
#                              labels = c('C', 'D'), label_size = 12,ncol=2)
# ggsave("spatial.species2.png",species.combined2,dpi=300,width=10,height=3.75)
# 
# species.combined3<-plot_grid(spatial.undosquamis, 
#                              labels = c('E'), label_size = 12,ncol=2)
# ggsave("spatial.species3.png",species.combined3,dpi=300,width=15,height=5.75)

  
# Habitat bubble plots  ----
glimpse(habitat)

hab<-habitat%>%
  dplyr::select(sample,bedform.bioturbated,bedform.none)%>%
  dplyr::rename(bioturbated=bedform.bioturbated,none=bedform.none)%>%
  left_join(metadata)%>%
  glimpse()

bedforms <- ggplot() +
  geom_polygon(data = shapefile_df, aes(x = long, y = lat, group = group),color = 'black', fill = 'grey90', size = .1)+
  geom_scatterpie(aes(x=longitude, y=latitude,r=150),
                    data=hab, cols=c("bioturbated","none"), color="black", alpha=.8,legend_name="Bedform")+
  xlab('Longitude')+
  ylab('Latitude')+
  labs(size = "Relative abundance")+
  theme_bw()+
  theme_collapse+
  theme.larger.text
bedforms

hab.sp<-habitat%>%
  left_join(metadata)

spatial.sand<-ggplot() +
  geom_polygon(data = shapefile_df, aes(x = long, y = lat, group = group),color = 'black', fill = 'grey90', size = .2)+
  geom_point(data=filter(hab.sp,broad.unconsolidated==0),aes(longitude,latitude,size=broad.unconsolidated),shape=21,colour="dodgerblue4",fill="white",alpha=0.75)+
  geom_point(data=filter(hab.sp,broad.unconsolidated>0),aes(longitude,latitude,size=broad.unconsolidated),shape=21,colour="dodgerblue4",fill="dodgerblue2",alpha=0.75)+
  xlab('Longitude')+
  ylab('Latitude')+
  labs(size = "Percent cover")+
  annotate("text",x=193000, y=7606500,label="Unconsolidated (fine)",color="Black",hjust=0,family="TN",cex=3.5,fontface="italic")+
  theme_bw()+
  theme_collapse+
  theme.larger.text

spatial.sand

hab.combined<-plot_grid(spatial.sand, bedforms,
                             labels = c('A','B'), label_size = 12,ncol=2)


hab.combined2<-plot_grid(spatial.sand, bedforms,
                        labels = c('A','B'), label_size = 12,ncol=1)

setwd(plots.dir)
ggsave("hab.combined.png",hab.combined,dpi=300,width=15,height=5.75)
ggsave("hab.combined.potrait.png",hab.combined2,dpi=300,width=12,height=17.5,unit="cm")



# Basic Map 
spatial.deployments<-ggplot() +
  geom_polygon(data = shapefile_df, aes(x = long, y = lat, group = group),color = 'black', fill = 'grey90', size = .2)+
  geom_point(data=metadata,aes(longitude,latitude),shape=21,colour="black",fill="white",size=4)+
  xlab('Longitude')+
  ylab('Latitude')+
  theme_bw()+
  theme_collapse+
  theme.species

spatial.deployments
setwd(plots.dir)
ggsave("deployment.map.png",spatial.deployments,dpi=300)

