# Set directories----
rm(list=ls())

# Study name ----
study <- "2020-01_Guardian-Ningaloo_stereoBRUVs" 

# Libraries required
install_github("UWAMEGFisheries/GlobalArchive") #to check for updates
library(GlobalArchive)

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

library(ggmap)
library(rgdal)
library(raster)
library(png)
library(cowplot)

## Set your working directory ----
working.dir<-dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

# Set sub directories----
plots.dir=paste(working.dir,"plots",sep="/")
tidy.dir=paste(working.dir,"tidy data",sep="/")
images.dir=paste(working.dir,"images",sep="/")


theme_collapse<-theme(      ## the commented values are from theme_grey
  panel.grid.major=element_line(colour = "white"), ## element_line(colour = "white")
  panel.grid.minor=element_line(colour = "white", size = 0.25), 
  plot.margin= grid::unit(c(0, 0, 0, 0), "in"))

theme.larger.text<-theme(
  strip.text.x = element_text(size = 5,angle = 0),
  strip.text.y = element_text(size = 5),
  axis.title.x=element_text(vjust=-0.0, size=10),
  axis.title.y=element_text(vjust=0.0,size=10),
  axis.text.x=element_text(size=8),
  axis.text.y=element_text(size=8),
  legend.title = element_text(family="TN",size=8),
  legend.text = element_text(family="TN",size=8))

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

c.p <- readPNG("Carcharinus plumbeus 5cmL_nb.png")
c.p <- as.raster(c.p)

l.l <- readPNG("Lagocephalus_lunaris_nb_BORNT.png")
l.l <- as.raster(l.l)

c.e <- readPNG("Carangoides_equula_nb_GIBBONS.png")
c.e <- as.raster(c.e)

d.c <- readPNG("Dentex_carpenteri_nb_GIBBONS.png")
d.c <- as.raster(d.c)

s.l <- readPNG("Sphyrna_lewini_nb_GIBBONS.png")
s.l <- as.raster(s.l)

# read in maxn
setwd(tidy.dir)
dir()

maxn <- read.csv("2020-01_Guardian-Ningaloo_stereoBRUVs.complete.maxn.csv")
metadata <- read.csv("2020-01_Guardian-Ningaloo_stereoBRUVs.checked.metadata.csv")

# workout total maxn for each species ---
maxn.sum<-maxn%>%
  mutate(scientific=paste(genus,species,sep=" "))%>%
  group_by(scientific)%>%
  dplyr::summarise(maxn=sum(maxn))%>%
  ungroup()%>%
  top_n(15)

## Total frequency of occurance
bar<-ggplot(maxn.sum, aes(x=reorder(scientific,maxn), y=maxn)) +   
  geom_bar(stat="identity",position=position_dodge())+
  coord_flip()+
  xlab("Species")+
  ylab(expression(Overall~abundance~(Sigma~MaxN)))+
  #scale_x_discrete(limits = rev(levels(scientific)))+
  #annotation_custom(lcpic, xmin=0.5, xmax=2.5, ymin=.75, ymax=1.5)+ 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_collapse
bar

## filter the top 10 species and re plot--
maxn.10<-maxn.sum%>%
  filter(maxn>2)%>%
  filter(!scientific=="Carcharhinus sp1")

## Make labels for x axis
allbarlabs<-c("Terapon jarbua","Carcharhinus sorrah","Saurida undosquamis","Netuma thalassina","Caranx ignobilis","Pristipomoides multidens","Decapterus spp","Lagocephalus lunaris","Nemipterus spp","Carangoides chrysophrys")

## Top ten plot ----
bar.top.10<-ggplot(maxn.10, aes(x=reorder(scientific,maxn), y=maxn)) +   
  geom_bar(stat="identity",colour="black",fill="lightgrey",position=position_dodge())+
  ylim (0, 250)+
  coord_flip()+
  xlab("Species")+
  ylab(expression(Overall~abundance~(Sigma~MaxN)))+
  #scale_x_discrete(labels=allbarlabs)+
  #annotation_custom(lcpic, xmin=0.5, xmax=2.5, ymin=.75, ymax=1.5)+ 
  #   scale_y_log10()+
  # Apperance
  theme_bw()+
  theme(axis.text.y = element_text(face="italic"))+
  theme_collapse+
  theme.larger.text+
  annotation_raster(n.b, xmin=9.75,xmax=10.25,ymin=210, ymax=250)+
  annotation_raster(d.c, xmin=8.6,xmax=9.4,ymin=175, ymax=214)+
  annotation_raster(d.spp, xmin=7.6, xmax=8.5, ymin=50, ymax=108)+
  annotation_raster(c.e, xmin=6.6,xmax=7.4,ymin=28, ymax=65)+
  annotation_raster(n.v, xmin=5.7,xmax=6.3,ymin=15, ymax=50)+
  annotation_raster(s.u, xmin=4.7,xmax=5.3,ymin=12, ymax=67)+
  annotation_raster(c.p, xmin=3.4,xmax=4.7,ymin=9, ymax=105)+
  annotation_raster(l.l, xmin=2.7,xmax=3.3,ymin=9, ymax=63)+
  annotation_raster(s.l, xmin=1.55,xmax=2.5,ymin=6, ymax=95)+
  annotation_raster(d.spp, xmin=0.7,xmax=1.3,ymin=5, ymax=50)
bar.top.10

setwd(plots.dir)
ggsave("stacked.bar.plot.png",bar.top.10,dpi=600,width=5)

