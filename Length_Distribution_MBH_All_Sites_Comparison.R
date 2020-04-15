# Set directories----
rm(list=ls())

# Study name ----
study <- "ningaloo" 

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
working.dir<-dirname(rstudioapi::getActiveDocumentContext()$path)

# Read in length data 
setwd(tidy.dir)
lengths <- read.csv('ningaloo.complete.length.csv') %>%
  dplyr::mutate(scientific=paste(family,genus,species,sep=" "))

# Filter top ten species 

top.10.lengths <- dplyr::filter(lengths, grepl('Argyrops spinifer|Carangoides chrysophrys|Carangoides gymnostethus|Decapterus spp|Gymnocranius grandoculis|Lagocephalus sceleratus|Lethrinus miniatus|Lethrinus rubrioperculatus|Pristipomoides multidens|Lutjanus sebae|Lethrinus punctulatus|Diagramma pictum labiosum|Lethrinus ravus', scientific))
glimpse(top.10.lengths)

# Create a new verion that is only MBH sites 
top.10.lengths.MBH <- filter(lengths, !sample %in% c('10.13','10.14','10.15','14.11',
                                   '14.12','14.13','23.09','30.01',
                                   '30.01','30.02','30.03','30.04',
                                   '30.05','30.06','31.01','31.02',
                                   '31.03','31.04','31.05','31.06'))

## Create length plot for each species 

theme_ga<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Decapterus spp
length.decapterus <- ggplot(data=filter(top.10.lengths, species%in%c('spp')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.decapterus

length.decapterus.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('spp')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.decapterus.MBH

# Pristipomoides multidens
length.multidens <- ggplot(data=filter(top.10.lengths, species%in%c('multidens')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.multidens

length.multidens.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('multidens')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.multidens.MBH

# Gymnocranius grandoculis
length.grandoculis <- ggplot(data=filter(top.10.lengths, species%in%c('grandoculis')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.grandoculis

length.grandoculis.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('grandoculis')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.grandoculis.MBH

# Carangoides chrysophrys
length.chrysophrys <- ggplot(data=filter(top.10.lengths, species%in%c('chrysophrys')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.chrysophrys

length.chrysophrys.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('chrysophrys')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.chrysophrys.MBH

# Lethrinus miniatus
length.miniatus <- ggplot(data=filter(top.10.lengths, species%in%c('miniatus')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.miniatus

length.miniatus.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('miniatus')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.miniatus.MBH

# Carangoides gymnostethus
length.gymnostethus <- ggplot(data=filter(top.10.lengths, species%in%c('gymnostethus')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.gymnostethus

length.gymnostethus.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('gymnostethus')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.gymnostethus.MBH

# Lethrinus rubrioperculatus
length.rubrioperculatus <- ggplot(data=filter(top.10.lengths, species%in%c('rubrioperculatus')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.rubrioperculatus

length.rubrioperculatus.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('rubrioperculatus')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.rubrioperculatus.MBH

# Argyrops spinifer
length.spinifer <- ggplot(data=filter(top.10.lengths, species%in%c('spinifer')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.spinifer

length.spinifer.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('spinifer')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.spinifer.MBH

# Lagocephalus sceleratus
length.sceleratus <- ggplot(data=filter(top.10.lengths, species%in%c('sceleratus')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.sceleratus

length.sceleratus.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('sceleratus')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.sceleratus.MBH

# # Pterocaesio diagramma  
# length.diagramma <- ggplot(data=filter(top.10.lengths, species%in%c('diagramma')), aes(as.numeric(length))) +
#   geom_histogram(aes(y =..density..),col="black",fill="grey",alpha = .5)+
#   scale_y_continuous(expand = expansion(mult = c(0, .1)))+
#   labs(x = "Length (mm)", y = "Density") +
#   theme_ga
# length.diagramma
# 
# length.diagramma.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('diagramma')), aes(as.numeric(length))) +
#   geom_histogram(aes(y =..density..),col="black",fill="grey",alpha = .5)+
#   scale_y_continuous(expand = expansion(mult = c(0, .1)))+
#   labs(x = "Length (mm)", y = "Density") +
#   theme_ga
# length.diagramma.MBH

# Lutjanus sebae
length.sebae <- ggplot(data=filter(top.10.lengths, species%in%c('sebae')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.sebae

length.sebae.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('sebae')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.sebae.MBH

# Lethrinus punctulatus
length.punctulatus <- ggplot(data=filter(top.10.lengths, species%in%c('punctulatus')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.punctulatus

length.punctulatus.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('punctulatus')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.punctulatus.MBH

# Diagramma pictum labiosum  
length.pictum.labiosum <- ggplot(data=filter(top.10.lengths, species%in%c('pictum labiosum')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.pictum.labiosum

length.pictum.labiosum.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('pictum labiosum')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.pictum.labiosum.MBH

# Lethrinus ravus
length.ravus <- ggplot(data=filter(top.10.lengths, species%in%c('ravus')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.ravus

length.ravus.MBH <- ggplot(data=filter(top.10.lengths.MBH, species%in%c('ravus')), aes(as.numeric(length))) +
  geom_histogram(aes(y =..ncount..),col="black",fill="grey",alpha = .5)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "Length (mm)", y = "Scaled Count") +
  theme_ga
length.ravus.MBH

# Save the pairs of plots

setwd(plots.dir)

comparison.decapterus<-plot_grid(length.decapterus, length.decapterus.MBH,
                             labels = c('All', 'MBH'), label_size = 12,ncol=2)
ggsave("comparison.decapterus.png",comparison.decapterus,dpi=300,width=10,height=3.75)

comparison.multidens<-plot_grid(length.multidens, length.multidens.MBH,
                                 labels = c('All', 'MBH'), label_size = 12,ncol=2)
ggsave("comparison.multidens.png",comparison.multidens,dpi=300,width=10,height=3.75)

comparison.grandoculis<-plot_grid(length.grandoculis, length.grandoculis.MBH,
                                labels = c('All', 'MBH'), label_size = 12,ncol=2)
ggsave("comparison.grandoculis.png",comparison.grandoculis,dpi=300,width=10,height=3.75)

comparison.chrysophrys<-plot_grid(length.chrysophrys, length.chrysophrys.MBH,
                                  labels = c('All', 'MBH'), label_size = 12,ncol=2)
ggsave("comparison.chrysophrys.png",comparison.chrysophrys,dpi=300,width=10,height=3.75)

comparison.miniatus<-plot_grid(length.miniatus, length.miniatus.MBH,
                                  labels = c('All', 'MBH'), label_size = 12,ncol=2)
ggsave("comparison.miniatus.png",comparison.miniatus,dpi=300,width=10,height=3.75)

comparison.gymnostethus<-plot_grid(length.gymnostethus, length.gymnostethus.MBH,
                                  labels = c('All', 'MBH'), label_size = 12,ncol=2)
ggsave("comparison.gymnostethus.png",comparison.gymnostethus,dpi=300,width=10,height=3.75)

comparison.rubrioperculatus<-plot_grid(length.rubrioperculatus, length.rubrioperculatus.MBH,
                                  labels = c('All', 'MBH'), label_size = 12,ncol=2)
ggsave("comparison.rubrioperculatus.png",comparison.rubrioperculatus,dpi=300,width=10,height=3.75)

comparison.spinifer<-plot_grid(length.spinifer, length.spinifer.MBH,
                                  labels = c('All', 'MBH'), label_size = 12,ncol=2)
ggsave("comparison.spinifer.png",comparison.spinifer,dpi=300,width=10,height=3.75)

comparison.sceleratus<-plot_grid(length.sceleratus, length.sceleratus.MBH,
                                  labels = c('All', 'MBH'), label_size = 12,ncol=2)
ggsave("comparison.sceleratus.png",comparison.sceleratus,dpi=300,width=10,height=3.75)

comparison.sebae<-plot_grid(length.sebae, length.sebae.MBH,
                                  labels = c('All', 'MBH'), label_size = 12,ncol=2)
ggsave("comparison.sebae.png",comparison.sebae,dpi=300,width=10,height=3.75)

comparison.punctulatus<-plot_grid(length.punctulatus, length.punctulatus.MBH,
                                  labels = c('All', 'MBH'), label_size = 12,ncol=2)
ggsave("comparison.punctulatus.png",comparison.punctulatus,dpi=300,width=10,height=3.75)

comparison.pictum.labiosum<-plot_grid(length.pictum.labiosum, length.pictum.labiosum.MBH,
                                  labels = c('All', 'MBH'), label_size = 12,ncol=2)
ggsave("comparison.pictum.labiosum.png",comparison.pictum.labiosum,dpi=300,width=10,height=3.75)

comparison.ravus<-plot_grid(length.ravus, length.ravus.MBH,
                                  labels = c('All', 'MBH'), label_size = 12,ncol=2)
ggsave("comparison.ravus.png",comparison.ravus,dpi=300,width=10,height=3.75)








  
  