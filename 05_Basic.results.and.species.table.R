# Set directories----
rm(list=ls())

# Study name ----
study <- "2020-01_Guardian-Ningaloo_stereoBRUVs" 

# Libraries required
install_github("UWAMEGFisheries/GlobalArchive") #to check for updates
library(GlobalArchive)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(googlesheets4)

## Set your working directory ----
working.dir<-dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

# Set sub directories----
plots.dir=paste(working.dir,"plots",sep="/")
tidy.dir=paste(working.dir,"tidy data",sep="/")

# Read in maxn ----
setwd(tidy.dir)
dir()

maxn <- read_csv("2020-01_Guardian-Ningaloo_stereoBRUVs.complete.maxn.csv")%>%
  mutate(scientific=paste(family,genus,species,sep=" "))%>%
  glimpse()

# Read in habitat ----
habitat<-read_csv("2020-01_Guardian-Ningaloo_stereoBRUVs._habitat.csv" )

# Read in metadata ----
metadata<-read_csv("2020-01_Guardian-Ningaloo_stereoBRUVs.checked.metadata.csv")%>%
  dplyr::select(sample,latitude,longitude,date,time,site,depth)

# Read in life history
url <- "https://docs.google.com/spreadsheets/u/1/d/1SMLvR9t8_F-gXapR2EemQMEPSw_bUbPLcXd3lJ5g5Bo/edit?usp=drive_web&ouid=100340010373917954123"

master<-read_sheet(url)%>%
  ga.clean.names()%>%
  filter(grepl('Australia', global.region))%>%
  filter(grepl('NW', marine.region))%>%
  dplyr::select(family,genus,species,iucn.ranking,fishing.mortality,fishing.type,australian.common.name)%>% 
  distinct()%>%
  glimpse()

names(master)

# Create Species list ----

species.table <- maxn%>%
  group_by(family,genus,species,scientific)%>%
  summarise_at(vars(matches("maxn")),funs(sum,mean,sd,se=sd(.)/sqrt(n())))%>%
  ungroup()%>%
  mutate(mean=round(mean,digits=2))%>%
  mutate(sd=round(sd,digits=2))%>%
  mutate(se=round(se,digits=2))%>%
  mutate(genus.species=paste(genus,species,sep=" "))%>%
  arrange(family)%>%
  left_join(master)%>%
  dplyr::select(-c(scientific))%>%
  dplyr::mutate(mean.relative.abundance.per.deployment.plus.minus.SE=paste(mean,"+/-",se,sep=" "))%>%
  dplyr::rename(total.relative.abundance = sum)%>%
  ungroup()

unique(species.table$fishing.type)

cleaned<-species.table%>%
  dplyr::select(family,genus.species,australian.common.name,fishing.type,iucn.ranking,mean.relative.abundance.per.deployment.plus.minus.SE,total.relative.abundance)%>%
  ## fix up variables
  mutate(fishing.type=ifelse(fishing.type%in%c("C/R","C","B/C"),"Commercial",""))
  ## Make names nicer for table

unique(cleaned$fishing.type)

# Descriptive stats

# total abundance
sum(maxn$maxn) # 504
length(unique(maxn$scientific)) # 16

length(unique(maxn$family)) # 13 families
length(unique(maxn$genus)) # 14 genus



# Make data for anita ----
summary<-maxn%>%
  #filter(species%in%c("bathybius","carpenteri","tabl","equula","virgatus","variegatus"))%>%
  group_by(sample,scientific)%>%
  dplyr::summarise_at(vars(matches("maxn")),funs(sum))%>%
  mutate(presence=as.integer(maxn != 0))
  
presence<-summary%>%
  dplyr::select(-c(maxn))%>%
  spread(scientific,value=presence)

maxn.summary<-summary%>%
  dplyr::select(-c(presence))%>%
  spread(scientific,value=maxn)

anita.presence<-left_join(metadata,habitat)%>%
  left_join(presence)

anita.maxn<-left_join(metadata,habitat)%>%
  left_join(maxn.summary)

setwd(tidy.dir)
write.csv(anita.presence,"presence.spatial.model.csv",row.names = FALSE)
write.csv(anita.maxn,"maxn.spatial.model.csv",row.names = FALSE)
