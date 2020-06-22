library(raster)
library(rstudioapi)
library(tidyverse)
library(sf)
library(dplyr)
library(sp)
library(rgdal)
library(rgeos)
library(spatialEco)

##Set working directory----
## Set work directory----
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # sets working directory to where this script is saved (DON't MOVE THE SCRIPT)

## Set sub directories----
d.dir = paste(working.dir,"Tidy data",sep="/") 
s.dir = paste(working.dir,"Spatial",sep="/") # spatial is where I keep spatial data files, rasters and shapefiles
p.dir <- paste(working.dir,"Plots",sep="/")

############# Extracting just targeted fish from dataset ##############
library(httpuv)
library(googlesheets4)


complete.species <- read.csv("ningaloo.checked.maxn.csv")
complete.species <- complete.species%>%
  mutate(scientific=paste(family,genus,species,sep=" "))%>%
  select("sample", "scientific", "family", "genus", "species")%>%
  filter(!grepl("NA NA NA", scientific))

# Just extract targeted species

url <- "https://docs.google.com/spreadsheets/d/1SMLvR9t8_F-gXapR2EemQMEPSw_bUbPLcXd3lJ5g5Bo/edit?ts=5e6f36e2#gid=825736197"

target<-googlesheets4::read_sheet(url)%>%
  filter(grepl('Australia', Global.Region))%>% # Change country here
  filter(grepl('NW', Marine.region))%>%
  dplyr::mutate(aLL=as.numeric(aLL))%>%
  dplyr::mutate(bLL=as.numeric(bLL))%>%
  dplyr::mutate(a=as.numeric(a))%>%
  dplyr::mutate(b=as.numeric(b))%>%
  dplyr::rename(family=Family)%>%
  dplyr::rename(genus=Genus)%>%
  dplyr::rename(species=Species)%>%
  mutate(scientific=paste(family,genus,species,sep=" "))%>%
  select(scientific, family, genus, species, FB.Length_MAX,Fishing.mortality, Fishing.type, MinLegal.WA)%>%
  glimpse()


target.species<-left_join(complete.species,target,by=c("scientific"))%>%
  select(!c("family.y","genus.y","species.y"))%>%
  dplyr::rename(family=family.x)%>%
  dplyr::rename(genus=genus.x)%>%
  dplyr::rename(species=species.x)%>%
  filter(!grepl("N", Fishing.mortality))%>% #Remove non-fished species
  filter(!grepl("NA", Fishing.mortality))%>% #Remove any individuals we couldn't ID or put as spp
  filter(grepl("R",Fishing.type)) #Only select those that are recreationally fished 

write.csv(target.species, 'target.species.csv')


################ Adding target classification and FB max length to all observations ##################
target.species<-read.csv('target.species.csv')
species.lengths<-read.csv('ningaloo.complete.length.csv')

# Create dataframe that just has one row for each targeted species
target.species<-target.species[!duplicated(target.species$scientific), ]

# Format species.lengths table
species.lengths<-species.lengths%>%
  mutate(scientific=paste(family,genus,species,sep=" "))%>%
  select(sample,scientific,family,genus,species,number,length,site)

#Add two dataframes together to give only targeted fish species at each sample location
target.lengths <- inner_join(species.lengths,target.species, by='scientific')%>%
  select(sample.x,scientific,family.x,genus.x,species.x,number,length,site,FB.Length_MAX,MinLegal.WA,Fishing.mortality,
         Fishing.type)%>%
  mutate(min.length=0.15*FB.Length_MAX)%>%
  rename(sample=sample.x)%>%
  rename(family=family.x)%>%
  rename(genus=genus.x)%>%
  rename(species=species.x)
  
target.lengths<-target.lengths[,c(1,2,3,4,5,6,7,8,9,10,13,11,12)]

################# Determining size groups ##################

target.size.groups <- target.lengths%>%
  mutate(Cut.off=(MinLegal.WA-min.length)/2)

target.size.groups$MinLegal.WA[is.na(target.size.groups$MinLegal.WA)] <- "None"
 
  
target.size.groups <- target.size.groups%>%
   mutate(Size.Class=case_when(
         length>MinLegal.WA~'Legal',
         length<MinLegal.WA & length>Cut.off ~ "Large.sub.legal",
         length>min.length & length<Cut.off ~ "Small.sub.legal",
         length>min.length & MinLegal.WA=='None' ~ "Legal"))


## Filtering out where no individuals of targeted species were seen at a site

target.size.groups <- target.size.groups%>%
  filter(!is.na(Size.Class))

## Split into legal and large sublegal (we have no small sublegal)

legal.target <- subset(target.size.groups, Size.Class=='Legal')
sublegal.target <- subset(target.size.groups, Size.Class=="Large.sub.legal")

setwd(d.dir)
write.csv(legal.target, "Legal.Target.csv")
write.csv(sublegal.target, "Sublegal.Target.csv")


################### Combining the response data with the covariate data ##################
covariates <- read.csv('covariates.csv')
legal.target <- read.csv('Legal.Target.csv')
sublegal.target <- read.csv('Sublegal.Target.csv')

Sum.legal <- aggregate(legal.target$number, by=list(sample=legal.target$sample), FUN=sum)%>%
  rename(target.fish=x)
 
Sum.sublegal <- aggregate(sublegal.target$number, by=list(sample=sublegal.target$sample), FUN=sum)%>%
  rename(target.fish=x)

final.data.legal <- full_join(covariates,Sum.legal, by='sample')%>%
  select(-c('X','ID'))
final.data.legal$target.fish[is.na(final.data.legal$target.fish)] <- 0

final.data.sublegal <- full_join(covariates,Sum.sublegal, by='sample')%>%
  select(-c('X','ID'))
final.data.sublegal$target.fish[is.na(final.data.sublegal$target.fish)] <- 0


write.csv(final.data.legal, "final.data.legal.csv")
write.csv(final.data.sublegal, "final.data.sublegal.csv")  






















