# Clear memory ----
rm(list=ls())

# Libraries required ----
# To connect to GlobalArchive
library(devtools)
install_github("UWAMEGFisheries/GlobalArchive") #to check for updates
library(GlobalArchive)
# To connect to life.history
library(httpuv)
# To tidy data
library(tidyr)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(fst)


# Study name---
study<-"ningaloo"  ## change for your project

## Set your working directory ----
working.dir<-dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

## Save these directory names to use later----
staging.dir<-paste(working.dir,"Staging",sep="/") 
download.dir<-paste(working.dir,"EM Export",sep="/")
tidy.dir<-paste(working.dir,"Tidy data",sep="/")
plots.dir=paste(working.dir,"Plots",sep="/")
error.dir=paste(working.dir,"Errors to check",sep="/")

# Read in the data----
setwd(staging.dir)
dir()

# Read in metadata----
metadata<-read_csv(file=paste(study,"metadata.csv",sep = "_"),na = c("", " "))%>%
  mutate(sample=as.character(sample))%>%
  glimpse()

# Read in habitat----
setwd(download.dir)
dir()

forwards <- read.delim("2019-08_ningaloo_forwards_Dot Point Measurements.txt",skip=4,header=T)%>%
  ga.clean.names()%>%
  mutate(sample=str_replace_all(filename,c(".jpg"="")))%>%
  dplyr::select(-c(filename,collector,date,depth,fishing.status,frame,latitude,longitude,location,radius.,rugosity,site.,spare,spare.1,time.mins,transect.))%>%
  mutate(direction="forwards")%>%
  glimpse()

backwards <- read.delim("2019-08_ningaloo_backwards_Dot Point Measurements.txt",skip=4,header=T)%>%
  ga.clean.names()%>%
  mutate(sample=str_replace_all(filename,c(".jpg"="",".JPG"="")))%>%
  dplyr::select(-c(filename,collector,date,depth,fishing.status,frame,latitude,longitude,location,radius.,rugosity,site.,spare,spare.1,time.mins,transect.))%>%
  mutate(direction="backwards")%>%
  mutate(image.col=image.col+10000)%>%
  mutate(image.row=image.row+10000)%>%
  glimpse()

test.forwards<-forwards%>%
  dplyr::group_by(sample)%>%
  dplyr::summarise(number=n())

missing <- anti_join(forwards,metadata, by = c("sample"))
missing <- anti_join(metadata,forwards, by = c("sample"))

test.backwards<-backwards%>%
  dplyr::group_by(sample)%>%
  dplyr::summarise(number=n())

missing <- anti_join(backwards,metadata, by = c("sample"))
missing <- anti_join(metadata,backwards, by = c("sample")) 

habitat<-bind_rows(forwards,backwards)

names(habitat)%>%sort()


fov.point.score <- habitat %>%
  dplyr::select(sample,starts_with("fieldofview"))%>%
  glimpse()

fov.point.score <- habitat %>%
  #dplyr::select(sample,starts_with("fieldofview"))%>%
  glimpse()



# FOV point score ----
fov.point.score <- habitat %>%
  dplyr::select(sample,starts_with("fieldofview"))%>%
  filter(!fieldofview%in%c(NA,""))%>% #check ones that don't have field of view
  dplyr::mutate(count = 1) %>%
  mutate(fieldofview=paste("fieldofview",fieldofview,sep="."))%>%
  dplyr::group_by(sample,fieldofview) %>%
  dplyr::summarise(count=sum(count))%>%
  spread(key = fieldofview, value = count, fill=0)%>%
  ungroup()

# FOV percent cover ----
fov.percent.cover <- fov.point.score%>%
  dplyr::mutate(total.sum=rowSums(.[,2:(ncol(.))],na.rm = TRUE ))%>%
  dplyr::group_by(sample) %>%
  mutate_at(vars(starts_with("fieldofview.")),funs(./total.sum*100))%>%
  mutate_at(vars(starts_with("fieldofview.")),funs(round(.,digits=2)))%>%
  dplyr::select(-total.sum) %>%
  ga.clean.names()%>%
  # dplyr::left_join(metadata) %>%
  glimpse()

# Create relief----
relief.mean.and.sd <- habitat %>%
  dplyr::select(sample,starts_with("relief"))%>%
  filter(!relief%in%c(NA,""))%>% #check ones that don't have field of view
  mutate(relief.rank=ifelse(relief=="0. Flat substrate, sandy, rubble with few features. ~0 substrate slope.",0,
                            ifelse(relief=="1. Some relief features amongst mostly flat substrate/sand/rubble. <45 degree substrate slope.",1,
                            ifelse(relief=="2. Mostly relief features amongst some flat substrate or rubble. ~45 substrate slope.",2,
                            ifelse(relief=="3. Good relief structure with some overhangs. >45 substrate slope.",3,
                            ifelse(relief=="4. High structural complexity, fissures and caves. Vertical wall. ~90 substrate slope.",4,
                            ifelse(relief=="5. Exceptional structural complexity, numerous large holes and caves. Vertical wall. ~90 substrate slope.",5,
                                   relief)))))))%>%
  mutate(relief.rank=as.numeric(relief.rank))%>%
  group_by(sample) %>%
  dplyr::summarise(mean.relief= mean (relief.rank), sd.relief= sd (relief.rank))%>%
  glimpse()

# CREATE catami_broad------
broad.point.score <- habitat %>%
  dplyr::select(sample,starts_with("broad"))%>%
  dplyr::filter(!broad%in%c("Open Water","Unknown"))%>%
  dplyr::mutate(count = 1) %>%
  mutate(broad=paste("broad",broad,sep="."))%>%
  dplyr::group_by(sample,broad) %>%
  dplyr::summarise(count=sum(count))%>%
  spread(key = broad, value = count, fill=0) %>%
  ga.clean.names()%>%
  ungroup()

unique(habitat$broad)

# broad percent cover ----
broad.percent.cover <- broad.point.score%>%
  dplyr::mutate(total.sum=rowSums(.[,2:(ncol(.))],na.rm = TRUE ))%>%
  dplyr::group_by(sample) %>%
  mutate_at(vars(starts_with("broad.")),funs(./total.sum*100))%>%
  mutate_at(vars(starts_with("broad.")),funs(round(.,digits=2)))%>%
  dplyr::select(-total.sum) %>%
  # dplyr::left_join(metadata) %>%
  glimpse()

# Write final habitat data----
setwd(tidy.dir)
dir()

habitat.relief.fov<-relief.mean.and.sd%>%
  full_join(fov.percent.cover,by=c("sample"))%>%
  full_join(broad.percent.cover,by=c("sample"))%>%
  # semi_join(metadata)%>%
  # left_join(metadata)%>%
  mutate(campaignid = "2019-08_Ningaloo_stereo-BRUVs")%>%
  glimpse()


habitat.fst<-habitat.relief.fov
names(habitat.fst)<-ga.capitalise(names(habitat.fst))

habitat.fst<-habitat.fst%>%
  mutate(CampaignID="test")

write.csv(habitat.relief.fov, file=paste(study,"complete.habitat.csv",sep = "."), row.names=FALSE)
write.fst(habitat.fst,paste(study,"complete.habitat.fst",sep = "."))
