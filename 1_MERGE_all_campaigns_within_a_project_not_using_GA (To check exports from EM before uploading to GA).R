### Merge EventMeasure database output tables into maxn and length files

### OBJECTIVES ###
# combine database tables into single Metadata, MaxN and Length files for subsequent validation and data analysis.

### Please forward any updates and improvements to tim.langlois@uwa.edu.au & brooke.gibbons@uwa.edu.au or raise an issue in the "globalarchive-query" GitHub repository

rm(list=ls()) # Clear memory

## Load Libraries ----
# To connect to GlobalArchive
library(devtools)
install_github("UWAMEGFisheries/GlobalArchive") #to check for updates
library(GlobalArchive)
# To connect to GitHub
library(RCurl)
library(R.utils)
# To tidy data
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
# to connect to googlesheets
library(googlesheets)

## Set Study Name ----
# Change this to suit your study name. This will also be the prefix on your final saved files.
study<-"ningaloo" 

## Folder Structure ----
# This script uses one main folder ('working directory')
# Three subfolders will be created within the 'working directory'. They are 'EM Export','Staging' and 'Tidy data'
# Save the database exports into the 'EM Export' folder
# The 'Staging' folder is used to save the combined files (e.g. metadata, maxn or length) NOTE: These initial outputs have not gone through any check (e.g. checks against the life-history sheet)

# Naming in the script is extremely important! Names need to be consistent or the script will break

# **The only folder you will need to create is your working directory**

## Set your working directory ----
working.dir<-dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

## Save these directory names to use later----
staging.dir<-paste(working.dir,"Staging",sep="/") 
download.dir<-paste(working.dir,"EM Export",sep="/")
tidy.dir<-paste(working.dir,"Tidy data",sep="/")

setwd(working.dir)

## Create EM Export, Staging and Tidy data folders ----
dir.create(file.path(working.dir, "EM Export"))
dir.create(file.path(working.dir, "Staging"))
dir.create(file.path(working.dir, "Tidy data"))

# BEFORE CONTINUING ----
# You should now copy your database tables into the EM Export folder 

# Combine all data----
# The below code will find all files that have the same ending (e.g. "_Metadata.csv") and bind them together.
# The end product is three data frames; metadata, maxn and length.

# Metadata ----
# You will need a metadata file. This can either be a .csv file or a google sheet 
# If using a .csv the file name MUST end in "_Metadata.csv"
# Both csv and googlesheet will need to match the global archive format
# See the user manual: https://globalarchivemanual.github.io/ for the correct format
# In this example we will use a csv file (you will need to create a csv file to upload to GlobalArchive anyway but can use this script to save the file to upload to globalarchive)

# For csv file ----
setwd(working.dir)
metadata <-read.csv("ningaloo_metadata.csv")%>% # list all files ending in "_Metadata.csv"
  ga.clean.names()%>%
  mutate(sample=as.factor(sample))%>%
  # dplyr::select(sample,latitude,longitude,date,time,location,status,site,depth,observer,successful.count,successful.length,comment)%>% # This line ONLY keep the 15 columns listed. Remove or turn this line off to keep all columns (Turn off with a # at the front).
  glimpse()

dir()

unique(metadata$campaignid) # check the number of campaigns in metadata, and the campaign name

setwd(staging.dir)
write.csv(metadata,paste(study,"metadata.csv",sep="_"),row.names = FALSE)

## Combine Points and Count files into maxn ----
maxn<-ga.create.em.maxn()%>%
  dplyr::inner_join(metadata)%>%
  dplyr::filter(successful.count=="Yes")%>%
  dplyr::filter(maxn>0)

# Save MaxN file ----
setwd(staging.dir)
write.csv(maxn,paste(study,"maxn.csv",sep="_"),row.names = FALSE)

## Combine Length, Lengths and 3D point files into length3dpoints----
length3dpoints<-ga.create.em.length3dpoints()%>%
  dplyr::select(-c(time,comment))%>% # take time out as there is also a time column in the metadata
  dplyr::inner_join(metadata)%>%
  # dplyr::filter(successful.length=="Yes")%>% charlotte to add later
  glimpse()

## Save length files ----
setwd(staging.dir)
write.csv(length3dpoints,paste(study,"length3dpoints.csv",sep="_"),row.names = FALSE)

