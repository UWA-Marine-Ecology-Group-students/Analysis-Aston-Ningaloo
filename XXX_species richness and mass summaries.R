library(dplyr)
library(readr)
library(tidyr)

setwd("C:/GitHub/Analysis-Aston-Ningaloo/Tidy data")

# Read in maxn ----
maxn <- read_csv("ningaloo.complete.maxn.csv", na = c("", " "))%>%
  mutate(sample=as.factor(sample))%>%
  glimpse()

metadata<-read_csv("ningaloo.checked.metadata.csv",na = c("", " "))%>%
  mutate(sample=as.factor(sample))%>%
  glimpse()

samples<-metadata%>%distinct(sample)

ta.sr <- maxn%>%
  group_by(scientific,sample)%>%
  dplyr::summarise(maxn = sum(maxn))%>%
  spread(scientific,maxn, fill = 0)%>%
  mutate(total.abundance=rowSums(.[,2:(ncol(.))],na.rm = TRUE ))%>% #Add in Totals
  mutate(species.richness=rowSums(.[,2:167] > 0))%>%
  dplyr::select(sample,total.abundance,species.richness)%>%
  gather(.,"scientific","maxn",2:3)%>%
  left_join(metadata)

mass<-read_csv("ningaloo.complete.mass.csv",na = c("", " "),col_types = cols(.default = "c"))%>%
  mutate(mass.g=as.numeric(mass.g))%>%
  mutate(length.cm=as.numeric(length.cm))%>%
  glimpse()

total.mass<-mass%>%
  replace_na(list(mass.g=0))%>%
  dplyr::group_by(sample)%>%
  dplyr::summarise(mass.total=sum(mass.g))%>%
  dplyr::ungroup()%>%
  tidyr::complete(nesting(sample)) %>%
  full_join(samples)%>%
  replace_na(list(mass.total = 0))

mass.30cm<-mass%>%
  mutate(length.cm=as.numeric(length.cm))%>%
  filter(length.cm>30)%>%
  dplyr::group_by(sample)%>%
  dplyr::summarise(mass.30cm=sum(mass.g))%>%
  dplyr::ungroup()%>%
  tidyr::complete(nesting(sample)) %>%
  full_join(samples)%>%
  replace_na(list(mass.30cm = 0))

mass.20cm<-mass%>%
  mutate(length.cm=as.numeric(length.cm))%>%
  filter(length.cm>20)%>%
  dplyr::group_by(sample)%>%
  dplyr::summarise(mass.20cm=sum(mass.g))%>%
  dplyr::ungroup()%>%
  tidyr::complete(nesting(sample)) %>%
  full_join(samples)%>%
  replace_na(list(mass.20cm = 0))

mass.summaries<-left_join(total.mass,mass.20cm)%>%
  left_join(.,mass.30cm)

write.csv(ta.sr, file=paste("ningaloo","total.abundance.and.species.richness.csv",sep = "."), row.names=FALSE)
write.csv(ta.sr, file=paste("ningaloo","mass.summaries.csv",sep = "."), row.names=FALSE)
