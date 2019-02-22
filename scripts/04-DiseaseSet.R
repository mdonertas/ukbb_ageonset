library(tidyverse)
library(igraph)

traits <- readRDS('./data/processed/traits_clean/traitData_baseline_additions.rds')
disCoding <- read_tsv('./data/raw/ukbb/datacoding/coding6.tsv')
disTree <- graph_from_data_frame(select(disCoding,parent_id,node_id),directed = T)
SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated.rds')

prevDF <- SRdisease %>%
  select(Disease,Age,eid) %>%
  right_join(select(traits,eid,Sex)) %>%
  na.omit()%>%
  group_by(Disease,Sex)%>%
  summarise(nCases = length(unique(eid))) %>%
  spread(Sex, nCases) %>%
  mutate(fPrev=Female/sum(traits$Sex=='Female'),
         mPrev=Male/sum(traits$Sex=='Male'),
         nCases=sum(Female,Male,na.rm=T)) %>%
  mutate(fCase=nCases/(sum(traits$Sex=='Female')+sum(traits$Sex=='Male')))%>%
  ungroup()
saveRDS(prevDF,file='./data/processed/traits_clean/SRdisease_prop_prevdf.rds')

disSet <- prevDF %>%
  na.omit() %>%
  filter(nCases>=2000 & fPrev>=0.001 & mPrev>=0.001)%>%
  filter(Disease!='unclassifiable')
saveRDS(disSet,file='./data/processed/traits_clean/SRdiseaseSet.rds')

sum(prevDF$fPrev>=0.001 & prevDF$mPrev>=0.001, na.rm=T)
