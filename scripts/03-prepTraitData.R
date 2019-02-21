library(tidyverse)
library(data.table)
library(RFRlib)

ukbfields <- read_tsv('./data/processed/ukbb/helperfiles/ukbfields.tsv')  

samplestouse <- read_tsv('./data/processed/ukbb/sampleQC/sampleIDs_passedQC.tsv')

## read data

fields <- ukbfields %>%
  filter(description %in% c('eid', 'Sex', 'Age at recruitment', 
                            'Age when attended assessment centre',
                            'Age at death', 'Standing height', 'Weight',
                            'Number of self-reported cancers',
                            'Number of self-reported non-cancer illnesses',
                            'Number of operations, self-reported', 
                            'Number of treatments/medications taken',
                            'Sleep duration', 'Facial ageing', 
                            'Maternal smoking around birth', 'Smoking status', 
                            'Alcohol drinker status', "Father's age at death",
                            "Mother's age at death", 'Overall health rating', 
                            'Health satisfaction', 
                            'Age when periods started (menarche)', 
                            'Age at menopause (last menstrual period)', 
                            'Non-accidental death in close genetic family')) %>%
  mutate(fields = as.character(fields)) %>%
  filter(visit == 0)

ukbb <- setDF(fread('./data/raw/ukbb/ukb10904.csv',
                    select = fields$fields,
                    col.names = fields$description)) %>%
  filter(eid %in% samplestouse$eid)
dim(ukbb)

system('mkdir -p data/processed/traits_clean')
saveRDS(ukbb, './data/processed/traits_clean/traitData_baseline.rds')
write_tsv(ukbb, './data/processed/traits_clean/traitData_baseline.tsv')

SRdisease_fields <- ukbfields %>%
  filter(description %in% c('eid', 'Non-cancer illness code, self-reported',
                            'Interpolated Age of participant when non-cancer illness first diagnosed')) %>%
  mutate(fields = as.character(fields)) %>%
  filter(visit == 0)

SRcancer_fields <- ukbfields %>%
  filter(description %in% c('eid', 'Cancer code, self-reported',
                            'Interpolated Age of participant when cancer first diagnosed')) %>%
  mutate(fields = as.character(fields)) %>%
  filter(visit == 0)

ukbb_srDisease <- setDF(fread('./data/raw/ukbb/ukb10904.csv',
                              select = SRdisease_fields$fields)) %>%
  filter(eid %in% samplestouse$eid)

diseaseCodes <- read_tsv('./data/raw/ukbb/datacoding/coding6.tsv') %>%
  rename(diseaseID = coding) 

ukbb_srDisease.diseases <- ukbb_srDisease %>%
  select(eid, grep('20002-', colnames(ukbb_srDisease))) %>%
  gather(key = order, value = diseaseID, -eid) %>%
  na.omit() %>%
  separate(order, into = c('field_id','visit','num')) %>%
  mutate(visit = as.numeric(visit),
         num = as.numeric(num)) %>%
  left_join(SRdisease_fields, by = c('field_id','visit','num')) %>%
  left_join(diseaseCodes, by = c('diseaseID'))%>%
  as.tibble()

ukbb_srDisease.age <- ukbb_srDisease %>%
  select(eid, grep('20009-', colnames(ukbb_srDisease))) %>%
  gather(key = order, value = Age, -eid) %>%
  na.omit() %>%
  separate(order, into = c('field_id','visit','num')) %>%
  mutate(visit = as.numeric(visit),
         num = as.numeric(num)) %>%
  left_join(SRdisease_fields, by = c('field_id','visit','num')) %>%
  as.tibble()

ukbb_srDisease <- left_join(ukbb_srDisease.diseases, ukbb_srDisease.age,
                            by = c('eid', 'visit', 'num')) %>%
  select(eid, diseaseID, meaning, node_id, parent_id, selectable, Age) %>%
  rename(Disease = meaning)

rm(ukbb_srDisease.age, ukbb_srDisease.diseases)

saveRDS(ukbb_srDisease, './data/processed/traits_clean/SRdisease_baseline.rds')
write_tsv(ukbb_srDisease, './data/processed/traits_clean/SRdisease_baseline.tsv')

##

ukbb_srCancer <- setDF(fread('./data/raw/ukbb/ukb10904.csv',
                             select = SRcancer_fields$fields)) %>%
  filter(eid %in% samplestouse$eid)

cancerCodes <- read_tsv('./data/raw/ukbb/datacoding/coding3.tsv') %>%
  rename(cancerID = coding) 

ukbb_srCancer.cancers <- ukbb_srCancer %>%
  select(eid, grep('20001-', colnames(ukbb_srCancer))) %>%
  gather(key = order, value = cancerID, -eid) %>%
  na.omit() %>%
  separate(order, into = c('field_id','visit','num')) %>%
  mutate(visit = as.numeric(visit),
         num = as.numeric(num)) %>%
  left_join(SRcancer_fields, by = c('field_id','visit','num')) %>%
  left_join(cancerCodes, by = c('cancerID'))%>%
  as.tibble()

ukbb_srCancer.age <- ukbb_srCancer %>%
  select(eid, grep('20007-', colnames(ukbb_srCancer))) %>%
  gather(key = order, value = Age, -eid) %>%
  na.omit() %>%
  separate(order, into = c('field_id','visit','num')) %>%
  mutate(visit = as.numeric(visit),
         num = as.numeric(num)) %>%
  left_join(SRcancer_fields, by = c('field_id','visit','num')) %>%
  as.tibble()

ukbb_srCancer <- left_join(ukbb_srCancer.cancers, ukbb_srCancer.age,
                           by = c('eid', 'visit', 'num')) %>%
  select(eid, cancerID, meaning, node_id, parent_id, selectable, Age) %>%
  rename(Cancer = meaning)

rm(ukbb_srCancer.cancers, ukbb_srCancer.age)

saveRDS(ukbb_srCancer, './data/processed/traits_clean/SRcancer_baseline.rds')
write_tsv(ukbb_srCancer, './data/processed/traits_clean/SRcancer_baseline.tsv')

traits <- readRDS('./data/processed/traits_clean/traitData_baseline.rds') %>%
  mutate(Sex = c('Female','Male')[Sex+1])%>%
  mutate(BMI = Weight/((`Standing height`/100)^2))
traits[traits==-1 | traits ==-3] <- NA

traits$`Parent Min Age at Death` <- traits %>%
  select(`Mother's age at death`, `Father's age at death`) %>%
  apply(.,1,function(x){
    xx=min(x[!is.na(x)])
    ifelse(xx==Inf,NA,xx)
  })
traits$`Parent Max Age at Death` <- traits %>%
  select(`Mother's age at death`, `Father's age at death`) %>%
  apply(.,1,function(x){
    xx=max(x[!is.na(x)])
    ifelse(xx==-Inf,NA,xx)
  })
traits$`Parent Avg Age at Death` <- traits %>%
  select(`Mother's age at death`, `Father's age at death`) %>%
  apply(.,1,function(x){
    xx=mean(x[!is.na(x)])
    ifelse(is.nan(xx),NA,xx)
  })

SRdisease_np <- readRDS('./data/processed/traits_clean/SRdisease_baseline.rds')
traits <- SRdisease_np%>%
  group_by(eid)%>%
  summarise(nSRdisease=length(unique(Disease)))%>%
  right_join(traits)
SRcancer <- readRDS('./data/processed/traits_clean/SRcancer_baseline.rds')
traits <- SRcancer%>%
  group_by(eid)%>%
  summarise(nSRcancer=length(unique(Cancer)))%>%
  right_join(traits)
SRdisease_np$Age[SRdisease_np$Age %in% c(-1,-3)]=NA
SRcancer$Age[SRcancer$Age %in% c(-1,-3)]=NA
disCoding <- read_tsv('./data/raw/ukbb/datacoding/coding6.tsv')

library(igraph)
disTree <- graph_from_data_frame(select(disCoding,parent_id,node_id),directed = T)

disDat <- lapply(setdiff(V(disTree)$name,'0'),function(nodex){
  print(nodex)
  childs <- subcomponent(disTree,v=nodex,mode='out')$name
  nodeinfo <- filter(disCoding,node_id==nodex)
  nodedat <- SRdisease_np %>%
    filter(node_id %in% childs) %>%
    group_by(eid) %>%
    summarise(diseaseID=nodeinfo$coding,
              Disease=nodeinfo$meaning,
              node_id=nodex,
              parent_id=nodeinfo$parent_id,
              selectable=nodeinfo$selectable,
              Age=min(Age,na.rm=T))
  nodedat$Age[nodedat$Age==Inf]=NA
  return(nodedat)
})
SRdisease <- reshape2::melt(disDat,id.vars=colnames(disDat[[1]]))%>%
  select(-L1)
saveRDS(SRdisease,'./data/processed/traits_clean/SRdisease_baseline_propagated.rds')

traits <- SRdisease%>%
  group_by(eid)%>%
  summarise(nSRdiseaseProp=length(unique(Disease)))%>%
  right_join(traits)

saveRDS(traits, './data/processed/traits_clean/traitData_baseline_additions.rds')

## after creating inplink_notin_bgen.fam using prepGWAS/02-removeFiles4bolt

traits <- readRDS('./data/processed/traits_clean/traitData_baseline_additions.rds')
SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline.rds')
SRdisease_prop <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated.rds')
SRcancer <- readRDS('./data/processed/traits_clean/SRcancer_baseline.rds')

notinbgen <- read_delim('./data/processed/ukbb/gwas/remove/inplink_notin_bgen.fam', col_names = F, delim = ' ')

mean(notinbgen$X1 %in% traits$eid)
# [1] 0
mean(notinbgen$X2 %in% traits$eid)
# [1] 0

mean(notinbgen$X1 %in% SRdisease$eid)
# [1] 0
mean(notinbgen$X2 %in% SRdisease$eid)
# [1] 0

mean(notinbgen$X1 %in% SRdisease_prop$eid)
# [1] 0
mean(notinbgen$X2 %in% SRdisease_prop$eid)
# [1] 0

mean(notinbgen$X1 %in% SRcancer$eid)
# [1] 0
mean(notinbgen$X2 %in% SRcancer$eid)
# [1] 0

