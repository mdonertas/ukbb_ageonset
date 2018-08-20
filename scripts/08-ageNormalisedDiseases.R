library(tidyverse)
library(RFRlib)
library(ggthemes)
library(ggridges)
library(igraph)
library(intergraph)
library(ggpubr)

theme_set(theme_rfr(base_family = 'Helvetica',x.text.angle = 90))

traits <- readRDS('./data/processed/traits_clean/traitData_baseline.rds')
traits <- traits %>%
  mutate(Sex = c('Female','Male')[Sex+1])%>%
  mutate(BMI = Weight/((`Standing height`/100)^2))
traits[traits==-1 | traits ==-3] <- NA


sexcolors <- setNames(c('lavenderblush3', 'lightgoldenrodyellow'), c('Female', 'Male'))

SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated.rds')

prevDF <- SRdisease %>%
  select(Disease,Age,eid) %>%
  right_join(select(traits,eid,Sex)) %>% 
  na.omit()%>%
  group_by(Disease,Sex)%>%
  summarise(nCases = length(unique(eid))) %>%
  spread(Sex, nCases) %>%
  mutate(fPrev=Female/262801,
         mPrev=Male/221865,
         nCases=sum(Female,Male,na.rm=T)) %>%
  mutate(fCase=nCases/(262801+221865))
disSet <- prevDF %>%
  filter(nCases>=2000 & fPrev>=0.001 & mPrev>=0.001)%>%
  filter(Disease!='unclassifiable')

for(nint in c(5,10)){
  disAgeMat=SRdisease %>%
    filter(Disease%in%disSet$Disease) %>%
    mutate(ageGr=cut(floor(Age),seq(0,80,by=nint),right = T))%>%
    group_by(Disease,ageGr) %>%
    summarise(n=length(unique(eid)))%>%
    spread(ageGr,n,fill=0)%>%
    as.data.frame()
  
  
  rownames(disAgeMat)=disAgeMat$Disease
  disAgeMat=disAgeMat[,colnames(disAgeMat)!='<NA>']
  disAgeMat$Disease=NULL
  disAgeMat=as.matrix(disAgeMat)
  
  totagegr5=traits %>%
    select(eid,'Age when attended assessment centre')%>%
    mutate(AgeGrTot=cut(floor(`Age when attended assessment centre`),seq(0,200,by=nint),right=T))%>%
    group_by(AgeGrTot)%>%
    summarise(n=length(unique(eid)))
  totagegr5 %>%
    spread(AgeGrTot,n,fill=0)
  totagegr5=setNames(totagegr5$n,totagegr5$AgeGrTot)[colnames(disAgeMat)]
  names(totagegr5)=colnames(disAgeMat)
  totagegr5[is.na(totagegr5)]=0
  totagegr5=rev(cumsum(rev(totagegr5)))
  totagegr5
  disAgeMat2=t(apply(disAgeMat,1,function(x){
    xx=(x/totagegr5)
    xx/sum(x)
  }))
  ul=max(c(abs(min(disAgeMat2)),abs(max(disAgeMat2))))
  pheatmap::pheatmap(disAgeMat2,cluster_cols = F, cellwidth=15,cellheight = 15,
                     breaks = seq(0,ul,length.out = 20),
                     color = c('white',(colorRampPalette(RColorBrewer::brewer.pal(5,'Blues'))(19))),
                     cutree_rows = 7,
                     filename=paste('./results/UKBB_disease_EDA/ageOnsetHeatmap_Norm_',nint,'.pdf',sep=''))
}
