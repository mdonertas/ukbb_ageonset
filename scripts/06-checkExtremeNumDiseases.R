library(tidyverse)
library(igraph)
library(ggnet)
library(intergraph)

SRdisease_prop <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated.rds')
SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline.rds')
SRdisease$Age[SRdisease$Age %in% c(-1,-3)]=NA
SRdisease_prop$Age[SRdisease_prop$Age %in% c(-1,-3)]=NA
disCoding <- read_tsv('./data/raw/ukbb/datacoding/coding6.tsv')
disTree <- graph_from_data_frame(select(disCoding,parent_id,node_id),directed = T)
igraph_options(plot.layout=layout_as_tree)

eidx=SRdisease %>%
        group_by(eid) %>%
        summarise(n=length(unique(diseaseID))) %>%
        arrange(-n) %>%
        head(5) %>%
        select(eid)

for( i in 1:5){
  selectedNodes1=unique((eidx[i,] %>%
                           left_join(SRdisease))$node_id)
  
  selectedNodes2=unique((eidx[i,] %>%
                           left_join(SRdisease_prop))$node_id)
  
  pdf(paste('./results/UKBB_disease_EDA/top5NumDis_',i,'_tree.pdf',sep=''),width=20,height = 5)
  plot(disTree,
       vertex.size=c(0.75,0.75)[1+(V(disTree)$name%in%selectedNodes1)],
       vertex.frame.color=NA,
       vertex.label=NA,
       asp=0.1, 
       main=paste('Disease tree for an individual with',length(selectedNodes1),'diseases'),
       edge.arrow.size=0.2,
       edge.color='gray80',
       vertex.color=c('gray70','gray25')[1+(V(disTree)$name%in%selectedNodes1)])
  
  plot(disTree,
       vertex.size=c(0.75,0.75)[1+(V(disTree)$name%in%selectedNodes2)],
       vertex.frame.color=NA,
       vertex.label=NA,
       asp=0.1,
       main=paste('Disease tree for an individual with',length(selectedNodes1),'diseases\nThis propagated version includes',length(selectedNodes2),'nodes'),
       edge.arrow.size=0.2,
       edge.color='gray80',
       vertex.color=c('gray70','gray25')[1+(V(disTree)$name%in%selectedNodes2)])
  dev.off()
} 
