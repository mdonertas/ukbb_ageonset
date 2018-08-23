library(tidyverse)
library(RFRlib)
library(igraph)
library(GGally)
library(ggrepel)
library(RColorBrewer)
theme_set(theme_rfr())

traits <- readRDS('./data/processed/traits_clean/traitData_baseline.rds')
traits <- traits %>%
  mutate(Sex = c('Female','Male')[Sex+1])%>%
  mutate(BMI = Weight/((`Standing height`/100)^2))
traits[traits==-1 | traits ==-3] <- NA
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

sexcolors <- setNames(c('lavenderblush3', 'lightgoldenrodyellow'), c('Female', 'Male'))

SRdisease$Age[SRdisease$Age %in% c(-1,-3)]=NA
disCoding <- read_tsv('./data/raw/ukbb/datacoding/coding6.tsv')
disTree <- graph_from_data_frame(select(disCoding,parent_id,node_id),directed = T)
igraph_options(plot.layout=layout_as_tree)

# disDat <- lapply(setdiff(V(disTree)$name,'0'),function(nodex){
#   print(nodex)
#   childs <- subcomponent(disTree,v=nodex,mode='out')$name
#   nodeinfo <- filter(disCoding,node_id==nodex)
#   nodedat <- SRdisease %>%
#     filter(node_id %in% childs) %>%
#     group_by(eid) %>%
#     summarise(diseaseID=nodeinfo$coding,
#               Disease=nodeinfo$meaning,
#               node_id=nodex,
#               parent_id=nodeinfo$parent_id,
#               selectable=nodeinfo$selectable,
#               Age=min(Age,na.rm=T))
#   nodedat$Age[nodedat$Age==Inf]=NA
#   return(nodedat)
# })
# SRdisease <- reshape2::melt(disDat,id.vars=colnames(disDat[[1]]))%>%
#   select(-L1)
# saveRDS(SRdisease,'./data/processed/traits_clean/SRdisease_baseline_propagated.rds')

SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated.rds')

traits <- SRdisease%>%
  group_by(eid)%>%
  summarise(nSRdiseaseProp=length(unique(Disease)))%>%
  right_join(traits)

pdf('./results/UKBB_disease_EDA/disTree.pdf',width=20,height = 5)
plot(disTree,
     vertex.size=0.75,
     vertex.frame.color=NA,
     vertex.label=NA,
     asp=0.1,
     edge.arrow.size=0.2,
     edge.color='gray80',
     vertex.color=c('gray70'))
dev.off()

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
saveRDS(prevDF,file='./data/processed/traits_clean/SRdisease_prevdf.rds')
nCase_dist <- prevDF %>%
  ggplot(aes(x=nCases))+
  geom_histogram(color='gray60',bins=50)+
  scale_x_continuous(trans='log10',labels = scales::comma,
                     breaks=c(1,10,100,1000,10000,100000))+
  geom_vline(xintercept = 2000,color='orange',linetype='dashed')+
  ylab('Frequency')+
  xlab('Number of Cases (log scale)')+
  annotate('text',x=2100,y=25,label='2000 cases', hjust=0, color='orange') +
  annotate('label',x=1900,y=30,label=paste('N=',sum(prevDF$nCases<2000),sep=''), hjust=1) +
  annotate('label',x=2100,y=30,label=paste('N=',sum(prevDF$nCases>=2000),sep=''), hjust=0)
# system('mkdir -p results/UKBB_disease_EDA')
ggsave('./results/UKBB_disease_EDA/nCase_dist.pdf', nCase_dist,  useDingbats=F, height=4)

sexPrev <- reshape2::melt(table(Female=c('<0.001','>=0.001')[1+(ifelse(is.na(prevDF$fPrev),0,prevDF$fPrev)>=0.001)],
                                Male=c('<0.001','>=0.001')[1+(ifelse(is.na(prevDF$mPrev),0,prevDF$mPrev)>=0.001)])) %>%
  mutate(type=c('-','-','-','Common Diseases')) %>%
  ggplot(aes(x=Female, y=Male))+
  geom_tile(aes(fill=type),color='gray25')+
  theme_minimal(base_size = 20)+
  scale_fill_manual(values = c('gray70','lightseagreen'))+
  guides(fill=F)+
  geom_label(aes(label=paste('N=',value,sep='')),size=12)
ggsave('./results/UKBB_disease_EDA/sexPrev.pdf', sexPrev, useDingbats= F, height=4,width=6)

disSet <- prevDF %>%
  na.omit() %>%
  filter(nCases>=2000 & fPrev>=0.001 & mPrev>=0.001)%>%
  filter(Disease!='unclassifiable')
saveRDS(disSet,file='./data/processed/traits_clean/SRdiseaseSet.rds')

selectedNodes <- filter(disCoding,meaning%in%disSet$Disease)$node_id
pdf('./results/UKBB_disease_EDA/disTree_selected.pdf',width=20,height = 5)
plot(disTree,
     vertex.size=c(0.75,0.75)[1+(V(disTree)$name%in%selectedNodes)],
     vertex.frame.color=NA,
     vertex.label=NA,
     asp=0.1,
     edge.arrow.size=0.2,
     edge.color='gray80',
     vertex.color=c('gray70','gray25')[1+(V(disTree)$name%in%selectedNodes)])
dev.off()

nCases_sex <- disSet %>%
  mutate(FemMinMale=(fPrev-mPrev)*1000,
         Male=mPrev*1000,
         Female=-fPrev*1000)%>%
  gather(Sex,numCases,-FemMinMale,-Disease,-fPrev,-mPrev,-nCases,-fCase) %>%
  ggplot(aes(x=reorder(Disease,FemMinMale),y=numCases, fill=Sex))+
  geom_bar(stat='identity', color='gray50')+
  theme_rfr() +
  scale_fill_manual(values=sexcolors[c('Female','Male')])+
  coord_flip()+
  xlab('')+
  ylab('Number of Cases in 1,000')+
  scale_y_continuous(breaks=c(-400,-300,-200,-100,0,100,200,300,400),
                     labels=c(400,300,200,100,0,100,200,300,400),limits = c(-450,450))
ggsave('./results/UKBB_disease_EDA/nCases_sex.pdf', nCases_sex, useDingbats=F, height=18, width = 12)

nCases <- disSet %>% 
  mutate(nCases=fCase*1000)%>%
  ggplot(aes(x=reorder(Disease,-nCases),y=nCases))+
  geom_bar(stat='identity')+
  theme_rfr(x.text.angle = 90)+
  xlab('')+ylab('Number of Cases in 1,000')
ggsave('./results/UKBB_disease_EDA/nCases.pdf', nCases, useDingbats=F,height=8,width=16)

nSRdisease_bySex <- traits %>%
  select(nSRdisease, Sex) %>%
  ggplot(aes(x = nSRdisease,fill=Sex)) +
  geom_bar(position = 'dodge',color='gray60') +
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  ylab('Number of participants') +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks=seq(1,30,by=1)) +
  xlab('Number of Self-Reported Diseases')+
  theme(legend.position = 'top')
ggsave('./results/UKBB_EDA/nSRdisease_bySex.pdf', nSRdisease_bySex, useDingbats=F, height = 4,width = 6)

nSRdiseaseProp_bySex <- traits %>%
  select(nSRdiseaseProp, Sex) %>%
  ggplot(aes(x = nSRdiseaseProp,fill=Sex)) +
  geom_bar(position = 'dodge',color='gray60') +
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  ylab('Number of participants') +
  scale_y_continuous(labels = scales::comma) +
  xlab('Number of Self-Reported Diseases (Propagated)')+
  theme(legend.position = 'top')
ggsave('./results/UKBB_EDA/nSRdiseaseProp_bySex.pdf', nSRdiseaseProp_bySex, useDingbats=F, height = 4,width = 6)

eidx=traits %>%
  arrange(-nSRdisease) %>%
  head(5) %>%
  select(eid)

for( i in 1:5){
  selectedNodes1=unique((eidx[i,] %>%
                           left_join(SRdisease_np))$node_id)
  
  selectedNodes2=unique((eidx[i,] %>%
                           left_join(SRdisease))$node_id)
  pdf(paste('./results/UKBB_disease_EDA/top5NumDis_',i,'_tree.pdf',sep=''),width=20,height = 12)
  par(mfcol=c(2,1))
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

SRcaseMat <- SRdisease %>%
  filter(Disease%in%disSet$Disease)%>%
  select(eid,Disease) %>%
  unique() %>%
  na.omit()%>%
  mutate(value=1) %>%
  spread(Disease,value,fill=0) %>%
  as.data.frame()

rownames(SRcaseMat) <- SRcaseMat$eid
SRcaseMat$eid <- NULL
SRcaseMat <- as.matrix(SRcaseMat)

SRComMat <- cor(SRcaseMat,use = 'pairwise')
diag(SRComMat)=NA
pheatmap::pheatmap(SRComMat,
                   cellwidth = 10, cellheight = 10, 
                   breaks = seq(-max(abs(SRComMat),na.rm=T),max(abs(SRComMat),na.rm=T),length.out = 20),
                   filename='./results/UKBB_disease_EDA/disComMat.pdf',
                   # display_numbers = T, number_format = '%.2f',
                   color = rev(colorRampPalette(RColorBrewer::brewer.pal(5,'RdBu'))(19)))

SRComMat2=SRComMat
SRComMat2[abs(SRComMat)<0.05]=0
SRComMat2=SRComMat2+1
SRComMat2[SRComMat2==1]=0
myg=graph_from_adjacency_matrix(SRComMat2,mode = 'undirected',weighted = T,diag = F)

nodelabels=setNames(disCoding$meaning,disCoding$node_id)

nms=sapply(strsplit(V(myg)$name,'/'),function(x)x[1])

disTreecl <- reshape::melt(sapply(neighbors(disTree,'0')$name,function(nm){
  subcomponent(disTree,nm,'out')$name
}))%>%
  rename(node=value,
         cluster=L1)%>%
  mutate(node=unname(nodelabels[as.character(node)]),
         cluster=unname(nodelabels[as.character(cluster)]))

disTreecl=setNames(disTreecl$cluster,disTreecl$node)
disTreecl=factor(unname(disTreecl[V(myg)$name]))


E(myg)$color=ifelse(E(myg)$weight>1,'darkred','midnightblue')
disCooccurNet <- ggnet2(myg,size=5, edge.color = "color",color=c(brewer.pal(8,'Dark2'),
                                                                 brewer.pal(8,'Set1'))[as.numeric(disTreecl)],
                        edge.alpha = 0.2,
                        edge.size = 0.5)+
  geom_text_repel(label=nms,
                  size=6,
                  color=c(brewer.pal(8,'Dark2'),
                          brewer.pal(8,'Set1'))[as.numeric(disTreecl)])

ggsave('./results/UKBB_disease_EDA/disCooccur.pdf',disCooccurNet,width = 20,height = 15)

disCooccurNet_nolabel <- ggnet2(myg,size=5, edge.color = "color",color=c(brewer.pal(8,'Dark2'),
                                                                         brewer.pal(8,'Set1'))[as.numeric(disTreecl)],
                                edge.alpha = 0.2,
                                edge.size = 0.5)

ggsave('./results/UKBB_disease_EDA/disCooccur_nolabel.pdf',disCooccurNet_nolabel,width = 20,height = 15)