respath='./'
library(tidyverse)
library(RFRlib)
library(ggridges)
library(igraph)

theme_set(theme_rfr(legend.pos = 'top'))

exclusions <- unique(c(read_delim('./data/processed/ukbb/gwas/remove/inplink_notin_bgen.fam',delim=' ',col_names = F)$X1,
                       read_delim('./data/processed/ukbb/gwas/remove/sampleQC_exc.fam',delim=' ',col_names = F)$X1,
                       read_delim('./data/processed/ukbb/gwas/remove/withdrawn.fam',delim=' ',col_names = F)$X1))

traits <- readRDS('./data/processed/traits_clean/traitData_baseline.rds') %>%
  filter(!eid %in% exclusions)
traits <- traits %>%
  mutate(Sex = c('Female','Male')[Sex+1])%>%
  mutate(BMI = Weight/((`Standing height`/100)^2))
traits[traits==-1 | traits ==-3] <- NA

sexcolors <- setNames(c('lavenderblush3', 'lightgoldenrodyellow'), c('Female', 'Male'))

diseaseCoding <- read_tsv('./data/raw/ukbb/datacoding/coding6.tsv')
disSet <- diseaseCoding %>%
  filter(coding %in% c(1278,1309,1075,1453,1111,1440,1065))
disTree <- graph_from_data_frame(data.frame(parent=diseaseCoding$parent_id,
                                            node=diseaseCoding$node_id), directed=T)

selectedNodes <- disSet$node_id
nodeLabels <- setNames(diseaseCoding$meaning,diseaseCoding$node_id)
V(disTree)$desc <- nodeLabels[as.character(V(disTree)$name)]

V(disTree)$desc[!V(disTree)$name %in% c(names(which(sapply(setdiff(neighbors(disTree,v='0')$name,'99999'),function(nm)any(subcomponent(disTree,nm,'out')$name%in%selectedNodes)))),selectedNodes)]=''
# V(disTree)$desc[!V(disTree)$name %in% selectedNodes]=''

pdf(paste(respath,'results/pilot/disTree_selected.pdf',sep=''),width=20,height = 5)
plot(disTree,
     layout = layout_as_tree,
     vertex.size=c(0.75,1)[1+(V(disTree)$name%in%selectedNodes)],
     vertex.frame.color=NA,
     vertex.label=V(disTree)$desc,
     vertex.label.color='gray25',
     vertex.label.size=0.5,
     vertex.label.dist=1,
     vertex.label.degree=sample(seq(0,360,by=30),1),
     asp=0.1,
     edge.arrow.size=0.2,
     edge.color='gray80',
     vertex.color=c('gray70','gray25')[1+(V(disTree)$name%in%selectedNodes)])
dev.off()

SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated.rds') %>%
  filter(!eid %in% exclusions)

SRdisease <- traits %>%
  select(eid,Sex)%>%
  right_join(SRdisease,by='eid')

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

saveRDS(prevDF,'./data/processed/traits_clean/SRdisease_prevdf.rds')

disSetUlt <- prevDF %>%
  na.omit() %>%
  filter(nCases>=2000 & fPrev>=0.001 & mPrev>=0.001)%>%
  filter(Disease!='unclassifiable')

saveRDS(disSetUlt,'./data/processed/traits_clean/SRdiseaseSet.rds')

SRdiseasesub <- SRdisease %>%
  filter(node_id %in% selectedNodes)

luniq(SRdiseasesub$Disease)
# 7

disStats <- SRdiseasesub %>%
  group_by(Disease) %>%
  summarise( n = luniq(eid),
             minAge = min(Age, na.rm=T),
             meanAge = mean(Age, na.rm=T),
             maxAge = max(Age, na.rm=T),
             medAge = median(Age, na.rm=T),
             sdAge = sd(Age, na.rm=T))

disStats %>%
  knitr::kable()

# |Disease                            |      n|    minAge|  meanAge|  maxAge|   medAge|     sdAge|
# |:----------------------------------|------:|---------:|--------:|-------:|--------:|---------:|
# |asthma                             |  56094| 0.0000000| 31.47984| 70.0548| 32.50000| 18.836419|
# |cataract                           |   7038| 0.0000000| 55.96208| 70.2437| 59.50000| 12.894936|
# |heart attack/myocardial infarction |  11082| 0.0136895| 53.07790| 69.5000| 54.50000|  8.713935|
# |hypertension                       | 129927| 0.4900850| 50.86008| 70.2821| 52.50000| 10.432344|
# |osteoporosis                       |   7751| 0.5000000| 56.22231| 70.0767| 57.43840|  7.736277|
# |psoriasis                          |   5429| 0.0000000| 31.34932| 69.7645| 29.50865| 17.320483|
# |tuberculosis (tb)                  |   2473| 0.0000000| 17.27097| 67.8480| 13.50000| 14.346507|

write_csv(x = disStats, path = paste(respath,'results/pilot/disStats.csv',sep=''))


disStats_bySex <- SRdiseasesub %>%
  group_by(Disease,Sex) %>%
  summarise( n = luniq(eid),
             minAge = min(Age, na.rm=T),
             meanAge = mean(Age, na.rm=T),
             maxAge = max(Age, na.rm=T),
             medAge = median(Age, na.rm=T),
             sdAge = sd(Age, na.rm=T))
disStats_bySex %>%
  knitr::kable()

# |Disease                            |Sex    |     n|    minAge|  meanAge|  maxAge|   medAge|     sdAge|
# |:----------------------------------|:------|-----:|---------:|--------:|-------:|--------:|---------:|
# |asthma                             |Female | 32369| 0.0000000| 33.80781| 69.5000| 35.50000| 17.901877|
# |asthma                             |Male   | 23725| 0.0000000| 28.40630| 70.0548| 27.50000| 19.584532|
# |cataract                           |Female |  3992| 0.0000000| 56.74727| 70.2437| 60.50000| 12.354857|
# |cataract                           |Male   |  3046| 0.0000000| 54.93588| 70.1561| 58.50000| 13.501827|
# |heart attack/myocardial infarction |Female |  2170| 0.5000000| 54.08720| 69.5000| 55.50000|  9.043412|
# |heart attack/myocardial infarction |Male   |  8912| 0.0136895| 52.83773| 69.5000| 53.50000|  8.616812|
# |hypertension                       |Female | 62591| 0.5000000| 50.11482| 70.1397| 51.50000| 11.231060|
# |hypertension                       |Male   | 67336| 0.4900850| 51.51712| 70.2821| 52.50000|  9.625983|
# |osteoporosis                       |Female |  6853| 1.5000000| 56.43824| 70.0767| 57.50000|  7.465782|
# |osteoporosis                       |Male   |   898| 0.5000000| 54.56548| 70.0384| 56.36520|  9.401777|
# |psoriasis                          |Female |  2527| 0.0000000| 29.37995| 69.3730| 25.08735| 18.210615|
# |psoriasis                          |Male   |  2902| 0.0000000| 33.07218| 69.7645| 32.30035| 16.311781|
# |tuberculosis (tb)                  |Female |  1350| 0.0000000| 16.69725| 67.8480| 12.50000| 13.892318|
# |tuberculosis (tb)                  |Male   |  1123| 0.0000000| 17.96443| 65.9643| 14.50000| 14.853779|

write_csv(x = disStats_bySex, path = paste(respath,'results/pilot/disStats_bySex.csv',sep=''))

numPeople <- disStats_bySex %>%
  ggplot(aes(x=reorder(Disease,-n),y=n,fill=Sex))+
  geom_bar(stat = 'identity',position = 'dodge',color='gray60') +
  ylab('Frequency')+
  xlab('')+
  scale_fill_manual(values = sexcolors)+
  guides(fill=guide_legend(''))+
  theme_rfr(x.text.angle = 90, legend.pos = 'right')

ggsave(filename = paste(respath,'results/pilot/numCases.pdf',sep=''),plot = numPeople,width = 5,height = 6)

ageDist <- SRdiseasesub %>%
  select(eid,Disease,Age,Sex) %>%
  unique() %>%
  na.omit()%>%
  left_join(disStats) %>%
  ggplot(aes(x=Age,y=reorder(Disease,-medAge),fill=Sex))+
  facet_grid(~Sex)+
  geom_density_ridges(bandwidth=2,color='gray50')+
  ylab('') +
  scale_fill_manual(values = sexcolors)+
  xlab('Age of onset')+
  guides(fill=guide_legend(''))

ggsave(filename = paste(respath,'results/pilot/ageDist.pdf',sep=''),plot = ageDist, height = 6,width = 8)

SRcaseMat <- SRdiseasesub %>%
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
SRComMat[SRComMat<0]=0
diag(SRComMat)=NA
pheatmap::pheatmap(SRComMat,
                   cellwidth = 10, cellheight = 10, 
                   breaks = seq(0,1,length.out = 20),
                   # display_numbers = T, number_format = '%.2f',
                   filename = paste(respath,'results/pilot/disComMat.pdf',sep=''),
                   color = c('white',(colorRampPalette(RColorBrewer::brewer.pal(5,'Reds'))(19))))

nSRdis <- SRdiseasesub %>%
  select(eid)%>%
  unique() %>%
  left_join(SRdisease)%>%
  group_by(eid)%>%
  summarise(nSRdis=luniq(Disease)) %>%
  ungroup()

nSRdis_byDisease <- nSRdis %>%
  left_join(SRdiseasesub) %>%
  left_join(disStats)%>%
  ggplot(aes(y=nSRdis,x=Disease,fill=Sex))+
  xlab('')+
  ylab('Number of SR diseases')+
  scale_fill_manual(values = sexcolors)+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  theme_rfr(x.text.angle = 90, legend.pos = 'top')

ggsave(filename = paste(respath,'results/pilot/nSRdis_byDisease.pdf',sep=''),plot = nSRdis_byDisease, width=6, height=7)

nSRcan_byDisease <- SRdiseasesub %>%
  select(eid,Sex,Disease) %>%
  left_join(select(traits,eid,'Number of self-reported cancers')) %>%
  ggplot(aes(y=`Number of self-reported cancers`, x= Disease, fill=Sex))+
  xlab('')+ ylab('Number of self-reported cancers')+
  scale_fill_manual(values=sexcolors)+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  theme_rfr(x.text.angle = 90, legend.pos = 'top')

ggsave(filename = paste(respath,'results/pilot/nSRcan_byDisease.pdf',sep=''),plot = nSRcan_byDisease, width=6, height=7)

dis_cancerStats <- SRdiseasesub %>%
  select(eid,Sex,Disease) %>%
  left_join(select(traits,eid,'Number of self-reported cancers')) %>%
  group_by(Disease,`Number of self-reported cancers`,Sex) %>%
  summarise(n=luniq(eid))

write_csv(x = dis_cancerStats,path = paste(respath,'results/pilot/dis_cancerStats.csv',sep=''))
