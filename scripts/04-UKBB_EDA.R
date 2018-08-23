library(tidyverse)
library(RFRlib)

theme_set(theme_rfr())

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

sexcolors <- setNames(c('lavenderblush3', 'lightgoldenrodyellow'), c('Female', 'Male'))

sexdistrib <- traits %>%
  group_by(Sex) %>%
  summarise(n = length(unique(eid)))%>%
  ggplot(aes(x = Sex, y = n)) +
  geom_bar(stat = 'identity', aes(fill = Sex)) + 
  geom_label(aes(label = n)) +
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  ylab('Number of participants') + 
  xlab('') +
  scale_y_continuous(labels = scales::comma)+
  guides(fill = F)
# system('mkdir -p results/UKBB_EDA')
ggsave('./results/UKBB_EDA/sex.pdf', sexdistrib, useDingbats=F, height = 4,width = 3)

traits %>%
  select(`Age at recruitment`, `Age when attended assessment centre`) %>%
  filter( `Age when attended assessment centre` != `Age at recruitment`)

agedistrib_bysex <- traits %>%
  select(`Age when attended assessment centre`, Sex) %>%
  ggplot(aes(x = `Age when attended assessment centre`, fill = Sex)) +
  facet_wrap(~Sex) + 
  geom_histogram(binwidth = 2, color = 'gray25') +
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  guides(fill = F) +
  ylab('Number of participants') +
  scale_y_continuous(labels = scales::comma)
ggsave('./results/UKBB_EDA/baselineAge_bySex.pdf', agedistrib_bysex, useDingbats=F, height = 4,width = 6)

agedistrib <- traits %>%
  select(`Age when attended assessment centre`) %>%
  ggplot(aes(x = `Age when attended assessment centre`)) +
  geom_histogram(binwidth = 2, color = 'gray60') +
  ylab('Number of participants') +
  scale_y_continuous(labels = scales::comma)
ggsave('./results/UKBB_EDA/baselineAge.pdf', agedistrib, useDingbats=F, height = 4,width = 6)

totagex <- data.frame(totNum=sapply(seq(0,90,by=5),function(agex){
  xx=select(traits,`Age when attended assessment centre`)%>%
    mutate(`Age when attended assessment centre`=`Age when attended assessment centre`+5)
  sum(xx[[1]]>=agex)
}),age=seq(0,90,by=5))

agebinner <- function(x,bin=5,min=0,max=100){
  sq <- seq(min,max,by=bin)
  sapply(x,function(i){sq[which.min(i>=sq)-1]})
}
annotdf <- traits %>%
  select(`Age at death`, Sex) %>%
  na.omit() %>%
  group_by(Sex) %>%
  summarise(n = n())
Age_at_death_bySex <- traits %>%
  select(`Age at death`, Sex) %>%
  na.omit()%>%
  mutate(`Age at death`=agebinner(`Age at death`,5))%>%
  group_by(`Age at death`,Sex)%>%
  summarise(n=n())%>%
  rename(age=`Age at death`)%>%
  left_join(totagex)%>%
  mutate(fracDead=n/totNum*1000)%>%
  mutate(fracDead=ifelse(Sex=='Male',fracDead,-fracDead))%>%
  ggplot(aes(x=age,y=fracDead,fill=Sex))+
  geom_bar(stat='identity',color='gray25',aes(fill=Sex))+
  scale_fill_manual(values = sexcolors)+
  geom_label(data = annotdf, x = 38, y= c(-185,185), aes( label = paste('N=',n,sep='')), hjust = c(0,1))+ 
  ylab('Number of Participants in 1000')+xlab('Age at death')+
  coord_flip()+
  theme(legend.position = 'top')+
  scale_y_continuous(breaks=seq(-150,150,by=50),labels=abs(seq(-150,150,by=50)),limits = c(-170,170))+
  scale_x_continuous(breaks=seq(40,75,5))
ggsave('./results/UKBB_EDA/Age_at_death_bySex.pdf', Age_at_death_bySex, useDingbats=F, height = 4,width = 5)

annotdf <- traits %>%
  select(`Age at death`) %>%
  na.omit() %>%
  summarise(n = n())
Age_at_death <- traits %>%
  select(`Age at death`) %>%
  na.omit()%>%
  mutate(`Age at death`=agebinner(`Age at death`,5))%>%
  group_by(`Age at death`)%>%
  summarise(n=n())%>%
  rename(age=`Age at death`)%>%
  left_join(totagex)%>%
  mutate(fracDead=n/totNum*1000)%>%
  ggplot(aes(x=age,y=fracDead))+
  geom_bar(stat='identity',color='gray25')+
  geom_label(data = annotdf, x = 38, y= 250, aes( label = paste('N=',n,sep='')), hjust = c(0))+
  ylab('Number of Participants in 1000')+xlab('Age at death')+
  scale_x_continuous(breaks=seq(40,75,5))

ggsave('./results/UKBB_EDA/Age_at_death.pdf', Age_at_death,useDingbats=F, height = 4,width = 4)

annotdf <- traits %>%
  select(`Parent Min Age at Death`) %>%
  na.omit() %>%
  summarise(n = n())
parentAgeatDeath <- traits %>%
  select(`Parent Min Age at Death`,`Parent Avg Age at Death`,`Parent Max Age at Death`)%>%
  gather(Type,`Age at death`)%>%
  mutate(Type=gsub('Parent | Age at Death','',Type))%>%
  na.omit()%>%
  ggplot(aes(x = Type, y= `Age at death`)) +
  geom_violin(fill='gray60')+
  geom_boxplot(width=0.05,outlier.shape = NA)+
  ggtitle(paste('Parent Age at Death, N=',annotdf$n,sep=''))+
  xlab('')
ggsave('./results/UKBB_EDA/Parent_Age_at_death.pdf', parentAgeatDeath, useDingbats=F, height = 4,width = 5)

cotraits <- traits[,c(2, 8:10, 12, 13, 15, 17:19, 21, 23:27)] %>%
  mutate(Sex=as.numeric(as.factor(Sex)))
cox=cor(cotraits,use='pairwise')
pheatmap::pheatmap(cox, cellwidth = 25, cellheight = 25, display_numbers = T, 
                   breaks = seq(-1,1,length.out = 20),
                   color = colorRampPalette(c('#2166AC','gray90','#B2182B'))(19),
                   filename = './results/UKBB_EDA/correlations.pdf',
                   number_format = '%.2f')
