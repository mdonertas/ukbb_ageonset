source('./scripts/00-setup.R')

traits <- readRDS('./data/processed/traits_clean/traitData_baseline_additions2.rds')

sexdistrib <- traits %>%
  group_by(Sex) %>%
  summarise(n = length(unique(eid)))%>%
  ggplot(aes(x = Sex, y = n)) +
  geom_bar(stat = 'identity', aes(fill = Sex)) +
  geom_label(aes(label = scales::comma(n)), size = 8 / pntnorm, vjust = 1, nudge_y = -2000) +
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  ylab('Number of participants') +
  xlab('') +
  scale_y_continuous(labels = scales::comma)+
  guides(fill = F) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5))
# system('mkdir -p results/UKBB_EDA')
ggsave('./results/UKBB_EDA/sex.pdf', sexdistrib, useDingbats=F, height = 6,width = 5, units = "cm")
ggsave('./results/UKBB_EDA/sex.png', sexdistrib, height = 6,width = 5, units = "cm")

traits %>%
  select(`Age at recruitment`, `Age when attended assessment centre`) %>%
  filter( `Age when attended assessment centre` != `Age at recruitment`)

traits %>%
  group_by(Sex) %>%
  select(`Age when attended assessment centre`) %>%
  mutate(age = `Age when attended assessment centre`) %>%
  summarise(min = min(age), max = max(age), median = median(age))

# # A tibble: 2 x 4
# Sex      min   max median
# <chr>  <dbl> <dbl>  <dbl>
# 1 Female    39    71     57
# 2 Male      37    73     58

agedistrib_bysex <- traits %>%
  select(`Age when attended assessment centre`, Sex) %>%
  ggplot(aes(x = `Age when attended assessment centre`, fill = Sex)) +
  facet_wrap(~Sex) +
  geom_histogram(binwidth = 2, color = 'gray25') +
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  guides(fill = F) +
  ylab('Number of participants') +
  xlab('Age') +
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5))
ggsave('./results/UKBB_EDA/baselineAge_bySex.pdf', agedistrib_bysex, useDingbats=F, height = 6,width = 8, units = "cm")
ggsave('./results/UKBB_EDA/baselineAge_bySex.png', agedistrib_bysex, height = 6,width = 8, units = "cm")

agedistrib <- traits %>%
  select(`Age when attended assessment centre`) %>%
  ggplot(aes(x = `Age when attended assessment centre`)) +
  geom_histogram(binwidth = 2, color = 'gray60') +
  ylab('Number of participants') + xlab('Age') + 
  scale_y_continuous(labels = scales::comma)
ggsave('./results/UKBB_EDA/baselineAge.pdf', agedistrib, useDingbats=F, height = 6,width = 8, units = "cm")
ggsave('./results/UKBB_EDA/baselineAge.png', agedistrib, height = 6,width = 8, units = "cm")

summary(traits$`Number of operations, self-reported`)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#   0.000   1.000   1.000   1.733   3.000  32.000     133

summary(traits$`Number of treatments/medications taken`)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#   0.000   0.000   2.000   2.456   4.000  48.000     133

numMed_Op <- traits %>%
  select(`Number of operations, self-reported`, 
         `Number of treatments/medications taken`, `Sex`) %>%
  rename(`Number of operations` = `Number of operations, self-reported`,
         `Number of medications` = `Number of treatments/medications taken`) %>% 
  gather(key = 'Type', value = 'number', `Number of operations`:`Number of medications`) %>%
  ggplot(aes(x = number, fill = Sex)) +
  facet_grid(Sex~Type, scales = 'free_y') + 
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  geom_bar() +
  scale_y_log10(labels = scales::comma) +
  xlab('') + ylab('Number of Participants\n(in log10 scale)')
ggsave('./results/UKBB_EDA/numMed_Op.pdf', numMed_Op, useDingbats=F, height = 6,width = 10, units = "cm")
ggsave('./results/UKBB_EDA/numMed_Op.png', numMed_Op, height = 6,width = 10, units = "cm")

facialAging = read_tsv("./data/raw/ukbb/datacoding/coding100435.tsv")
facialAging = setNames(c('Prefer not to answer', 'Do not know', 'Younger', 'Older', 'About same'), facialAging$coding)

healthrate = read_tsv("./data/raw/ukbb/datacoding/coding100508.tsv")
healthrate = setNames(healthrate$meaning, healthrate$coding)
healthsatisfaction = read_tsv("./data/raw/ukbb/datacoding/coding100478.tsv")
healthsatisfaction = setNames(healthsatisfaction$meaning, healthsatisfaction$coding)
smokealcohol = read_tsv("./data/raw/ukbb/datacoding/coding90.tsv")
smokealcohol = setNames(smokealcohol$meaning, smokealcohol$coding)

facialagingPlot = traits %>%
  select(Sex, `Facial ageing`) %>%
  na.omit() %>%
  ggplot(aes(fill = Sex, x = as.character(`Facial ageing`))) +
  geom_bar(position = 'dodge') + 
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  scale_x_discrete(labels = facialAging[as.character(na.omit(unique(traits$`Facial ageing`)))]) +
  scale_y_continuous(labels = scales::unit_format('',1/1000)) +
  ylab('Number of participants\n(in thousands)') +
  xlab('') +
  coord_flip() + ggtitle('Facial age')
ggsave('./results/UKBB_EDA/facialaging.pdf', facialagingPlot, useDingbats=F, height = 6,width = 7, units = "cm")
ggsave('./results/UKBB_EDA/facialaging.png', facialagingPlot, height = 6,width = 7, units = "cm")

healthratePlot = traits %>%
  select(Sex, `Overall health rating`) %>%
  na.omit() %>%
  ggplot(aes(fill = Sex, x = as.character(`Overall health rating`))) +
  geom_bar(position = 'dodge') + 
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  scale_x_discrete(labels = healthrate[as.character(na.omit(unique(traits$`Overall health rating`)))]) +
  scale_y_continuous(labels = scales::unit_format('',1/1000)) +
  ylab('Number of participants\n(in thousands)') +
  xlab('') +
  coord_flip() + ggtitle('Overall health rating')
ggsave('./results/UKBB_EDA/healthratePlot.pdf', healthratePlot, useDingbats=F, height = 5,width = 7, units = "cm")
ggsave('./results/UKBB_EDA/healthratePlot.png', healthratePlot, height = 5,width = 7, units = "cm")

healthsatisfactionPlot = traits %>%
  select(Sex, `Health satisfaction`) %>%
  na.omit() %>%
  ggplot(aes(fill = Sex, x = as.character(`Health satisfaction`))) +
  geom_bar(position = 'dodge') + 
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  scale_x_discrete(labels = healthsatisfaction[as.character(na.omit(unique(traits$`Health satisfaction`)))]) +
  scale_y_continuous(labels = scales::unit_format('',1/1000)) +
  ylab('Number of participants\n(in thousands)') +
  xlab('') +
  coord_flip() + ggtitle('Health satisfaction')
ggsave('./results/UKBB_EDA/healthsatisPlot.pdf', healthsatisfactionPlot, useDingbats=F, height = 6,width = 8, units = "cm")
ggsave('./results/UKBB_EDA/healthsatisPlot.png', healthsatisfactionPlot, height = 6, width = 8, units = "cm")

smokingPlot = traits %>%
  select(Sex, `Smoking status`) %>%
  na.omit() %>%
  ggplot(aes(fill = Sex, x = as.character(`Smoking status`))) +
  geom_bar(position = 'dodge') + 
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  scale_x_discrete(labels = smokealcohol[as.character(na.omit(unique(traits$`Smoking status`)))]) +
  scale_y_continuous(labels = scales::unit_format('',1/1000)) +
  ylab('Number of participants\n(in thousands)') +
  xlab('') +
  coord_flip() + ggtitle("Smoking Status")
ggsave('./results/UKBB_EDA/smoking.pdf', smokingPlot, useDingbats=F, height = 5,width = 7, units = "cm")
ggsave('./results/UKBB_EDA/smoking.png', smokingPlot, height = 5,width = 7, units = "cm")

alcoholPlot = traits %>%
  select(Sex, `Alcohol drinker status`) %>%
  na.omit() %>%
  ggplot(aes(fill = Sex, x = as.character(`Alcohol drinker status`))) +
  geom_bar(position = 'dodge') + 
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  scale_x_discrete(labels = smokealcohol[as.character(na.omit(unique(traits$`Alcohol drinker status`)))]) +
  scale_y_continuous(labels = scales::unit_format('',1/1000)) +
  ylab('Number of participants\n(in thousands)') +
  xlab('') +
  coord_flip() + ggtitle("Alcohol drinker status")
ggsave('./results/UKBB_EDA/alcohol.pdf', alcoholPlot, useDingbats=F, height = 5,width = 7, units = "cm")
ggsave('./results/UKBB_EDA/alcohol.png', alcoholPlot, height = 5,width = 7, units = "cm")


traits %>%
  select(`Age at death`) %>%
  na.omit() %>%
  nrow()
# 13697


totagex.f <- data.frame(totNum=sapply(seq(0,90,by=5),function(agex){
  xx=select(filter(traits, Sex == "Female"), `Age when last person died`)
  sum(xx[[1]]>=agex)
}),age=seq(0,90,by=5), Sex = "Female")
totagex.m <- data.frame(totNum=sapply(seq(0,90,by=5),function(agex){
  xx=select(filter(traits, Sex == "Male"), `Age when last person died`)
  sum(xx[[1]]>=agex)
}),age=seq(0,90,by=5), Sex = "Male")
totagex = rbind(totagex.f, totagex.m)
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
  geom_label(data = annotdf, x = 38, y= c(-45,45), aes( label = paste('N=',n,sep='')), hjust = c(0,1),
             color = "white", size = 8 / pntnorm)+
  ylab('Number of Participants in 1000')+xlab('Age at death')+
  coord_flip()+
  theme(legend.position = 'top')+
  scale_y_continuous(breaks=seq(-40,40,by=10),labels=abs(seq(-40,40,by=10)),limits = c(-45,45))+
  scale_x_continuous(breaks=seq(40,75,5))
ggsave('./results/UKBB_EDA/Age_at_death_bySex.pdf', Age_at_death_bySex, useDingbats=F, height = 6,width = 6, units = "cm")
ggsave('./results/UKBB_EDA/Age_at_death_bySex.png', Age_at_death_bySex, height = 6,width = 6, units = "cm")

totagex<- data.frame(totNum=sapply(seq(0,90,by=5),function(agex){
  xx=select(traits, `Age when last person died`)
  sum(xx[[1]]>=agex)
}),age=seq(0,90,by=5))
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
ggsave('./results/UKBB_EDA/Age_at_death.pdf', Age_at_death,useDingbats=F, height = 6,width = 6, units ="cm")
ggsave('./results/UKBB_EDA/Age_at_death.png', Age_at_death, height = 6,width = 6, units ="cm")

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
  ggtitle(paste('Parent Age at Death\nN=',scales::comma(annotdf$n),sep=''))+
  xlab('')
ggsave('./results/UKBB_EDA/Parent_Age_at_death.pdf', parentAgeatDeath, useDingbats=F, height = 6,width = 6, units = "cm")
ggsave('./results/UKBB_EDA/Parent_Age_at_death.png', parentAgeatDeath, height = 6,width = 6, units = "cm")

traits %>%
  select(`Parent Min Age at Death`,`Parent Avg Age at Death`,`Parent Max Age at Death`) %>%
  summary()
# Parent Min Age at Death Parent Avg Age at Death Parent Max Age at Death
# Min.   : 10.00          Min.   : 10.00          Min.   : 10.00         
# 1st Qu.: 60.00          1st Qu.: 65.00          1st Qu.: 69.00         
# Median : 69.00          Median : 73.00          Median : 78.00         
# Mean   : 67.33          Mean   : 71.56          Mean   : 75.79         
# 3rd Qu.: 77.00          3rd Qu.: 79.50          3rd Qu.: 84.00         
# Max.   :115.00          Max.   :115.00          Max.   :117.00         
# NA's   :92756           NA's   :92756           NA's   :92756
NSRdisease <- traits %>%
  select(Sex, nSRdiseaseProp) %>%
  na.omit() %>%
  ggplot(aes(x = Sex, y= nSRdiseaseProp)) +
  geom_violin(aes(fill = Sex)) +
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) + 
  geom_boxplot(width = 0.05, outlier.shape = NA) + 
  xlab('') + ylab('Number of Self-reported\nDiseases (Propagated)')
ggsave('./results/UKBB_EDA/NSRdisease.pdf', NSRdisease, useDingbats=F, height = 7,width = 6, units = "cm")
ggsave('./results/UKBB_EDA/NSRdisease.png', NSRdisease, height = 7,width = 6, units = "cm")

sum(!is.na(traits$nSRdiseaseProp))
# [1] 362772
mean(!is.na(traits$nSRdiseaseProp))
# [1] 0.748604
sum(is.na(traits$nSRdiseaseProp))
# [1] 121826

NSRcancer <- traits %>%
  select(Sex, nSRcancer) %>%
  na.omit() %>%
  ggplot(aes(x= nSRcancer)) +
  geom_bar(aes(fill = Sex), position = 'dodge') +
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  ylab('Number of Participants') + xlab('Number of Self-reported Cancers')
ggsave('./results/UKBB_EDA/nSRcancer.pdf', NSRcancer, useDingbats=F, height = 6,width = 7, units = "cm")
ggsave('./results/UKBB_EDA/nSRcancer.png', NSRcancer, height = 6,width = 7, units = "cm")

sum(!is.na(traits$nSRcancer))
# [1] 39910

height <- traits %>%
  select(Sex, `Standing height`) %>%
  na.omit() %>%
  ggplot(aes(x = Sex, y= `Standing height`)) +
  geom_violin(aes(fill = Sex)) +
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) + 
  geom_boxplot(width = 0.05, outlier.shape = NA) + 
  xlab('') + ylab('Height (cm)')
ggsave('./results/UKBB_EDA/height.pdf',  height, useDingbats=F, height = 8,width = 6, units = "cm")
ggsave('./results/UKBB_EDA/height.png',  height, height = 8,width = 6, units = "cm")

weight <- traits %>%
  select(Sex, Weight) %>%
  na.omit() %>%
  ggplot(aes(x = Sex, y= Weight)) +
  geom_violin(aes(fill = Sex)) +
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) + 
  geom_boxplot(width = 0.05, outlier.shape = NA) + 
  xlab('') + ylab('Weight (kg)')
ggsave('./results/UKBB_EDA/weight.pdf',  weight, useDingbats=F, height = 8,width = 6, units = "cm")
ggsave('./results/UKBB_EDA/weight.png',  weight, height = 8,width = 6, units = "cm")

BMI <- traits %>%
  select(Sex, BMI) %>%
  na.omit() %>%
  ggplot(aes(x = Sex, y= BMI)) +
  geom_violin(aes(fill = Sex)) +
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) + 
  geom_boxplot(width = 0.05, outlier.shape = NA) + 
  xlab('') + ylab('BMI')
ggsave('./results/UKBB_EDA/BMI.pdf',  BMI, useDingbats=F, height = 8,width = 6, units = "cm")
ggsave('./results/UKBB_EDA/BMI.png',  BMI, height = 8,width = 6, units = "cm")

menarche <- traits %>%
  select(Sex, `Age when periods started (menarche)`) %>%
  na.omit() %>%
  ggplot(aes(x= `Age when periods started (menarche)`)) +
  geom_bar(aes(fill = Sex)) + 
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) + 
  scale_y_continuous(labels = scales::comma) +
  ylab('Number of participants') + xlab('Age when periods started') + guides(fill=F)
ggsave('./results/UKBB_EDA/menarche.pdf', menarche  , useDingbats=F, height = 6,width = 6, units = "cm")
ggsave('./results/UKBB_EDA/menarche.png', menarche  , height = 6, width = 6, units = "cm")

traits %>%
  select( `Age when periods started (menarche)`) %>%
  summary()

# Age when periods started (menarche)
# Min.   : 5.00                      
# 1st Qu.:12.00                      
# Median :13.00                      
# Mean   :12.97                      
# 3rd Qu.:14.00                      
# Max.   :25.00                      
# NA's   :229996 

menopause <- traits %>%
  select(Sex, `Age at menopause (last menstrual period)`) %>%
  na.omit() %>%
  ggplot(aes(x= `Age at menopause (last menstrual period)`)) +
  geom_bar(aes(fill = Sex)) + 
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) + 
  scale_y_continuous(labels = scales::comma) +
  ylab('Number of participants') + xlab('Age at menopause') + guides(fill=F)
ggsave('./results/UKBB_EDA/menopause.pdf', menopause  , useDingbats=F, height = 6,width = 6, units = "cm")
ggsave('./results/UKBB_EDA/menopause.png', menopause  , height = 6, width = 6, units = "cm")

traits %>%
  select(`Age at menopause (last menstrual period)`) %>%
  summary()
# Age at menopause (last menstrual period)
# Min.   :18.0                            
# 1st Qu.:48.0                            
# Median :50.0                            
# Mean   :49.7                            
# 3rd Qu.:53.0                            
# Max.   :68.0                            
# NA's   :335673
deathsinfamily = traits %>%
  select(Sex, `Non-accidental death in close genetic family`) %>%
  na.omit() %>%
  ggplot(aes(fill = Sex, x = as.character(`Non-accidental death in close genetic family`))) +
  geom_bar(position = 'dodge') + 
  scale_fill_manual(values = sexcolors[unique(traits$Sex)]) +
  scale_x_discrete(labels = c('No','Yes')[na.omit(unique(traits$`Non-accidental death in close genetic family`)) + 1]) +
  scale_y_continuous(labels = scales::unit_format('',1/1000)) +
  ylab('Number of participants\n(in thousands)') +
  xlab('') +
  coord_flip() + ggtitle('Non-accidental death\nin family')
ggsave('./results/UKBB_EDA/deathsinfamily.pdf', deathsinfamily, useDingbats=F, height = 6,width = 8, units = "cm")
ggsave('./results/UKBB_EDA/deathsinfamily.png', deathsinfamily, height = 6,width = 8, units = "cm")

p1 <- ggarrange(sexdistrib, 
          agedistrib_bysex,
          Age_at_death_bySex, 
          height, weight, BMI,
          legend = 'bottom',
          nrow = 2, ncol=3, labels = 'auto', common.legend = T, align = 'hv')
ggsave('./results/UKBB_EDA/figure1.pdf', p1, useDingbats=F, height = 13,width = 18, units = "cm")
ggsave('./results/UKBB_EDA/figure1.png', p1, height = 13,width = 18, units = "cm")

p2 <- ggarrange(healthratePlot, healthsatisfactionPlot, smokingPlot, alcoholPlot,
          facialagingPlot, deathsinfamily,common.legend = T, legend = 'bottom',
          ncol = 2, nrow = 3, align = 'hv', labels = 'auto')
ggsave('./results/UKBB_EDA/figure2.pdf', p2, useDingbats=F, height = 16,width = 18, units = "cm")
ggsave('./results/UKBB_EDA/figure2.png', p2, height = 16,width = 18, units = "cm")

p3 <- ggarrange(parentAgeatDeath, menarche, menopause, ncol = 3, nrow = 1, align = 'h', labels = 'auto')
ggsave('./results/UKBB_EDA/figure3.pdf', p3, useDingbats=F, height = 6,width = 18, units = "cm")
ggsave('./results/UKBB_EDA/figure3.png', p3, height =6,width = 18, units = "cm")

p4 <- ggarrange(numMed_Op, NSRcancer, nrow=1,ncol=2,widths = c(2,1.2), common.legend = T, legend = 'bottom', labels = c('a','b'))

ggsave('./results/UKBB_EDA/figure4.pdf', p4, useDingbats=F, height = 7,width = 18, units = "cm")
ggsave('./results/UKBB_EDA/figure4.png', p4, height =7,width = 18, units = "cm")

cotraits <- traits[,c(2:3, 5:6, 10, 12, 15, 16, 18:24, 26:30)] %>%
  mutate(Sex=as.numeric(as.factor(Sex))) %>%
  rename("Number of non-cancer disease" = nSRdiseaseProp,
         "Number of cancers" = nSRcancer) %>%
  mutate(`Facial ageing` = -1 * ifelse(`Facial ageing` == 2, 4, `Facial ageing`),
         `Overall health rating` = -`Overall health rating`,
         `Health satisfaction` = -`Health satisfaction`) 
cox=cor(cotraits,use='pairwise', method = 's')
pheatmap::pheatmap(cox, cellwidth = 22, cellheight = 22, display_numbers = T,
                   breaks = seq(-1,1,length.out = 20),
                   color = colorRampPalette(c('#2166AC','gray90','#B2182B'))(19),
                   filename = './results/UKBB_EDA/correlations.pdf',
                   number_format = '%.2f')
pheatmap::pheatmap(cox, cellwidth = 22, cellheight = 22, display_numbers = T,
                   breaks = seq(-1,1,length.out = 20),
                   color = colorRampPalette(c('#2166AC','gray90','#B2182B'))(19),
                   filename = './results/UKBB_EDA/correlations.png',
                   number_format = '%.2f')
cox_trimmed=cox
cox_trimmed[abs(cox)<0.2]=0
pheatmap::pheatmap(cox_trimmed, cellwidth = 22, cellheight = 22, display_numbers = T,
                   breaks = seq(-1,1,length.out = 20),
                   color = colorRampPalette(c('#2166AC','gray90','#B2182B'))(19),
                   filename = './results/UKBB_EDA/correlations_trimmed_below02.pdf',
                   number_format = '%.2f')
pheatmap::pheatmap(cox_trimmed, cellwidth = 22, cellheight = 22, display_numbers = T,
                   breaks = seq(-1,1,length.out = 20),
                   color = colorRampPalette(c('#2166AC','gray90','#B2182B'))(19),
                   filename = './results/UKBB_EDA/correlations_trimmed_below02.png',
                   number_format = '%.2f')
