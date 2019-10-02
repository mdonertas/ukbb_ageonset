source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('../ukbb_ageonset/results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)
disCoding = disCoding[as.character(disIDs)]
signifSNPs=readRDS('./data/processed/evoAnalysis/UKBBRAF_Pleiotropy.rds')
LDblocks = read_tsv('./data/raw/LDblocks1000G_EUR.bed')%>%
  mutate(CHR=as.numeric(gsub('chr','',chr))) %>%
  select(-chr)
signifSNPs = signifSNPs %>%
  filter(cluster !=4) 
signifSNPs = signifSNPs %>%
  mutate(gr_agonist = grepl('agonist',cl1_cl1) |grepl('agonist',cl2_cl2)|grepl('agonist',cl3_cl3),
         gr_antagonist = grepl('antagonist',cl1_cl1) |grepl('antagonist',cl2_cl2)|grepl('antagonist',cl3_cl3),
         ac_agonist = grepl('agonist',cl1_cl2) | grepl('agonist',cl1_cl3) | grepl('agonist',cl2_cl3),
         ac_antagonist = grepl('antagonist',cl1_cl2) | grepl('antagonist',cl1_cl3) | grepl('antagonist',cl2_cl3))
LDblocks$block=1:nrow(LDblocks)
ld_onedis = apply(LDblocks,1,function(x){
  st=x['start']
  en=x['stop']
  ch=x['CHR']
  xx = signifSNPs %>%
    filter(numDis == 1) %>%
    filter(CHR==ch & BP<=en & BP>=st) %>%
    select(-disID,-BETA,-P_BOLT_LMM_INF) %>%
    unique() %>%
    mutate(cluster = factor(cluster)) %>% 
    select(SNP,cluster,UKBB_RAF,ALL_RAF,AFR_RAF,EAS_RAF,AMR_RAF,SAS_RAF,EUR_RAF) %>% 
    gather(key = 'pop', value ='raf',-SNP,-cluster) %>%
    mutate(pop = factor(gsub('_RAF','',pop),levels = c('UKBB','ALL','AFR','AMR','EAS','EUR','SAS'))) %>%
    na.omit()
  if(nrow(xx)>0){
    xx %>%
      group_by(cluster,pop) %>% 
      summarise(RAF_mean = mean(raf,na.rm=T),
                RAF_min = min(raf,na.rm=T),
                RAF_max = max(raf,na.rm=T),
                RAF_sd = sd(raf,na.rm=T),
                RAF_q1 = quantile(raf,probs = 0.25, na.rm=T),
                RAF_q3 = quantile(raf,probs = 0.75, na.rm=T),
                RAF_median = quantile(raf,probs = 0.75, na.rm=T)) %>%
      mutate(LDblock=x['block'])
  } else {NULL}
}) 
ld_onedis=ld_onedis[!sapply(ld_onedis,is.null)]
ld_onedis = reshape2::melt(ld_onedis,id.vars = colnames(ld_onedis[[1]])) 
ukkb_ld_onedis = ld_onedis %>% 
  mutate(cluster=as.factor(cluster)) %>%
  filter(pop == 'UKBB') %>%
  ggplot(aes(y = RAF_median, x= cluster, fill = cluster)) +
  # facet_wrap(~pop)+
  geom_violin() +
  geom_boxplot(width = 0.25, fill ='white', outlier.shape = NA) +
  geom_sina(size=0.3,alpha=0.5)+
  scale_fill_manual(values = ageonsetcolors) +
  stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label",fill='white')  +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency') +
  stat_compare_means(comparisons = list(c('1','2'),c('2','3'),c('1','3')), label = 'p.signif')

onekg_ld_onedis = ld_onedis %>% 
  mutate(cluster=as.factor(cluster)) %>%
  filter(pop != 'UKBB') %>%
  ggplot(aes(y = RAF_median, x= cluster, fill = cluster)) +
  facet_wrap(~pop)+
  geom_violin() +
  geom_boxplot(width = 0.25, fill ='white', outlier.shape = NA) +
  geom_sina(size=0.3,alpha=0.5)+
  scale_fill_manual(values = ageonsetcolors) +
  stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label",fill='white')  +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency') +
  stat_compare_means(comparisons = list(c('1','2'),c('2','3'),c('1','3')), label = 'p.signif')

all_onedisease = ggarrange(ukkb_ld_onedis+ggtitle('UK Biobank'), onekg_ld_onedis, labels = 'auto', nrow=1,ncol=2, widths = c(1,2))

ggsave('./results/evoAnalysis/RAF_all_onedisease_ldblock.pdf',all_onedisease, units = 'cm', width = 24, height = 12, useDingbats = F)
ggsave('./results/evoAnalysis/RAF_all_onedisease_ldblock.png',all_onedisease, units = 'cm', width = 24, height = 12)

xx = ld_onedis %>%
  filter(pop=='UKBB')%>%
  select(RAF_median, cluster) %>%
  unique() 

xx %>% group_by(cluster) %>% summarise(n = length(RAF_median))
# # A tibble: 3 x 2
# cluster     n
# <fct>   <int>
# 1 1       485
# 2 2       232
# 3 3       169

cl1perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==1)%>%na.omit())$RAF_median,100)))
cl2perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==2)%>%na.omit())$RAF_median,100)))
cl3perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==3)%>%na.omit())$RAF_median,100)))

onedisease_ukbb_sampling = data.frame(cluster = factor(rep(1:3,each=1000),levels = 1:3), median = c(cl1perm,cl2perm,cl3perm)) %>%
  ggplot(aes(x = cluster, y= median,color=cluster)) +
  geom_sina(alpha=0.5,size=0.5) + 
  stat_summary(fun.y = 'mean', size = 0.5, geom = 'hline',
               aes(yintercept = ..y..), linetype = 'dashed') +
  xlab('Age of onset cluster') + ylab('Median RAF (for 100 Blocks)') +
  scale_color_manual(values = ageonsetcolors) +
  guides(color=F)

ld_onecl = apply(LDblocks,1,function(x){
  st=x['start']
  en=x['stop']
  ch=x['CHR']
  xx = signifSNPs %>%
    filter(numCl == 1) %>%
    filter(CHR==ch & BP<=en & BP>=st) %>%
    select(-disID,-BETA,-P_BOLT_LMM_INF) %>%
    unique() %>%
    mutate(cluster = factor(cluster)) %>% 
    select(SNP,cluster,UKBB_RAF,ALL_RAF,AFR_RAF,EAS_RAF,AMR_RAF,SAS_RAF,EUR_RAF) %>% 
    gather(key = 'pop', value ='raf',-SNP,-cluster) %>%
    mutate(pop = factor(gsub('_RAF','',pop),levels = c('UKBB','ALL','AFR','AMR','EAS','EUR','SAS'))) %>%
    na.omit()
  if(nrow(xx)>0){
    xx %>%
      group_by(cluster,pop) %>% 
      summarise(RAF_mean = mean(raf,na.rm=T),
                RAF_min = min(raf,na.rm=T),
                RAF_max = max(raf,na.rm=T),
                RAF_sd = sd(raf,na.rm=T),
                RAF_q1 = quantile(raf,probs = 0.25, na.rm=T),
                RAF_q3 = quantile(raf,probs = 0.75, na.rm=T),
                RAF_median = quantile(raf,probs = 0.75, na.rm=T)) %>%
      mutate(LDblock=x['block'])
  } else {NULL}
}) 
ld_onecl=ld_onecl[!sapply(ld_onecl,is.null)]
ld_onecl = reshape2::melt(ld_onecl,id.vars = colnames(ld_onecl[[1]])) 
ukkb_ld_onecl = ld_onecl %>% 
  mutate(cluster=as.factor(cluster)) %>%
  filter(pop == 'UKBB') %>%
  ggplot(aes(y = RAF_median, x= cluster, fill = cluster)) +
  # facet_wrap(~pop)+
  geom_violin() +
  geom_boxplot(width = 0.25, fill ='white', outlier.shape = NA) +
  geom_sina(size=0.3,alpha=0.5)+
  scale_fill_manual(values = ageonsetcolors) +
  stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label",fill='white')  +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency') +
  stat_compare_means(comparisons = list(c('1','2'),c('2','3'),c('1','3')), label = 'p.signif')

onekg_ld_onecl = ld_onecl %>% 
  mutate(cluster=as.factor(cluster)) %>%
  filter(pop != 'UKBB') %>%
  ggplot(aes(y = RAF_median, x= cluster, fill = cluster)) +
  facet_wrap(~pop)+
  geom_violin() +
  geom_boxplot(width = 0.25, fill ='white', outlier.shape = NA) +
  geom_sina(size=0.3,alpha=0.5)+
  scale_fill_manual(values = ageonsetcolors) +
  stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label",fill='white')  +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency') +
  stat_compare_means(comparisons = list(c('1','2'),c('2','3'),c('1','3')), label = 'p.signif')

all_onecluster = ggarrange(ukkb_ld_onecl+ggtitle('UK Biobank'), onekg_ld_onecl, labels = 'auto', nrow=1,ncol=2, widths = c(1,2))

ggsave('./results/evoAnalysis/RAF_all_onecluster_ldblock.pdf',all_onecluster, units = 'cm', width = 24, height = 12, useDingbats = F)
ggsave('./results/evoAnalysis/RAF_all_onecluster_ldblock.png',all_onecluster, units = 'cm', width = 24, height = 12)

xx = ld_onecl %>%
  filter(pop=='UKBB')%>%
  select(RAF_median, cluster) %>%
  unique() 

xx %>% group_by(cluster) %>% summarise(n = length(RAF_median))
# # A tibble: 3 x 2
# cluster     n
# <fct>   <int>
# 1 1       489
# 2 2       248
# 3 3       174

cl1perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==1)%>%na.omit())$RAF_median,100)))
cl2perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==2)%>%na.omit())$RAF_median,100)))
cl3perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==3)%>%na.omit())$RAF_median,100)))

onecluster_ukbb_sampling = data.frame(cluster = factor(rep(1:3,each=1000),levels = 1:3), median = c(cl1perm,cl2perm,cl3perm)) %>%
  ggplot(aes(x = cluster, y= median,color=cluster)) +
  geom_sina(alpha=0.5,size=0.5) + 
  stat_summary(fun.y = 'mean', size = 0.5, geom = 'hline',
               aes(yintercept = ..y..), linetype = 'dashed') +
  xlab('Age of onset cluster') + ylab('Median RAF (for 100 Blocks)') +
  scale_color_manual(values = ageonsetcolors) +
  guides(color=F)


####

ld_cl1cl2_anta = apply(LDblocks,1,function(x){
  st=x['start']
  en=x['stop']
  ch=x['CHR']
  xx = signifSNPs %>%
    filter(cl1_cl2 == ' antagonist' & !grepl('agonist',cl1_cl3) & !grepl('agonist', cl2_cl3)) %>%
    filter(CHR==ch & BP<=en & BP>=st) %>%
    select(-disID,-BETA,-P_BOLT_LMM_INF) %>%
    unique() %>%
    mutate(cluster = factor(cluster)) %>% 
    select(SNP,cluster,UKBB_RAF,ALL_RAF,AFR_RAF,EAS_RAF,AMR_RAF,SAS_RAF,EUR_RAF) %>% 
    gather(key = 'pop', value ='raf',-SNP,-cluster) %>%
    mutate(pop = factor(gsub('_RAF','',pop),levels = c('UKBB','ALL','AFR','AMR','EAS','EUR','SAS'))) %>%
    na.omit()
  if(nrow(xx)>0){
    xx %>%
      group_by(cluster,pop) %>% 
      summarise(RAF_mean = mean(raf,na.rm=T),
                RAF_min = min(raf,na.rm=T),
                RAF_max = max(raf,na.rm=T),
                RAF_sd = sd(raf,na.rm=T),
                RAF_q1 = quantile(raf,probs = 0.25, na.rm=T),
                RAF_q3 = quantile(raf,probs = 0.75, na.rm=T),
                RAF_median = quantile(raf,probs = 0.75, na.rm=T)) %>%
      mutate(LDblock=x['block'])
  } else {NULL}
}) 
ld_cl1cl2_anta=ld_cl1cl2_anta[!sapply(ld_cl1cl2_anta,is.null)]
ld_cl1cl2_anta = reshape2::melt(ld_cl1cl2_anta,id.vars = colnames(ld_cl1cl2_anta[[1]])) 
ukkb_ld_anta= ld_cl1cl2_anta %>% 
  mutate(cluster=as.factor(cluster)) %>%
  filter(pop == 'UKBB') %>%
  ggplot(aes(y = RAF_median, x= cluster, fill = cluster)) +
  geom_violin()+
  geom_sina(size=1,alpha=1)+
  geom_boxplot(width=0.1,fill='white') +
  scale_fill_manual(values = ageonsetcolors) +
  stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label",fill='white')  +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency')+
  ylim(0,1) +
  stat_compare_means(label = 'p.format',label.y=0.01,label.x = 0.5,hjust=0)

onekg_ld_anta = ld_cl1cl2_anta %>% 
  mutate(cluster=as.factor(cluster)) %>%
  filter(pop != 'UKBB') %>%
  ggplot(aes(y = RAF_median, x= cluster, fill = cluster)) +
  facet_wrap(~pop) +
  geom_violin()+
  geom_sina(size=1,alpha=1)+
  geom_boxplot(width=0.1,fill='white') +
  scale_fill_manual(values = ageonsetcolors) +
  stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label",fill='white')  +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency')+
  ylim(0,1) +
  stat_compare_means(label = 'p.format',label.y=0.01,label.x = 0.5,hjust=0)

all_anta = ggarrange(ukkb_ld_anta+ggtitle('UK Biobank'), onekg_ld_anta, labels = 'auto', nrow=1,ncol=2, widths = c(1,2))

ggsave('./results/evoAnalysis/RAF_all_antagonist_ldblock.pdf',all_anta, units = 'cm', width = 24, height = 12, useDingbats = F)
ggsave('./results/evoAnalysis/RAF_all_antagonist_ldblock.png',all_anta, units = 'cm', width = 24, height = 12)

xx = ld_cl1cl2_anta %>%
  filter(pop=='UKBB')%>%
  select(RAF_median, cluster) %>%
  unique() 

xx %>% group_by(cluster) %>% summarise(n = length(RAF_median))
# # A tibble: 2 x 2
# cluster     n
# <fct>   <int>
# 1 1          20
# 2 2          20

cl1perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==1)%>%na.omit())$RAF_median,10)))
cl2perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==2)%>%na.omit())$RAF_median,10)))

antagonist_ukbb_sampling = data.frame(cluster = factor(rep(1:2,each=1000),levels = 1:2), median = c(cl1perm,cl2perm)) %>%
  ggplot(aes(x = cluster, y= median,color=cluster)) +
  geom_sina(alpha=1,size=1) + 
  stat_summary(fun.y = 'mean', size = 0.5, geom = 'hline',
               aes(yintercept = ..y..), linetype = 'dashed') +
  xlab('Age of onset cluster') + ylab('Median RAF (for 100 Blocks)') +
  scale_color_manual(values = ageonsetcolors) +
  guides(color=F)

allukbb=ggarrange(ukkb_ld_onedis+ggtitle('SNPs associated with\none disease'), ukkb_ld_onecl+ggtitle('SNPs associated with\none cluster'), ukkb_ld_anta+ggtitle('SNPs Antagonistic Between\nCluster1 and Cluster3'), ncol=3,labels = 'auto',nrow=1) 

allukbb_sampling=ggarrange(onedisease_ukbb_sampling+ggtitle('SNPs associated with\none disease'), onecluster_ukbb_sampling+ggtitle('SNPs associated with\none cluster'), antagonist_ukbb_sampling+ggtitle('SNPs Antagonistic Between\nCluster1 and Cluster3'), ncol=3,labels = 'auto',nrow=1) 

ggsave('./results/evoAnalysis/RAF_ukbb_ldblock.pdf',allukbb,units = 'cm',width = 20 ,height = 7,useDingbats =F)
ggsave('./results/evoAnalysis/RAF_ukbb_ldblock.png',allukbb,units = 'cm',width = 20,height = 7)

ggsave('./results/evoAnalysis/RAF_ukbb_sampling_ldblock.pdf',allukbb_sampling,units = 'cm',width = 20 ,height = 7,useDingbats =F)
ggsave('./results/evoAnalysis/RAF_ukbb_sampling_ldblock.png',allukbb_sampling,units = 'cm',width = 20,height = 7)

onekg_ld_anta=ld_cl1cl2_anta %>% 
  mutate(cluster=as.factor(cluster)) %>%
  filter(pop != 'UKBB') %>%
  # filter(cluster==1) %>%
  ggplot(aes(y = RAF_median, x= cluster, fill = cluster)) +
  facet_wrap(~pop,ncol = 6) +
  geom_violin()+
  geom_sina(size=1,alpha=1)+
  geom_boxplot(width=0.1,fill='white') +
  scale_fill_manual(values = ageonsetcolors) +
  stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label",fill='white')  +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency')+
  ylim(0,1) +
  stat_compare_means(label = 'p.format',label.y=0.01,label.x = 0.5,hjust=0)

fig4=ggarrange(allukbb,onekg_ld_anta,nrow=2,ncol=1,labels=c('','d'),widths=c(1.1,1))

ggsave('./results/evoAnalysis/fig4.pdf',fig4,units = 'cm',width = 22 ,height = 15,useDingbats =F)
ggsave('./results/evoAnalysis/fig4.png',fig4,units = 'cm',width = 22,height = 15)

