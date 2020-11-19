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
    select(SNP,cluster,disease,UKBB_RAF,ALL_RAF,AFR_RAF,EAS_RAF,AMR_RAF,SAS_RAF,EUR_RAF) %>% 
    gather(key = 'pop', value ='raf',-SNP,-cluster,-disease) %>%
    mutate(pop = factor(gsub('_RAF','',pop),levels = c('UKBB','ALL','AFR','AMR','EAS','EUR','SAS'))) %>%
    na.omit()
  if(nrow(xx)>0){
    xx %>%
      group_by(cluster,disease, pop) %>% 
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

xx = readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated_selected.rds')
ld_onedis = xx %>%
  group_by(Disease) %>%
  summarise(n = length(unique(eid))) %>%
  rename(disease = Disease) %>%
  right_join(ld_onedis)

p1 = ld_onedis %>%
  mutate(maf = ifelse(RAF_median>0.5,1-RAF_median,RAF_median)) %>%
  group_by(disease,cluster,pop,n) %>%
  summarise(MAF = median(maf)) %>%
  filter(pop == 'UKBB') %>%
  ggplot(aes(x = n, y= MAF, color = cluster)) +
  geom_smooth(method = 'lm', alpha= 0.1)+
  geom_point(size =3) +
  scale_color_manual(values =ageonsetcolors) +
  xlab('Number of cases') +
  ylab('Median MAF') +
  guides(color = guide_legend('Age of onset cluster'))+
  scale_x_log10(labels = scales::comma)

ggsave('./results/evoAnalysis/MAF_cases_onedis.pdf', p1, units = 'cm', width = 16.8, height = 10, useDingbats =F)
ggsave('./results/evoAnalysis/MAF_cases_onedis.png', p1, units = 'cm', width = 16.8, height = 10)

######

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
    select(SNP,cluster,disease,UKBB_RAF,ALL_RAF,AFR_RAF,EAS_RAF,AMR_RAF,SAS_RAF,EUR_RAF) %>% 
    gather(key = 'pop', value ='raf',-SNP,-cluster,-disease) %>%
    mutate(pop = factor(gsub('_RAF','',pop),levels = c('UKBB','ALL','AFR','AMR','EAS','EUR','SAS'))) %>%
    na.omit()
  if(nrow(xx)>0){
    xx %>%
      group_by(disease,cluster,pop) %>% 
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



ld_onecl = xx %>%
  group_by(Disease) %>%
  summarise(n = length(unique(eid))) %>%
  rename(disease = Disease) %>%
  right_join(ld_onecl)

p2 = ld_onecl %>%
  mutate(maf = ifelse(RAF_median>0.5,1-RAF_median,RAF_median)) %>%
  group_by(disease,cluster,pop,n) %>%
  summarise(MAF = median(maf)) %>%
  filter(pop == 'UKBB') %>%
  ggplot(aes(x = n, y= MAF, color = cluster)) +
  geom_smooth(method = 'lm', alpha = 0.1) +
  geom_point(size =3) +
  scale_color_manual(values =ageonsetcolors) +
  xlab('Number of cases') +
  ylab('Median MAF') +
  guides(color = guide_legend('Age of onset cluster')) +
  scale_x_log10(labels = scales::comma) 

ggsave('./results/evoAnalysis/MAF_cases_onecl.pdf',p2, units = 'cm', width = 16.8, height = 10, useDingbats =F)
ggsave('./results/evoAnalysis/MAF_cases_onecl.png',p2, units = 'cm', width = 16.8, height = 10)
