source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('../ukbb_ageonset/results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)
disCoding = disCoding[as.character(disIDs)]
# signifSNPs <- lapply(paste('../ukbb_ageonset/data/processed/caseControl/a',disIDs,'/signif_gwasRes_proxyGenes.rds',sep=''),function(x){
#   x=readRDS(x)
#   x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend) & SNP == 'rs4480762')
#   select(x,SNP,CHR,BP,Ref,Alt,A1FREQ,BETA,P_BOLT_LMM_INF) %>% rename(UKBB_RefAF = A1FREQ) %>%unique()
# })
# names(signifSNPs)=disIDs
# clusters = readRDS('../ukbb_ageonset/data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')
# signifSNPs <- reshape2::melt(signifSNPs, id.vars = colnames(signifSNPs[[1]])) %>%
#   mutate(RA = ifelse(BETA>0, Ref, Alt)) %>%
#   rename(disID = L1) %>%
#   left_join(data.frame(disease = unname(disCoding), disID = names(disCoding), cluster = clusters$clustering[unname(disCoding)])) %>%
#   unique()
# 
# Tgenome = readRDS('../melike/projects/shared_data/1kG/data/processed/1kgMergedUKBB.rds')
# Tgenome=Tgenome[,colnames(Tgenome)[c(6,1,2,4,5,3,7,9,13:18)]]
# Tgenome$AF = as.numeric(as.character(Tgenome$AF))
# Tgenome$EAS_AF = as.numeric(as.character(Tgenome$EAS_AF))
# Tgenome$AMR_AF = as.numeric(as.character(Tgenome$AMR_AF))
# Tgenome$AFR_AF = as.numeric(as.character(Tgenome$AFR_AF))
# Tgenome$EUR_AF = as.numeric(as.character(Tgenome$EUR_AF))
# Tgenome$SAS_AF = as.numeric(as.character(Tgenome$SAS_AF))
# Tgenome$AA = as.character(Tgenome$AA)
# aa = sapply(strsplit(Tgenome$AA,'[|]'),function(x)x[1])
# Tgenome$AA = aa
# rm(aa)
# signifSNPs = signifSNPs %>%
#   mutate(CHR = as.numeric(CHR)) %>%
#   left_join(Tgenome) %>% select(-A1FREQ) %>%
#   mutate(UKBB_RAF = ifelse(RA == Ref, UKBB_RefAF, 1-UKBB_RefAF),
#          ALL_RAF = ifelse(RA == Ref, 1-AF, AF),
#          EAS_RAF = ifelse(RA == Ref, 1-EAS_AF, EAS_AF),
#          AMR_RAF = ifelse(RA == Ref, 1-AMR_AF, AMR_AF),
#          AFR_RAF = ifelse(RA == Ref, 1-AFR_AF, AFR_AF),
#          EUR_RAF = ifelse(RA == Ref, 1-EUR_AF, EUR_AF),
#          SAS_RAF = ifelse(RA == Ref, 1-SAS_AF, SAS_AF)) %>% 
#   mutate(AA = ifelse(toupper(AA)%in%c(Ref,Alt),AA,NA)) %>% 
#   mutate(AA_LC = toupper(AA)) %>%
#   mutate(Ancestral = AA == RA,
#          Ancestral_LC = AA_LC == RA)
# rm(Tgenome)
# 
# signifSNPs = signifSNPs %>%
#   select(SNP,RA,disID,cluster) %>%
#   group_by(SNP) %>%
#   summarise(numDis = length(unique(disID)),
#             numCl = length(unique(cluster))) %>%
#   right_join(signifSNPs) 

# xx = signifSNPs %>%
#   filter(numDis>1) %>%
#   select(SNP,RA,disID, cluster) %>%
#   group_by(SNP,RA) %>%
#   summarise(diseases = list(sort(unique(disID))),
#             clusters = list(sort(unique(cluster))),
#             cl1_num = sum(cluster==1),
#             cl2_num = sum(cluster==2),
#             cl3_num = sum(cluster==3),
#             cl4_num = sum(cluster==4))
# 
# xx = t(sapply(unique(xx$SNP),function(snp){
#   xx2=filter(xx,SNP==snp)
#   xx3 = apply(xx2[,5:8],2,function(x){
#     apply(xx2[,5:8],2,function(y){
#       agonist = any(((x>0) + (y>0)) > 1)
#       antagonist = ((!all(x==0)) & (!all(y==0)) & (any((x>0) & (y==0)) | any((x==0) & y>0)))
#       data.frame(agonist = agonist, antagonist = antagonist)
#     })
#   }) %>% reshape2::melt(id.vars = c('agonist','antagonist')) %>%
#     mutate(cluster1 = gsub('_num','',L1),
#            cluster2 = gsub('_num','',L2)) %>% select(cluster1,cluster2,agonist,antagonist) %>% filter(cluster1!=cluster2)
#   xx3 = rbind(xx3,data.frame(cluster1 = paste('cl',1:4,sep=''),
#              cluster2 = paste('cl',1:4,sep=''),
#              agonist =colSums(xx2[,5:8]>1)>0,
#              antagonist = colSums(xx2[,5:8])!=apply(xx2[,5:8],2,max)))
#   rownames(xx3)=NULL
#   xx3 %>%
#     mutate(SNP=snp,clr = paste(cluster1,cluster2,sep='-')) %>%
#     mutate(agonist = c(NA,'agonist')[agonist+1]) %>%
#     mutate(antagonist = c(NA,'antagonist')[antagonist+1]) %>%
#     mutate(res = gsub('NA','',paste(agonist,antagonist,sep=' '))) %>%
#     select(SNP,clr,res) %>%
#     spread(clr,res)
#   })) %>% as.data.frame()
# rownames(xx) = NULL  

# xx = readRDS('./data/processed/evoAnalysis.rds')
# all(unlist(xx$`cl1-cl2`) == unlist(xx$`cl2-cl1`))
# xx$`cl2-cl1`=NULL
# all(unlist(xx$`cl1-cl3`) == unlist(xx$`cl3-cl1`))
# xx$`cl3-cl1`=NULL
# all(unlist(xx$`cl2-cl3`) == unlist(xx$`cl3-cl2`))
# xx$`cl3-cl2`=NULL
# all(unlist(xx$`cl1-cl4`) == unlist(xx$`cl4-cl1`))
# xx$`cl4-cl1`=NULL
# all(unlist(xx$`cl2-cl4`) == unlist(xx$`cl4-cl2`))
# xx$`cl4-cl2`=NULL
# all(unlist(xx$`cl3-cl4`) == unlist(xx$`cl4-cl3`))
# xx$`cl4-cl3`=NULL
# 
# signifSNPs = data.frame(SNP = unlist(xx$SNP),
#            cl1_cl1 = unlist(xx$`cl1-cl1`),
#            cl1_cl2 = unlist(xx$`cl1-cl2`),
#            cl1_cl3 = unlist(xx$`cl1-cl3`),
#            cl1_cl4 = unlist(xx$`cl1-cl4`),
#            cl2_cl2 = unlist(xx$`cl2-cl2`),
#            cl2_cl3 = unlist(xx$`cl2-cl3`),
#            cl2_cl4 = unlist(xx$`cl2-cl4`),
#            cl3_cl3 = unlist(xx$`cl3-cl3`),
#            cl3_cl4 = unlist(xx$`cl3-cl4`),
#            cl4_cl4 = unlist(xx$`cl4-cl4`)) %>%
#   right_join(signifSNPs)
# 
# saveRDS(signifSNPs,'./data/processed/evoAnalysis/UKBBRAF_Pleiotropy.rds')
signifSNPs=readRDS('./data/processed/evoAnalysis/UKBBRAF_Pleiotropy.rds')
library(ggthemes)
numsnps = reshape2::melt((table(unique(select(signifSNPs,SNP,numDis,numCl))%>%select(-SNP)))) %>%
  mutate(value = ifelse(value ==0, NA, value)) %>%
  mutate(numDis = paste('# Diseases = ', numDis, sep='')) %>%
  ggplot(aes(x = as.factor(numCl), fill = as.factor(numCl), y = value)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_label(aes(label = scales::comma(value)), fill = 'white', y=0.5) +
  facet_wrap(~numDis) + 
  scale_y_log10(label = scales::comma) +
  guides(fill = F) +
  xlab('Number of Clusters') + ylab('Number of SNPs (in log10 scale)') +
  scale_fill_wsj()

# numsnps = reshape2::melt((table(unique(select(signifSNPs,SNP,numDis,numCl))%>%select(-SNP)))) %>%
#   ggplot(aes(x = factor(numCl,levels = 1:3), y= factor(numDis,levels = 9:1))) +
#   geom_tile(aes(fill = log10(value))) +
#   scale_fill_gradient(low = 'gray90',high = 'gray25',na.value = 'white') +
#   geom_label(aes(label = scales::comma(value))) +
#   xlab('Number of Clusters') + ylab('Number of Diseases') +
#   theme_pubr(base_size = 14) +
#   guides(fill = F) +
#   theme(panel.grid.major = element_blank(), 
#         axis.line = element_blank(),
#         axis.ticks = element_blank())

ggsave('./results/evoAnalysis/pleiotropy.pdf', numsnps, units = 'cm', width = 16, height = 12, useDingbats = F)
ggsave('./results/evoAnalysis/pleiotropy.png', numsnps, units = 'cm', width = 16, height = 12)

select(signifSNPs, SNP,cluster) %>%
  group_by(cluster) %>%
  summarise(num = length(unique(SNP)))
# cluster   num
# <int> <int>
# 1       1 59568
# 2       2 26458
# 3       3 17417
# 4       4     7

# filter out cluster 4 diseases
signifSNPs = signifSNPs %>%
  filter(cluster !=4) 

# define within and across group agonists and antagonists

signifSNPs = signifSNPs %>%
  mutate(gr_agonist = grepl('agonist',cl1_cl1) |grepl('agonist',cl2_cl2)|grepl('agonist',cl3_cl3),
         gr_antagonist = grepl('antagonist',cl1_cl1) |grepl('antagonist',cl2_cl2)|grepl('antagonist',cl3_cl3),
         ac_agonist = grepl('agonist',cl1_cl2) | grepl('agonist',cl1_cl3) | grepl('agonist',cl2_cl3),
         ac_antagonist = grepl('antagonist',cl1_cl2) | grepl('antagonist',cl1_cl3) | grepl('antagonist',cl2_cl3))

signifSNP_summary = signifSNPs %>%
  select(1,12:17,19:24,30:46)

onecluster_noantagonist_1kG = signifSNPs %>%
  filter((!gr_antagonist) & (!ac_antagonist)) %>%
  filter(numCl == 1) %>% 
  select(-disID,-disease,-BETA,-P_BOLT_LMM_INF) %>%
  unique() %>%
  mutate(cluster = factor(cluster)) %>% 
  select(SNP,cluster,numDis,ALL_RAF,AFR_RAF,EAS_RAF,AMR_RAF,SAS_RAF,EUR_RAF) %>% unique() %>%
  gather(key = 'pop', value ='raf',-SNP,-cluster,-numDis) %>%
  mutate(pop = factor(gsub('_RAF','',pop),levels = c('ALL','AFR','AMR','EAS','EUR','SAS'))) %>%
  # filter(cluster == 2) %>%
  ggplot(aes(y = raf, x= cluster, fill = cluster)) +
  # geom_boxplot() +
  geom_violin() +
  geom_boxplot(width = 0.1, fill ='white', outlier.shape = NA) +
  scale_fill_manual(values = ageonsetcolors) +
  stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label", fill = 'white')  +
  facet_wrap(~pop, ncol = 3, nrow = 2) +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency') +
  ggtitle('1000 Genomes')

onecluster_noantagonist_ukbb = signifSNPs %>%
  filter((!gr_antagonist) & (!ac_antagonist)) %>%
  filter(numCl == 1) %>% 
  select(-disID,-disease,-BETA,-P_BOLT_LMM_INF) %>%
  unique() %>%
  mutate(cluster = factor(cluster)) %>% 
  ggplot(aes(y = UKBB_RAF, x= cluster, fill = cluster)) +
  # geom_boxplot() +
  geom_violin() +
  geom_boxplot(width = 0.1, fill ='white', outlier.shape = NA) +
  scale_fill_manual(values = ageonsetcolors) +
  stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label", fill = 'white')  +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency') +
  ggtitle('UK Biobank')

all_onecluster_noantagonist = ggarrange(onecluster_noantagonist_ukbb,onecluster_noantagonist_1kG,ncol = 2, nrow = 1, widths = c(1,2), labels = 'auto')

xx = signifSNPs %>%
  filter((!gr_antagonist) & (!ac_antagonist)) %>%
  filter(numCl == 1) %>% 
  select(SNP, RA, UKBB_RAF, cluster) %>%
  unique() 

xx %>% group_by(cluster) %>% summarise(n = length(unique(SNP)))
# # A tibble: 3 x 2
# cluster     n
# <fct>   <int>
# 1 1       51098
# 2 2       18075
# 3 3       15772

cl1perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==1)%>%na.omit())$UKBB_RAF,5000)))
cl2perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==2)%>%na.omit())$UKBB_RAF,5000)))
cl3perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==3)%>%na.omit())$UKBB_RAF,5000)))

onecluster_ukbb_sampling = data.frame(cluster = factor(rep(1:3,each=1000),levels = 1:3), median = c(cl1perm,cl2perm,cl3perm)) %>%
  ggplot(aes(x = cluster, y= median)) +
  geom_sina(alpha=0.3,size=0.5) + 
  stat_summary(fun.y = 'mean', colour = "gray60", size = 0.5, geom = 'hline',
               aes(yintercept = ..y..), linetype = 'dashed') +
  xlab('Age of onset cluster') + ylab('Median RAF (for 5,000 SNPs)') 
  

onedisease_1kG = signifSNPs %>%
  filter((!gr_antagonist) & (!ac_antagonist)) %>%
  filter(numCl == 1) %>%
  filter(numDis == 1) %>%
  select(-disID,-disease,-BETA,-P_BOLT_LMM_INF) %>%
  unique() %>%
  mutate(cluster = factor(cluster)) %>% 
  select(SNP,cluster,numDis,ALL_RAF,AFR_RAF,EAS_RAF,AMR_RAF,SAS_RAF,EUR_RAF) %>% unique() %>%
  gather(key = 'pop', value ='raf',-SNP,-cluster,-numDis) %>%
  mutate(pop = factor(gsub('_RAF','',pop),levels = c('ALL','AFR','AMR','EAS','EUR','SAS'))) %>%
  # filter(cluster == 2) %>%
  ggplot(aes(y = raf, x= cluster, fill = cluster)) +
  # geom_boxplot() +
  geom_violin() +
  geom_boxplot(width = 0.1, fill ='white', outlier.shape = NA) +
  scale_fill_manual(values = ageonsetcolors) +
  stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label", fill = 'white')  +
  facet_wrap(~pop, ncol = 3, nrow = 2) +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency') +
  ggtitle('1000 Genomes')

onedisease_ukbb = signifSNPs %>%
  filter((!gr_antagonist) & (!ac_antagonist)) %>%
  filter(numCl == 1) %>%
  filter(numDis == 1) %>%
  select(-disID,-disease,-BETA,-P_BOLT_LMM_INF) %>%
  unique() %>%
  mutate(cluster = factor(cluster)) %>% 
  # filter(cluster == 2) %>%
  ggplot(aes(y = UKBB_RAF, x= cluster, fill = cluster)) +
  # geom_boxplot() +
  geom_violin() +
  geom_boxplot(width = 0.1, fill ='white', outlier.shape = NA) +
  scale_fill_manual(values = ageonsetcolors) +
  stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label", fill = 'white')  +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency') +
  ggtitle('UK Biobank')

all_onedisease = ggarrange(onedisease_ukbb,onedisease_1kG,ncol = 2, nrow = 1, widths = c(1,2), labels = 'auto')

xx = signifSNPs %>%
  filter((!gr_antagonist) & (!ac_antagonist)) %>%
  filter(numCl == 1) %>%
  filter(numDis == 1) %>%
  select(SNP, RA, UKBB_RAF, cluster) %>%
  unique() 

xx %>% group_by(cluster) %>% summarise(n = length(unique(SNP)))
# # A tibble: 3 x 2
# cluster     n
# <int> <int>
# 1       1 31641
# 2       2 10313
# 3       3  8059

cl1perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==1)%>%na.omit())$UKBB_RAF,5000)))
cl2perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==2)%>%na.omit())$UKBB_RAF,5000)))
cl3perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==3)%>%na.omit())$UKBB_RAF,5000)))

onedisease_ukbb_sampling = data.frame(cluster = factor(rep(1:3,each=1000),levels = 1:3), median = c(cl1perm,cl2perm,cl3perm)) %>%
  ggplot(aes(x = cluster, y= median)) +
  geom_sina(alpha=0.3,size=0.5) + 
  stat_summary(fun.y = 'mean', colour = "gray60", size = 0.5, geom = 'hline',
               aes(yintercept = ..y..), linetype = 'dashed') +
  xlab('Age of onset cluster') + ylab('Median RAF (for 5,000 SNPs)') 

cl13_antagonist_1kg = signifSNPs %>%
  filter(cl1_cl3 == ' antagonist' & !grepl('agonist',cl1_cl2)) %>%
  # filter(!gr_antagonist) %>%
  # filter((!gr_antagonist) & (ac_antagonist)) %>%
  # filter(numCl == 1) %>% 
  # filter(numDis == 1) %>%
  select(-disID,-disease,-BETA,-P_BOLT_LMM_INF) %>%
  unique() %>%
  mutate(cluster = factor(cluster)) %>% 
  select(SNP,cluster,numDis,ALL_RAF,AFR_RAF,EAS_RAF,AMR_RAF,SAS_RAF,EUR_RAF) %>% unique() %>%
  gather(key = 'pop', value ='raf',-SNP,-cluster,-numDis) %>%
  mutate(pop = factor(gsub('_RAF','',pop),levels = c('ALL','AFR','AMR','EAS','EUR','SAS'))) %>%
  # filter(cluster == 2) %>%
  ggplot(aes(y = raf, x= cluster, fill = cluster)) +
  # geom_boxplot() +
  geom_violin() +
  geom_boxplot(width = 0.1, fill ='white') +
  scale_fill_manual(values = ageonsetcolors) +
  stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label", fill = 'white')  +
  facet_wrap(~pop, ncol = 3, nrow = 2) +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency') +
  ggtitle('1000 Genomes')

cl13_antagonist_ukbb = signifSNPs %>%
  filter(cl1_cl3 == ' antagonist' & !grepl('agonist',cl1_cl2)) %>%
  # filter(!gr_antagonist) %>%
  # filter((!gr_antagonist) & (ac_antagonist)) %>%
  # filter(numCl == 1) %>% 
  # filter(numDis == 1) %>%
  select(-disID,-disease,-BETA,-P_BOLT_LMM_INF) %>%
  unique() %>%
  mutate(cluster = factor(cluster)) %>% 
  # filter(cluster == 2) %>%
  ggplot(aes(y = UKBB_RAF, x= cluster, fill = cluster)) +
  # geom_boxplot() +
  geom_violin() +
  geom_boxplot(width = 0.1, fill ='white') +
  scale_fill_manual(values = ageonsetcolors) +
  stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label", fill = 'white')  +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency') +
  ggtitle('UK Biobank')

all_cl13_antagonist = ggarrange(cl13_antagonist_ukbb,cl13_antagonist_1kg,ncol = 2, nrow = 1, widths = c(1,2), labels = 'auto')

xx = signifSNPs %>%
  filter(cl1_cl3 == ' antagonist' & !grepl('agonist',cl1_cl2)) %>%
  select(SNP, RA, UKBB_RAF, cluster) %>%
  unique() 

xx %>% group_by(cluster) %>% summarise(n = length(unique(SNP)))
# # A tibble: 2 x 2
# cluster     n
# <int> <int>
# 1       1   226
# 2       3   226

cl1perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==1)%>%na.omit())$UKBB_RAF,100)))
cl3perm = sapply(1:1000,function(i)median(sample((filter(xx,cluster==3)%>%na.omit())$UKBB_RAF,100)))

antagonist_ukbb_sampling = data.frame(cluster = factor(rep(c(1,3),each=1000),levels = c(1,3)), median = c(cl1perm,cl3perm)) %>%
  ggplot(aes(x = cluster, y= median)) +
  geom_sina(alpha=0.3,size=0.5) + 
  stat_summary(fun.y = 'mean', colour = "gray60", size = 0.5, geom = 'hline',
               aes(yintercept = ..y..), linetype = 'dashed') +
  xlab('Age of onset cluster') + ylab('Median RAF (for 100 SNPs)')

allukbb=ggarrange(onedisease_ukbb+ggtitle('SNPs associated with\none disease'), onecluster_noantagonist_ukbb+ggtitle('SNPs associated with\none cluster'), cl13_antagonist_ukbb+ggtitle('SNPs Antagonistic Between\nCluster1 and Cluster3'), ncol=3,labels = 'auto',nrow=1) 

allukbb_sampling=ggarrange(onedisease_ukbb_sampling+ggtitle('SNPs associated with\none disease'), onecluster_ukbb_sampling+ggtitle('SNPs associated with\none cluster'), antagonist_ukbb_sampling+ggtitle('SNPs Antagonistic Between\nCluster1 and Cluster3'), ncol=3,labels = 'auto',nrow=1) 

ggsave('./results/evoAnalysis/RAF_ukbb.pdf',allukbb,units = 'cm',width = 20 ,height = 7,useDingbats =F)
ggsave('./results/evoAnalysis/RAF_ukbb.png',allukbb,units = 'cm',width = 20,height = 7)

ggsave('./results/evoAnalysis/RAF_ukbb_sampling.pdf',allukbb_sampling,units = 'cm',width = 20 ,height = 7,useDingbats =F)
ggsave('./results/evoAnalysis/RAF_ukbb_sampling.png',allukbb_sampling,units = 'cm',width = 20,height = 7)


ggsave('./results/evoAnalysis/RAF_all_onedisease.pdf',all_onedisease, units = 'cm', width = 24, height = 12, useDingbats = F)
ggsave('./results/evoAnalysis/RAF_all_onedisease.png',all_onedisease, units = 'cm', width = 24, height = 12)


ggsave('./results/evoAnalysis/RAF_all_onecluster_noantagonist.pdf',all_onecluster_noantagonist, units = 'cm', width = 24, height = 12, useDingbats = F)
ggsave('./results/evoAnalysis/RAF_all_onecluster_noantagonist.png',all_onecluster_noantagonist, units = 'cm', width = 24, height = 12)


ggsave('./results/evoAnalysis/all_cl13_antagonist.pdf',all_cl13_antagonist, units = 'cm', width = 24, height = 12, useDingbats = F)
ggsave('./results/evoAnalysis/all_cl13_antagonist.png',all_cl13_antagonist, units = 'cm', width = 24, height = 12)



a = unique(filter(signifSNPs,cl1_cl1 == 'agonist ' | cl2_cl2 == 'agonist ')$SNP)
b = unique(filter(signifSNPs,cl1_cl1 == ' antagonist' | cl2_cl2 == ' antagonist')$SNP)
c = unique(filter(signifSNPs,cl1_cl2 == 'agonist ')$SNP)
d = unique(filter(signifSNPs,cl1_cl2 == ' antagonist')$SNP)
cl1cl2_mat = matrix(c(length(a),length(b),length(c),length(d)),byrow = T,ncol=2)
cl1cl2_fi = fisher.test(cl1cl2_mat)
# Fisher's Exact Test for Count Data
# 
# data:  cl1cl2_mat
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  36.69307 55.21408
# sample estimates:
# odds ratio 
#   44.88721 
log2(cl1cl2_fi$estimate)
# odds ratio 
# 5.488233
excx = names(which(table(c(a,b,c,d))!=1))
a = setdiff(a,excx)
b = setdiff(b,excx)
c = setdiff(c,excx)
d = setdiff(d,excx)
cl1cl2_mat_exc = matrix(c(length(a),length(b),length(c),length(d)),byrow = T,ncol=2)
cl1cl2_fi_exc = fisher.test(cl1cl2_mat_exc)
# Fisher's Exact Test for Count Data
# 
# data:  cl1cl2_mat_exc
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   777.9397 1863.9628
# sample estimates:
# odds ratio 
#   1169.171
log2(cl1cl2_fi_exc$estimate)
# odds ratio 
# 10.19127
a = unique(filter(signifSNPs,cl1_cl1 == 'agonist ' | cl3_cl3 == 'agonist ')$SNP)
b = unique(filter(signifSNPs,cl1_cl1 == ' antagonist' | cl3_cl3 == ' antagonist')$SNP)
c = unique(filter(signifSNPs,cl1_cl3 == 'agonist ')$SNP)
d = unique(filter(signifSNPs,cl1_cl3 == ' antagonist')$SNP)
cl1cl3_mat = matrix(c(length(a),length(b),length(c),length(d)),byrow = T,ncol=2)
cl1cl3_fi = fisher.test(cl1cl3_mat)
# Fisher's Exact Test for Count Data
# 
# data:  cl1cl3_mat
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  142.8577 224.6479
# sample estimates:
# odds ratio 
#   177.9271
log2(cl1cl3_fi$estimate)
# odds ratio 
# 7.475142
excx = names(which(table(c(a,b,c,d))!=1))
a = setdiff(a,excx)
b = setdiff(b,excx)
c = setdiff(c,excx)
d = setdiff(d,excx)
cl1cl3_mat_exc = matrix(c(length(a),length(b),length(c),length(d)),byrow = T,ncol=2)
cl1cl3_fi_exc = fisher.test(cl1cl3_mat_exc)
# Fisher's Exact Test for Count Data
# 
# data:  cl1cl3_mat_exc
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  186.5824 341.8882
# sample estimates:
# odds ratio 
#   251.7659
log2(cl1cl3_fi_exc$estimate)
# odds ratio 
# 7.975939


a = unique(filter(signifSNPs,cl2_cl2 == 'agonist ' | cl3_cl3 == 'agonist ')$SNP)
b = unique(filter(signifSNPs,cl2_cl2 == ' antagonist' | cl3_cl3 == ' antagonist')$SNP)
c = unique(filter(signifSNPs,cl2_cl3 == 'agonist ')$SNP)
d = unique(filter(signifSNPs,cl2_cl3 == ' antagonist')$SNP)
cl2cl3_mat = matrix(c(length(a),length(b),length(c),length(d)),byrow = T,ncol=2)
cl2cl3_fi = fisher.test(cl2cl3_mat)
# Fisher's Exact Test for Count Data
# 
# data:  cl2cl3_mat
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   525.1282 1251.2094
# sample estimates:
# odds ratio 
#   821.6346
log2(cl2cl3_fi$estimate)
# odds ratio 
# 9.682353
excx = names(which(table(c(a,b,c,d))!=1))
a = setdiff(a,excx)
b = setdiff(b,excx)
c = setdiff(c,excx)
d = setdiff(d,excx)
cl2cl3_mat_exc = matrix(c(length(a),length(b),length(c),length(d)),byrow = T,ncol=2)
cl2cl3_fi_exc = fisher.test(cl2cl3_mat_exc)
# Fisher's Exact Test for Count Data
# 
# data:  cl2cl3_mat_exc
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1204.706 4300.262
# sample estimates:
# odds ratio 
#   2002.389 
log2(cl2cl3_fi_exc$estimate)
# odds ratio 
# 10.96751 

xx = signifSNPs %>%
  filter((!gr_antagonist) & (!ac_antagonist)) %>%
  filter(numCl == 1) %>% 
  select(-disID,-disease,-BETA,-P_BOLT_LMM_INF) %>%
  unique() %>%
  mutate(cluster = factor(cluster)) %>%
  mutate(rafgr = cut(UKBB_RAF,breaks =seq(0,1,by=0.1), include.lowest = T)) %>%
  group_by(cluster,rafgr) %>%
  summarise(ancestral = sum(Ancestral,na.rm=T),
            ancestral_lc = sum(Ancestral_LC,na.rm=T),
            nancestral = sum(!Ancestral,na.rm=T),
            nancestral_lc = sum(!Ancestral_LC,na.rm=T),
            mean_ancestral = mean(Ancestral,na.rm=T),
            mean_ancestral_lc = mean(Ancestral_LC,na.rm=T))

onecluster_noantagonist_aa = xx %>%
  ggplot(aes(x = rafgr, y = mean_ancestral, fill = cluster, group  = cluster)) +
  geom_bar(stat='identity', position = 'dodge') +
  scale_fill_manual(values = ageonsetcolors) +
  xlab('Risk Allele Frequency (UK Biobank)') +
  ylab('Ancestral Allele %') +
  guides(fill = guide_legend('Age of Onset Cluster')) +
  theme(legend.position = 'top')

xx = signifSNPs %>%
  # filter((!gr_antagonist) & (!ac_antagonist)) %>%
  filter(numDis == 1) %>% 
  select(-disID,-disease,-BETA,-P_BOLT_LMM_INF) %>%
  unique() %>%
  mutate(cluster = factor(cluster)) %>%
  mutate(rafgr = cut(UKBB_RAF,breaks =seq(0,1,by=0.1), include.lowest = T)) %>%
  group_by(cluster,rafgr) %>%
  summarise(ancestral = sum(Ancestral,na.rm=T),
            ancestral_lc = sum(Ancestral_LC,na.rm=T),
            nancestral = sum(!Ancestral,na.rm=T),
            nancestral_lc = sum(!Ancestral_LC,na.rm=T),
            mean_ancestral = mean(Ancestral,na.rm=T),
            mean_ancestral_lc = mean(Ancestral_LC,na.rm=T))

onedisease_aa = xx %>%
  ggplot(aes(x = rafgr, y = mean_ancestral, fill = cluster, group  = cluster)) +
  geom_bar(stat='identity', position = 'dodge') +
  scale_fill_manual(values = ageonsetcolors) +
  xlab('Risk Allele Frequency (UK Biobank)') +
  ylab('Ancestral Allele %') +
  guides(fill = guide_legend('Age of Onset Cluster')) +
  theme(legend.position = 'top')

snplist = names(which(table((signifSNPs %>%
  filter(cl1_cl3 == ' antagonist') %>%
  select(-disID,-disease,-BETA,-P_BOLT_LMM_INF) %>%
  unique() %>%
  mutate(cluster = factor(cluster)) %>% 
  filter(cluster%in%c(1,3))%>%
  select(SNP,Ancestral_LC) %>% unique())$SNP)==2))
snplist = (signifSNPs %>%
  filter(SNP %in% snplist) %>%
  filter(cluster==1) %>%
  filter(Ancestral_LC))$SNP
cl1raaa = signifSNPs %>%
  filter(SNP %in%snplist) %>%
  filter(cl1_cl3 == ' antagonist' & grepl('agonist',cl1_cl2)) %>%
  select(-disID,-disease,-BETA,-P_BOLT_LMM_INF) %>%
  unique() %>%
  mutate(cluster = factor(cluster)) %>% 
  # filter(cluster == 2) %>%
  ggplot(aes(y = UKBB_RAF, x= cluster, fill = cluster)) +
  # geom_boxplot() +
# geom_violin() +
geom_boxplot() +
  stat_summary(fun.y = 'mean', aes(label = round(..y..,2)), geom = "label", fill = 'white')  +
  geom_jitter(width = 0.2, alpha = 0.7, size = 0.3)+
  scale_fill_manual(values = ageonsetcolors)  +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency (UKBB)') +
  ggtitle('Cluster 1 Risk Allele = Ancestral') + ylim( 0,1)

snplist = names(which(table((signifSNPs %>%
                               filter(cl1_cl3 == ' antagonist') %>%
                               select(-disID,-disease,-BETA,-P_BOLT_LMM_INF) %>%
                               unique() %>%
                               mutate(cluster = factor(cluster)) %>% 
                               filter(cluster%in%c(1,3))%>%
                               select(SNP,Ancestral_LC) %>% unique())$SNP)==2))
snplist = (signifSNPs %>%
             filter(SNP %in% snplist) %>%
             filter(cluster==3) %>%
             filter(Ancestral_LC))$SNP
cl3raaa = signifSNPs %>%
  filter(SNP %in%snplist) %>%
  filter(cl1_cl3 == ' antagonist' & grepl('agonist',cl1_cl2)) %>%
  select(-disID,-disease,-BETA,-P_BOLT_LMM_INF) %>%
  unique() %>%
  mutate(cluster = factor(cluster)) %>% 
  # filter(cluster == 2) %>%
  ggplot(aes(y = UKBB_RAF, x= cluster, fill = cluster)) +
  # geom_boxplot() +
  # geom_violin() +
  geom_boxplot() +
  stat_summary(fun.y = 'mean', aes(label = round(..y..,2)), geom = "label", fill = 'white')  +
  geom_jitter(width = 0.2, alpha = 0.7, size = 0.3)+
  scale_fill_manual(values = ageonsetcolors) +
  guides(fill = F) +
  xlab('Age of onset cluster') + ylab('Risk Allele Frequency (UKBB)') +
  ggtitle('Cluster 3 Risk Allele = Ancestral') + ylim(0,1)

ggarrange(cl1raaa,cl3raaa,labels = 'auto')

#####

xx = signifSNPs %>%
  filter((!gr_antagonist) & (!ac_antagonist)) %>%
  filter(numDis == 1) %>% 
  select(SNP, RA, UKBB_RAF, cluster, Ancestral_LC) %>%
  unique()
library(ggthemes)
xx  = xx %>%
  group_by(Ancestral_LC,cluster) %>%
  summarise(n = length(unique(SNP))) %>%
  mutate(cluster = factor(cluster)) %>% 
  spread(cluster,n) %>% as.data.frame()
rownames(xx) = c('Derived','Ancestral','-')
xx$Ancestral_LC= NULL
colnames(xx) = paste('Cluster',1:3)
xx = as.matrix(xx)
ancestral_st = reshape2::melt(t(t(xx)/colSums(xx))) %>%
  ggplot(aes(x = as.factor(Var1), y = value, fill = as.factor(Var1))) +
  geom_bar(stat = 'identity', position='dodge') +
  scale_fill_ptol()+
  geom_label(aes(label = paste(round(value*100,1),'%',sep='')),fill='white') +
  facet_wrap(~Var2) +
  xlab('') + ylab('Proportion') +
  guides(fill = guide_legend('Ancestral State'))

ggsave('./results/evoAnalysis/ancestral_prop.pdf', ancestral_st,units = 'cm', width = 16.7, height = 9, useDingbats = F)
ggsave('./results/evoAnalysis/ancestral_prop.png', ancestral_st,units = 'cm', width = 16.7, height = 9)



xx = signifSNPs %>%
  filter((!gr_antagonist) & (!ac_antagonist)) %>%
  filter(numDis == 1) %>% 
  select(SNP, RA, UKBB_RAF, cluster, Ancestral_LC,disID) %>%
  unique()
library(ggthemes)
xx  = xx %>%
  group_by(Ancestral_LC,disID) %>%
  summarise(n = length(unique(SNP))) %>%
  spread(disID,n,fill=0) %>% as.data.frame()
rownames(xx) = c('Derived','Ancestral','-')
xx$Ancestral_LC= NULL
xx = as.matrix(xx)
xx = xx[1:2,]

ancestral_st_dis = reshape2::melt(t(t(xx)/colSums(xx,na.rm=T))) %>%
  mutate(cluster=factor(unname(readRDS('../ukbb_ageonset/data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$clustering[disCoding[as.character(Var2)]]))) %>% 
  filter(Var1=='Ancestral') %>%
  ggplot(aes(x = cluster, y = value, fill = cluster)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = ageonsetcolors) +
  xlab('Age of Onset Cluster') +
  guides(fill =F) +
  ylab('Proportion of Ancestral Alleles')
ggsave('./results/evoAnalysis/ancestral_prop_dis.pdf', ancestral_st_dis,units = 'cm', width = 8, height = 9, useDingbats = F)
ggsave('./results/evoAnalysis/ancestral_prop_dis.png', ancestral_st_dis,units = 'cm', width = 8, height = 9)

ancestralxx = ggarrange(ancestral_st + theme(legend.position='top'),
                        ancestral_st_dis,labels = 'auto',
                        nrow = 2,ncol=1, heights = c(1.1,1))

ggsave('./results/evoAnalysis/ancestral_props.png', ancestralxx,units = 'cm', width = 18, height = 15)
ggsave('./results/evoAnalysis/ancestral_props.pdf', ancestralxx,units = 'cm', width = 18, height = 15, useDingbats =F)

antagonistic_table = data.frame(clusters = c('cl1-cl2','cl1-cl3','cl2-cl3'),
           odds_ratio = c(cl1cl2_fi$estimate,
                          cl1cl3_fi$estimate,
                          cl2cl3_fi$estimate),
           p_value = c(cl1cl2_fi$p.value,
                       cl1cl3_fi$p.value,
                       cl2cl3_fi$p.value),
           odds_ratio_ind = c(cl1cl2_fi_exc$estimate,
                          cl1cl3_fi_exc$estimate,
                          cl2cl3_fi_exc$estimate),
           p_value_ind = c(cl1cl2_fi_exc$p.value,
                       cl1cl3_fi_exc$p.value,
                       cl2cl3_fi_exc$p.value)) %>%
  mutate(log2_odds = log2(odds_ratio),
         log2_odds_ind = log2(odds_ratio_ind))
library(gridExtra)
xx = tableGrob(antagonistic_table,rows=NULL)
ggsave('./results/evoAnalysis/antagonistic.pdf',xx,units = 'cm',width=22,height = 4, useDingbats = F)
ggsave('./results/evoAnalysis/antagonistic.png',xx,units = 'cm',width=22,height = 4)
