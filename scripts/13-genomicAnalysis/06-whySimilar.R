source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)

snpsim_odds = readRDS('./data/processed/genomicAnalysis/diseaseSim_odds_basedonSNPs.rds')
snpsim_p = readRDS('./data/processed/genomicAnalysis/diseaseSim_p_basedonSNPs.rds')

genesim_odds = readRDS('./data/processed/genomicAnalysis/diseaseSim_odds_basedonGenes.rds')
genesim_p = readRDS('./data/processed/genomicAnalysis/diseaseSim_p_basedonGenes.rds')

proxy_genesim_odds = readRDS('./data/processed/genomicAnalysis/diseaseSim_odds_basedonProxyGenes.rds')
proxy_genesim_p = readRDS('./data/processed/genomicAnalysis/diseaseSim_p_basedonProxyGenes.rds')

eqtl_genesim_odds = readRDS('./data/processed/genomicAnalysis/diseaseSim_odds_basedoneQTLGenes.rds')
eqtl_genesim_p = readRDS('./data/processed/genomicAnalysis/diseaseSim_p_basedoneQTLGenes.rds')

disnm = unname(disCoding[as.character(disIDs)])

catsim = sapply(disnm, function(x){
  x = unname(disTreecl[x])
  sapply(disnm, function(y){
    y = unname(disTreecl[y])
    ifelse(x==y,return(1),return(0))
  })
})

RRtable <- readRDS('./data/processed/diseaseCooccur/RRtable.rds')
cormat <- readRDS('./data/processed/diseaseCooccur/disCorrelations.rds')

mydis = readRDS('./data/processed/ageonset/distmat_cort.rds')
ageonsetclusters = readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')

snpsim_odds = reshape2::melt(snpsim_odds) %>%
  set_names(c('disA','disB','snp_odds'))

snpsim_p = reshape2::melt(snpsim_p) %>%
  set_names(c('disA','disB','snp_p'))

genesim_odds = reshape2::melt(genesim_odds) %>%
  set_names(c('disA','disB','gene_odds'))

genesim_p = reshape2::melt(genesim_p) %>%
  set_names(c('disA','disB','gene_p'))

proxy_genesim_odds = reshape2::melt(proxy_genesim_odds) %>%
  set_names(c('disA','disB','proxy_gene_odds'))

proxy_genesim_p = reshape2::melt(proxy_genesim_p) %>%
  set_names(c('disA','disB','proxy_gene_p'))

eqtl_genesim_odds = reshape2::melt(eqtl_genesim_odds) %>%
  set_names(c('disA','disB','eqtl_gene_odds'))

eqtl_genesim_p = reshape2::melt(eqtl_genesim_p) %>%
  set_names(c('disA','disB','eqtl_gene_p'))

catsim = reshape2::melt(catsim) %>%
  set_names(c('disA','disB','sameCat'))

RRtable = RRtable %>%
  select(disA,disB,RR,sublevel,uplevel)

cormat = reshape2::melt(cormat) %>%
  set_names(c('disA','disB','cooccur'))

mydis = reshape2::melt(as.matrix(mydis)) %>%
  set_names(c('disA','disB','age_distance'))

agesim = sapply(disnm, function(x){
  x = unname(ageonsetclusters$clustering[x])
  sapply(disnm, function(y){
    y = unname(ageonsetclusters$clustering[y])
    # ifelse(x==y,paste(x,y,sep='-'),paste(min(x,y),'other',sep='-'))
    # paste(x,y,sep='-')
    ifelse(x==y,1,0)
  })
})

agesim = reshape2::melt(agesim) %>%
  set_names(c('disA','disB','sameAge'))

alldat = snpsim_odds %>%
  full_join(snpsim_p) %>%
  full_join(genesim_odds) %>%
  full_join(genesim_p) %>%
  full_join(proxy_genesim_odds) %>%
  full_join(proxy_genesim_p) %>%
  full_join(eqtl_genesim_odds) %>%
  full_join(eqtl_genesim_p) %>%
  full_join(catsim) %>%
  full_join(RRtable) %>%
  full_join(cormat) %>%
  full_join(mydis) %>%
  full_join(agesim)

alldat = alldat %>%
  filter(disA != disB)

alldat = alldat %>%
  filter(!(sublevel | uplevel))

mydat = alldat %>%
  mutate(snp_odds = log2(snp_odds),
         RR = log2(RR),
         gene_odds = log2(gene_odds),
         proxy_gene_odds = log2(proxy_gene_odds),
         eqtl_gene_odds = log2(eqtl_gene_odds)) %>%
  # mutate(snpodds = ifelse(snp_p<0.01,snp_odds,1)) %>%
  filter(snp_p<0.01) %>%
  select(disA,disB,snp_odds,sameCat,cooccur,sameAge) 

a=ifelse(mydat$disA>mydat$disB,mydat$disA,mydat$disB)
b=ifelse(mydat$disA>mydat$disB,mydat$disB,mydat$disA)
mydat$disA = a
mydat$disB = b
nrow(mydat)
nrow(unique(mydat))
mydat = unique(mydat)

library(igraph)
library(ggnetwork)
library(ggnet)
library(ggforce)

mynet = mydat %>%
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(snp_odds ~ sameCat + cooccur, data=mydat)$resid) %>% 
  select(disA,disB,snp_odds) %>%
  na.omit() %>%
  unique() %>% graph_from_data_frame()

V(mynet)$ageonset_cl = ageonsetclusters$clustering[V(mynet)$name]
V(mynet)$discat = as.character(disTreecl[V(mynet)$name])
E(mynet)$edgewidth = summary(E(mynet)$snp_odds/30)
net0 = ggnet2(mynet, size = 0, edge.size = 'edgewidth', edge.color = 'gray80')+
  theme_void()
net1 = net0 +
  geom_point(size = 2, color = ageonsetcolors[as.character(V(mynet)$ageonset_cl)])
net2 = net0 +
  geom_point(size = 2, color = discatcolors[V(mynet)$discat])

# ggsave('./results/genomicAnalysis/SNPnet_ageonsetcl.pdf', net1, units = 'cm', width = 8, height = 8, useDingbats = F)
# ggsave('./results/genomicAnalysis/SNPnet_ageonsetcl.png', net1, units = 'cm', width = 8, height = 8)

# ggsave('./results/genomicAnalysis/SNPnet_discat.pdf', net2, units = 'cm', width = 8, height = 8, useDingbats = F)
# ggsave('./results/genomicAnalysis/SNPnet_discat.png', net2, units = 'cm', width = 8, height = 8)

net3 = ggarrange(net1,net2,labels = 'auto')
# ggsave('./results/genomicAnalysis/SNPnet_discat_onset_together.pdf', net3, units = 'cm', width = 16, height = 8, useDingbats = F)
# ggsave('./results/genomicAnalysis/SNPnet_discat_onset_together.png', net3, units = 'cm', width = 16, height = 8)



mynet = mydat %>%
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(snp_odds ~ sameCat + cooccur, data=mydat)$resid) %>% 
  select(disA,disB,corrected_val) %>%
  na.omit() %>%
  unique() %>% graph_from_data_frame()

V(mynet)$ageonset_cl = ageonsetclusters$clustering[V(mynet)$name]
V(mynet)$discat = as.character(disTreecl[V(mynet)$name])
E(mynet)$edgewidth = summary((6+E(mynet)$corrected_val)/30)
net0 = ggnet2(mynet, size = 0, edge.size = 'edgewidth', edge.color = 'gray80')+
  theme_void()
net1 = net0 +
  geom_point(size = 2, color = ageonsetcolors[as.character(V(mynet)$ageonset_cl)])
net2 = net0 +
  geom_point(size = 2, color = discatcolors[V(mynet)$discat])

# ggsave('./results/genomicAnalysis/SNPnet_ageonsetcl_corr.pdf', net1, units = 'cm', width = 8, height = 8, useDingbats = F)
# ggsave('./results/genomicAnalysis/SNPnet_ageonsetcl_corr.png', net1, units = 'cm', width = 8, height = 8)

# ggsave('./results/genomicAnalysis/SNPnet_discat_corr.pdf', net2, units = 'cm', width = 8, height = 8, useDingbats = F)
# ggsave('./results/genomicAnalysis/SNPnet_discat_corr.png', net2, units = 'cm', width = 8, height = 8)

net3 = ggarrange(net1,net2,labels = 'auto')
# ggsave('./results/genomicAnalysis/SNPnet_corr_discat_onset_together.pdf', net3, units = 'cm', width = 16, height = 8, useDingbats = F)
# ggsave('./results/genomicAnalysis/SNPnet_corr_discat_onset_together.png', net3, units = 'cm', width = 16, height = 8)

myodds = mynet[]
myodds[myodds!=0]=E(mynet)$edgewidth
myodds = as.matrix(myodds)
myodds = sapply(rownames(myodds),function(a){
  sapply(rownames(myodds),function(b){
    aa=c(myodds[a,b],myodds[b,a])
    if(all(aa==0)){
      0
    } else{
        setdiff(aa,0)
      }
  })
})


disannot = data.frame(diseaseCategories = unname(disTreecl[rownames(myodds)]), ageonset_clusters = as.character(unname(readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$clustering[rownames(myodds)])))
rownames(disannot) =rownames(myodds)
annotcolors = list(diseaseCategories = discatcolors,ageonset_clusters = ageonsetcolors)
library(pheatmap)
pheatmap(myodds, color = colorRampPalette(brewer.pal(8,'Reds'))(11),cellwidth = 5, cellheight = 5, annotation_row = disannot, annotation_col = disannot[,-3], annotation_colors = annotcolors, fontsize = 6, filename = './results/genomicAnalysis/signifSNPoverlap_heatmap_all.pdf')


lmmod = lm(snp_odds ~ sameCat + cooccur + sameAge, data=mydat) 
aovobj = aov(lmmod)
summary(aovobj)
par(mfrow=c(1,2))
boxplot(snp_odds~sameAge,data=mydat,main = 'uncorrected')
boxplot(lm(snp_odds ~ sameCat + cooccur, data=mydat)$resid~sameAge,data=mydat,main = 'corrected')

rawp = mydat %>%
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(snp_odds ~ sameCat + cooccur, data=mydat)$resid) %>%
  ggplot(aes(x = sameAge, y= snp_odds)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, color = 'gray40') +
  ylab('Genetic Similarity (on log2)') + xlab('') +
  ggtitle('\nRaw values')


correctedp = mydat %>%
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(snp_odds ~ sameCat + cooccur, data=mydat)$resid) %>%
  ggplot(aes(x = sameAge, y= corrected_val)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = 'gray40', size= 0.5, width = 0.2) +
  ylab('Genetic Similarity (on log2)') + xlab('') +
  ggtitle('Corrected values\n(by category and co-occurrence)') 


snpp = ggarrange(rawp,correctedp,labels='auto')
ggsave('./results/genomicAnalysis/SNPsimilarity_bysameAge.pdf',snpp, units = 'cm', width = 16, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPsimilarity_bysameAge.png',snpp, units = 'cm', width = 16, height = 7)

netandbox = ggarrange(net1,correctedp,labels='auto')
ggsave('./results/genomicAnalysis/SNPsimilarity_age_netbox.pdf',netandbox, units = 'cm', width = 16, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPsimilarity_age_netbox.png',netandbox, units = 'cm', width = 16, height = 7)

mydat = alldat %>%
  mutate(snp_odds = log2(snp_odds),
         RR = log2(RR),
         gene_odds = log2(gene_odds),
         proxy_gene_odds = log2(proxy_gene_odds),
         eqtl_gene_odds = log2(eqtl_gene_odds)) %>%
  filter(snp_p<0.01) %>%
  select(disA,disB,snp_odds,sameCat,cooccur,sameAge) %>% 
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(snp_odds ~ sameCat + cooccur, data=.)$resid) %>%
  mutate(disA_ageonset = factor(ageonsetclusters$clustering[as.character(disA)], levels = 1:4),
         disB_ageonset = factor(ageonsetclusters$clustering[as.character(disB)], levels = 1:4),
         disA_category = disTreecl[as.character(disA)],
         disB_category = disTreecl[as.character(disB)])

rawp1 = mydat %>%
  filter(disA_ageonset ==1) %>%
  ggplot(aes(x = disB_ageonset, y= snp_odds)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(size = 0.5) +
  ylab('Genetic Similarity (on log2)') + xlab('') + scale_x_discrete(drop=F,breaks=factor(1:4)) +
  ggtitle('Similarity between age-of-onset cluster 1\nand other clusters (with raw values)') 

correctedp1 = mydat %>%
  filter(disA_ageonset ==1) %>%
  mutate(corrected_val = lm(snp_odds ~ sameCat + cooccur, data=.)$resid) %>%
  ggplot(aes(x = disB_ageonset, y= corrected_val)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter( size = 0.5) +
  ylab('Genetic Similarity (on log2)') + xlab('') + scale_x_discrete(drop=F,breaks=factor(1:4)) +
  ggtitle('Similarity between age-of-onset cluster 1 and others\nwith corrected values (by category and co-occurrence)') 

snpp1 = ggarrange(rawp1,correctedp1,labels='auto', nrow = 2, ncol=1, common.legend = T, legend = 'right')
ggsave('./results/genomicAnalysis/SNPsimilarity_cl1_bysameAge.pdf',snpp1, units = 'cm', width = 16, height = 14, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPsimilarity_cl1_bysameAge.png',snpp1, units = 'cm', width = 16, height = 14)


rawp2 = mydat %>%
  filter(disA_ageonset ==2) %>%
  ggplot(aes(x = disB_ageonset, y= snp_odds)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter( size = 0.5) +
  ylab('Genetic Similarity (on log2)') + xlab('') + scale_x_discrete(drop=F,breaks=factor(1:4)) +
  ggtitle('Similarity between age-of-onset cluster 2\nand other clusters (with raw values)')

correctedp2 = mydat %>%
  filter(disA_ageonset ==2) %>%
  mutate(corrected_val = lm(snp_odds ~ sameCat + cooccur, data=.)$resid) %>%
  ggplot(aes(x = disB_ageonset, y= corrected_val)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter( size = 0.5) +
  ylab('Genetic Similarity (on log2)') + xlab('') + scale_x_discrete(drop=F,breaks=factor(1:4)) +
  ggtitle('Similarity between age-of-onset cluster 2 and others\nwith corrected values (by category and co-occurrence)')

snpp2 = ggarrange(rawp2,correctedp2,labels='auto', nrow = 2, ncol=1, common.legend = T, legend = 'right')
ggsave('./results/genomicAnalysis/SNPsimilarity_cl2_bysameAge.pdf',snpp2, units = 'cm', width = 16, height = 14, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPsimilarity_cl2_bysameAge.png',snpp2, units = 'cm', width = 16, height = 14)

rawp3 = mydat %>%
  filter(disA_ageonset ==3) %>%
  ggplot(aes(x = disB_ageonset, y= snp_odds)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.5) +
  ylab('Genetic Similarity (on log2)') + xlab('') + scale_x_discrete(drop=F,breaks=factor(1:4)) +
  ggtitle('Similarity between age-of-onset cluster 3\nand others (with raw values)')

correctedp3 = mydat %>%
  filter(disA_ageonset ==3) %>%
  mutate(corrected_val = lm(snp_odds ~ sameCat + cooccur, data=.)$resid) %>%
  ggplot(aes(x = disB_ageonset, y= corrected_val)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter( size = 0.5) +
  ylab('Genetic Similarity (on log2)') + xlab('') + scale_x_discrete(drop=F,breaks=factor(1:4)) +
  ggtitle('Similarity between age-of-onset cluster 3 and others \nwith corrected values (by category and co-occurrence)')

snpp3 = ggarrange(rawp3,correctedp3,labels='auto', nrow = 2, ncol=1, common.legend = T, legend = 'right')
ggsave('./results/genomicAnalysis/SNPsimilarity_cl3_bysameAge.pdf',snpp3, units = 'cm', width = 16, height = 14, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPsimilarity_cl3_bysameAge.png',snpp3, units = 'cm', width = 16, height = 14)

allcltgthr = ggarrange(rawp1,correctedp1,rawp2,correctedp2,rawp3,correctedp3, nrow = 3, ncol =2,labels = 'auto', align = 'hv')
ggsave('./results/genomicAnalysis/SNPsimilarity_clusters_bysameAge.pdf',allcltgthr, units = 'cm', width = 24, height = 24, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPsimilarity_clusters_bysameAge.png',allcltgthr, units = 'cm', width = 24, height = 24)



####

mydat = alldat %>%
  mutate(snp_odds = log2(snp_odds),
         RR = log2(RR),
         gene_odds = log2(gene_odds),
         proxy_gene_odds = log2(proxy_gene_odds),
         eqtl_gene_odds = log2(eqtl_gene_odds)) %>%
  # mutate(snpodds = ifelse(snp_p<0.01,snp_odds,1)) %>%
  filter(gene_p<0.01) %>%
  select(disA,disB,gene_odds,sameCat,cooccur,sameAge) 

a=ifelse(mydat$disA>mydat$disB,mydat$disA,mydat$disB)
b=ifelse(mydat$disA>mydat$disB,mydat$disB,mydat$disA)
mydat$disA = a
mydat$disB = b
nrow(mydat)
nrow(unique(mydat))
mydat = unique(mydat)

lmmod = lm(gene_odds ~ sameCat + cooccur + sameAge, data=mydat) 
aovobj = aov(lmmod)
summary(aovobj)
par(mfrow=c(1,2))
boxplot(gene_odds~sameAge,data=mydat,main = 'uncorrected')
boxplot(lm(gene_odds ~ sameCat + cooccur, data=mydat)$resid~sameAge,data=mydat,main = 'corrected')

rawp = mydat %>%
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(gene_odds ~ sameCat + cooccur, data=mydat)$resid) %>%
  ggplot(aes(x = sameAge, y= gene_odds)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_sina(color = 'gray40', size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('\nRaw values')

correctedp = mydat %>%
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(gene_odds ~ sameCat + cooccur, data=mydat)$resid) %>%
  ggplot(aes(x = sameAge, y= corrected_val)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_sina(color = 'gray40', size= 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('Corrected values\n(by category and co-occurrence)')

snpp = ggarrange(rawp,correctedp,labels='auto')
ggsave('./results/genomicAnalysis/Genesimilarity_bysameAge.pdf',snpp, units = 'cm', width = 18, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/Genesimilarity_bysameAge.png',snpp, units = 'cm', width = 18, height = 7)

mydat = alldat %>%
  mutate(snp_odds = log2(snp_odds),
         RR = log2(RR),
         gene_odds = log2(gene_odds),
         proxy_gene_odds = log2(proxy_gene_odds),
         eqtl_gene_odds = log2(eqtl_gene_odds)) %>%
  filter(gene_p<0.01) %>%
  select(disA,disB,gene_odds,sameCat,cooccur,sameAge) %>% 
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(gene_odds ~ sameCat + cooccur, data=.)$resid) %>%
  mutate(disA_ageonset = factor(ageonsetclusters$clustering[as.character(disA)]),
         disB_ageonset = factor(ageonsetclusters$clustering[as.character(disB)]),
         disA_category = disTreecl[as.character(disA)],
         disB_category = disTreecl[as.character(disB)])

rawp = mydat %>%
  filter(disA_ageonset ==1) %>%
  ggplot(aes(x = disB_ageonset, y= gene_odds)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('Similarity between age-of-onset cl1\nand other clusters (with raw values)')

correctedp = mydat %>%
  filter(disA_ageonset ==1) %>%
  mutate(corrected_val = lm(gene_odds ~ sameCat + cooccur, data=.)$resid) %>%
  ggplot(aes(x = disB_ageonset, y= corrected_val)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('\nwith corrected values (by category and co-occurrence)')

snpp = ggarrange(rawp,correctedp,labels='auto', nrow = 2, ncol=1, common.legend = T, legend = 'right')
ggsave('./results/genomicAnalysis/Genesimilarity_cl1_bysameAge.pdf',snpp, units = 'cm', width = 16, height = 14, useDingbats = F)
ggsave('./results/genomicAnalysis/Genesimilarity_cl1_bysameAge.png',snpp, units = 'cm', width = 16, height = 14)


rawp = mydat %>%
  filter(disA_ageonset ==2) %>%
  ggplot(aes(x = disB_ageonset, y= gene_odds)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('Similarity between age-of-onset cl2\nand other clusters (with raw values)')

correctedp = mydat %>%
  filter(disA_ageonset ==2) %>%
  mutate(corrected_val = lm(gene_odds ~ sameCat + cooccur, data=.)$resid) %>%
  ggplot(aes(x = disB_ageonset, y= corrected_val)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('\nwith corrected values (by category and co-occurrence)')

snpp = ggarrange(rawp,correctedp,labels='auto', nrow = 2, ncol=1, common.legend = T, legend = 'right')
ggsave('./results/genomicAnalysis/Genesimilarity_cl2_bysameAge.pdf',snpp, units = 'cm', width = 16, height = 14, useDingbats = F)
ggsave('./results/genomicAnalysis/Genesimilarity_cl2_bysameAge.png',snpp, units = 'cm', width = 16, height = 14)

rawp = mydat %>%
  filter(disA_ageonset ==3) %>%
  ggplot(aes(x = disB_ageonset, y= gene_odds)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('Similarity between age-of-onset cl3\nand other clusters (with raw values)')

correctedp = mydat %>%
  filter(disA_ageonset ==3) %>%
  mutate(corrected_val = lm(gene_odds ~ sameCat + cooccur, data=.)$resid) %>%
  ggplot(aes(x = disB_ageonset, y= corrected_val)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('\nwith corrected values (by category and co-occurrence)')

snpp = ggarrange(rawp,correctedp,labels='auto', nrow = 2, ncol=1, common.legend = T, legend = 'right')
ggsave('./results/genomicAnalysis/Genesimilarity_cl3_bysameAge.pdf',snpp, units = 'cm', width = 16, height = 14, useDingbats = F)
ggsave('./results/genomicAnalysis/Genesimilarity_cl3_bysameAge.png',snpp, units = 'cm', width = 16, height = 14)

##


mydat = alldat %>%
  mutate(snp_odds = log2(snp_odds),
         RR = log2(RR),
         gene_odds = log2(gene_odds),
         proxy_gene_odds = log2(proxy_gene_odds),
         eqtl_gene_odds = log2(eqtl_gene_odds)) %>%
  # mutate(snpodds = ifelse(snp_p<0.01,snp_odds,1)) %>%
  filter(proxy_gene_p<0.01) %>%
  select(disA,disB,proxy_gene_odds,sameCat,cooccur,sameAge) 

a=ifelse(mydat$disA>mydat$disB,mydat$disA,mydat$disB)
b=ifelse(mydat$disA>mydat$disB,mydat$disB,mydat$disA)
mydat$disA = a
mydat$disB = b
nrow(mydat)
nrow(unique(mydat))
mydat = unique(mydat)

lmmod = lm(proxy_gene_odds ~ sameCat + cooccur + sameAge, data=mydat) 
aovobj = aov(lmmod)
summary(aovobj)
par(mfrow=c(1,2))
boxplot(proxy_gene_odds~sameAge,data=mydat,main = 'uncorrected')
boxplot(lm(proxy_gene_odds ~ sameCat + cooccur, data=mydat)$resid~sameAge,data=mydat,main = 'corrected')

rawp = mydat %>%
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(proxy_gene_odds ~ sameCat + cooccur, data=mydat)$resid) %>%
  ggplot(aes(x = sameAge, y= proxy_gene_odds)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_sina(color = 'gray40', size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('\nRaw values')

correctedp = mydat %>%
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(proxy_gene_odds ~ sameCat + cooccur, data=mydat)$resid) %>%
  ggplot(aes(x = sameAge, y= corrected_val)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_sina(color = 'gray40', size= 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('Corrected values\n(by category and co-occurrence)')

snpp = ggarrange(rawp,correctedp,labels='auto')
ggsave('./results/genomicAnalysis/ProxyGenesimilarity_bysameAge.pdf',snpp, units = 'cm', width = 18, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/ProxyGenesimilarity_bysameAge.png',snpp, units = 'cm', width = 18, height = 7)

mydat = alldat %>%
  mutate(snp_odds = log2(snp_odds),
         RR = log2(RR),
         gene_odds = log2(gene_odds),
         proxy_gene_odds = log2(proxy_gene_odds),
         eqtl_gene_odds = log2(eqtl_gene_odds)) %>%
  filter(proxy_gene_p<0.01) %>%
  select(disA,disB,proxy_gene_odds,sameCat,cooccur,sameAge) %>% 
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(proxy_gene_odds ~ sameCat + cooccur, data=.)$resid) %>%
  mutate(disA_ageonset = factor(ageonsetclusters$clustering[as.character(disA)]),
         disB_ageonset = factor(ageonsetclusters$clustering[as.character(disB)]),
         disA_category = disTreecl[as.character(disA)],
         disB_category = disTreecl[as.character(disB)])

rawp = mydat %>%
  filter(disA_ageonset ==1) %>%
  ggplot(aes(x = disB_ageonset, y= proxy_gene_odds)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('Similarity between age-of-onset cl1\nand other clusters (with raw values)')

correctedp = mydat %>%
  filter(disA_ageonset ==1) %>%
  mutate(corrected_val = lm(proxy_gene_odds ~ sameCat + cooccur, data=.)$resid) %>%
  ggplot(aes(x = disB_ageonset, y= corrected_val)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('\nwith corrected values (by category and co-occurrence)')

snpp = ggarrange(rawp,correctedp,labels='auto', nrow = 2, ncol=1, common.legend = T, legend = 'right')
ggsave('./results/genomicAnalysis/ProxyGenesimilarity_cl1_bysameAge.pdf',snpp, units = 'cm', width = 16, height = 14, useDingbats = F)
ggsave('./results/genomicAnalysis/ProxyGenesimilarity_cl1_bysameAge.png',snpp, units = 'cm', width = 16, height = 14)


rawp = mydat %>%
  filter(disA_ageonset ==2) %>%
  ggplot(aes(x = disB_ageonset, y= proxy_gene_odds)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('Similarity between age-of-onset cl2\nand other clusters (with raw values)')

correctedp = mydat %>%
  filter(disA_ageonset ==2) %>%
  mutate(corrected_val = lm(proxy_gene_odds ~ sameCat + cooccur, data=.)$resid) %>%
  ggplot(aes(x = disB_ageonset, y= corrected_val)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('\nwith corrected values (by category and co-occurrence)')

snpp = ggarrange(rawp,correctedp,labels='auto', nrow = 2, ncol=1, common.legend = T, legend = 'right')
ggsave('./results/genomicAnalysis/ProxyGenesimilarity_cl2_bysameAge.pdf',snpp, units = 'cm', width = 16, height = 14, useDingbats = F)
ggsave('./results/genomicAnalysis/ProxyGenesimilarity_cl2_bysameAge.png',snpp, units = 'cm', width = 16, height = 14)

rawp = mydat %>%
  filter(disA_ageonset ==3) %>%
  ggplot(aes(x = disB_ageonset, y= proxy_gene_odds)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('Similarity between age-of-onset cl3\nand other clusters (with raw values)')

correctedp = mydat %>%
  filter(disA_ageonset ==3) %>%
  mutate(corrected_val = lm(proxy_gene_odds ~ sameCat + cooccur, data=.)$resid) %>%
  ggplot(aes(x = disB_ageonset, y= corrected_val)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('\nwith corrected values (by category and co-occurrence)')

snpp = ggarrange(rawp,correctedp,labels='auto', nrow = 2, ncol=1, common.legend = T, legend = 'right')
ggsave('./results/genomicAnalysis/ProxyGenesimilarity_cl3_bysameAge.pdf',snpp, units = 'cm', width = 16, height = 14, useDingbats = F)
ggsave('./results/genomicAnalysis/ProxyGenesimilarity_cl3_bysameAge.png',snpp, units = 'cm', width = 16, height = 14)

##


mydat = alldat %>%
  mutate(snp_odds = log2(snp_odds),
         RR = log2(RR),
         gene_odds = log2(gene_odds),
         proxy_gene_odds = log2(proxy_gene_odds),
         eqtl_gene_odds = log2(eqtl_gene_odds)) %>%
  # mutate(snpodds = ifelse(snp_p<0.01,snp_odds,1)) %>%
  filter(eqtl_gene_p<0.01) %>%
  select(disA,disB,eqtl_gene_odds,sameCat,cooccur,sameAge) 

a=ifelse(mydat$disA>mydat$disB,mydat$disA,mydat$disB)
b=ifelse(mydat$disA>mydat$disB,mydat$disB,mydat$disA)
mydat$disA = a
mydat$disB = b
nrow(mydat)
nrow(unique(mydat))
mydat = unique(mydat)

lmmod = lm(eqtl_gene_odds ~ sameCat + cooccur + sameAge, data=mydat) 
aovobj = aov(lmmod)
summary(aovobj)
par(mfrow=c(1,2))
boxplot(eqtl_gene_odds~sameAge,data=mydat,main = 'uncorrected')
boxplot(lm(eqtl_gene_odds ~ sameCat + cooccur, data=mydat)$resid~sameAge,data=mydat,main = 'corrected')

rawp = mydat %>%
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(eqtl_gene_odds ~ sameCat + cooccur, data=mydat)$resid) %>%
  ggplot(aes(x = sameAge, y= eqtl_gene_odds)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_sina(color = 'gray40', size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('\nRaw values')

correctedp = mydat %>%
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(eqtl_gene_odds ~ sameCat + cooccur, data=mydat)$resid) %>%
  ggplot(aes(x = sameAge, y= corrected_val)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_sina(color = 'gray40', size= 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('Corrected values\n(by category and co-occurrence)')

snpp = ggarrange(rawp,correctedp,labels='auto')
ggsave('./results/genomicAnalysis/eQTLGenesimilarity_bysameAge.pdf',snpp, units = 'cm', width = 18, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/eQTLGenesimilarity_bysameAge.png',snpp, units = 'cm', width = 18, height = 7)

mydat = alldat %>%
  mutate(snp_odds = log2(snp_odds),
         RR = log2(RR),
         gene_odds = log2(gene_odds),
         proxy_gene_odds = log2(proxy_gene_odds),
         eqtl_gene_odds = log2(eqtl_gene_odds)) %>%
  filter(eqtl_gene_p<0.01) %>%
  select(disA,disB,eqtl_gene_odds,sameCat,cooccur,sameAge) %>% 
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(eqtl_gene_odds ~ sameCat + cooccur, data=.)$resid) %>%
  mutate(disA_ageonset = factor(ageonsetclusters$clustering[as.character(disA)]),
         disB_ageonset = factor(ageonsetclusters$clustering[as.character(disB)]),
         disA_category = disTreecl[as.character(disA)],
         disB_category = disTreecl[as.character(disB)])

rawp = mydat %>%
  filter(disA_ageonset ==1) %>%
  ggplot(aes(x = disB_ageonset, y= eqtl_gene_odds)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('Similarity between age-of-onset cl1\nand other clusters (with raw values)')

correctedp = mydat %>%
  filter(disA_ageonset ==1) %>%
  mutate(corrected_val = lm(eqtl_gene_odds ~ sameCat + cooccur, data=.)$resid) %>%
  ggplot(aes(x = disB_ageonset, y= corrected_val)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('\nwith corrected values (by category and co-occurrence)')

snpp = ggarrange(rawp,correctedp,labels='auto', nrow = 2, ncol=1, common.legend = T, legend = 'right')
ggsave('./results/genomicAnalysis/eQTLGenesimilarity_cl1_bysameAge.pdf',snpp, units = 'cm', width = 16, height = 14, useDingbats = F)
ggsave('./results/genomicAnalysis/eQTLGenesimilarity_cl1_bysameAge.png',snpp, units = 'cm', width = 16, height = 14)


rawp = mydat %>%
  filter(disA_ageonset ==2) %>%
  ggplot(aes(x = disB_ageonset, y= eqtl_gene_odds)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('Similarity between age-of-onset cl2\nand other clusters (with raw values)')

correctedp = mydat %>%
  filter(disA_ageonset ==2) %>%
  mutate(corrected_val = lm(eqtl_gene_odds ~ sameCat + cooccur, data=.)$resid) %>%
  ggplot(aes(x = disB_ageonset, y= corrected_val)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('\nwith corrected values (by category and co-occurrence)')

snpp = ggarrange(rawp,correctedp,labels='auto', nrow = 2, ncol=1, common.legend = T, legend = 'right')
ggsave('./results/genomicAnalysis/eQTLGenesimilarity_cl2_bysameAge.pdf',snpp, units = 'cm', width = 16, height = 14, useDingbats = F)
ggsave('./results/genomicAnalysis/eQTLGenesimilarity_cl2_bysameAge.png',snpp, units = 'cm', width = 16, height = 14)

rawp = mydat %>%
  filter(disA_ageonset ==3) %>%
  ggplot(aes(x = disB_ageonset, y= eqtl_gene_odds)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('Similarity between age-of-onset cl3\nand other clusters (with raw values)')

correctedp = mydat %>%
  filter(disA_ageonset ==3) %>%
  mutate(corrected_val = lm(eqtl_gene_odds ~ sameCat + cooccur, data=.)$resid) %>%
  ggplot(aes(x = disB_ageonset, y= corrected_val)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), size = 1) +
  geom_jitter(aes(color = cooccur), size = 0.5) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('\nwith corrected values (by category and co-occurrence)')

snpp = ggarrange(rawp,correctedp,labels='auto', nrow = 2, ncol=1, common.legend = T, legend = 'right')
ggsave('./results/genomicAnalysis/eQTLGenesimilarity_cl3_bysameAge.pdf',snpp, units = 'cm', width = 16, height = 14, useDingbats = F)
ggsave('./results/genomicAnalysis/eQTLGenesimilarity_cl3_bysameAge.png',snpp, units = 'cm', width = 16, height = 14)
