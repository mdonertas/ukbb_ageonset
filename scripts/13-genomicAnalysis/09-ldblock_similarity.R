#run 03-SNPanalysis.R before line 183
LDblocks = read_tsv('./data/raw/LDblocks1000G_EUR.bed')%>%
  mutate(CHR=as.numeric(gsub('chr','',chr))) %>%
  select(-chr)
xx=t(sapply(strsplit(rownames(snpmat),'_'),function(x)c(paste(x,collapse='_'),x[1],x[2])))
xx = apply(LDblocks,1,function(x){
  st=x['start']
  en=x['stop']
  ch=x['CHR']
  chr = as.numeric(xx[,2])
  bp= as.numeric(xx[,3])
  xx[which(chr==ch & bp<=en & bp>=st),1]
})
xx=xx[which(sapply(xx,length)>0)]
ldmat = t(sapply(xx,function(i){
  if(length(i)>1){
    as.numeric(colSums(snpmat[i,])>0) 
  } else if(length(i)==1){
    snpmat[i,]
  } 
}))
colnames(ldmat)=colnames(snpmat)
l=nrow(ldmat)
odds = apply(ldmat,2,function(x){
  apply(ldmat,2,function(y){
    obs = sum((x == y) & (x != 0))
    p = (sum(x!=0)/ l) * (sum(y!=0)/ l) 
    bit = binom.test(x = obs, n = l, p = p, alternative = 'g')
    bit$estimate/(p)
  })
})
diag(odds) = NA
p = apply(ldmat,2,function(x){
  apply(ldmat,2,function(y){
    obs = sum((x == y) & (x != 0))
    p = (sum(x!=0)/ l) * (sum(y!=0)/ l) 
    bit = binom.test(x = obs, n = l, p = p, alternative = 'greater')
    bit$p.value
  })
})

saveRDS(odds, './data/processed/genomicAnalysis/ldblock_odds.rds')
saveRDS(p, './data/processed/genomicAnalysis/ldblock_p.rds')
odds=readRDS('./data/processed/genomicAnalysis/ldblock_odds.rds')
p=readRDS('./data/processed/genomicAnalysis/ldblock_p.rds')
odds[p>0.01] = 1
l2odds = log2(odds)
i = which(rowMeans(l2odds==0, na.rm =T) != 1)
l2odds = l2odds[i,i]
# disannot = data.frame(diseaseCategories = unname(disTreecl[numSignif$disease]),
                      # ageonset_clusters = as.character(unname(readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$clustering[numSignif$disease])))
# rownames(disannot) = numSignif$disease
# disannot = disannot %>%
  # mutate(disease = rownames(disannot))
# annotcolors = list(diseaseCategories = discatcolors,ageonset_clusters = ageonsetcolors)

rownames(disannot) = disannot$disease

pheatmap::pheatmap(l2odds, color = colorRampPalette(brewer.pal(8,'Reds'))(20),cellwidth = 5, cellheight = 5, cutree_rows = 10, cutree_cols = 10, annotation_row = disannot[,-3], annotation_col = disannot[,-3], annotation_colors = annotcolors, fontsize = 6, filename = './results/genomicAnalysis/signifSNPoverlap_heatmap_ldblock.pdf')

pheatmap::pheatmap(l2odds, color = colorRampPalette(brewer.pal(8,'Reds'))(20),cellwidth = 5, cellheight = 5, cutree_rows = 10, cutree_cols = 10, annotation_row = disannot[,-3], annotation_col = disannot[,-3], annotation_colors = annotcolors, fontsize = 6, filename = './results/genomicAnalysis/signifSNPoverlap_heatmap_ldblock.png')

l2odds = log2(odds)
# i = which(rowMeans(l2odds==0, na.rm =T) != 1)
# l2odds = l2odds[i,i]
rownames(disannot) = disannot$disease
pheatmap::pheatmap(l2odds, color = colorRampPalette(brewer.pal(8,'Reds'))(20),cellwidth = 5, cellheight = 5, annotation_row = disannot[,-3], annotation_col = disannot[,-3], annotation_colors = annotcolors, fontsize = 6, filename = './results/genomicAnalysis/signifSNPoverlap_heatmap_all_ldblock.pdf')

pheatmap::pheatmap(l2odds, color = colorRampPalette(brewer.pal(8,'Reds'))(20),cellwidth = 5, cellheight = 5, annotation_row = disannot[,-3], annotation_col = disannot[,-3], annotation_colors = annotcolors, fontsize = 6, filename = './results/genomicAnalysis/signifSNPoverlap_heatmap_all_ldblock.png')

xx = lapply(as.character(1:4),function(i){
  gr1 = rownames(disannot)[which(disannot$ageonset_clusters==i)]
  gr1 = intersect(gr1,rownames(l2odds))
  xx = lapply(as.character(1:4),function(k){
    gr2 = rownames(disannot)[which(disannot$ageonset_clusters==k)]
    gr2 = intersect(gr2,rownames(l2odds))
    l2odds[gr1,gr2]
  })
  names(xx) = paste('cl',1:4,sep='')
  xx
})
names(xx) = paste('cl',1:4,sep='')
xx = reshape2::melt(xx) %>%
  set_names(c('disease1','disease2','value','cluster1','cluster2')) %>%
  na.omit()
# a=ifelse(xx$cluster1<=xx$cluster2,xx$cluster1,xx$cluster2)
# b=ifelse(xx$cluster1<=xx$cluster2,xx$cluster2,xx$cluster1)
# 
# xx$cluster1 = a
# xx$cluster2 = b
# xx$disease1=as.character(xx$disease1)
# xx$disease2=as.character(xx$disease2)
# a=ifelse(xx$disease1<=xx$disease2,xx$disease1,xx$disease2)
# b=ifelse(xx$disease1<=xx$disease2,xx$disease2,xx$disease1)
# 
# xx$disease1 = a
# xx$disease2 = b

xx = xx %>%
  mutate(val = ifelse(value==0,NA,value))

ggplot(xx,aes(x=cluster1,y=val)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(cluster2~.) +
  geom_jitter(aes(y=value), width = 0.3, size =0.3, alpha=0.3) +
  theme_bw()

library(igraph)
library(ggnetwork)
library(ggnet)
library(ggforce)
mynet = graph_from_data_frame(select(xx,disease1,disease2,val) %>%na.omit() %>%unique())
V(mynet)$ageonset_cl = setNames(disannot$ageonset_clusters,rownames(disannot))[V(mynet)$name]
V(mynet)$discat = as.character(setNames(disannot$diseaseCategories,rownames(disannot))[V(mynet)$name])
E(mynet)$edgewidth = summary(E(mynet)$val/30)
net0 = ggnet2(mynet, size = 0, edge.size = 'edgewidth', edge.color = 'gray80')+
  theme_void()

net1 = net0 +
  geom_point(size = 2, color = ageonsetcolors[V(mynet)$ageonset_cl])+
  theme_void()
net2 = net0 +
  geom_point(size = 2, color = discatcolors[V(mynet)$discat])+
  theme_void()

netp = ggarrange(net1,net2, labels = c('a','b'))

ggsave('./results/genomicAnalysis/SNPnet_ldblock.pdf', netp, units = 'cm', width = 18, height = 8, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPnet_ldblock.png', netp, units = 'cm', width = 18, height = 8)

###

source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)

snpsim_odds = readRDS('./data/processed/genomicAnalysis/ldblock_odds.rds')
snpsim_p = readRDS('./data/processed/genomicAnalysis/ldblock_p.rds')

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
         log10p=-log10(snp_p)) %>%
  # mutate(snpodds = ifelse(snp_p<0.01,snp_odds,1)) %>%
  filter(snp_p<=0.01) %>%
  select(disA,disB,snp_odds,sameCat,cooccur,sameAge,log10p) 

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
net0 = ggnet2(mynet, size = 0, edge.color = 'gray80') +
  theme_void()
net1 = net0 +
  geom_point(size = 2, color = ageonsetcolors[V(mynet)$ageonset_cl]) +
  theme_void()
net2 = net0 +
  geom_point(size = 2, color = discatcolors[V(mynet)$discat])+
  theme_void()

ggsave('./results/genomicAnalysis/SNPnet_ageonsetcl_ldblock.pdf', net1, units = 'cm', width = 18, height = 8, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPnet_ageonsetcl_ldblock.png', net1, units = 'cm', width = 18, height = 8)

ggsave('./results/genomicAnalysis/SNPnet_discat_ldblock.pdf', net2, units = 'cm', width = 18, height = 8, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPnet_discat_ldblock.png', net2, units = 'cm', width = 18, height = 8)

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
  geom_point(size = 2, color = ageonsetcolors[V(mynet)$ageonset_cl])+
  theme_void()
net2 = net0 +
  geom_point(size = 2, color = discatcolors[V(mynet)$discat])+
  theme_void()

ggsave('./results/genomicAnalysis/SNPnet_ageonsetcl_corr_ldblock.pdf', net1, units = 'cm', width = 18, height = 8, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPnet_ageonsetcl_corr_ldblock.png', net1, units = 'cm', width = 18, height = 8)

ggsave('./results/genomicAnalysis/SNPnet_discat_corr_ldblock.pdf', net2, units = 'cm', width = 18, height = 8, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPnet_discat_corr_ldblock.png', net2, units = 'cm', width = 18, height = 8)


lmmod = lm(snp_odds ~ sameCat + cooccur + sameAge, data=mydat) 
aovobj = aov(lmmod)
summary(aovobj)
par(mfrow=c(1,2))
boxplot(log10p~sameAge,data=mydat,main = 'uncorrected')
boxplot(lm(log10p ~ sameCat + cooccur, data=mydat)$resid~sameAge,data=mydat,main = 'corrected')

rawp = mydat %>%
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(snp_odds ~ sameCat + cooccur, data=mydat)$resid) %>%
  ggplot(aes(x = sameAge, y= snp_odds)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, color = 'gray40') +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('\nRaw values')


correctedp = mydat %>%
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) %>%
  mutate(corrected_val = lm(log10p ~ sameCat + cooccur, data=mydat)$resid) %>%
  ggplot(aes(x = sameAge, y= corrected_val)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = 'gray40', size= 0.5, width = 0.2) +
  ylab('Genetic Similarity\nlog2(odds ratio)') + xlab('') +
  ggtitle('Corrected values\n(by category and co-occurrence)') 


snpp = ggarrange(rawp,correctedp,labels='auto')
ggsave('./results/genomicAnalysis/SNPsimilarity_bysameAge_ldblock.pdf',snpp, units = 'cm', width = 18, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPsimilarity_bysameAge_ldblock.png',snpp, units = 'cm', width = 18, height = 7)

netandbox = ggarrange(net1,correctedp,labels='auto')
ggsave('./results/genomicAnalysis/SNPsimilarity_age_netbox_ldblock.pdf',netandbox, units = 'cm', width = 18, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPsimilarity_age_netbox_ldblock.png',netandbox, units = 'cm', width = 18, height = 7)
