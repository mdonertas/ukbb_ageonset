source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)

#check if all significant snps are in signif_gwasRes_proxyGenes.rds even if they cannot be mapped:
proxyall = readRDS('./data/processed/caseControl/a1071/gwasRes_proxyGenes.rds')
proxysignif = readRDS('./data/processed/caseControl/a1071/signif_gwasRes_proxyGenes.rds')
length(unique(filter(proxyall,P_BOLT_LMM_INF<=5e-8)$SNP))
length(unique(proxysignif$SNP))
# yes so continue with only signif_gwasRes_proxyGenes.rds
rm(proxyall,proxysignif)

allSNPs = read_tsv('./data/processed/ukbb/gwas/bolt/a1071.imp.stats')
totSNPnum = nrow(allSNPs)

signifSNPs <- sapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_proxyGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  setdiff(unique(x$SNP),c('',NA))
})
names(signifSNPs)=disIDs

numSignif = sapply(signifSNPs,length)

numSignif = data.frame(numSignif) %>%
  mutate(disease = disCoding[names(signifSNPs)])

numSignif_perDisease = numSignif %>%
  ggplot(aes(x = numSignif)) +
  geom_histogram(bins=50) +
  xlab('# of significant (p<=5e-8) SNPs')

xx = numSignif %>%
  filter(numSignif >=10000) %>%
  select(disease, numSignif) %>%
  arrange(-numSignif) %>%
  rename(`#` = numSignif, Disease=disease)
text.p= ggtexttable(xx, rows = NULL,theme = ttheme('mOrange',base_size = 6,padding = unit(c(2,2),'mm')))
numSignif_perDisease = numSignif_perDisease + annotation_custom(ggplotGrob(text.p),
                                                                xmin = 5000, ymin = 20,
                                                                xmax = 35000)

disannot = data.frame(diseaseCategories = unname(disTreecl[numSignif$disease]),
                      ageonset_clusters = as.character(unname(readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$clustering[numSignif$disease])))
rownames(disannot) = numSignif$disease
disannot = disannot %>%
  mutate(disease = rownames(disannot))
ageonsetcolors = setNames(rev(brewer.pal(4,'Oranges')),1:4)
annotcolors = list(diseaseCategories = discatcolors,ageonset_clusters = setNames(rev(brewer.pal(4,'Oranges')),1:4))

numSignif_by_ageonset = numSignif %>%
  left_join(disannot) %>%
  ggplot(aes(x = ageonset_clusters, y = numSignif)) +
  geom_boxplot(outlier.shape = NA, aes(fill = ageonset_clusters)) +
  geom_sina( size = 0.5) +
  scale_y_log10(breaks = c(10,100,1000), labels = c(10,100,1000)) +
  xlab('') +
  ylab('# of significant (p<=5e-8) SNPs') +
  scale_fill_manual(values = ageonsetcolors, guide = F) +
  ggtitle('# of SNPs by\nAge of Onset Clusters')

numSignif_by_cat = numSignif %>%
  left_join(disannot) %>%
  mutate(val = (numSignif)) %>%
  mutate(diseaseCategories = fct_reorder(diseaseCategories,val,.fun = function(x)median(x[which(x!=0)]))) %>%
  mutate(diseaseCategories = factor(diseaseCategories, levels = rev(levels(diseaseCategories)))) %>%
  ggplot(aes(x = diseaseCategories, y = numSignif)) +
  geom_boxplot(outlier.shape = NA, aes(fill = diseaseCategories)) +
  geom_sina(size = 0.5) +
  scale_y_log10(breaks = c(10,100,1000), labels = c(10,100,1000)) +
  xlab('') +
  ylab('# of significant (p<=5e-8) SNPs') +
  scale_fill_manual(values = discatcolors, guide = guide_legend('category')) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = 'right') +
  ggtitle('# of SNPs by Disease Categories')

p_numbersnps = ggarrange(ggarrange(numSignif_perDisease,numSignif_by_ageonset,labels = c('a','b'), nrow = 1, ncol =2, align = 'hv'),
                         numSignif_by_cat,labels = c('','c'), nrow = 2, ncol =1, heights = c(1,1))
ggsave('./results/genomicAnalysis/numberSNPs.pdf',p_numbersnps,units = 'cm',width = 18,height = 14,useDingbats = F)
ggsave('./results/genomicAnalysis/numberSNPs.png',p_numbersnps,units = 'cm',width = 18,height = 14)

ggsave('./results/genomicAnalysis/numSNP_perDisease.pdf', numSignif_perDisease, units = 'cm', width = 8, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/numSNP_perDisease.png', numSignif_perDisease, units = 'cm', width = 8, height = 7)

ggsave('./results/genomicAnalysis/numSNP_perDisease_byAgeOnset.pdf', numSignif_by_ageonset, units = 'cm', width = 8, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/numSNP_perDisease_byAgeOnset.png', numSignif_by_ageonset, units = 'cm', width = 8, height = 7)

ggsave('./results/genomicAnalysis/numSNP_perDisease_byCategory.pdf', numSignif_by_cat, units = 'cm', width = 16, height = 8, useDingbats = F)
ggsave('./results/genomicAnalysis/numSNP_perDisease_byCategory.png', numSignif_by_cat, units = 'cm', width = 16, height = 8)

summary(numSignif$numSignif)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     0.0    13.5  1389.8   737.0 35008.0 
sum(numSignif$numSignif == 0)
# [1] 36

anova(lm(numSignif ~ ageonset_clusters, data = left_join(numSignif,disannot)))
tapply(X = left_join(numSignif,disannot)$numSignif, INDEX = left_join(numSignif,disannot)$diseaseCategories,summary)
### 

signifSNPs <- lapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_proxyGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  unique(select(x,CHR,BP,Ref,Alt,BETA)) %>%
    mutate(SNP = paste(CHR,BP,Ref,Alt,sep='_')) %>%
    select(SNP,BETA)
})
names(signifSNPs)=disIDs

signifSNPs = reshape2::melt(signifSNPs,id.vars = colnames(signifSNPs[[1]]))  %>% 
  rename(disCode = L1) %>%
  mutate(disease = disCoding[disCode]) %>%
  mutate(value = sign(BETA))

SNPmat = signifSNPs %>% 
  select(SNP,disease,value) %>%
  spread(disease,value,fill=0)
rownames(SNPmat) = SNPmat$SNP
SNPmat$SNP =NULL
snpmat = as.matrix(SNPmat)
dim(snpmat)
# [1] 93879    80
table(rowSums(snpmat!=0)>1)
# FALSE  TRUE 
# 50043 43836
mean((rowSums(snpmat!=0)>1))
# [1] 0.4669415

overview = t(apply(snpmat[which(rowSums(snpmat!=0)>1),],1,function(x){
  pos = sum(x==1)
  neg = sum(x==-1)
  c(pos,neg)
}))
nrow(overview[which((overview[,1] != 0) & (overview[,2] != 0)),] )
# 2240
nrow(overview[which((overview[,1] != 0) & (overview[,2] != 0)),] )/ sum(rowSums(snpmat!=0)>1)
# 0.05109955

cat_per_snp = apply(snpmat,1,function(x){
  disTreecl[names(which(x!=0))]
})

cat_per_snp_l=sapply(cat_per_snp,function(x)length(unique(unname(x))))

summary(cat_per_snp_l)
maxsnp = strsplit(names(which.max(rowSums(snpmat!=0))),'_')[[1]]
(allSNPs %>% filter(CHR == maxsnp[1] & BP == maxsnp[2] & ALLELE1 == maxsnp[3] & ALLELE0 == maxsnp[4]))$SNP
which(snpmat[which.max(rowSums(snpmat!=0)),]!=0)

dispersnp = data.frame(num = rowSums(snpmat!=0)) %>%
  group_by(num) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = num, y = count)) +
  geom_bar(stat = 'identity') +
  geom_label(aes(label = comma(count)), size = 6/pntnorm) +
  xlab('Number of diseases per SNP (p<=5e-8)') +
  ylab('Count')

ggsave('./results/genomicAnalysis/diseaseperSNP.pdf',dispersnp,units = 'cm', width = 8, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/diseaseperSNP.png',dispersnp,units = 'cm', width = 8, height = 7)

catpersnp = data.frame(num = cat_per_snp_l) %>%
  group_by(num) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = num, y = count)) +
  geom_bar(stat = 'identity') +
  geom_label(aes(label = comma(count)), size = 6/pntnorm) +
  xlab('Number of categories per SNP (p<=5e-8)') +
  ylab('Count')

ggsave('./results/genomicAnalysis/categoryperSNP.pdf',catpersnp,units = 'cm', width = 8, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/categoryperSNP.png',catpersnp,units = 'cm', width = 8, height = 7)


persnpplots = ggarrange(dispersnp,catpersnp,labels = c('a','b'),align = 'hv')
ggsave('./results/genomicAnalysis/perSNP.pdf',persnpplots,units = 'cm', width = 16, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/perSNP.png',persnpplots,units = 'cm', width = 16, height = 7)

mean(cat_per_snp_l==1)

odds = apply(snpmat,2,function(x){
  apply(snpmat,2,function(y){
    obs = sum((x == y) & (x != 0))
    p = (sum(x!=0)/ totSNPnum) * (sum(y!=0)/ totSNPnum) 
    bit = binom.test(x = obs, n = totSNPnum, p = p, alternative = 'g')
    bit$estimate/(p)
  })
})
diag(odds) = NA
p = apply(snpmat,2,function(x){
  apply(snpmat,2,function(y){
    obs = sum((x == y) & (x != 0))
    p = (sum(x!=0)/ totSNPnum) * (sum(y!=0)/ totSNPnum) 
    bit = binom.test(x = obs, n = totSNPnum, p = p, alternative = 'greater')
    bit$p.value
  })
})
odds[p>0.01] = 1
l2odds = log2(odds)
i = which(rowMeans(l2odds==0, na.rm =T) != 1)
l2odds = l2odds[i,i]
rownames(disannot) = disannot$disease
pheatmap::pheatmap(l2odds, color = colorRampPalette(brewer.pal(8,'Reds'))(20),cellwidth = 5, cellheight = 5, cutree_rows = 4, cutree_cols = 4, annotation_row = disannot[,-3], annotation_col = disannot[,-3], annotation_colors = annotcolors, fontsize = 6, filename = './results/genomicAnalysis/signifSNPoverlap_heatmap.pdf')

pheatmap::pheatmap(l2odds, color = colorRampPalette(brewer.pal(8,'Reds'))(20),cellwidth = 5, cellheight = 5, cutree_rows = 4, cutree_cols = 4, annotation_row = disannot[,-3], annotation_col = disannot[,-3], annotation_colors = annotcolors, fontsize = 6, filename = './results/genomicAnalysis/signifSNPoverlap_heatmap.png')

l2odds = log2(odds)
# i = which(rowMeans(l2odds==0, na.rm =T) != 1)
# l2odds = l2odds[i,i]
rownames(disannot) = disannot$disease
pheatmap::pheatmap(l2odds, color = colorRampPalette(brewer.pal(8,'Reds'))(20),cellwidth = 5, cellheight = 5, cutree_rows = 4, cutree_cols = 4, annotation_row = disannot[,-3], annotation_col = disannot[,-3], annotation_colors = annotcolors, fontsize = 6, filename = './results/genomicAnalysis/signifSNPoverlap_heatmap_all.pdf')

pheatmap::pheatmap(l2odds, color = colorRampPalette(brewer.pal(8,'Reds'))(20),cellwidth = 5, cellheight = 5, cutree_rows = 4, cutree_cols = 4, annotation_row = disannot[,-3], annotation_col = disannot[,-3], annotation_colors = annotcolors, fontsize = 6, filename = './results/genomicAnalysis/signifSNPoverlap_heatmap_all.png')

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
net0 = ggnet2(mynet, size = 0, edge.size = 'edgewidth', edge.color = 'gray80')

net1 = net0 +
  geom_point(size = 2, color = ageonsetcolors[V(mynet)$ageonset_cl])
net2 = net0 +
  geom_point(size = 2, color = discatcolors[V(mynet)$discat])

netp = ggarrange(net1,net2, labels = c('a','b'))

ggsave('./results/genomicAnalysis/SNPnet.pdf', netp, units = 'cm', width = 18, height = 8, useDingbats = F)
ggsave('./results/genomicAnalysis/SNPnet.png', netp, units = 'cm', width = 18, height = 8)


