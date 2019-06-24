source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)

signifGenes <- sapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_eQTLGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  setdiff(unique(x$eQTL_hgnc),c('',NA))
})
names(signifGenes)=disIDs

numSignif = sapply(signifGenes,length)

numSignif = data.frame(numSignif) %>%
  mutate(disease = disCoding[names(signifGenes)])

numSignif_perDisease = numSignif %>%
  ggplot(aes(x = numSignif)) +
  geom_histogram(bins=50) +
  xlab('# of Genes')

xx = numSignif %>%
  filter(numSignif >=500) %>%
  select(disease, numSignif) %>%
  arrange(-numSignif) %>%
  rename(`#` = numSignif, Disease=disease)
text.p= ggtexttable(xx, rows = NULL,theme = ttheme('mOrange',base_size = 6,padding = unit(c(2,2),'mm')))
numSignif_perDisease = numSignif_perDisease + annotation_custom(ggplotGrob(text.p),
                                                                xmin = 100, ymin = 20,
                                                                xmax = 1000)

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
  ylab('# of Genes') +
  scale_fill_manual(values = ageonsetcolors, guide = F) +
  ggtitle('# of Genes by\nAge of Onset Clusters')

numSignif_by_cat = numSignif %>%
  left_join(disannot) %>%
  mutate(val = (numSignif)) %>%
  mutate(diseaseCategories = fct_reorder(diseaseCategories,val,.fun = function(x)median(x[which(x!=0)]))) %>%
  mutate(diseaseCategories = factor(diseaseCategories, levels = rev(levels(diseaseCategories)))) %>%
  na.omit() %>%
  filter(numSignif !=0)%>%
  ggplot(aes(x = diseaseCategories, y = numSignif)) +
  geom_boxplot(outlier.shape = NA, aes(fill = diseaseCategories)) +
  geom_sina(size = 0.5) +
  scale_y_log10(breaks = c(10,100,1000), labels = c(10,100,1000)) +
  xlab('') +
  ylab('# of Genes') +
  scale_fill_manual(values = discatcolors, guide = guide_legend('category')) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = 'right') +
  ggtitle('# of Genes by Disease Categories')

p_numbergenes = ggarrange(ggarrange(numSignif_perDisease,numSignif_by_ageonset,labels = c('a','b'), nrow = 1, ncol =2, align = 'hv'),
                          numSignif_by_cat,labels = c('','c'), nrow = 2, ncol =1, heights = c(1,1))
ggsave('./results/genomicAnalysis/numbereQTLGenes.pdf',p_numbergenes,units = 'cm',width = 18,height = 14,useDingbats = F)
ggsave('./results/genomicAnalysis/numbereQTLGenes.png',p_numbergenes,units = 'cm',width = 18,height = 14)

ggsave('./results/genomicAnalysis/numeQTLGene_perDisease.pdf', numSignif_perDisease, units = 'cm', width = 8, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/numeQTLGene_perDisease.png', numSignif_perDisease, units = 'cm', width = 8, height = 7)

ggsave('./results/genomicAnalysis/numeQTLGene_perDisease_byAgeOnset.pdf', numSignif_by_ageonset, units = 'cm', width = 8, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/numeQTLGene_perDisease_byAgeOnset.png', numSignif_by_ageonset, units = 'cm', width = 8, height = 7)

ggsave('./results/genomicAnalysis/numeQTLGene_perDisease_byCategory.pdf', numSignif_by_cat, units = 'cm', width = 16, height = 8, useDingbats = F)
ggsave('./results/genomicAnalysis/numeQTLGene_perDisease_byCategory.png', numSignif_by_cat, units = 'cm', width = 16, height = 8)

summary(numSignif$numSignif)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    1.00   35.91   21.50  743.00
sum(numSignif$numSignif == 0)
# [1] 53

anova(lm(numSignif ~ ageonset_clusters, data = left_join(numSignif,disannot)))
tapply(X = left_join(numSignif,disannot)$numSignif, INDEX = left_join(numSignif,disannot)$diseaseCategories,summary)
### 
genemat = reshape2::melt(signifGenes) %>%
  rename(gene = value, disID = L1) %>%
  mutate(val = 1) %>%
  spread(disID, val, fill = 0)
rownames(genemat) = genemat$gene
genemat$gene= NULL
genemat = as.matrix(genemat)
dim(genemat)
# [1] 1763   63
table(rowSums(genemat!=0)>1)
# FALSE  TRUE 
# 648  1115
mean((rowSums(genemat!=0)>1))
# [1] 0.6324447

colnames(genemat) = disCoding[as.character(colnames(genemat))]
cat_per_gene = apply(genemat,1,function(x){
  disTreecl[names(which(x!=0))]
})

cat_per_gene_l=sapply(cat_per_gene,function(x)length(unique(unname(x))))

summary(cat_per_gene_l)
which(genemat[which.max(rowSums(genemat!=0)),]!=0)

dispergene = data.frame(num = rowSums(genemat!=0)) %>%
  group_by(num) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = num, y = count)) +
  geom_bar(stat = 'identity') +
  geom_label(aes(label = comma(count)), size = 6/pntnorm) +
  xlab('Number of diseases per Gene') +
  ylab('Count')

ggsave('./results/genomicAnalysis/diseasepereQTLGene.pdf',dispergene,units = 'cm', width = 8, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/diseasepereQTLGene.png',dispergene,units = 'cm', width = 8, height = 7)

catpergene = data.frame(num = cat_per_gene_l) %>%
  group_by(num) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = num, y = count)) +
  geom_bar(stat = 'identity') +
  geom_label(aes(label = comma(count)), size = 6/pntnorm) +
  xlab('Number of categories per Gene') +
  ylab('Count')

ggsave('./results/genomicAnalysis/categorypereQTLGene.pdf',catpergene,units = 'cm', width = 8, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/categorypereQTLGene.png',catpergene,units = 'cm', width = 8, height = 7)


pergeneplots = ggarrange(dispergene,catpergene,labels = c('a','b'),align = 'hv')
ggsave('./results/genomicAnalysis/pereQTLGene.pdf',pergeneplots,units = 'cm', width = 16, height = 7, useDingbats = F)
ggsave('./results/genomicAnalysis/pereQTLGene.png',pergeneplots,units = 'cm', width = 16, height = 7)

mean(cat_per_gene_l==1)
totgenes = readRDS('./data/processed/caseControl/a1071/gwasRes_eQTLGenes.rds')
totgenes= setdiff(unique(totgenes$eQTL_hgnc),c('',NA))
totGeneNum = length(totgenes)
odds = apply(genemat,2,function(x){
  apply(genemat,2,function(y){
    obs = sum((x == y) & (x != 0))
    p = (sum(x!=0)/  totGeneNum) * (sum(y!=0)/ totGeneNum) 
    bit = binom.test(x = obs, n = totGeneNum, p = p, alternative = 'g')
    bit$estimate/(p)
  })
})
saveRDS(odds,'./data/processed/genomicAnalysis/diseaseSim_odds_basedoneQTLGenes.rds')
diag(odds) = NA
p = apply(genemat,2,function(x){
  apply(genemat,2,function(y){
    obs = sum((x == y) & (x != 0))
    p = (sum(x!=0)/ totGeneNum) * (sum(y!=0)/ totGeneNum) 
    bit = binom.test(x = obs, n = totGeneNum, p = p, alternative = 'greater')
    bit$p.value
  })
})
saveRDS(p,'./data/processed/genomicAnalysis/diseaseSim_p_basedoneQTLGenes.rds')
odds[p>0.01] = 1
l2odds = log2(odds)
i = which(rowMeans(l2odds==0, na.rm =T) != 1)
l2odds = l2odds[i,i]
rownames(disannot) = disannot$disease
pheatmap::pheatmap(l2odds, color = colorRampPalette(brewer.pal(8,'Reds'))(20),cellwidth = 5, cellheight = 5, cutree_rows = 4, cutree_cols = 4, annotation_row = disannot[,-3], annotation_col = disannot[,-3], annotation_colors = annotcolors, fontsize = 6, filename = './results/genomicAnalysis/signifeqtlGeneoverlap_heatmap.pdf')

pheatmap::pheatmap(l2odds, color = colorRampPalette(brewer.pal(8,'Reds'))(20),cellwidth = 5, cellheight = 5, cutree_rows = 4, cutree_cols = 4, annotation_row = disannot[,-3], annotation_col = disannot[,-3], annotation_colors = annotcolors, fontsize = 6, filename = './results/genomicAnalysis/signifeqtlGeneoverlap_heatmap.png')

l2odds = log2(odds)
# i = which(rowMeans(l2odds==0, na.rm =T) != 1)
# l2odds = l2odds[i,i]
rownames(disannot) = disannot$disease
pheatmap::pheatmap(l2odds, color = colorRampPalette(brewer.pal(8,'Reds'))(20),cellwidth = 5, cellheight = 5, cutree_rows = 4, cutree_cols = 4, annotation_row = disannot[,-3], annotation_col = disannot[,-3], annotation_colors = annotcolors, fontsize = 6, filename = './results/genomicAnalysis/signifeqtlGeneoverlap_heatmap_all.pdf')

pheatmap::pheatmap(l2odds, color = colorRampPalette(brewer.pal(8,'Reds'))(20),cellwidth = 5, cellheight = 5, cutree_rows = 4, cutree_cols = 4, annotation_row = disannot[,-3], annotation_col = disannot[,-3], annotation_colors = annotcolors, fontsize = 6, filename = './results/genomicAnalysis/signifeqtlGeneoverlap_heatmap_all.png')

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
a=cluster_optimal(mynet,weights = E(mynet)$val)
V(mynet)$clusters = a$membership
net0 = ggnet2(mynet, size = 0, edge.size = 'edgewidth', edge.color = 'gray80')

net1 = net0 +
  geom_point(size = 1.5, color = ageonsetcolors[V(mynet)$ageonset_cl])+
  geom_mark_ellipse(group=V(mynet)$clusters,expand = 0.015)
net2 = net0 +
  geom_point(size = 1.5, color = discatcolors[V(mynet)$discat])+
  geom_mark_ellipse(group=V(mynet)$clusters,expand = 0.015)

netp = ggarrange(net1,net2, labels = c('a','b'))
netp
ggsave('./results/genomicAnalysis/eqtlGENEnet.pdf', netp, units = 'cm', width = 18, height = 8, useDingbats = F)
ggsave('./results/genomicAnalysis/eqtlGENEnet.png', netp, units = 'cm', width = 18, height = 8)
