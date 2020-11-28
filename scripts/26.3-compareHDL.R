source('./scripts/00-setup.R')
library(ggthemes)
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
mydat = alldat %>%
  mutate(snp_odds = log2(snp_odds),
         RR = log2(RR),
         gene_odds = log2(gene_odds),
         proxy_gene_odds = log2(proxy_gene_odds),
         eqtl_gene_odds = log2(eqtl_gene_odds)) %>%
  # mutate(snpodds = ifelse(snp_p<0.01,snp_odds,1)) %>%
  # filter(snp_p<0.01) %>%
  select(disA,disB,snp_odds,sameCat,cooccur,sameAge,snp_p) %>%
  mutate(sameCat = factor(c('different category','same category')[1+sameCat]),
         sameAge = factor(c('different age clusters','same age cluster')[1+sameAge])) 

# a=ifelse(mydat$disA>mydat$disB,mydat$disA,mydat$disB)
# b=ifelse(mydat$disA>mydat$disB,mydat$disB,mydat$disA)
# mydat$disA = a
# mydat$disB = b
# nrow(mydat)
# nrow(unique(mydat))
# mydat = unique(mydat)

xx = readRDS('./data/processed/HDL_summary/HDLsummaryTable.rds') %>%
  rename(disA = disease1,
         disB = disease2)

p1=left_join(mydat, xx) %>%
  select(disA,disB,snp_odds,sameCat,cooccur,sameAge,rg,rg.se,p,cor,cluster1, cluster2,cat1,cat2,samecat,samecl,uplevel,downlevel) %>%
  na.omit() %>%
  filter(snp_odds!=(-Inf)) %>%
  ggplot(aes(x = snp_odds, y= rg)) +
  geom_point(aes( color = sameAge),size = 0.5) +
  geom_smooth(method = 'lm') +
  stat_cor(method='spearman', cor.coef.name = 'rho') +
  scale_color_gdocs()+
  guides(color = guide_legend(title=NULL)) +
  xlab('SNP Overlap-based Odds Ratio') +
  ylab('HDL rg') +
  theme(legend.direction = 'vertical',
        legend.position = c(0.3,0.77),
        legend.background = element_rect(fill='gray80',color = 'gray35'),
        legend.spacing.y = ggplot2::unit(0.01,'mm'),
        legend.margin = margin(r=0.5,0,0,0,unit = 'mm'))

ggsave('./results/HDL_overlap_compare.pdf',p1,units = 'cm', width = 8, height = 8, useDingbats =F)
ggsave('./results/HDL_overlap_compare.png',p1,units = 'cm', width = 8, height = 8)

p2=left_join(mydat, xx) %>%
  select(disA,disB,snp_odds,sameCat,cooccur,sameAge,rg,rg.se,p,cor,cluster1, cluster2,cat1,cat2,samecat,samecl,uplevel,downlevel) %>%
  na.omit() %>%
  filter(!uplevel  & !downlevel) %>%
  mutate(padj = p.adjust(p,method = 'fdr')) %>%
  filter(padj<0.05) %>%
  ggplot(aes(x = sameAge, y= rg, fill = sameAge)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.5,size = 0.3, color = 'gray35',width = 0.20) +
  stat_compare_means() +
  scale_fill_gdocs() +
  guides(fill = F,color =F) +
  xlab(NULL) + ylab('HDL rg')

ggsave('./results/HDLres.pdf',p2,units = 'cm', width = 8, height = 8, useDingbats =F)
ggsave('./results/HDLres.png',p2,units = 'cm', width = 8, height = 8)

p = ggarrange(p1,p2,align = 'hv', legend=F, labels = 'auto')
p

ggsave('./results/HDLfig.pdf',p,units = 'cm', width = 16.7, height = 8, useDingbats =F)
ggsave('./results/HDLfig.png',p,units = 'cm', width = 16.7, height = 8)

