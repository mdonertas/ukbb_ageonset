source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)

proxygenes = readRDS('./data/processed/genomicAnalysis/signif_proxygenes.rds')
signifSNPs = lapply(proxygenes,function(proxy){
  proxy %>%
    mutate(SNPid = paste(CHR,BP,Ref,Alt,sep='_')) %>%
    filter(!(CHR==mhcchr & BP>= mhcstart & BP<=mhcend)) %>%
    select(SNPid, BETA) %>%
    unique()
})

signifSNPs = reshape2::melt(signifSNPs, id.vars = c('SNPid','BETA')) %>%
  setNames(c('SNPid', 'BETA', 'disID')) %>%
  mutate(disease = disCoding[as.character(disID)]) %>%
  mutate(disCat = disTreecl[disease],
         ageonset = (readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$cluster)[disease]) %>%
  filter(ageonset != 4)

eqtldata = readRDS('./data/processed/genomicAnalysis/snp2gene_eQTL.rds') %>%
  mutate(SNPid = paste(CHR, BP, Ref, Alt, sep = '_')) %>%
  select(-SNP)

eqtl = inner_join(signifSNPs, eqtldata)
rm(eqtldata)

pos = eqtl %>%
  filter(BETA>=0)

neg = eqtl %>%
  filter(BETA<0)

neg = neg %>%
  mutate(ref = Alt,
         alt = Ref,
         eQTL_slope = -eQTL_slope,
         BETA = -BETA) %>%
  mutate(Ref = ref,
         Alt = alt) %>%
  select(-ref,-alt) %>%
  mutate(SNPid = paste(CHR, BP, Ref, Alt, sep ='_')) 

eqtl = rbind(pos,neg) %>%
  mutate(posID = paste(CHR,BP,sep="_"))
rm(pos,neg)

signifSNPinfo = group_by(eqtl, SNPid) %>%
  summarise(numDis = length(unique(disease)),
            numCat = length(unique(disCat)),
            numAC = length(unique(ageonset)),
            aooclusters = paste(sort(unique(ageonset)), collapse = ' & ')) 

eqtl = left_join(eqtl, signifSNPinfo)
rm(signifSNPinfo)

eqtl = eqtl %>%
  group_by(posID) %>%
  summarise(antagonistic = length(unique(SNPid))>1) %>%
  right_join(eqtl)

eqtl = eqtl %>%
  mutate(type = ifelse(numDis == 1, 'Unique', ifelse(numCat > 1, 'Multicategory','Multidisease'))) %>%
  mutate(type = factor(type, levels = c('Unique','Multidisease','Multicategory')))

dissum = group_by(eqtl, disease) %>%
  summarise(numSNP = length(unique(SNPid)))
dissum = setNames(dissum$numSNP,dissum$disease)

tissum = group_by(eqtl, eQTL_tissue) %>%
  summarise(numSNP = length(unique(SNPid)))
tissum = setNames(tissum$numSNP,tissum$eQTL_tissue)

tot = length(unique(eqtl$SNPid))

dis_tissue_mat = group_by(filter(eqtl,type == 'Multicategory'), disease, eQTL_tissue) %>%
  summarise(numSNP = length(unique(SNPid))) %>%
  spread(eQTL_tissue,numSNP,fill = 0) %>%
  as.data.frame()

rownames(dis_tissue_mat) = dis_tissue_mat$disease
dis_tissue_mat$`<NA>` = NULL
dis_tissue_mat$disease = NULL
dis_tissue_mat = as.matrix(dis_tissue_mat)
# dis_tissue_mat = dis_tissue_mat[rowMeans(dis_tissue_mat== 0 )<1,]

x_multicategory = lapply(rownames(dis_tissue_mat),function(disease){
  x = as.data.frame(t(sapply(colnames(dis_tissue_mat),function(tissue){
    a = dis_tissue_mat[disease,tissue]
    b = dissum[disease] - a
    c = tissum[tissue] - a
    d = tot - a - b - c
    fi = fisher.test(matrix(c(a,b,c,d),byrow = T, ncol = 2),alternative = 'greater')
    c(a,b,c,d,fi$est, fi$p.val)
  }))) %>%
    set_names(c('a','b','c','d','odds','p')) 
  x =  x %>% mutate(tissue = rownames(x))
}) %>%
  reshape2::melt(id.vars = c('a','b','c','d','odds','p', 'tissue')) %>%
  rename(disease = L1) %>%
  mutate(disease = rownames(dis_tissue_mat)[disease]) %>%
  mutate(type = 'Multicategory')

dis_tissue_mat = group_by(filter(eqtl,type == 'Multidisease'), disease, eQTL_tissue) %>%
  summarise(numSNP = length(unique(SNPid))) %>%
  spread(eQTL_tissue,numSNP,fill = 0) %>%
  as.data.frame()

rownames(dis_tissue_mat) = dis_tissue_mat$disease
dis_tissue_mat$`<NA>` = NULL
dis_tissue_mat$disease = NULL
dis_tissue_mat = as.matrix(dis_tissue_mat)
# dis_tissue_mat = dis_tissue_mat[rowMeans(dis_tissue_mat== 0 )<1,]

x_multidisease = lapply(rownames(dis_tissue_mat),function(disease){
  x = as.data.frame(t(sapply(colnames(dis_tissue_mat),function(tissue){
    a = dis_tissue_mat[disease,tissue]
    b = dissum[disease] - a
    c = tissum[tissue] - a
    d = tot - a - b - c
    fi = fisher.test(matrix(c(a,b,c,d),byrow = T, ncol = 2),alternative = 'greater')
    c(a,b,c,d,fi$est, fi$p.val)
  }))) %>%
    set_names(c('a','b','c','d','odds','p')) 
  x =  x %>% mutate(tissue = rownames(x))
}) %>%
  reshape2::melt(id.vars = c('a','b','c','d','odds','p', 'tissue')) %>%
  rename(disease = L1) %>%
  mutate(disease = rownames(dis_tissue_mat)[disease]) %>%
  mutate(type = 'Multidisease')

dis_tissue_mat = group_by(filter(eqtl,type == 'Unique'), disease, eQTL_tissue) %>%
  summarise(numSNP = length(unique(SNPid))) %>%
  spread(eQTL_tissue,numSNP,fill = 0) %>%
  as.data.frame()

rownames(dis_tissue_mat) = dis_tissue_mat$disease
dis_tissue_mat$`<NA>` = NULL
dis_tissue_mat$disease = NULL
dis_tissue_mat = as.matrix(dis_tissue_mat)
# dis_tissue_mat = dis_tissue_mat[rowMeans(dis_tissue_mat== 0 )<1,]

x_unique = lapply(rownames(dis_tissue_mat),function(disease){
  x = as.data.frame(t(sapply(colnames(dis_tissue_mat),function(tissue){
    a = dis_tissue_mat[disease,tissue]
    b = dissum[disease] - a
    c = tissum[tissue] - a
    d = tot - a - b - c
    fi = fisher.test(matrix(c(a,b,c,d),byrow = T, ncol = 2),alternative = 'greater')
    c(a,b,c,d,fi$est, fi$p.val)
  }))) %>%
    set_names(c('a','b','c','d','odds','p')) 
  x =  x %>% mutate(tissue = rownames(x))
}) %>%
  reshape2::melt(id.vars = c('a','b','c','d','odds','p', 'tissue')) %>%
  rename(disease = L1) %>%
  mutate(disease = rownames(dis_tissue_mat)[disease]) %>%
  mutate(type = 'Unique')

x = rbind(x_multicategory, x_multidisease, x_unique)


x = x %>%
  # filter(a>10) %>%
  mutate(log2odds = log2(odds),
         padj = p.adjust(p, method = 'fdr')) 

p = x %>%
  select(tissue,disease,type,log2odds) %>%
  unique() %>%
  spread(type, log2odds) %>%
  mutate(xx = Multicategory - Unique) %>%
  ggplot(aes(x = tissue, y = xx)) +
  stat_summary(fun.min = function(x){quantile(x,0.25)},
               fun.max = function(x){quantile(x,0.75)},
               fun = function(x)median(x)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'darkred') +
  # stat_summary() +
  # geom_boxplot() +
  coord_flip() +
  theme_bw() +
  ylab('Odds Ratio Difference\n(Multicategory - Private SNPs)') +
  xlab(NULL)

ggsave('./results/eQTL/multicat_uniq_oddsdiff.pdf', p, units = 'cm',width = 16, height = 16, useDingbats = F)
ggsave('./results/eQTL/multicat_uniq_oddsdiff.png', p, units = 'cm',width = 16, height = 16)

dis_tissue_odds = x %>%
  # filter(a>0 & b>0) %>%
  filter(type == 'Unique') %>%
  filter(padj < 0.1) %>%
  select(disease,log2odds,tissue) %>%
  spread(tissue,log2odds,fill = 0)

rownames(dis_tissue_odds) = dis_tissue_odds$disease
dis_tissue_odds$disease = NULL
dis_tissue_odds = as.matrix(dis_tissue_odds)
dis_tissue_odds = dis_tissue_odds[rowMeans(dis_tissue_odds==0)<1,]
# dis_tissue_odds[(dis_tissue_odds==Inf)] = 7
pheatmap::pheatmap(dis_tissue_odds, color = c('gray',colorRampPalette(c('white','#B2182B'))(24)),
                   breaks = c(0,seq(0.01,7,length.out = 25)), main = 'eQTL Enrichment for Private SNPs',
                   filename = './results/eQTL/odds_uniq.pdf', cellwidth = 10, cellheight = 10)
pheatmap::pheatmap(dis_tissue_odds, color = c('gray',colorRampPalette(c('white','#B2182B'))(24)),
                   breaks = c(0,seq(0.01,7,length.out = 25)), main = 'eQTL Enrichment for Private SNPs',
                   filename = './results/eQTL/odds_uniq.png', cellwidth = 10, cellheight = 10)


dis_tissue_odds = x %>%
  # filter(a>0 & b>0) %>%
  filter(type == 'Multicategory') %>%
  filter(padj < 0.1) %>%
  select(disease,log2odds,tissue) %>%
  spread(tissue,log2odds,fill = 0)

rownames(dis_tissue_odds) = dis_tissue_odds$disease
dis_tissue_odds$disease = NULL
dis_tissue_odds = as.matrix(dis_tissue_odds)
dis_tissue_odds = dis_tissue_odds[rowMeans(dis_tissue_odds==0)<1,]
dis_tissue_odds[(dis_tissue_odds==Inf)] = 7
pheatmap::pheatmap(dis_tissue_odds, color = c('gray',colorRampPalette(c('white','#B2182B'))(24)),
                   breaks = c(0,seq(0.01,7,length.out = 25)), main = 'eQTL Enrichment for Multicategory SNPs',
                   filename = './results/eQTL/odds_multicat.pdf', cellwidth = 10, cellheight = 10)
pheatmap::pheatmap(dis_tissue_odds, color = c('gray',colorRampPalette(c('white','#B2182B'))(24)),
                   breaks = c(0,seq(0.01,7,length.out = 25)), main = 'eQTL Enrichment for Multicategory SNPs',
                   filename = './results/eQTL/odds_multicat.png', cellwidth = 10, cellheight = 10)

n_fun <- function(x){
  return(data.frame(y = -3, label = paste("n = ",scales::comma(length(x)))))
}

p2 = group_by(eqtl, SNPid, type, aooclusters) %>%
  summarise(numTissue = length(unique(eQTL_tissue))) %>%
  ggplot(aes(x = type, y = numTissue)) +
  geom_sina(size =0.1, alpha = 0.1, color = 'gray') +
  stat_summary(fun.min = function(x){quantile(x,0.25)},
               fun.max = function(x){quantile(x,0.75)},
               fun = function(x)median(x),
               color = 'darkred') +
  stat_summary(fun.data = n_fun, geom = "text", size = 8/pntnorm) +
  # geom_boxplot(width = 0.05, color = 'darkred', outlier.shape = NA) +
  stat_compare_means(comparisons = list(c('Unique','Multidisease'),c('Unique','Multicategory'))) +
  xlab(NULL) +
  ylab('Number of significant tissues')+ 
  ylim(-5,61)

ggsave('./results/eQTL/numtissue.pdf', p2, units = 'cm',width = 10, height = 8, useDingbats = F)
ggsave('./results/eQTL/numtissue.png', p2, units = 'cm',width = 10, height = 8)

xx = group_by(eqtl, SNPid, type, aooclusters) %>%
  summarise(numTissue = length(unique(eQTL_tissue))) %>%
  ungroup()

bootstrapxx = sapply(as.character(unique(xx$type)),function(tip){
  sapply(1:1000,function(i){
    mean(sample(filter(xx,type == tip)$numTissue, 3000))
  })
})

bootstrapxx = reshape2::melt(bootstrapxx) 

p3 = ggplot(bootstrapxx, aes(x = Var2, y = value)) +
  geom_sina(size =0.1, alpha = 0.1, color = 'gray') +
  stat_summary(fun.min = function(x){quantile(x,0.25)},
               fun.max = function(x){quantile(x,0.75)},
               fun = function(x)median(x),
               color = 'darkred') +
  stat_compare_means(comparisons = list(c('Unique','Multidisease'),c('Unique','Multicategory'))) +
  xlab(NULL) +
  ylab('Number of significant tissues')+ 
  ylim(12,23) +
  ggtitle('1,000 samples of 3,000 SNPs from each category')

ggsave('./results/eQTL/numtissue_bootstrap.pdf', p3, units = 'cm',width = 12, height = 8, useDingbats = F)
ggsave('./results/eQTL/numtissue_bootstrap.png', p3, units = 'cm',width = 12, height = 8)


##############################################
##direction of change in expression with age##
##############################################

exprtissues = list.files('../aging_in_GTEx_v8/data/processed/expression/change_w_age/', full.names = T)
exprtissuename = gsub('.rds','',sapply(strsplit(exprtissues,'/'),function(x)x[8]))
expr = lapply(exprtissues,function(x){
  y = as.data.frame(readRDS(x))
  y$eQTL_ensembl = rownames(y)
  y
  })
names(expr) = gsub('-','_',exprtissuename)
expr = reshape2::melt(expr, id.vars = c('eQTL_ensembl','rho','pval','padj')) %>%
  rename(eQTL_tissue = L1,
         aging_rho = rho,
         aging_p = pval,
         aging_padj = padj) 

# remove antagonistic SNPs
eqtl = eqtl %>%
  filter(!antagonistic)

eqtl = left_join(eqtl,expr)

eqtl = eqtl %>%
  mutate(eQTLage = c('opposite','same')[((eQTL_slope * aging_rho)>0)+1]) 

gene_info = eqtl %>%
  select(disease,disCat,ageonset,eQTL_ensembl) %>%
  unique() %>%
  group_by(eQTL_ensembl) %>%
  summarise(geneDis = length(unique(disease)),
            geneCat = length(unique(disCat)),
            gene_aoocluster = paste(sort(unique(ageonset)), collapse = ' & '))

eqtl = left_join(eqtl, gene_info)

eqtl = eqtl %>%
  mutate(genetype = ifelse(geneDis == 1, 'Unique', ifelse(geneCat > 1, 'Multicategory','Multidisease'))) %>%
  mutate(genetype = factor(genetype, levels = c('Unique','Multidisease','Multicategory')))

eqtl %>%
  # filter(!is.na(eQTLage)) %>%
  filter(abs(aging_rho) >= 0.1) %>%
  select(type,aooclusters,SNPid,eQTLage,disease,eQTL_tissue) %>%
  group_by(type,aooclusters,SNPid,disease,eQTL_tissue) %>%
  summarise(consistent = length(unique(eQTLage))==1,
            eQTLage = paste(eQTLage, collapse=', ')) %>%
  filter(consistent) %>%
  ungroup() %>%
  group_by(type,aooclusters,SNPid, disease) %>%
  summarise( perc = mean(eQTLage == 'same')) %>%
  ungroup() %>%
  group_by(SNPid, type, aooclusters) %>%
  summarise(perc = mean(perc)) %>%
  filter(aooclusters %in% c('1','2','3')) %>%
  ggplot(aes( x = type , y = perc)) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'darkred') +
  geom_violin(fill = 'gray70') +
  geom_boxplot(width = 0.1) +
  facet_wrap(~aooclusters) +
  xlab(NULL) + ylab('% change in the same direction as ageing')


p3 = eqtl %>% filter(abs(aging_rho)>=0.1) %>%
  group_by(SNPid, type, aooclusters, disease, disCat, ageonset, eQTL_tissue) %>%
  summarise(dir = mean(eQTLage == 'same', na.rm = T)) %>%
  na.omit() %>%
  ungroup() %>%
  group_by(type, aooclusters, disease, disCat, ageonset, eQTL_tissue) %>%
  summarise(dir = mean(dir)) %>% # hastalik basina dokuda tipine gore SNP ortalamasi
  ungroup(eQTL_tissue) %>%
  summarise(dir = median(dir)) %>% #per disease
  filter(aooclusters %in% c('1','2','3')) %>%
  ggplot(aes(x = type , y= dir)) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color  ='darkred') +
  facet_wrap(~aooclusters) +
  geom_boxplot(outlier.shape = NA, aes(fill = aooclusters), size = 0.25) +
  scale_fill_manual(values = ageonsetcolors) +
  geom_jitter(width = 0.1, size = 0.5, color = 'gray25') +
  xlab(NULL) + ylab('% change in the same direction as ageing') +
  guides(fill = F) +
  theme_pubr(base_size = 8) +
  coord_flip() +
  theme(panel.border = element_rect(color = 'black', size = 0.5, fill = NA)) 

ggsave('./results/eQTL/samedir_age.pdf', p3, units = 'cm',width = 16.8, height = 5, useDingbats = F)
ggsave('./results/eQTL/samedir_age.png', p3, units = 'cm',width = 16.8, height = 5)

