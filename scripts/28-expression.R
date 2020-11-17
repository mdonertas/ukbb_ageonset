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
  mutate(type = ifelse(numDis == 1, 'Unique', ifelse(numCat > 1, 'Multicategory','Multidisease'))) %>%
  mutate(type = factor(type, levels = c('Unique','Multidisease','Multicategory')))

eqtl = eqtl %>%
  group_by(posID) %>%
  summarise(antagonistic = length(unique(SNPid))>1) %>%
  right_join(eqtl)

eqtl = eqtl %>%
  filter(aooclusters %in% c('1','2','3')) %>%
  filter(!antagonistic)

eqtl %>%
  group_by(aooclusters) %>%
  summarise(n = length(unique(SNPid)))
# # A tibble: 3 x 2
# aooclusters     n
# <chr>       <int>
# 1 1           34675
# 2 2           13624
# 3 3           10162


eqtl_summary = group_by(eqtl, type, aooclusters, eQTL_tissue, eQTL_ensembl) %>%
  summarise(snp_direction = mean(eQTL_slope >= 0))

exprtissues = list.files('../aging_in_GTEx_v8/data/processed/expression/tpm_lm_l2_qn_sc/', full.names = T)
exprtissuename = gsub('.rds','',sapply(strsplit(exprtissues,'/'),function(x)x[8]))
expr = lapply(exprtissues,function(x){
  y = as.data.frame(readRDS(x))
  y$eQTL_ensembl = rownames(y)
  y
})
names(expr) = gsub('-','_',exprtissuename)
expr = reshape2::melt(expr, id.vars = c('eQTL_ensembl')) %>%
  rename(sample_id = variable,
         expression = value,
         eQTL_tissue = L1)

attr = readRDS('../aging_in_GTEx_v8/data/processed/attr.rds')
expr = left_join(expr, attr)
eqtl_summary = left_join(eqtl_summary, expr) 
head(eqtl_summary)

dircut = c(1,0.9,0.5)

for(dircut in c(1,0.9,0.75,0.6)){
  xx = eqtl_summary %>%
    ungroup() %>%
    filter(!is.na(expression)) %>%
    filter(snp_direction >= dircut) %>%
    select(type,aooclusters, eQTL_tissue,eQTL_ensembl, snp_direction, sample_id, expression, age) %>%
    group_by(type,aooclusters, eQTL_tissue,eQTL_ensembl,age) %>%
    summarise(exp = mean(expression)) %>%
    spread(age,exp) %>% 
    mutate(`30` = `30-39`-`20-29`,
           `40` = `40-49`-`30-39`,
           `50` = `50-59`-`40-49`,
           `60` = `60-69`-`50-59`,
           `70` = `70-79`-`60-69`) %>%
    select(type,aooclusters,eQTL_tissue,eQTL_ensembl, `30`,`40`,`50`,`60`,`70`) %>%
    gather('age','diff',-type,-aooclusters,-eQTL_tissue,-eQTL_ensembl) %>%
    ungroup() %>% 
    group_by(type, aooclusters, eQTL_tissue, age) %>%
    summarise(meandiff = mean(diff),
              mediandiff = median(diff),
              num_gene = length(unique(eQTL_ensembl))) 
  xa = lapply(unique(xx$type),function(tip){
    xa = sapply(unique(xx$aooclusters),function(agex){
      xa = filter(xx, type == tip & aooclusters == agex)
      p30_40 = wilcox.test(filter(xa,age=='30')$mediandiff,filter(xa,age=='40')$mediandiff, alternative = 'less')$p.val
      p40_50 = wilcox.test(filter(xa,age=='40')$mediandiff,filter(xa,age=='50')$mediandiff, alternative = 'greater')$p.val
      # p50_60 = wilcox.test(filter(xa,age=='50')$mediandiff,filter(xa,age=='60')$mediandiff)$p.val
      # p60_70 = wilcox.test(filter(xa,age=='60')$mediandiff,filter(xa,age=='70')$mediandiff)$p.val
      c(p30_40,p40_50)
    })
    rownames(xa) = c('p30_40','p40_50')
    xa
  })
  names(xa) = unique(xx$type)
  xa = reshape2::melt(xa) %>%
    set_names(c('agegr','aoocluster','p','type'))
  xa = filter(xa,p<=0.05)
  print(list(dircut, xa))
  p = xx %>%
    ggplot(aes(x = age, y = mediandiff, fill = aooclusters)) +
    facet_grid(aooclusters~type) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, aes(size = num_gene), alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'darkred') +
    xlab(NULL) + ylab('Median difference betweeen the mean expression levels between ages per tissue') +
    theme_pubr(base_size = 10, legend = 'right') +
    scale_fill_manual(values = ageonsetcolors) +
    guides(fill = F, size= guide_legend('# Genes')) +
    scale_size_continuous(range = c(0.5,2)) +
    ggtitle(paste('eQTL direction Up (eqtldir >= ',dircut,')',sep=''))
  ggsave(paste('./results/eQTL/diff_',gsub('[.]','',dircut),'_median.pdf', sep = ''), p, units = 'cm',width = 16, height =18, useDingbats = F)
  ggsave(paste('./results/eQTL/diff_',gsub('[.]','',dircut),'_median.png', sep = ''), p, units = 'cm',width = 16, height =18)
  p = xx %>%
    ggplot(aes(x = age, y = meandiff, fill = aooclusters)) +
    facet_grid(aooclusters~type) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, aes(size = num_gene), alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'darkred') +
    xlab(NULL) + ylab('Mean difference betweeen the mean expression levels between ages per tissue') +
    theme_pubr(base_size = 10, legend = 'right') +
    scale_fill_manual(values = ageonsetcolors) +
    guides(fill = F, size= guide_legend('# Genes')) +
    scale_size_continuous(range = c(0.5,2)) +
    ggtitle(paste('eQTL direction Up (eqtldir >= ',dircut,')',sep=''))
  ggsave(paste('./results/eQTL/diff_',gsub('[.]','',dircut),'_mean.pdf', sep = ''), p, units = 'cm',width = 16, height =18, useDingbats = F)
  ggsave(paste('./results/eQTL/diff_',gsub('[.]','',dircut),'_mean.png', sep = ''), p, units = 'cm',width = 16, height =18)
}

for(dircut in c(0,0.1,0.25,0.4)){
  xx = eqtl_summary %>%
    ungroup() %>%
    filter(!is.na(expression)) %>%
    filter(snp_direction <= dircut) %>%
    select(type,aooclusters, eQTL_tissue,eQTL_ensembl, snp_direction, sample_id, expression, age) %>%
    group_by(type,aooclusters, eQTL_tissue,eQTL_ensembl,age) %>%
    summarise(exp = mean(expression)) %>%
    spread(age,exp) %>% 
    mutate(`30` = `30-39`-`20-29`,
           `40` = `40-49`-`30-39`,
           `50` = `50-59`-`40-49`,
           `60` = `60-69`-`50-59`,
           `70` = `70-79`-`60-69`) %>%
    select(type,aooclusters,eQTL_tissue,eQTL_ensembl, `30`,`40`,`50`,`60`,`70`) %>%
    gather('age','diff',-type,-aooclusters,-eQTL_tissue,-eQTL_ensembl) %>%
    ungroup() %>% 
    group_by(type, aooclusters, eQTL_tissue, age) %>%
    summarise(meandiff = mean(diff),
              mediandiff = median(diff),
              num_gene = length(unique(eQTL_ensembl))) 
  xa = lapply(unique(xx$type),function(tip){
    xa = sapply(unique(xx$aooclusters),function(agex){
      xa = filter(xx, type == tip & aooclusters == agex)
      p30_40 = wilcox.test(filter(xa,age=='30')$mediandiff,filter(xa,age=='40')$mediandiff, alternative = 'greater')$p.val
      p40_50 = wilcox.test(filter(xa,age=='40')$mediandiff,filter(xa,age=='50')$mediandiff, alternative = 'less')$p.val
      # p50_60 = wilcox.test(filter(xa,age=='50')$mediandiff,filter(xa,age=='60')$mediandiff)$p.val
      # p60_70 = wilcox.test(filter(xa,age=='60')$mediandiff,filter(xa,age=='70')$mediandiff)$p.val
      c(p30_40,p40_50)
    })
    rownames(xa) = c('p30_40','p40_50')
    xa
  })
  names(xa) = unique(xx$type)
  xa = reshape2::melt(xa) %>%
    set_names(c('agegr','aoocluster','p','type'))
  xa = filter(xa,p<=0.05)
  print(list(dircut, xa))
  p = xx %>%
    ggplot(aes(x = age, y = mediandiff, fill = aooclusters)) +
    facet_grid(aooclusters~type) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, aes(size = num_gene), alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'darkred') +
    xlab(NULL) + ylab('Median difference betweeen the mean expression levels between ages per tissue') +
    theme_pubr(base_size = 10, legend = 'right') +
    scale_fill_manual(values = ageonsetcolors) +
    guides(fill = F, size= guide_legend('# Genes')) +
    scale_size_continuous(range = c(0.5,2)) +
    ggtitle(paste('eQTL direction Down (eqtldir >= ',dircut,')',sep=''))
  ggsave(paste('./results/eQTL/diff_',gsub('[.]','',dircut),'_median.pdf', sep = ''), p, units = 'cm',width = 16, height =18, useDingbats = F)
  ggsave(paste('./results/eQTL/diff_',gsub('[.]','',dircut),'_median.png', sep = ''), p, units = 'cm',width = 16, height =18)
  p = xx %>%
    ggplot(aes(x = age, y = meandiff, fill = aooclusters)) +
    facet_grid(aooclusters~type) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, aes(size = num_gene), alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'darkred') +
    xlab(NULL) + ylab('Mean difference betweeen the mean expression levels between ages per tissue') +
    theme_pubr(base_size = 10, legend = 'right') +
    scale_fill_manual(values = ageonsetcolors) +
    guides(fill = F, size= guide_legend('# Genes')) +
    scale_size_continuous(range = c(0.5,2)) +
    ggtitle(paste('eQTL direction Up (eqtldir >= ',dircut,')',sep=''))
  ggsave(paste('./results/eQTL/diff_',gsub('[.]','',dircut),'_mean.pdf', sep = ''), p, units = 'cm',width = 16, height =18, useDingbats = F)
  ggsave(paste('./results/eQTL/diff_',gsub('[.]','',dircut),'_mean.png', sep = ''), p, units = 'cm',width = 16, height =18)
}

# - dircut 1
# - agegr aoocluster p type
# 1 p30_40 2 0.001042861 Unique
# 2 p40_50 2 0.005182969 Unique
# 3 p30_40 2 0.043375945 Multidisease
# 4 p30_40 1 0.012140784 Multicategory
# 5 p40_50 1 0.032992904 Multicategory
# 6 p40_50 2 0.043054248 Multicategory
# - dircut 0.9
# - agegr aoocluster p type
# 1 p30_40 2 0.003603334 Unique
# 2 p40_50 2 0.002573135 Unique
# 3 p30_40 2 0.049699335 Multidisease
# 4 p30_40 1 0.012140784 Multicategory
# 5 p40_50 1 0.032992904 Multicategory
# 6 p40_50 2 0.043054248 Multicategory
# - dircut 0.75
# - agegr aoocluster p type
# 1 p30_40 2 0.006915503 Unique
# 2 p40_50 2 0.002407209 Unique
# 3 p30_40 2 0.045408846 Multidisease
# 4 p30_40 1 0.012140784 Multicategory
# 5 p40_50 1 0.032992904 Multicategory
# 6 p40_50 2 0.043054248 Multicategory
# - dircut 0.6
# - agegr aoocluster p type
# 1 p30_40 2 0.006494900 Unique
# 2 p40_50 2 0.001711611 Unique
# 3 p30_40 1 0.012140784 Multicategory
# 4 p40_50 1 0.032992904 Multicategory
# 5 p40_50 2 0.043054248 Multicategory
# - dircut 0.4
# - agegr aoocluster p type
# 1 p30_40 1 0.0089036559 Unique
# 2 p40_50 3 0.0038376263 Unique
# 3 p30_40 2 0.0027332432 Multidisease
# 4 p30_40 2 0.0004622526 Multicategory
# 5 p40_50 2 0.0031331887 Multicategory
# 6 p30_40 3 0.0075343934 Multicategory
# 7 p40_50 3 0.0007439967 Multicategory
# - dircut 0.25
# - agegr aoocluster p type
# 1 p30_40 1 0.0066716514 Unique
# 2 p40_50 3 0.0057595867 Unique
# 3 p30_40 2 0.0016477228 Multidisease
# 4 p30_40 2 0.0004622526 Multicategory
# 5 p40_50 2 0.0031331887 Multicategory
# 6 p30_40 3 0.0075343934 Multicategory
# 7 p40_50 3 0.0007439967 Multicategory
# - dircut 0.1
# - agegr aoocluster p type
# 1 p30_40 1 0.0062887512 Unique
# 2 p40_50 3 0.0019423069 Unique
# 3 p30_40 2 0.0023722626 Multidisease
# 4 p30_40 2 0.0002295567 Multicategory
# 5 p40_50 2 0.0017116107 Multicategory
# 6 p30_40 3 0.0075343934 Multicategory
# 7 p40_50 3 0.0007439967 Multicategory
# - dircut 0
# - agegr aoocluster p type
# 1 p30_40 1 0.0188242780 Unique
# 2 p40_50 3 0.0071970999 Unique
# 3 p30_40 2 0.0014182072 Multidisease
# 4 p30_40 2 0.0002295567 Multicategory
# 5 p40_50 2 0.0017116107 Multicategory
# 6 p30_40 3 0.0075343934 Multicategory
# 7 p40_50 3 0.0007439967 Multicategory