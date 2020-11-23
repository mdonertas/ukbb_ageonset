source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)

proxygenes = readRDS('./data/processed/genomicAnalysis/signif_proxygenes.rds')
signifSNPs = lapply(proxygenes,function(proxy){
  proxy %>%
    mutate(SNPid = paste(CHR,BP,Ref,Alt,sep='_'),
           posID = paste(CHR,BP,sep="_")) %>%
    filter(!(CHR==mhcchr & BP>= mhcstart & BP<=mhcend)) %>%
    select(posID) %>%
    unique()
})
signifSNPs = reshape2::melt(signifSNPs) %>%
  setNames(c('posID','disID')) %>%
  mutate(disease = disCoding[as.character(disID)]) %>%
  mutate(disCat = disTreecl[disease],
         ageonset = (readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$cluster)[disease]) %>%
  filter(ageonset != 4)
signifSNPinfo = group_by(signifSNPs, posID) %>%
  summarise(numDis = length(unique(disease)),
            numCat = length(unique(disCat)),
            numAC = length(unique(ageonset)),
            aooclusters = paste(sort(unique(ageonset)), collapse = ' & ')) 

signifSNPs = left_join(signifSNPs, signifSNPinfo)

signifSNPs = signifSNPs %>%
  mutate(type = ifelse(numDis == 1, 'Unique', ifelse(numCat > 1, 'Multicategory','Multidisease'))) %>%
  mutate(type = factor(type, levels = c('Unique','Multidisease','Multicategory')))

edudietfn = grep('bgz',list.files('./data/raw/edu_diet_gwas/',full.names = T),invert=T,v=T)[-6]
edudietn = sapply(strsplit(edudietfn,'/'),function(x)sapply(strsplit(x[length(x)],'[.]'),function(x)x[1]))
edudiet = lapply(edudietfn,function(x){
  read_tsv(x) %>%
    select(variant,pval) %>%
    separate(variant, into = c('CHR','BP','Ref','Alt'), remove = F) %>%
    mutate(posID = paste(CHR,BP,sep='_')) %>%
    filter(posID %in% signifSNPs$posID) %>%
    # filter(pval<=0.01) %>%
    unique()
})
names(edudiet) = edudietn

snpgrs = group_by(signifSNPs,aooclusters,type) %>%
  summarise(snplist = list(sort(unique(posID)))) %>%
  mutate(snp_type = paste('cluster',aooclusters,'-',type))

snpgrs = setNames(snpgrs$snplist,snpgrs$snp_type)
ukbbsnps = unique((read_tsv('./data/processed/ukbb/gwas/bolt/a1071.imp.stats') %>%
                     mutate(posID = paste(CHR,BP,sep="_")))$posID)
nealesnps = unique((read_tsv(edudietfn[1]) %>%
                      separate(variant, into = c('CHR','BP','Ref','Alt'), remove = F) %>%
                      mutate(posID = paste(CHR,BP,sep='_')))$posID)
testsnps = intersect(ukbbsnps,nealesnps)
rm(ukbbsnps, nealesnps)
x = lapply(c(0.05,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,5e-8),function(pcut){
  edudietsnps = sapply(edudiet,function(x)unique(filter(x,CHR%in%1:22 & pval<=pcut)$posID))
  print(pcut)
  x = lapply(snpgrs[c(2,3,4,5,6,7,9,10,11,12,14,15)],function(x){
    y = t(sapply(edudietsnps,function(y){
      a = length(unique(intersect(x,y)))
      b = length(unique(setdiff(x,y)))
      c = length(unique(setdiff(y,x)))
      d = length(unique(setdiff(testsnps,union(x,y))))
      mat = matrix(c(a,b,c,d),ncol=2,byrow=T)
      fi = fisher.test(mat)
      c(a,b,c,d,fi$est,fi$p.val)
    }))
    colnames(y) = c('a','b','c','d','odds','pval')
    print(y)
    y
  })
  x
})
names(x) = c(0.05,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,5e-8)

names(ageonsetcolors) = paste('cluster ',gsub('-',' & ',names(ageonsetcolors)),sep='')
p = reshape2::melt(x) %>%
  spread(Var2,value) %>%
  rename(edudiet = Var1,
         nm = L2,
         pcut = L1) %>%
  separate(nm, into = c('aoo','type'), sep = ' - ',remove = F) %>%
  mutate(log2odds = log2(odds),
         p = -log10(pval)) %>%
  mutate(edudiet = factor(setNames(c('Gluten free','Lactose free','Low Calorie','Other Diet','Townsend Deprivation Ind','Vegan','Vegetarian'),names(edudiet))[edudiet],levels = c('Townsend Deprivation Ind','Gluten free','Lactose free','Low Calorie','Vegan','Vegetarian','Other Diet'))) %>%
  # filter(type == 'Multicategory') %>% 
  ggplot(aes(x = pcut, y = log2odds, color = aoo, group = aoo)) +
  facet_grid(type~edudiet) +
  geom_point(aes(size = a)) +
  geom_line() +
  scale_color_manual(values = ageonsetcolors) +
  scale_size_continuous(range =c(0.1,3)) +
  theme_pubr(base_size = 6, legend = 'bottom') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  guides(size = guide_legend('Number of overlapping SNPs'),
        color = guide_legend(title = NULL)) +
  ylab('log2 Odds Ratio') +
  xlab('p cutoff') 

ggsave('./results/edudiet.pdf', units = 'cm', width = 16.8, height = 8, useDingbats = F)
ggsave('./results/edudiet.png', units = 'cm', width = 16.8, height = 8)
