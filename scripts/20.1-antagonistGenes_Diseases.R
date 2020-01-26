source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('../ukbb_ageonset/results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)
disCoding = disCoding[as.character(disIDs)]
proxyGenes <- lapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_proxyGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  x %>% 
    select(CHR,BP,Ref,Alt,proxy_hgnc) %>%
    na.omit() %>%
    unique() %>%
    mutate(snpid = paste(CHR,BP,Alt,Ref,sep='_'))
})
names(proxyGenes)=disIDs
proxyGenes = reshape2::melt(proxyGenes,id.vars = colnames(proxyGenes[[1]])) %>% 
  rename(disID = L1, geneid = proxy_hgnc) %>%
  mutate(proxy = TRUE)
eqtlGenes <- lapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_eQTLGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  x %>% 
    select(CHR,BP,Ref,Alt,eQTL_hgnc) %>%
    na.omit() %>%
    unique() %>%
    mutate(snpid = paste(CHR,BP,Alt,Ref,sep='_'))
})
names(eqtlGenes)=disIDs
eqtlGenes = reshape2::melt(eqtlGenes,id.vars = colnames(eqtlGenes[[1]])) %>% 
  rename(disID = L1, geneid = eQTL_hgnc) %>%
  mutate(eQTL = TRUE)
signifGenes <- full_join(proxyGenes, eqtlGenes) 
signifGenes <- signifGenes %>%
  mutate(disease = disCoding[as.character(disID)]) %>%
  mutate(disCat = disTreecl[disease],
         ageonset = (readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$cluster)[disease]) 
signifSNPs=readRDS('./data/processed/evoAnalysis/UKBBRAF_Pleiotropy.rds')
signifSNPs = signifSNPs %>%
  mutate(gr_agonist = grepl('agonist',cl1_cl1) |grepl('agonist',cl2_cl2)|grepl('agonist',cl3_cl3),
         gr_antagonist = grepl('antagonist',cl1_cl1) |grepl('antagonist',cl2_cl2)|grepl('antagonist',cl3_cl3),
         ac_agonist = grepl('agonist',cl1_cl2) | grepl('agonist',cl1_cl3) | grepl('agonist',cl2_cl3),
         ac_antagonist = grepl('antagonist',cl1_cl2) | grepl('antagonist',cl1_cl3) | grepl('antagonist',cl2_cl3))
signifSNPs = signifSNPs %>%
  mutate(snpid = paste(CHR,BP,Alt,Ref,sep='_'),
         snppos = paste(CHR,BP,Alt,Ref,sep='_'))
mafxx = signifSNPs %>%
  mutate(betagr = c(1:20/20)[cut(abs(signifSNPs$BETA),breaks = quantile(abs(signifSNPs$BETA),probs = seq(0,1,length.out = 21)),include.lowest = T)])
mafxxx = mafxx %>% select(snpid,betagr) %>% group_by(snpid) %>% summarise(betagr = min(betagr))
# snpexc = unique(filter(mafxx,betagr<0)$snpid)

signifSNPs = signifSNPs %>%
  filter(snppos %in% unique((signifSNPs %>% group_by(snppos) %>% summarise(num = length(unique(snpid))) %>% filter(num==1))$snppos))
signifGenes = unique(select(signifGenes,snpid,geneid)) %>% 
  filter(!(geneid=='' | is.na(geneid)))
signifSNPs = signifSNPs %>%
  filter(!snpid %in% unique((signifSNPs %>% filter(cluster %in% c(3,4)))$snpid))


xx = signifSNPs %>%
  # filter(!snpid %in% snpexc) %>%
  filter(cl1_cl2 == ' antagonist' & !grepl('agonist',cl1_cl3) & !grepl('agonist', cl2_cl3)) %>%
  select(snpid,disease,cluster,UKBB_RAF) %>%
  na.omit() %>%
  unique()
xxx = xx %>% 
  group_by(snpid,cluster) %>%
  summarise(diseases = paste(sort(unique(disease)),collapse = ', '),
            RAF=unique(UKBB_RAF)) %>%
  ungroup()
xxx1 = xxx %>% select(-RAF) %>%
  spread(cluster,diseases) %>%
  na.omit() 

xxx11 = xxx %>% select(-cluster) %>%
  rename(`1`=diseases,
         RAF1 = RAF) 
xxx12 = xxx %>% select(-cluster) %>%
  rename(`2`=diseases,
         RAF2 = RAF) 
xx = left_join(left_join(xxx1,xxx11),
left_join(xxx1,xxx12)) %>%
  mutate(RAF_difference = RAF1 - RAF2) %>%
  select(-RAF1,-RAF2) %>%
  set_names(c('snpid','cl1','cl2','RAF_diff'))

left_join(xx,signifGenes) %>%
  left_join(mafxxx) %>%
  separate(snpid,into=c('chr','bp','ref','alt'),remove = F) %>%
  # group_by(geneid,cl1,cl2) %>%
  # summarise(nsnp = length(unique(snpid)),
  #           median_RAF_diff = median(RAF_diff),
  #           chr = unique(chr),
  #           bp_st = min(bp),
  #           bp_e = max(bp)) %>%
  # ungroup() %>%
  select(geneid,cl1,cl2,RAF_diff,chr,bp,snpid,betagr) %>%
  # arrange(-abs(median_RAF_diff),-nsnp) %>%
  # mutate(geneid = ifelse(is.na(geneid),'-',geneid)) %>%
  mutate(geneid = fct_reorder(geneid,as.numeric(bp))) %>%
  mutate(chr = as.numeric(chr)) %>%
  na.omit() %>%
  ggplot(aes( x= geneid, y= RAF_diff))+
  geom_hline(color = 'red', linetype = 'dashed', size = 0.2, yintercept = 0) + 
  facet_wrap(~chr, scales = 'free_x', nrow = 2,ncol=6) +
  # geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.5, alpha = 1, aes(color = betagr>=0.3), width = 0.1) +
  scale_color_manual(values = c('gray90','gray25')) +
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) +
  xlab('') + ylab('Difference in Risk Allele Frequencies') +
  guides(color = guide_legend('>= 30 % BETA'))

ggsave('./results/evoAnalysis/antagonisticGenes.pdf', units = 'cm', width = 35, height = 15, useDingbats = F)
ggsave('./results/evoAnalysis/antagonisticGenes.png', units = 'cm', width = 35, height = 15)

left_join(xx,signifGenes) %>%
  left_join(mafxxx) %>%
  separate(snpid,into=c('chr','bp','ref','alt'),remove = F) %>%
  # group_by(geneid,cl1,cl2) %>%
  # summarise(nsnp = length(unique(snpid)),
  #           median_RAF_diff = median(RAF_diff),
  #           chr = unique(chr),
  #           bp_st = min(bp),
  #           bp_e = max(bp)) %>%
  # ungroup() %>%
  select(geneid,cl1,cl2,RAF_diff,chr,bp,snpid,betagr) %>%
  # arrange(-abs(median_RAF_diff),-nsnp) %>%
  # mutate(geneid = ifelse(is.na(geneid),'-',geneid)) %>%
  mutate(geneid = fct_reorder(geneid,as.numeric(bp))) %>%
  mutate(chr = as.numeric(chr)) %>%
  na.omit() %>%
  filter(betagr>=0.3) %>%
  mutate(cl1 = gsub(', ','\n',cl1),
         cl2 = gsub(', ','\n',cl2)) %>%
  group_by(cl1,cl2) %>%
  summarise(RAF= median(RAF_diff),
            geneid = paste(unique(sort(geneid)),collapse ='\n')) %>%
  ggplot(aes(x = cl2, y= cl1, fill =RAF)) +
  geom_tile() +
  scale_fill_gradient2(low = 'midnightblue',mid = 'white',high = 'firebrick3',midpoint = 0) +
  geom_label(aes(label = geneid), size = 4/pntnorm, nudge_y = 0.5, vjust = 1, fill = 'white') + 
  xlab('') + ylab('') + guides(fill = guide_colorbar('Difference in RAF')) +
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        legend.position = 'right') 

ggsave('./results/evoAnalysis/antagonisticDiseases.pdf', units = 'cm', width = 24, height = 12, useDingbats = F)
ggsave('./results/evoAnalysis/antagonisticDiseases.png', units = 'cm', width = 24, height = 12)



left_join(xx,signifGenes) %>%
  left_join(mafxxx) %>%
  separate(snpid,into=c('chr','bp','ref','alt'),remove = F) %>%
  # group_by(geneid,cl1,cl2) %>%
  # summarise(nsnp = length(unique(snpid)),
  #           median_RAF_diff = median(RAF_diff),
  #           chr = unique(chr),
  #           bp_st = min(bp),
  #           bp_e = max(bp)) %>%
  # ungroup() %>%
  select(geneid,cl1,cl2,RAF_diff,chr,bp,snpid,betagr) %>%
  # arrange(-abs(median_RAF_diff),-nsnp) %>%
  # mutate(geneid = ifelse(is.na(geneid),'-',geneid)) %>%
  mutate(geneid = fct_reorder(geneid,as.numeric(bp))) %>%
  mutate(chr = as.numeric(chr)) %>%
  na.omit() %>%
  # filter(betagr>=0.3) %>%
  group_by(cl1,cl2) %>%
  summarise(RAF= median(RAF_diff),
            geneid = paste(unique(sort(geneid)),collapse =', ')) %>%
  select(cl1, cl2, geneid, RAF) %>%
  arrange(-RAF) %>%
  set_names(c('Cluster1 Diseases', 'Cluster2 Diseases', 'Gene ID', 'Median RAF Difference')) %>%
  write_tsv('./results/evoAnalysis/antagonisticDiseaseGenes.tsv')
