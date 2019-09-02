source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)
disCoding = disCoding[as.character(disIDs)]
signifSNPs <- lapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_proxyGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  select(x,SNP,CHR,BP,Ref,Alt,BETA,P_BOLT_LMM_INF) %>% unique()
})
names(signifSNPs)=disIDs
signifSNPs = signifSNPs[which(sapply(signifSNPs,nrow) >0)]
clusters = readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')

varmap <- read_tsv('./data/processed/codingVariants/varmap_output.tsv') %>%
  filter(AA_CHANGE !='-')
saveRDS(varmap,'./data/processed/codingVariants/varmap_codingChange.rds')

signifs = lapply(signifSNPs,function(x){
  mutate(x, RA = ifelse(BETA<0,Alt,Ref))
})

names(signifs)=unname(disCoding[names(signifs)])

signifs = reshape2::melt(signifs,id.vars = colnames(signifs[[1]])) %>%
  rename(disease = L1) %>%
  mutate(CHR = as.numeric(CHR)) %>%
  mutate(aoocl = as.factor(clusters$clustering[as.character(disease)]))

xx = signifs %>%
  rename(CHROMOSOME = CHR, COORDS = BP, USER_BASE = Ref, USER_VARIANT = Alt) %>%
  inner_join(varmap) %>%
  separate(AA_CHANGE,into=c('aa1','aa2'),remove =F,sep='/')

disdat = xx %>% 
  select(CHROMOSOME,COORDS,USER_BASE,USER_VARIANT,RA,disease,aoocl, AA_CHANGE,aa1,aa2,SYNONYMOUS) %>% 
  group_by(disease, SYNONYMOUS) %>%
  summarise(cnt = n()) %>%
  spread(SYNONYMOUS,cnt,fill=0) %>%
  set_names(c('disease','missense','syn')) %>%
  mutate(sum = sum(missense, syn)) %>%
  mutate(misperc = missense/syn) %>%
  mutate(aoocl = as.factor(clusters$clustering[as.character(disease)]))

missyn_cl = disdat %>%
  filter(sum>=5) %>%
  ggplot(aes(x = aoocl, y = misperc, fill = aoocl)) +
  geom_hline(yintercept = 1,  color = 'gray60', linetype = 'dashed') +
  geom_boxplot()+
  geom_jitter() +
  scale_fill_manual(values = ageonsetcolors) +
  scale_y_continuous(trans = 'log2') +
  xlab('Age of Onset Clusters') + ylab('Missense / Synonymous Variant') +
  guides(fill =F) 
  

dismissense = disdat %>%
  filter(sum>=5) %>%
  ggplot(aes(x = reorder(disease,-misperc), y = misperc, fill = aoocl)) +
  geom_hline(yintercept = 1,  color = 'gray60', linetype = 'dashed') +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = ageonsetcolors) +
  coord_flip() +
  scale_y_continuous(trans = 'log2') +
  xlab('') + ylab('Missense / Synonymous Variant') +
  guides(fill = guide_legend('Age of Onset Cluster'))

cldat = xx %>% 
  select(CHROMOSOME,COORDS,USER_BASE,USER_VARIANT,RA,aoocl, AA_CHANGE,aa1,aa2,SYNONYMOUS) %>% 
  group_by(aoocl, SYNONYMOUS) %>%
  summarise(cnt = n()) %>%
  spread(SYNONYMOUS,cnt,fill=0) %>%
  set_names(c('aoocl','missense','syn')) %>%
  mutate(sum = sum(missense, syn)) %>%
  mutate(misperc = missense/syn) 

clmissense = cldat %>%
  ggplot(aes(x = aoocl, y = misperc, fill = aoocl)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = ageonsetcolors) +
  xlab('Age of Onset Cluster') + ylab('Missense / Synonymous Variant') + guides(fill =F) +
  scale_y_continuous(trans ='log2', labels = c(1,1.15),breaks = c(1,1.15)) +
  geom_hline(yintercept = 1,  color = 'gray60', linetype = 'dashed') 

clsift = xx %>% 
  select(CHROMOSOME,COORDS,USER_BASE,USER_VARIANT,RA,aoocl, AA_CHANGE,aa1,aa2,SYNONYMOUS, SIFTS_SCORE) %>% 
  unique() %>%
  select(aoocl,SIFTS_SCORE) %>%
  mutate(SIFTS_SCORE = as.numeric(as.character(SIFTS_SCORE))) %>%
  ggplot(aes(x = aoocl, y= SIFTS_SCORE,fill=aoocl)) +
  geom_boxplot() +
  geom_sina(size =0.5) +
  scale_fill_manual(values = ageonsetcolors) + guides(fill=F) +
  xlab('Age of Onset Cluster') + ylab('SIFTS Score')

clpolyphen = xx %>% 
  select(CHROMOSOME,COORDS,USER_BASE,USER_VARIANT,RA,aoocl, AA_CHANGE,aa1,aa2,SYNONYMOUS, POLYPHEN_SCORE) %>% 
  unique() %>%
  select(aoocl,POLYPHEN_SCORE) %>%
  mutate(POLYPHEN_SCORE = as.numeric(as.character(POLYPHEN_SCORE))) %>%
  ggplot(aes(x = aoocl, y= POLYPHEN_SCORE,fill=aoocl)) +
  geom_boxplot(outlier.shape = NA) +
  geom_sina(size =0.5) +
  scale_fill_manual(values = ageonsetcolors) + guides(fill=F) +
  xlab('Age of Onset Cluster') + ylab('Polyphen Score') +
  scale_y_sqrt()

clcath = xx %>%
  filter(!SYNONYMOUS) %>%
  select(aoocl, SNP, CATH_NAME) %>%
  unique() %>%
  group_by(aoocl,CATH_NAME) %>%
  summarise(cnt = length(unique(SNP))) %>%
  na.omit() %>%
  filter(CATH_NAME!='-') %>%
  spread(aoocl,cnt,fill=0)%>%
  as.data.frame
rownames(clcath) = clcath$CATH_NAME
clcath$CATH_NAME = NULL
clcath = as.matrix(clcath)

clcath = t(apply(clcath[which(rowSums(clcath)>1),],1,function(x)x/sum(x)))
pheatmap::pheatmap(clcath, color = brewer.pal(8,'Reds'))

clpfam = xx %>%
  filter(!SYNONYMOUS) %>%
  select(aoocl, SNP, PFAM_NAME) %>%
  unique() %>%
  group_by(aoocl,PFAM_NAME) %>%
  summarise(cnt = length(unique(SNP))) %>%
  na.omit() %>%
  filter(PFAM_NAME!='-') %>%
  spread(aoocl,cnt,fill=0)%>%
  as.data.frame
rownames(clpfam) = clpfam$PFAM_NAME
clpfam$PFAM_NAME = NULL
clpfam = as.matrix(clpfam)

clpfam = t(apply(clpfam[which(rowSums(clpfam)>1),],1,function(x)x/sum(x)))
pheatmap::pheatmap(clpfam, color = brewer.pal(8,'Reds'))

xx %>%
  select(SNP,CLOSEST_PDB_CODE,RES_NUM, aoocl) %>%
  unique() %>% filter(RES_NUM!='-') %>%
  group_by(aoocl,CLOSEST_PDB_CODE) %>%     
  summarise(cnt = n(  )) %>% filter(CLOSEST_PDB_CODE %in% c('2n4i','2ziy','2nbi')) %>%
  spread(aoocl,cnt,fill=0)
  # filter(aoocl == 3) %>% summary()

# figures
missyn_cl
dismissense
clmissense
clsift
clpolyphen
pheatmap::pheatmap(clcath, color = brewer.pal(8,'Reds'))
pheatmap::pheatmap(clpfam, color = brewer.pal(8,'Reds'))
