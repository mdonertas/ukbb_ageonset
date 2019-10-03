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

# varmap <- read_tsv('./data/processed/codingVariants/varmap_output.tsv') %>%
#   filter(AA_CHANGE !='-')
# saveRDS(varmap,'./data/processed/codingVariants/varmap_codingChange.rds')
varmap = readRDS('./data/processed/codingVariants/varmap_codingChange.rds')
signifs = lapply(signifSNPs,function(x){
  mutate(x, RA = ifelse(BETA<0,Alt,Ref))
})
names(signifs)=unname(disCoding[names(signifs)])
signifs = reshape2::melt(signifs,id.vars = colnames(signifs[[1]])) %>%
  rename(disease = L1) %>%
  mutate(CHR = as.numeric(CHR)) %>%
  mutate(cluster = as.factor(clusters$clustering[as.character(disease)]))
signifs = signifs %>%
  rename(CHROMOSOME = CHR, COORDS = BP, USER_BASE = Ref, USER_VARIANT = Alt) %>%
  inner_join(varmap) %>%
  separate(AA_CHANGE,into=c('aa1','aa2'),remove =F,sep='/')
signifs = signifs %>%
  mutate(SNP2 = paste(CHROMOSOME,COORDS,USER_BASE,USER_VARIANT,sep='_')) 

xx = signifs %>%
  filter((!SYNONYMOUS) & UNIPROT_ACCESSION!='-') %>%
  select(SNP2,UNIPROT_ACCESSION) %>%
  unique() 
# How many SNP - protein association: (only missense)
nrow(xx)
# 449
# How many unique SNPs:
length(unique(xx$SNP2))
# 449
# How many unique proteins:
length(unique(xx$UNIPROT_ACCESSION))
# 313

xx = signifs %>%
  filter((!SYNONYMOUS) & UNIPROT_ACCESSION!='-') %>%
  select(SNP2,UNIPROT_ACCESSION,cluster) %>%
  unique() 
xx1 = (xx %>% filter(cluster == 1))
xx2 = (xx %>% filter(cluster == 2))
xx3 = (xx %>% filter(cluster == 3))
# How many SNP - protein association: (only missense)
sapply(list(xx1,xx2,xx3),function(xx)nrow(xx))
# [1] 291 141  87
# How many unique SNPs:
sapply(list(xx1,xx2,xx3),function(xx)length(unique(xx$SNP2)))
# [1] 291 141  87
# How many unique proteins:
sapply(list(xx1,xx2,xx3),function(xx)length(unique(xx$UNIPROT_ACCESSION)))
# [1] 210 103  55

missenseprops = signifs %>% 
  select(SNP2,cluster,SYNONYMOUS) %>% 
  unique() %>% 
  group_by(cluster, SYNONYMOUS) %>%
  summarise(cnt = length(unique(SNP2))) %>%
  ungroup() %>%
  mutate(conseq = ifelse(SYNONYMOUS,'synonymous','missense')) %>%
  select(-SYNONYMOUS) %>%
  spread(conseq,cnt) %>%
  mutate(sum = missense+synonymous) %>%
  gather(key='conseq',value='num',-cluster,-sum) %>%
  mutate(prop = num/sum) %>%
  mutate(cluster=paste('Cluster',cluster,sep=' ')) %>%
  ggplot(aes(x = conseq, y= prop, fill = conseq)) +
  facet_wrap(~cluster)+
  geom_bar(stat='identity') +
  scale_fill_brewer(type='qual', palette= 4) +
  geom_label(aes(label= num),fill='white',vjust=1) +
  xlab('') + ylab('Proportions') +
  guides(fill = F) +
  theme_pubr(x.text.angle = 90, legend = 'right',base_size = 10)

ggsave('./results/codingVar/cl_missenseprops.pdf', missenseprops, units = 'cm', width = 8,height = 8, useDingbats=F)
ggsave('./results/codingVar/cl_missenseprops.png', missenseprops, units = 'cm', width = 8,height = 8)

checkx = signifs %>% 
  select(SNP2,cluster,SYNONYMOUS) %>% 
  unique() %>% 
  group_by(cluster, SYNONYMOUS) %>% summarise(snps = list(unique(SNP2)))
sapply(checkx$snps,function(x){
  sapply(checkx$snps,function(y){
    mean(x%in%y)
  })
})
#      [,1]       [,2]      [,3]       [,4]      [,5]       [,6]
# [1,] 1.00000000 0.00000000 0.3829787 0.00000000 0.1704545 0.00000000
# [2,] 0.00000000 1.00000000 0.0000000 0.31914894 0.0000000 0.08139535
# [3,] 0.17821782 0.00000000 1.0000000 0.00000000 0.1818182 0.00000000
# [4,] 0.00000000 0.14950166 0.0000000 1.00000000 0.0000000 0.11627907
# [5,] 0.04950495 0.00000000 0.1134752 0.00000000 1.0000000 0.00000000
# [6,] 0.00000000 0.02325581 0.0000000 0.07092199 0.0000000 1.00000000

siftpred = signifs %>%
  filter(!SYNONYMOUS) %>%
  select(SNP2, cluster, SIFTS_SCORE) %>%
  filter(SIFTS_SCORE!='-') %>%
  unique() %>%
  mutate(SIFTS_SCORE = as.numeric(as.character(SIFTS_SCORE))) %>%
  mutate(tolerated = c('deleterious','tolerated')[1+(SIFTS_SCORE>0.05)]) %>%
  group_by(cluster,tolerated) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  spread(tolerated,n) %>%
  mutate(sum = deleterious + tolerated) %>%
  gather(key = 'type', value = 'num',-cluster,-sum) %>%
  mutate(prop = num/sum) %>%
  mutate(cluster=paste('Cluster',cluster,sep=' ')) %>%
  ggplot(aes(x = type, y= prop, fill = type)) +
  facet_wrap(~cluster)+
  geom_bar(stat='identity') +
  scale_fill_brewer(type='qual', palette= 4) +
  geom_label(aes(label= num),fill='white',vjust=1) +
  xlab('') + ylab('Proportions') +
  guides(fill =F) +
  theme_pubr(x.text.angle = 90, legend = 'right',base_size = 10)

ggsave('./results/codingVar/cl_siftpred.pdf', siftpred, units = 'cm', width = 8,height = 8, useDingbats=F)
ggsave('./results/codingVar/cl_siftpred.png', siftpred, units = 'cm', width = 8,height = 8)

conseq = ggarrange(missenseprops,siftpred + ggtitle('SIFT Predictions'),ncol=2,labels='auto',align = 'hv')

ggsave('./results/codingVar/cl_conseq.pdf', conseq, units = 'cm', width = 16.7,height = 9, useDingbats=F)
ggsave('./results/codingVar/cl_conseq.png', conseq, units = 'cm', width = 16.7,height = 9)

fisher.test(matrix(c(96,197,23,117),byrow = T, nrow=2))
# Fisher's Exact Test for Count Data
# 
# data:  matrix(c(96, 197, 23, 117), byrow = T, nrow = 2)
# p-value = 0.0003387
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.460374 4.324188
# sample estimates:
# odds ratio 
#     2.4741
fisher.test(matrix(c(96,197,19,64),byrow = T, nrow=2))
# Fisher's Exact Test for Count Data
# 
# data:  matrix(c(96, 197, 19, 64), byrow = T, nrow = 2)
# p-value = 0.105
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.9081521 3.0686504
# sample estimates:
# odds ratio 
#   1.639416 
fisher.test(matrix(c(23,117,19,64),byrow = T, nrow=2))
# Fisher's Exact Test for Count Data
# 
# data:  matrix(c(23, 117, 19, 64), byrow = T, nrow = 2)
# p-value = 0.2878
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.3183666 1.3934549
# sample estimates:
# odds ratio 
#  0.6634622 

signif_summary = signifs %>%
  filter(UNIPROT_ACCESSION!='-') %>%
  select(UNIPROT_ACCESSION, disease, cluster) %>%
  mutate(category = disTreecl[as.character(disease)]) %>%
  group_by(UNIPROT_ACCESSION) %>%
  summarise(n_cluster = length(unique(cluster)),
            n_disease = length(unique(disease)),
            n_category = length(unique(category)),
            clusters = paste(sort(unique(cluster)),collapse = '&'),
            diseases = paste(sort(unique(disease)),collapse = '&'),
            categories = paste(sort(unique(category)),collapse = '&'))

signif_summary %>%
  filter(clusters == '1' & n_category > 1)
# # A tibble: 4 x 7
# UNIPROT_ACCESSION n_cluster n_disease n_category clusters diseases     categories  
# <chr>                 <int>     <int>      <int> <chr>    <chr>        <chr>       
# 1 Q01970                    1         2          2 1        cardiovascu… cardiovascu…
# 2 Q09428                    1         2          2 1        diabetes&hy… cardiovascu…
# 3 Q6YHU6                    1         3          2 1        diabetes&en… cardiovascu…
# 4 Q9BZW4                    1         2          2 1        diabetes&hi… cardiovascu…

cl1uniprot = (signif_summary %>%
                filter(clusters == '1' & n_category > 1))$UNIPROT_ACCESSION

signifs %>%
  filter(UNIPROT_ACCESSION %in% cl1uniprot & CLOSEST_PDB_CODE!='-') %>%
  filter(!SYNONYMOUS) %>%
  select(UNIPROT_ACCESSION, CLOSEST_PDB_CODE, GENE, USER_BASE, USER_VARIANT, RA) %>%
  unique()

xx = signifs %>%
  filter(!SYNONYMOUS & CATH_NAME!='-') %>%
  select(SNP2, cluster,CATH_NAME) %>%
  unique() %>%
  group_by(cluster,CATH_NAME) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  spread(CATH_NAME,n,fill=0) 
xx$sum = rowSums(xx[,-1])
cathpl=xx %>%
  gather(key = 'cath', value='num',-cluster,-sum) %>%
  mutate(prop = num/sum) %>%
  ggplot(aes(x = reorder(cath,prop), y= prop, fill=cluster)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=ageonsetcolors) +
  xlab('') + ylab('Proportions') +
  guides(fill =F) +
  coord_flip() +
  ggtitle('CATH Domains')
rownames(xx) = xx$cluster
xx$cluster=NULL
xx = as.matrix(xx)
xx = rbind(xx,colSums(xx))
rownames(xx)[4]='sum'
cath_enrich = t(apply(xx,2,function(x){
  a = x[1]
  b = xx[1,'sum']-a
  c = x[4]-a
  d = xx['sum','sum'] - a - b - c
  mat = matrix(c(a,b,c,d),byrow = T, nrow = 2)
  fi=fisher.test(mat)
  c(a,b,c,d,fi$est,fi$p.val)
}))
ggsave('./results/codingVar/cl_cath.pdf',cathpl,units='cm',width = 16.7,height = 8,useDingbats=F)
ggsave('./results/codingVar/cl_cath.png',cathpl,units='cm',width = 16.7,height = 8)

xx = signifs %>%
  filter(!SYNONYMOUS & PFAM_DOMAIN!='-') %>%
  select(SNP2, cluster,PFAM_DOMAIN) %>%
  unique() %>%
  group_by(cluster,PFAM_DOMAIN) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  spread(PFAM_DOMAIN,n,fill=0) 
xx$sum = rowSums(xx[,-1])
pfampl=xx %>%
  gather(key = 'pfam', value='num',-cluster,-sum) %>%
  mutate(prop = num/sum) %>%
  ggplot(aes(x = reorder(pfam,prop), y= prop, fill=cluster)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=ageonsetcolors) +
  xlab('') + ylab('Proportions') +
  guides(fill =F) +
  coord_flip()+
  # theme_pubr(base_size = 10, x.text.angle = 90) +
  ggtitle('PFAM Domains') +
  theme(axis.text.y=element_blank())
rownames(xx) = xx$cluster
xx$cluster=NULL
xx = as.matrix(xx)
xx = rbind(xx,colSums(xx))
rownames(xx)[4]='sum'
pfam_enrich = t(apply(xx,2,function(x){
  a = x[1]
  b = xx[1,'sum']-a
  c = x[4]-a
  d = xx['sum','sum'] - a - b - c
  mat = matrix(c(a,b,c,d),byrow = T, nrow = 2)
  fi=fisher.test(mat)
  c(a,b,c,d,fi$est,fi$p.val)
}))
any(p.adjust(pfam_enrich[,6],method='fdr')<0.1)
ggsave('./results/codingVar/cl_pfam.pdf',pfampl,units='cm',width = 8,height = 8,useDingbats=F)
ggsave('./results/codingVar/cl_pfam.png',pfampl,units='cm',width = 8,height = 8)

