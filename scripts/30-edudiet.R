source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)

edudietfn = list.files('./data/raw/edu_diet_gwas/',full.names = T)
edudietn = sapply(strsplit(edudietfn,'/'),function(x)sapply(strsplit(x[length(x)],'[.]'),function(x)x[1]))
edudiet = lapply(edudietfn,function(x){
  read_tsv(x) %>%
    select(variant,pval) %>%
    filter(pval<=5e-10) %>%
    unique()
})
names(edudiet) = edudietn

edudiet = lapply(edudiet,function(x){
  x %>%
    separate(variant, into = c('CHR','BP','Ref','Alt'), remove = F) %>%
    mutate(posID = paste(CHR,BP,sep='_'))
})

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

snpgrs = group_by(signifSNPs,aooclusters,type) %>%
  summarise(snplist = list(sort(unique(posID)))) %>%
  mutate(snp_type = paste('cluster',aooclusters,'-',type))

snpgrs = setNames(snpgrs$snplist,snpgrs$snp_type)
edudietsnps = sapply(edudiet,function(x)unique(filter(x,CHR%in%1:22)$posID))
ukbbsnps = unique((read_tsv('./data/processed/ukbb/gwas/bolt/a1071.imp.stats') %>%
  mutate(posID = paste(CHR,BP,sep="_")))$posID)
nealesnps = unique((read_tsv(edudietfn[1]) %>%
  separate(variant, into = c('CHR','BP','Ref','Alt'), remove = F) %>%
  mutate(posID = paste(CHR,BP,sep='_')))$posID)
testsnps = intersect(ukbbsnps,nealesnps)
x = lapply(snpgrs,function(x){
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
  y
})
reshape2::melt(x) %>%
  spread(Var2,value) %>%
  filter(a>0)

# [1] Var1 L1   a    b    c    d    odds pval
# <0 rows> (or 0-length row.names)

# There is no trait with at least one overlap with our SNP sets. (not only multicat or multidis but in general with any of the SNP sets for any of the diseases)
