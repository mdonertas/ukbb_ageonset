source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)

proxyGenes <- sapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_proxyGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  setdiff(unique(x$proxy_hgnc),c('',NA))
})
names(proxyGenes)=disIDs

eqtlGenes <- sapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_eQTLGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  setdiff(unique(x$eQTL_hgnc),c('',NA))
})
names(eqtlGenes)=disIDs

proxyGenes <- reshape2::melt(proxyGenes) %>%
  set_names(c('geneid','disID')) %>%
  mutate(proxy = TRUE)

eqtlGenes <- reshape2::melt(eqtlGenes) %>%
  set_names(c('geneid','disID')) %>%
  mutate(eqtl = TRUE)

signifGenes <- full_join(proxyGenes, eqtlGenes) 

signifGenes <- signifGenes %>%
  mutate(disease = disCoding[as.character(disID)]) %>%
  mutate(disCat = disTreecl[disease],
         ageonset = (readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$cluster)[disease]) 

rm(proxyGenes,eqtlGenes)

ageonsetsum = signifGenes %>%
  group_by(geneid) %>%
  summarise(ageonsetclusters = paste(sort(unique(ageonset)),collapse = '-',sep='-'),
            numagecluster = length(unique(ageonset))) %>%
  ungroup()

discatsum = signifGenes %>%
  group_by(geneid) %>%
  summarise(discategories = paste(sort(unique(disCat)),collapse = ', ',sep=', '),
            numdiscat = length(unique(disCat))) %>%
  ungroup()

dissum = signifGenes %>%
  group_by(geneid) %>%
  summarise(all_diseases = paste(sort(unique(disease)),collapse = ', ',sep=', '),
            numdiseases = length(unique(disease))) %>%
  ungroup()

genedat = full_join(dissum,full_join(ageonsetsum,discatsum)) %>%
  select(geneid, numdiseases,numdiscat,numagecluster,ageonsetclusters, everything()) 

rm(ageonsetsum,discatsum,dissum)

cl1genes= (genedat %>%
                    filter(ageonsetclusters == '1') %>%
                    filter(numdiscat > 0))$geneid

cl2genes = (genedat %>%
                    filter(ageonsetclusters == '2') %>%
                    filter(numdiscat > 0))$geneid

cl3genes= (genedat %>%
                    filter(ageonsetclusters == '3') %>%
                    filter(numdiscat > 0))$geneid

#### Read GWASCat

gwascat <-  read_tsv('../melike/projects/shared_data/GWASCatalog/20190730/data/raw/gwas_catalog_v1.0.2-associations_e96_r2019-07-30.tsv') %>%
  filter(`P-VALUE`<=5e-8)
gwasassocs <- lapply(unique(unlist(strsplit(gwascat$MAPPED_TRAIT,', '))),function(tr){
  rowx=grep(tr,gwascat$MAPPED_TRAIT,ignore.case = T)
  if(length(rowx)>0){
    unique(unlist(strsplit(unlist(strsplit(gwascat$MAPPED_GENE[rowx],', ')),' - ')))
  } else{
    NA
  }
})
names(gwasassocs) <- unique(unlist(strsplit(gwascat$MAPPED_TRAIT,', ')))
saveRDS(gwasassocs,'./data/processed/GWASCat/assocs_20190730.rds')

xx = read_tsv('../melike/projects/shared_data/GWASCatalog/20190607/data/raw/trait_mappings') %>%
  select(-`Disease trait`) %>%
  unique()
parents= unique(xx$`Parent term`)[c(1,2,5,6,7,17)]
# [1] "Other disease"             "Neurological disorder"    
# [3] "Digestive system disorder" "Immune system disorder"   
# [5] "Cardiovascular disease"    "Metabolic disorder"
xx = filter(xx, `Parent term`%in%parents)

ardannots = read_csv('/nfs/research1/thornton/melike/projects/ARDNet/data/processed/ARDannots/GWAScat_diseases.csv') 

ardannots = full_join(xx,ardannots,by= "EFO URI") %>%
  select(1,3,8) %>%
  set_names(c('gwascat_disease','parent_disease','ARD'))
rm(xx)

filter(ardannots,is.na(ARD)) %>% write_tsv('./data/processed/GWASCat/toannotate.tsv')

diseaselist = ardannots$gwascat_disease
gwasdisassocs = gwasassocs[which(names(gwasassocs) %in% diseaselist)]
gwasdisassocs = reshape2::melt(gwasdisassocs) %>%
  set_names(c('geneid','gwascat_disease'))

gwasdisassocs = left_join(gwasdisassocs,ardannots) %>% 
  select(geneid,gwascat_disease,parent_disease,ARD) %>%
  unique()

calcprop = function(genelist,maybe=F){
  xx = filter(gwasdisassocs,geneid %in% genelist) %>%
    select(gwascat_disease,ARD) %>%
    unique()
  ifelse(maybe,mean(xx$ARD %in% c('Maybe','Yes'), na.rm=T), mean(xx$ARD == 'Yes', na.rm=T))
}

cl1prop = calcprop(cl1genes)
cl1prop
# 0.3813559
cl1prop2 = calcprop(cl1genes,maybe = T)
cl1prop2
# 0.519084
cl1perm = sapply(1:1000, function(i){calcprop(sample(unique(gwasdisassocs$geneid),length(intersect(gwasdisassocs$geneid,cl1genes))))})
cl1perm2 = sapply(1:1000, function(i){calcprop(sample(unique(gwasdisassocs$geneid),length(intersect(gwasdisassocs$geneid,cl1genes))),maybe=T)})
mean(cl1perm>=cl1prop)
# 0
mean(cl1perm2>=cl1prop2)
# 0

cl2prop = calcprop(cl2genes)
cl2prop
# 0.3181818
cl2prop2 = calcprop(cl2genes,maybe = T)
cl2prop2
# 0.5368421
cl2perm = sapply(1:1000, function(i){calcprop(sample(unique(gwasdisassocs$geneid),length(intersect(gwasdisassocs$geneid,cl2genes))))})
cl2perm2 = sapply(1:1000, function(i){calcprop(sample(unique(gwasdisassocs$geneid),length(intersect(gwasdisassocs$geneid,cl2genes))),maybe=T)})
mean(cl2perm>=cl2prop)
# 0.296
mean(cl2perm2>=cl2prop2)
# 0.029

cl3prop = calcprop(cl3genes)
cl3prop
# 0.25
cl3prop2 = calcprop(cl3genes,maybe = T)
cl3prop2
# 0.4337349
cl3perm = sapply(1:1000, function(i){calcprop(sample(unique(gwasdisassocs$geneid),length(intersect(gwasdisassocs$geneid,cl3genes))))})
cl3perm2 = sapply(1:1000, function(i){calcprop(sample(unique(gwasdisassocs$geneid),length(intersect(gwasdisassocs$geneid,cl3genes))),maybe=T)})
mean(cl3perm>=cl3prop)
# 0.963
mean(cl3perm2>=cl3prop2)
# 0.949

## GEnes with more than one category
cl1genes_h1cat = (genedat %>%
                    filter(ageonsetclusters == '1') %>%
                    filter(numdiscat > 1))$geneid

cl2genes_h1cat = (genedat %>%
                    filter(ageonsetclusters == '2') %>%
                    filter(numdiscat > 1))$geneid

cl3genes_h1cat = (genedat %>%
                    filter(ageonsetclusters == '3') %>%
                    filter(numdiscat > 1))$geneid

cl1prop_h1 = calcprop(cl1genes_h1cat)
cl1prop_h1
# 0.375
cl1prop2_h1 = calcprop(cl1genes_h1cat,maybe = T)
cl1prop2_h1
# 0.5853659
cl1perm_h1 = sapply(1:1000, function(i){calcprop(sample(unique(gwasdisassocs$geneid),length(intersect(gwasdisassocs$geneid,cl1genes_h1cat))))})
cl1perm2_h1 = sapply(1:1000, function(i){calcprop(sample(unique(gwasdisassocs$geneid),length(intersect(gwasdisassocs$geneid,cl1genes_h1cat))),maybe=T)})
mean(cl1perm_h1>=cl1prop_h1)
# 0.085
mean(cl1perm2_h1>=cl1prop2_h1)
# 0.085

cl2prop_h1 = calcprop(cl2genes_h1cat)
cl2prop_h1
# 0.1904762
cl2prop2_h1 = calcprop(cl2genes_h1cat,maybe = T)
cl2prop2_h1
# 0.4545455
cl2perm_h1 = sapply(1:1000, function(i){calcprop(sample(unique(gwasdisassocs$geneid),length(intersect(gwasdisassocs$geneid,cl2genes_h1cat))))})
cl2perm2_h1 = sapply(1:1000, function(i){calcprop(sample(unique(gwasdisassocs$geneid),length(intersect(gwasdisassocs$geneid,cl2genes_h1cat))),maybe=T)})
mean(cl2perm_h1>=cl2prop_h1)
# 0.563
mean(cl2perm2_h1>=cl2prop2_h1)
# 0.538

cl3prop_h1 = calcprop(cl3genes_h1cat)
cl3prop_h1
# 0.1956522
cl3prop2_h1 = calcprop(cl3genes_h1cat,maybe = T)
cl3prop2_h1
# 0.4313725
cl3perm_h1 = sapply(1:1000, function(i){calcprop(sample(unique(gwasdisassocs$geneid),length(intersect(gwasdisassocs$geneid,cl3genes_h1cat))))})
cl3perm2_h1 = sapply(1:1000, function(i){calcprop(sample(unique(gwasdisassocs$geneid),length(intersect(gwasdisassocs$geneid,cl3genes_h1cat))),maybe=T)})
mean(cl3perm_h1>=cl3prop_h1)
# 0.923
mean(cl3perm2_h1>=cl3prop2_h1)
# 0.895

save(list = ls(), file = './data/processed/agingrelevance/gwascatresults.RData')
