source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)
disCoding = disCoding[as.character(disIDs)]
aooclusters = readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$cluster
allerrx = list.files('./errorfiles_hdl/')
allfilesx = grep('.rds',list.files('./data/processed/hdl/'),v=T)
allpairs = gsub('errorfile.|.txt','',allerrx)
length((setdiff(allpairs ,gsub('.rds','',allfilesx))))
#331 
# there are 331 jobs which did not finish. after a quick check with the errorfiles it seems to be caused by the memory limit. 
# run these pairs again

rerunlist = setdiff(allpairs ,gsub('.rds','',allfilesx))
for(pair in rerunlist){
  pair = strsplit(pair,'_')[[1]]
  i = pair[1]
  k = pair[2]
  print(c(i,k))
  system(paste('bsub -o ./errorfiles_hdl/newerrorfile.',i,'_',k,'.txt -M 24000 -R "rusage[mem=10000]" Rscript ./scripts/26.1-runHDL.R ',i,' ',k,sep=''))
}

# check again
allfilesx = grep('.rds',list.files('./data/processed/hdl/'),v=T)
allpairs = gsub('errorfile.|.txt','',allerrx)
length((setdiff(allpairs ,gsub('.rds','',allfilesx))))
# there are still 11 pairs that cannot be run. 
# program gives the following error:
# "Algorithm failed to converge after trying different initial values." 
# And thus cannot be analysed. 

allresx = lapply(paste('./data/processed/hdl/',allfilesx,sep=''),readRDS)
names(allresx) = gsub('.rds','',allfilesx)
a1 = t(sapply(allresx,function(x)c(x$rg,x$rg.se,x$P))) %>% as.data.frame()
colnames(a1) = c('rg','rg.se','p')
a1$pairs = rownames(a1)
a2 = lapply(allresx,function(x)x$estimates.df) %>%
  reshape2::melt() %>% 
  spread(key = Var1, value = value) 
a3 = filter(a2, Var2 == 'Estimate') %>%
  select(-1) %>%
  set_names(c('pairs','h1','h2','cov','cor'))
a2 = filter(a2, Var2 == 'se') %>%
  select(-1) %>%
  set_names(c('pairs','h1.se','h2.se','cov.se','cor.se'))

restable = full_join(a1,a2) %>% full_join(a3) %>%
  select(pairs, rg, rg.se, everything()) %>%
  separate(pairs, into= c('dis1','dis2'),remove = F)%>%
  mutate(disease1 = disCoding[as.character(dis1)]) %>%
  mutate(disease2 = disCoding[as.character(dis2)]) %>%
  mutate(cluster1 = aooclusters[as.character(disease1)]) %>%
  mutate(cluster2 = aooclusters[as.character(disease2)]) %>%
  mutate(cat1 = disTreecl[as.character(disease1)]) %>%
  mutate(cat2 = disTreecl[as.character(disease2)]) %>%
  mutate(samecat = cat1==cat2) %>%
  mutate(samecl = cluster1==cluster2) %>%
  mutate(padj = p.adjust(p,method='fdr'))

restable$uplevel = sapply(1:nrow(restable),function(i){
  (restable$disease2[i] %in% subcomponent(disTree,as.character(restable$disease1[i]),mode = 'out')$name)
})
restable$downlevel = sapply(1:nrow(restable),function(i){
  (restable$disease1[i] %in% subcomponent(disTree,as.character(restable$disease2[i]),mode = 'out')$name)
})


rm(a1,a2,a3,allresx)

restable %>%
  filter(!uplevel & !downlevel & padj<=0.05 & rg!=Inf) %>%
  # filter(rg>1) %>% 
  ggplot(aes(x = samecl, y= rg))+
  geom_boxplot() +
  stat_compare_means()

saveRDS(restable,'./data/processed/HDL_summary/HDLsummaryTable.rds')
