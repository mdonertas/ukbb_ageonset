# system('mkdir ./data/processed/LCV')
source('./scripts/00-setup.R')
source('./scripts/LCV/RunLCV.R')

disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)
disCoding = disCoding[as.character(disIDs)]
signifSNPs <- lapply(paste('../ukbb_ageonset/data/processed/caseControl/a',
                           disIDs, '/signif_gwasRes_proxyGenes.rds', sep = ''),
                     function(x) {
  x = readRDS(x)
  x = filter(x, !(CHR == mhcchr & BP >= mhcstart & BP <= mhcend))
  select(x,SNP,CHR,BP,Ref,Alt) %>% unique()
})
names(signifSNPs) = disIDs
signifSNPs = sapply(signifSNPs,function(x){
  filter(x, !SNP %in% unique(x$SNP[which(duplicated(x$SNP))])) %>%
    unique() %>%
    nrow()
})
traitids = names(which(signifSNPs >= 10))

for (i in 1:(length(traitids) - 1)) {
  for (k in (i + 1):length(traitids)) {
    system(paste('bsub -o ./errorfiles/errorfile.',traitids[i],'_',traitids[k],'.txt -M 16000 -R "rusage[mem=16000]" RScript ./scripts/25.1-runLCV.R ',traitids[i],' ',traitids[k],sep=''))
  }
}


