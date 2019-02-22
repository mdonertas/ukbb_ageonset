library(tidyverse)

hwe.exc=c()
for(i in 1:22){
  hwe.exc <- unique(c(hwe.exc,(read_tsv(paste('./data/processed/ukbb/gwas/hwe/chr',i,'.hardy',sep='')) %>%
                                 filter(P<=1e-6))$ID))
}

length(hwe.exc)
# [1] 202473

system('mkdir -p ./data/processed/ukbb/gwas/exclude')

data.frame(SNPID=hwe.exc)%>%
  write_delim('./data/processed/ukbb/gwas/exclude/hwe.txt',delim=' ', col_names = F)

maf.exc=c()
for(i in 1:22){
  maf.exc <- unique(c(maf.exc,(read_tsv(paste('./data/processed/ukbb/gwas/freq/chr',i,'.afreq',sep='')) %>%
                                 mutate(maf=ifelse(ALT_FREQS<0.5,ALT_FREQS,1-ALT_FREQS))%>%
                                 filter(maf<0.01))$ID))
}

length(maf.exc)
# [1] 127969

length(union(maf.exc,hwe.exc))
# [1] 314697

data.frame(SNPID=maf.exc)%>%
  write_delim('./data/processed/ukbb/gwas/exclude/maf.txt',delim=' ', col_names = F)
