library(tidyverse)
disname = commandArgs(trailingOnly = T)
geneatlasfile = paste('../melike/projects/ARD-GWAS/data/raw/assocs/selfReported_n_',disname,'.csv',sep='')
if(file.exists(geneatlasfile)){
  gn1 = read_delim(geneatlasfile, delim=' ') %>%
    set_names('SNP','allele','score','beta','se','p') %>%
    select(SNP,p)
  gn2 = read_tsv(paste('./data/processed/ukbb/gwas/bolt/a',disname,'.imp.stats',sep='')) %>%
    rename(pgwas = P_BOLT_LMM_INF) %>%
    select(SNP, pgwas)
  
  gn1 = left_join(gn2,gn1) %>%
    na.omit() %>%
    mutate(gwassig = pgwas<=5e-8,
           geneatsig = p<=5e-8) %>%
    mutate(tot = gwassig + geneatsig)
  res = table(gn1$tot)
  as.data.frame(res) %>%
    set_names(c('Cat','Freq')) %>%
    write_tsv(paste('/nfs/research1/thornton/ukbb_ageonset/data/processed/GWASvsGeneATLAS/',disname,'.tsv',sep=''))
} else {
  write_tsv(NULL,paste('/nfs/research1/thornton/ukbb_ageonset/data/processed/GWASvsGeneATLAS/',disname,'noGeneATLAS.tsv',sep=''))
}

