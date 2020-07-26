library(HDL)
source('./scripts/00-setup.R')
tr1id = commandArgs(trailingOnly = T)[1]
tr2id = commandArgs(trailingOnly = T)[2]
n=484598 
g1 = read_tsv(paste('./data/processed/ukbb/gwas/bolt/a', tr1id, '.imp.stats',
                    sep = '')) %>%
  mutate(CHR=as.numeric(CHR))%>%
  arrange(CHR,BP) %>%
  select(SNP,ALLELE1,ALLELE0,BETA,SE) %>%
  unique() %>%
  set_names(c('SNP','A1','A2','b','se')) %>%
  mutate(N = n) %>%
  select(SNP,A1,A2,N,b,se) %>%
  as.data.frame()

g2 = read_tsv(paste('./data/processed/ukbb/gwas/bolt/a', tr2id, '.imp.stats',
                    sep = '')) %>%
  mutate(CHR=as.numeric(CHR))%>%
  arrange(CHR,BP) %>%
  select(SNP,ALLELE1,ALLELE0,BETA,SE) %>%
  unique() %>%
  set_names(c('SNP','A1','A2','b','se')) %>%
  mutate(N = n) %>%
  select(SNP,A1,A2,N,b,se) %>%
  as.data.frame()
LD.path <- "./data/raw/hdl_refpanel/UKB_imputed_SVD_eigen99_extraction"
outx = paste('./data/processed/hdl/',tr1id,'_',tr2id,sep='')
res.HDL <- HDL.rg(g1, g2, LD.path, Nref = 335265, N0= n, eigen.cut = 0.99, output.file = paste(outx,'.txt',sep=''))
saveRDS(res.HDL,file = paste(outx,'.rds',sep=''))

