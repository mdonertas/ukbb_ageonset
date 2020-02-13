source('./scripts/00-setup.R')
source('./scripts/LCV/RunLCV.R')

ldscfiles = paste('./data/raw/ldsc/eur_w_ld_chr/', 1:22, '.l2.ldscore', 
                  sep = '')
ldsc = lapply(ldscfiles,function(f)read_tsv(f)) 
ldsc = reshape2::melt(ldsc, id.vars = colnames(ldsc[[1]])) %>%
  filter(MAF > 0.05) %>%
  filter(!(CHR == mhcchr & BP <= mhcend & BP >= mhcstart)) %>%
  select(CHR, BP, SNP, L2) 

tr1id = commandArgs(trailingOnly = T)[1]
tr2id = commandArgs(trailingOnly = T)[2]

tr1 = read_tsv(paste('./data/processed/ukbb/gwas/bolt/a', tr1id, '.imp.stats',
                       sep = '')) %>%
    mutate(tr1 = BETA/SE) %>%
    select(SNP, CHR, BP, ALLELE1, ALLELE0, tr1) 

tr2 = read_tsv(paste('./data/processed/ukbb/gwas/bolt/a', tr2id, '.imp.stats',
                     sep = '')) %>%
  mutate(tr2 = BETA/SE) %>%
  select(SNP, CHR, BP, ALLELE1, ALLELE0, tr2) 

mydat = full_join(tr1,tr2) 

mydat = filter(mydat, 
               !SNP %in% unique(mydat$SNP[which(duplicated(mydat$SNP))])) %>%
  left_join(ldsc) %>%
  na.omit() %>%
  arrange(CHR, BP)

LCV = RunLCV(mydat$L2, mydat$tr1, mydat$tr2)
saveRDS(LCV,file = paste('./data/processed/LCV/', tr1id, '_', tr2id, 
                         '.rds', sep = ''))


