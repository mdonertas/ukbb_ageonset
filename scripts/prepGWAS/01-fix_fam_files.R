library(tidyverse)

system('mkdir -p ./data/processed/ukbb/gwas/fam4bolt/')

for( i in 1:22){
  print(i)
  oldfile <- read_delim(paste('./data/raw/ukbb/fam/ukb30688_cal_chr',i,'_v2_s488346.fam',sep=''),delim = ' ',col_names = F)
  oldfile %>%
    mutate(X6 = as.numeric(as.factor(X6))) %>%
    write_delim(paste('./data/processed/ukbb/gwas/fam4bolt/ukb30688_cal_chr',i,'_v2_s488346.fam',sep=''), col_names = F, delim = ' ')
}
