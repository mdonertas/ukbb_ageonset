library('tidyverse')
vars = read_tsv('./data/processed/ukbb/gwas/bolt/a1071.imp.stats')
vars = vars %>%
  select(CHR,BP,SNP,ALLELE1,ALLELE0) %>%
  mutate(SNP = ifelse(substr(SNP,1,2)=='rs',SNP,'.')) %>%
  set_names(c('Chromosome','Coords','Identifier',	'Original base',	'Variant base'))
write_tsv(vars, './data/processed/codingVariants/varmap_input.vcf')