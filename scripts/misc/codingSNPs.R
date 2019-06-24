library(tidyverse)
xx = read_tsv('./data/processed/codingVariants/varmap_output.tsv')
which(xx$UNIPROT_ACCESSION!='-')[1:10]
xx[which(xx$UNIPROT_ACCESSION!='-')[1:10],] %>%
  View()
mean(xx$CODON_CHANGE!='-')
codingvar = xx %>%
  filter(CODON_CHANGE != '-')

data.frame(coding_variants = unique((codingvar %>%filter(SYNONYMOUS==F))$USER_ID)) %>%
  filter(coding_variants !='.') %>% 
  write_tsv('./data/processed/codingVariants/codingvars_synremoved.tsv')
