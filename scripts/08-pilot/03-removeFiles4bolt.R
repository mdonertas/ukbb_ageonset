savepath <- './data/processed/ukbb/gwas/remove/'

library(tidyverse)
# prep remove file for withdrawn:

wth <- read_tsv('./data/raw/ukbb/withdrawn',col_names = 'withdrawn')

famx <- data.frame(X1=NA,X2=NA,X3=NA,X4=NA,X5=NA,X6=NA)
for(x in list.files('./data/processed/ukbb/gwas/fam4bolt/',full.names = T)){
  famx <- unique(rbind(famx,read_delim(x,col_names = F, delim = ' ')))
}

withdrawn <- famx %>%
  filter(X1 %in% wth$withdrawn)

withdrawn %>% write_delim(paste(savepath,'withdrawn.fam',sep=''),col_names = F)

traits <- readRDS('./data/processed/traits_clean/traitData_baseline.rds')

touse <- c(intersect(traits$eid,famx$X1),NA)

excluded <- famx %>%
  filter(!X1 %in%touse)

excluded %>% write_delim(paste(savepath,'sampleQC_exc.fam',sep=''),col_names = F)

samplex <- data.frame(ID_1=NA,ID_2=NA,missing=NA,sex=NA)
for(x in list.files('./data/raw/ukbb/sample/',full.names = T)){
  samplex <- unique(rbind(samplex,read_delim(x, delim = ' ')))
}

inplink_notin_bgen <- famx %>%
  filter((!X1 %in% samplex$ID_1) | (!X1 %in% samplex$ID_2))


inplink_notin_bgen %>% write_delim(paste(savepath,'inplink_notin_bgen.fam',sep=''),col_names = F)
