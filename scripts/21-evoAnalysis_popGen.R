source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('../ukbb_ageonset/results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)
disCoding = disCoding[as.character(disIDs)]
signifSNPs=readRDS('./data/processed/evoAnalysis/UKBBRAF_Pleiotropy.rds')
signifSNPs = signifSNPs %>%
  filter(cluster !=4) 
signifSNPs = signifSNPs %>%
  mutate(gr_agonist = grepl('agonist',cl1_cl1) |grepl('agonist',cl2_cl2)|grepl('agonist',cl3_cl3),
         gr_antagonist = grepl('antagonist',cl1_cl1) |grepl('antagonist',cl2_cl2)|grepl('antagonist',cl3_cl3),
         ac_agonist = grepl('agonist',cl1_cl2) | grepl('agonist',cl1_cl3) | grepl('agonist',cl2_cl3),
         ac_antagonist = grepl('antagonist',cl1_cl2) | grepl('antagonist',cl1_cl3) | grepl('antagonist',cl2_cl3))
signifSNPs$SNP2 = paste(signifSNPs$CHR,signifSNPs$BP,signifSNPs$Ref,signifSNPs$Alt,sep='_')
snplist = signifSNPs %>%
  select(SNP2,cluster) %>%
  unique() %>%
  group_by(cluster) %>%
  summarise(SNPs = list(unique(sort(SNP2))))
unqcl1 = setdiff(snplist$SNPs[[1]],union(snplist$SNPs[[2]],snplist$SNPs[[3]]))
unqcl2 = setdiff(snplist$SNPs[[2]],union(snplist$SNPs[[1]],snplist$SNPs[[3]]))
unqcl3 = setdiff(snplist$SNPs[[3]],union(snplist$SNPs[[1]],snplist$SNPs[[2]]))
unqcl1 = tapply(unqcl1, sapply(strsplit(unqcl1,'_'),function(x)x[1]), c)
unqcl1 = sort(names(unqcl1))[rank(as.character(1:22))]
unqcl2 = tapply(unqcl2, sapply(strsplit(unqcl2,'_'),function(x)x[1]), c)
unqcl2 = sort(names(unqcl2))[rank(as.character(1:22))]
unqcl3 = tapply(unqcl3, sapply(strsplit(unqcl3,'_'),function(x)x[1]), c)
unqcl3 = sort(names(unqcl3))[rank(as.character(1:22))]
