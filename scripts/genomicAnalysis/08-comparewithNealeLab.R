source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))

neale = list.files('../melike/projects/shared_data/NealeLab_UKBB/data/processed/selfreported_bothsexes/',full.names = T)
nealecodes = sapply(strsplit(sapply(strsplit(sapply(strsplit(neale,'/'),function(x)x[10]),'_'),function(x)x[2]),'[.]'),function(x)x[1])
nealediseases = unname(setNames(disCoding$meaning,disCoding$coding)[nealecodes])
diseasenames = unname(setNames(disCoding$meaning,disCoding$node_id)[disIDs])

nealeres = lapply(neale[which(nealediseases %in% diseasenames)],function(x)
  {read_tsv(x) %>%
    mutate(CHR = as.character(CHR))})
names(nealeres) = nealediseases[which(nealediseases %in% diseasenames)]
nealeres = nealeres[which(sapply(nealeres,nrow)>5)]

signifSNPs <- lapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_proxyGenes.rds',sep=''),function(x){
  x=readRDS(x)
  select(x, SNP,CHR,BP,Ref,Alt,BETA,P_BOLT_LMM_INF) %>%
    mutate(CHR = as.character(CHR))
})
names(signifSNPs)=diseasenames
signifSNPs = signifSNPs[which(sapply(signifSNPs,nrow)>5)]
diseases=intersect(names(nealeres),names(signifSNPs))

allres = lapply(diseases,function(dis){
  full_join(nealeres[[dis]],signifSNPs[[dis]]) %>%
    mutate(beta = -beta) %>%
    rename(neale_beta = beta,
           neale_pval = pval,
           bolt_beta = BETA,
           bolt_pval = P_BOLT_LMM_INF) %>%
    filter((CHR==mhcchr & (BP<=mhcend & BP>=mhcstart)))
})
names(allres) = diseases
allres = allres[which(sapply(allres,nrow)>5)]
rm(neale,nealecodes,nealediseases,diseasenames,nealeres,signifSNPs)

allres_summary = lapply(allres,function(x){reshape2::melt(table(neale = sign(x$neale_beta),bolt = sign(x$bolt_beta),useNA='always'))%>%filter(value>0)})

allres_summary = reshape2::melt(allres_summary, id.vars = colnames(allres_summary[[1]])) %>%
  rename(disease = L1)

allres_summary = allres_summary %>%
  mutate(sign = sign(neale) * sign(bolt),
         inneale = !is.na(neale),
         inbolt = !is.na(bolt))

allres_summary %>%
  filter(sign == -1)
# nothing is found in opposite signs

jaccards = allres_summary %>%
  mutate(inboth = inneale * inbolt) %>%
  group_by(disease, inboth) %>%
  summarise(val = sum(value)) %>%
  spread(inboth,val,fill=0) %>%
  mutate(jaccard = `1`/(`1`+`0`)) %>%
  mutate(disCat = unname(disTreecl[disease]))

nealeresp = jaccards %>% mutate(totsig = `1`+`0`) %>%
  ggplot(aes(x = totsig, y=jaccard, color = disCat)) +
  geom_point() +
  ggtitle('Jaccard Similarity Index between BOLT-LMM and Neale Lab results') +
  ylim(0,1) +
  geom_point(size = 3, alpha = 0.7) +
  scale_x_log10() +
  scale_color_manual(values = discatcolors) +
  guides(color = guide_legend('')) +
  xlab('total # of significant SNPs') + ylab('Jaccard similarity') +
  theme_pubr(base_size = 8)

# this comes from comparison with GeneAtlas
geneatlasp = ggplot(allres2, aes(x = (OneSig + BothSig), y=jaccard, color = disCat)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_x_log10() +
  scale_color_manual(values = discatcolors) +
  guides(color = guide_legend(''),label = F) +
  geom_label_repel(data = filter(allres2,jaccard<=0.5),aes(label = Disease),size = 6/pntnorm) +
  ggtitle('Jaccard Similarity Index between BOLT-LMM and GeneATLAS results') +
  xlab('total # of significant SNPs') + ylab('Jaccard similarity') +
  ylim(0,1)+
  theme_pubr(base_size = 8)

p = ggarrange(geneatlasp, nealeresp, ncol=1,nrow = 2,common.legend = T,legend = 'right')
ggsave('./results/genomicAnalysis/compare_with_Neale_GeneATLAS.pdf',p,units = 'cm',width = 18, height = 16, useDingbats =F)
ggsave('./results/genomicAnalysis/compare_with_Neale_GeneATLAS.png',p,units = 'cm',width = 18, height = 16)
