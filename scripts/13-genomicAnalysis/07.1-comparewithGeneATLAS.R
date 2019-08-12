source('./scripts/00-setup.R')
filesx = grep('noGeneATLAS.tsv',list.files('./data/processed/GWASvsGeneATLAS/'),v=T,invert = T)
disNames = gsub('.tsv','',filesx)

library(tidyverse)
allres = lapply(disNames,function(dis){
  read_tsv(paste('./data/processed/GWASvsGeneATLAS/',dis,'.tsv',sep=''))
})
allres = reshape2::melt(allres, id.vars = c('Cat','Freq')) %>%
  mutate(Disease = setNames(disCoding$meaning,disCoding$node_id)[as.character(disNames[L1])]) %>%
  select(Disease, Cat, Freq) %>%
  spread(Cat, Freq, fill = 0) %>%
  set_names(c('Disease','NoSig','OneSig','BothSig')) %>%
  mutate(jaccard = (BothSig/(OneSig+BothSig))) 

allres2 = allres %>%
  mutate(disCat = disTreecl[Disease]) %>%
  filter((OneSig + BothSig)>=10) 

p = ggplot(allres2, aes(x = (OneSig + BothSig), y=jaccard, color = disCat)) +
  geom_point(size = 5, alpha = 0.7) +
  scale_x_log10() +
  scale_color_manual(values = discatcolors) +
  guides(color = guide_legend('')) +
  geom_label_repel(data = filter(allres2,jaccard<=0.5),aes(label = Disease)) +
  xlab('total # of significant SNPs') + ylab('Jaccard similarity')
p

ggplot(allres2, aes(x = (OneSig + BothSig), y=jaccard, color = disCat)) +
  geom_text(size = 2, aes(label = Disease)) +
  scale_x_log10() +
  scale_color_manual(values = discatcolors) +
  guides(color = guide_legend('')) +
  # geom_label_repel(data = filter(allres2,jaccard<=0.5),aes(label = Disease)) +
  xlab('total # of significant SNPs') + ylab('Jaccard similarity')
