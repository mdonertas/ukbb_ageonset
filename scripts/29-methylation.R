source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)

## methylation data
# adelman2019 - tableS3 (hg19)
adelman2019 = setdiff(unique(read_csv('./data/raw/methylation/adelman2019.csv')$Gene_Symbol),NA)
# vandiver2015 - tableS4 (hg19)
vandiver2015 = read_csv('./data/raw/methylation/vandiver2015.csv')
# johansson2013 - tableS1 (hg19)
johansson2013 = read_csv('./data/raw/methylation/johansson2013.csv') 
# marttila2015 - tableS1
marttila2015 = setdiff(unique(unlist(strsplit(read_csv('./data/raw/methylation/marttila2015.csv')$UCSC_REFGENE_ACCESSION,';'))),NA)

martx = biomaRt::useEnsembl('ensembl','hsapiens_gene_ensembl')
marttila2015_hgnc = biomaRt::getBM(c('hgnc_symbol'),filters = 'refseq_mrna', values = marttila2015, mart = martx)
marttila2015_hgnc = setdiff(unique(marttila2015_hgnc$hgnc_symbol),c(NA,''))

aginggenes = list(adelman = adelman2019, marttila = marttila2015_hgnc, combined = unique(c(adelman2019,marttila2015_hgnc)))
sapply(aginggenes,length)
# adelman marttila combined 
# 748     3278     3893 
# saveRDS(aginggenes,'./data/processed/methylationgenes_hgnc.rds')
## disease data
genegrs = readRDS('./data/processed/clustergenes_hgnc.rds')
eqtlgenes = readRDS('./data/processed/caseControl/a1071/gwasRes_eQTLGenes.rds')
eqtlgenes = setdiff(unique(eqtlgenes$eQTL_hgnc),NA)
proxygenes = readRDS('./data/processed/caseControl/a1071/gwasRes_proxyGenes.rds')
proxygenes = setdiff(unique(proxygenes$proxy_hgnc),NA)
allgenes = unique(c(eqtlgenes,proxygenes))
genegr = setNames(genegrs$genelist, genegrs$name)

permres = lapply(genegr,function(x){
  xx = t(sapply(aginggenes,function(y){
    xx = sapply(1:10000,function(i){100*mean(sample(allgenes,length(x))%in%y)})
    c(perm=mean(xx),val=100*mean(x%in%y),p=mean(xx>=(100*mean(x%in%y))),num = sum(x%in%y),genesize = length(x), agsize = length(intersect(y,allgenes)))
  }))
  rownames(xx) = c('Adelman2019', 'Marttila2015', 'Combined')
  xx
})

xx = reshape2::melt(permres) %>% 
  spread(Var2,value) %>% 
  # mutate(cluster = sapply(strsplit(L1,'_'),function(x)x[1])) %>%
  # mutate(type = sapply(strsplit(L1,'_'),function(x)paste(x[2],x[3],sep="_"))) %>%
  # separate(L1,into=c('cluster','more than one disease', 'more than one category'),sep = '_',remove = F) %>% 
  # mutate(type=ifelse(is.na(type),'all',type)) %>%
  separate(L1,into=c('cluster','type'),sep='_',remove=F) %>%
  rename(Aging = Var1) %>%
  mutate(prop = log2(val/perm))%>%
  left_join(rename(select(ungroup(genegrs),name,numgenes),L1=name))

xx = xx %>% 
  mutate(cluster=gsub('cl','',cluster),
         type = c('Multidisease','Multicategory')[as.numeric(as.factor(type))]) %>%
  mutate(cluster = factor(cluster, levels = c('1','2','3','1-2','1-3','2-3','1-2-3'))) %>%
  mutate(type = factor(type,levels = c('Multidisease','Multicategory'))) %>%
  mutate(prop = ifelse(prop==-Inf,NA,prop)) %>%
  mutate(cluster = factor(gsub('-',' & ',cluster),levels = rev(c('1 & 2 & 3', '2 & 3', '1 & 3','1 & 2','3','2','1'))))

library(ggthemes)
oddsplot = ggplot(xx, aes(x=cluster, y = prop, color = Aging, size = p<=0.05)) +
  geom_hline(yintercept = 0, color = 'darkred',linetype = 'dashed', size = 0.2, alpha = 0.5) +
  geom_text(data=dplyr::select(xx,cluster,numgenes,type),aes(label=paste(numgenes,'genes')), y = -1.5, color = 'gray25', size = 6/pntnorm,hjust=0,nudge_x = 0.3) +
  geom_jitter(alpha=0.7,aes(shape=p<0.05),width = 0.3) + 
  geom_text(data = filter(xx,p<=0.05),aes(label = num), size = 6/pntnorm) +
  scale_size_manual(values=setNames(c(0.5,2),c(F,T))) +
  scale_shape_manual(values=setNames(c(4,19),c(F,T))) +
  scale_color_gdocs() +
  geom_vline(xintercept = seq(1.5,7.5,by=1),color='gray70',linetype='dotted',size=0.3)+
  facet_wrap(~type, ncol = 2, nrow=1) +
  coord_flip() +
  xlab('') +
  ylab('Log2 Enrichment Score') +
  theme_pubr(base_size=6)+
  theme(legend.direction = 'horizontal', legend.position = 'bottom',
        panel.border = element_rect(fill=NA)) +
  guides(size = guide_legend('Significance'), shape=F,
         color = guide_legend('Methylation data')) +
  ylim(-1.5,1.2)

ggsave(filename = './results/methylation/overlap.pdf',oddsplot,units = 'cm',width = 16,height = 6,useDingbats=F)
ggsave(filename = './results/methylation/overlap.png',oddsplot,units = 'cm',width = 16,height = 6)


