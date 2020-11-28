source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)
load('./data/processed/gwascatcomparison.RData')
gwasassocs = gwasassocs[sapply(gwasassocs,length)>=15]
gwasmap = read_tsv('./data/raw/gwascat/gwas_catalog_trait-mappings_r2020-11-20.tsv')
cats = c('longevity',"Alzheimer's disease", "Parkinson's disease", "age-related macular degeneration","aging")
cancer =  unique(filter(gwasmap, `Parent term` == 'Cancer')$`EFO term`)
aginggenes = gwasassocs[names(gwasassocs) %in% c(cats,cancer)]
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
  rownames(xx) = names(aginggenes)
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
xxx = xx %>%
  filter(Aging %in% cats)
oddsplot = ggplot(xxx, aes(x=cluster, y = prop, color = Aging, size = p<=0.05)) +
  geom_hline(yintercept = 0, color = 'darkred',linetype = 'dashed', size = 0.2, alpha = 0.5) +
  geom_text(data=dplyr::select(xxx,cluster,numgenes,type),aes(label=paste(numgenes,'genes')), y = -0.5, color = 'gray25', size = 6/pntnorm,hjust=0,nudge_x = 0.3) +
  geom_jitter(alpha=0.7,aes(shape=p<=0.05),width = 0.3) + 
  geom_text(data = filter(xxx,p<=0.05),aes(label = num), size = 6/pntnorm) +
  scale_size_manual(values=setNames(c(0.5,2),c(F,T))) +
  scale_shape_manual(values=setNames(c(4,19),c(F,T))) +
  scale_color_gdocs() +
  # scale_color_manual(values = rainbow(5)) + 
  geom_vline(xintercept = seq(1.5,7.5,by=1),color='gray70',linetype='dotted',size=0.3)+
  facet_wrap(~type, ncol = 2, nrow=1) +
  coord_flip() +
  xlab('') +
  ylab('Log2 Enrichment Score') +
  theme_pubr(base_size=6)+
  theme(legend.direction = 'vertical', legend.position = 'right',
        panel.border = element_rect(fill=NA)) +
  guides(size = guide_legend('Significance'), shape=F,
         color = guide_legend('GWAS Catalog')) 

oddsplot

ggsave(filename = './results/agingGWASCat/aginggwascat.pdf',oddsplot,units = 'cm',width = 16,height = 6,useDingbats=F)
ggsave(filename = './results/agingGWASCat/aginggwascat.png',oddsplot,units = 'cm',width = 16,height = 6)

xxx = xx %>%
  filter(Aging %in% cancer)

xxx = xxx %>%
  mutate(signif = p<=0.05) %>%
  select(Aging,L1,signif) %>%
  spread(L1,signif) %>%
  as.data.frame()

rownames(xxx) = xxx$Aging
xxx$Aging = NULL
xxx = as.matrix(xxx + 1)
colnames(xxx) = c('Cluster 1, Multidisease',
                  "Cluster 1, Multicategory",
                  "Cluster 1 & 2, Multidisease",
                  "Cluster 1 & 2, Multicategory",
                  "Cluster 1 & 2 & 3, Multicategory",
                  "Cluster 1 & 3, Multicategory",
                  'Cluster 2, Multidisease',
                  "Cluster 2, Multicategory",
                  "Cluster 2 & 3, Multidisease",
                  "Cluster 2 & 3, Multicategory",
                  'Cluster 3, Multidisease',
                  "Cluster 3, Multicategory")
library(pheatmap)
pheatmap(t(xxx), cellwidth = 10,cellheight = 10, legend = F, filename = './results/agingGWASCat/cancer.pdf')
pheatmap(t(xxx), cellwidth = 10,cellheight = 10, legend = F, filename = './results/agingGWASCat/cancer.png')

saveRDS(permres, './data/processed/aginggwascat_permres.rds')
