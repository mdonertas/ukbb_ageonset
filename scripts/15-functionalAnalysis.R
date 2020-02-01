source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)

proxyGenes <- sapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_proxyGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  setdiff(unique(x$proxy_ensembl),c('',' ',NA))
})
names(proxyGenes)=disIDs

eqtlGenes <- sapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_eQTLGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  setdiff(unique(x$eQTL_ensembl),c('',' ', NA))
})
names(eqtlGenes)=disIDs

proxyGenes <- reshape2::melt(proxyGenes) %>%
  set_names(c('geneid','disID')) %>%
  mutate(proxy = TRUE)

eqtlGenes <- reshape2::melt(eqtlGenes) %>%
  set_names(c('geneid','disID')) %>%
  mutate(eqtl = TRUE)

signifGenes <- full_join(proxyGenes, eqtlGenes) 

signifGenes <- signifGenes %>%
  mutate(disease = disCoding[as.character(disID)]) %>%
  mutate(disCat = disTreecl[disease],
         ageonset = (readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$cluster)[disease]) 

rm(proxyGenes,eqtlGenes)

ageonsetsum = signifGenes %>%
  filter(ageonset!='4') %>%
  group_by(geneid) %>%
  summarise(ageonsetclusters = paste(sort(unique(ageonset)),collapse = '-',sep='-'),
            numagecluster = length(unique(ageonset))) %>%
  ungroup()

discatsum = signifGenes %>%
  group_by(geneid) %>%
  summarise(discategories = paste(sort(unique(disCat)),collapse = ', ',sep=', '),
            numdiscat = length(unique(disCat))) %>%
  ungroup()

dissum = signifGenes %>%
  group_by(geneid) %>%
  summarise(all_diseases = paste(sort(unique(disease)),collapse = ', ',sep=', '),
            numdiseases = length(unique(disease))) %>%
  ungroup()

genedat = full_join(dissum,full_join(ageonsetsum,discatsum)) %>%
  select(geneid, numdiseases,numdiscat,numagecluster,ageonsetclusters, everything()) 

rm(ageonsetsum,discatsum,dissum)

genegrs = genedat %>%
  mutate(monedis = numdiseases>1,
         monecat = numdiscat>1) %>%
  group_by(ageonsetclusters,monedis,monecat) %>%
  summarise(genelist = list(sort(unique(geneid))),
            numgenes = length(unique(geneid))) %>%
  filter(ageonsetclusters %in% c('1','2','3','1-2','1-2-3','1-3','2-3'))%>%
  filter(monedis | monecat) %>%
  mutate(name = paste('cl',ageonsetclusters,'_',c('Multidisease','Multicategory')[1+monecat],sep=''))

genegrs[1,]$genelist = list(unique(unlist(c(genegrs[1,]$genelist,genegrs[2,]$genelist))))
genegrs[1,]$numgenes = genegrs[1,]$numgenes + genegrs[2,]$numgenes

genegrs[3,]$genelist = list(unique(unlist(c(genegrs[3,]$genelist,genegrs[4,]$genelist))))
genegrs[3,]$numgenes = genegrs[3,]$numgenes + genegrs[4,]$numgenes

genegrs[7,]$genelist = list(unique(unlist(c(genegrs[7,]$genelist,genegrs[8,]$genelist))))
genegrs[7,]$numgenes = genegrs[7,]$numgenes + genegrs[8,]$numgenes

genegrs[9,]$genelist = list(unique(unlist(c(genegrs[9,]$genelist,genegrs[10,]$genelist))))
genegrs[9,]$numgenes = genegrs[9,]$numgenes + genegrs[10,]$numgenes

genegrs[11,]$genelist = list(unique(unlist(c(genegrs[11,]$genelist,genegrs[12,]$genelist))))
genegrs[11,]$numgenes = genegrs[11,]$numgenes + genegrs[12,]$numgenes

saveRDS(genegrs,file='./data/processed/clustergenes_ensembl.rds')

eqtlgenes = readRDS('./data/processed/caseControl/a1071/gwasRes_eQTLGenes.rds')
eqtlgenes = setdiff(unique(eqtlgenes$eQTL_ensembl),NA)
proxygenes = readRDS('./data/processed/caseControl/a1071/gwasRes_proxyGenes.rds')
proxygenes = setdiff(unique(proxygenes$proxy_ensembl),NA)

allgenes = unique(c(eqtlgenes,proxygenes))
rm(eqtlgenes,proxygenes)
allgenes = setdiff(allgenes,c(""," ",NA))

genegr = setNames(genegrs$genelist,genegrs$name)
library(goseq)
allgo=getgo(allgenes,"hg19","ensGene")
allgo = allgo[!sapply(allgo,is.null)]
allgo = reshape2::melt(allgo) %>%
  set_names(c('category','gene'))

GOres = lapply(genegr,function(x){
  gnsl = setNames(rep(0,length(allgenes)),allgenes)
  gnsl[names(gnsl)%in%x]=1
  pwf=nullp(gnsl,"hg19","ensGene")
  goseq(pwf,"hg19","ensGene")
})

GOres2 = reshape2::melt(GOres,id.vars = colnames(GOres[[1]])) %>%
  rename(L1='name') %>%
  left_join(dplyr::select(genegrs,numgenes,name)) %>%
  rename(numgenes= 'numGenes') %>%
  separate(name,into=c('cluster','type'),sep='_',remove = F) %>%
  mutate(allgenenum = length(allgenes)) %>%
  filter(numInCat >= 10 & numInCat <=500) %>%
  mutate(over_rep_padj = p.adjust(over_represented_pvalue, method='BY'),
         under_rep_padj = p.adjust(under_represented_pvalue, method='BY')) %>%
  mutate(a = numDEInCat) %>%
  mutate(b = numGenes - a) %>% 
  mutate(c = numInCat-a) %>%
  mutate(d = allgenenum -a - b - c) %>%
  mutate(enrichScore = (a/b)/(c/d)) %>%
  mutate(log2enrich = log2(enrichScore)) %>%
  mutate(cluster=gsub('cl','',cluster)) %>%
  mutate(cluster = factor(cluster, levels = c('1','2','3','1-2','1-3','2-3','1-2-3'))) %>%
  mutate(type = factor(type,levels = c('Multidisease','Multicategory')))

signifGO = GOres2 %>% 
  filter(over_rep_padj<=0.05 | under_rep_padj<=0.05)
signifGO %>%filter(cluster=='1' ) %>% View()
signifGO_genes = left_join(signifGO,allgo)

signifGO_genes$gene_cluster = sapply(signifGO_genes$gene,function(x){
  nm = unique(sapply(strsplit(names(which(sapply(genegr,function(y)x%in%y))),'_'),function(x)x[1]))
  ifelse(length(nm)==0,'-',gsub('cl','',nm))
})

gogenemat = signifGO_genes %>% 
  filter(gene_cluster!='-') %>%
  dplyr::select(category, gene) %>%
  unique() %>%
  mutate(value = 1) %>%
  spread(category, value, fill = 0)

gogenemat = as.data.frame(gogenemat)
rownames(gogenemat) = gogenemat$gene
gogenemat$gene = NULL
gogenemat = as.matrix(gogenemat)
library(pheatmap)
jaccardsim = apply(gogenemat,2,function(x){
  apply(gogenemat,2,function(y){
    sum(x==1 & y==1) / sum(x==1 | y==1)
  })
})
newdf = data.frame()
for(ont in c('BP','MF','CC')){
  for(cl in c('1','2','3','1-2','1-3','2-3','1-2-3')) {
    print(c(ont,cl))
    gocatx = (signifGO %>% filter(ontology == ont & cluster == cl))
    gocat= unique(gocatx$category)
    if(length(gocat)==0){
      
    } else if( length(gocat)==1){
      newdf = rbind(newdf, mutate(gocatx,rep=gocat))
    } else {
      simmat = jaccardsim[gocat,gocat]
      k = min(which(sapply(1:(nrow(simmat)),function(i){
        mycl = cutree(hclust(dist(simmat)),i)
        all(sapply(1:i,function(k){
          xx = names(which(mycl==k))
          if(length(xx)==1){
            return(1)
          } else {
            xx = simmat[xx,xx]
            median(xx[upper.tri(xx)]) 
          }
        })>=0.5)
      })))
      k=ifelse(k==1,nrow(simmat),k)
      treex =hclust(dist(simmat))
      treecl = cutree(treex,k)
      reps = sapply(1:k,function(i){
        xx=names(which(treecl==i))
        if(length(xx)>1){
          xx = simmat[xx,xx]
          names(which.max(rowMeans(xx)))
        } else{
          xx
        }
      })
      repclus=setNames(lapply(1:k,function(i)names(which(treecl==i))),reps)
      newdf = rbind(newdf, reshape2::melt(repclus) %>% 
                      set_names(c('category','rep')) %>%
                      arrange(rep) %>% 
                      mutate(ontology = ont, cluster = cl) %>%
                      left_join(signifGO) %>%
                      unique() %>%
                      dplyr::select(1,3:23,2))
    }}}

d1 = newdf %>%
  filter(rep==category) %>% 
  filter(cluster=='1') %>%
  mutate(termname = ifelse(nchar(term)>40,paste(substr(term,1,37),'...',sep=''),term)) %>%
  unique()
d2 = newdf %>%
  filter(rep==category) %>% 
  filter(cluster=='2') %>%
  mutate(termname = ifelse(nchar(term)>40,paste(substr(term,1,37),'...',sep=''),term)) %>%
  unique()
d3 = newdf %>%
  filter(rep==category) %>% 
  filter(cluster=='3') %>%
  mutate(termname = ifelse(nchar(term)>40,paste(substr(term,1,37),'...',sep=''),term)) %>%
  unique()
d12 = newdf %>%
  filter(rep==category) %>% 
  filter(cluster=='1-2') %>%
  mutate(termname = ifelse(nchar(term)>40,paste(substr(term,1,37),'...',sep=''),term)) %>%
  unique()
d123 = newdf %>%
  filter(rep==category) %>% 
  filter(cluster=='1-2-3') %>%
  mutate(termname = ifelse(nchar(term)>40,paste(substr(term,1,37),'...',sep=''),term)) %>%
  unique()

p1 = ggplot(d1, aes(x = reorder(termname,log2enrich), y = log2enrich, fill=type)) +
  geom_bar(stat='identity',position=position_dodge2(preserve = "single",padding = 0)) +
  coord_flip() +
  ylab('Log2 Enrichment Score') +
  xlab(NULL) +
  scale_fill_brewer(type='qual',palette=3) +
  ggtitle('Cluster 1') +
  ylim(0,7.1) +
  theme_pubr(base_size = 6)
p2 = ggplot(d2, aes(x = reorder(termname,log2enrich), y = log2enrich, fill=type)) +
  geom_bar(stat='identity',position=position_dodge2(preserve = "single",padding = 0)) +
  coord_flip() +
  ylab('Log2 Enrichment Score') +
  xlab(NULL) +
  scale_fill_brewer(type='qual',palette=3) +
  ggtitle('Cluster 2')+
  ylim(0,7.1)+
  theme_pubr(base_size = 6)
p3 = ggplot(d3, aes(x = reorder(termname,log2enrich), y = log2enrich, fill=type)) +
  geom_bar(stat='identity',position=position_dodge2(preserve = "single",padding = 0)) +
  coord_flip() +
  ylab('Log2 Enrichment Score') +
  xlab(NULL) +
  scale_fill_brewer(type='qual',palette=3) +
  ggtitle('Cluster 3')+
  ylim(0,7.1)+
  theme_pubr(base_size = 6)
p12 = ggplot(d12, aes(x = reorder(termname,log2enrich), y = log2enrich, fill=type)) +
  geom_bar(stat='identity',position=position_dodge2(preserve = "single",padding = 0)) +
  coord_flip() +
  ylab('Log2 Enrichment Score') +
  xlab(NULL) +
  scale_fill_brewer(type='qual',palette=3) +
  ggtitle('Cluster 1 & 2')+
  ylim(0,7.1)+
  theme_pubr(base_size = 6)
p123 = ggplot(d123, aes(x = reorder(termname,log2enrich), y = log2enrich, fill=type)) +
  geom_bar(stat='identity',position=position_dodge2(preserve = "single",padding = 0)) +
  coord_flip() +
  ylab('Log2 Enrichment Score') +
  xlab(NULL) +
  scale_fill_brewer(type='qual',palette=3) +
  ggtitle('Cluster 1 & 2 & 3')+
  ylim(0,7.1)+
  theme_pubr(base_size = 6)

px1 = ggarrange(p1,p2,ncol=1,nrow=2, heights = c(4,1),align='v',common.legend = T, legend='none')
px2 = p3
px12 = ggarrange(px1,px2,ncol=2,nrow=1,common.legend = T,legend='none')
px3=ggarrange(p12,p123,ncol=1,nrow=2,align='v',heights = c(1.2,1),legend='none')
pxfin=ggarrange(px12,px3,ncol=2,nrow=1,widths = c(2,1))
ggsave('./results/functionalAnalysis/gores_final.pdf',units='cm',width=20,height = 10,useDingbats = F)
ggsave('./results/functionalAnalysis/gores_final.png',units='cm',width=20,height = 10)

newdf %>%
  write_csv('./results/functionalAnalysis/gores_signif_reps.csv')

GOres2 %>%
  left_join(newdf)%>%
  write_csv('./results/functionalAnalysis/gores_all.csv')
