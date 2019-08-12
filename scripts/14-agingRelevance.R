source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)

proxyGenes <- sapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_proxyGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  setdiff(unique(x$proxy_ensembl),c('',NA))
})
names(proxyGenes)=disIDs

eqtlGenes <- sapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_eQTLGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  setdiff(unique(x$eQTL_ensembl),c('',NA))
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

length(unique(signifGenes$geneid))
# 3359

ageonsetsum = signifGenes %>%
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

cl1genes_h1cat = (genedat %>%
  filter(ageonsetclusters == '1') %>%
  filter(numdiscat > 1))$geneid

cl2genes_h1cat = (genedat %>%
                    filter(ageonsetclusters == '2') %>%
                    filter(numdiscat > 1))$geneid

cl3genes_h1cat = (genedat %>%
                    filter(ageonsetclusters == '3') %>%
                    filter(numdiscat > 1))$geneid


cl1genes= (genedat %>%
                    filter(ageonsetclusters == '1'))$geneid

gtexdat = readRDS('./data/processed/GTEx/rhomat.rds')
gtexp = readRDS('./data/processed/GTEx/padjmat.rds')

gtexsum = (abs(gtexdat)>=0.3 | gtexp<=0.1)
majtissues = sapply(strsplit(colnames( gtexsum),'-'),function(x)x[1])

gtex_sumsum = sapply(unique(majtissues),function(mt){
  if(sum(majtissues == mt)>1){
    xx = rowMeans(gtexsum[,which(majtissues == mt)],na.rm=T)>=0.2 
  } else {
    xx = gtexsum[,which(majtissues == mt)]
  }
})

gtex_sumsum[is.na(gtex_sumsum)] = F

martx = biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
genesx = intersect(rownames(gtex_sumsum), cl1genes_h1cat)
# genesx = sample(rownames(gtex_sumsum),length(genesx))
genesx = names(sort(rowSums(gtex_sumsum[genesx,]),dec=T))
gtex_sumsum = gtex_sumsum[,names(sort(colSums(gtex_sumsum[genesx,]),dec=T))]
idmap = biomaRt::getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),filters = 'ensembl_gene_id',values = genesx,mart = martx)
idmap = idmap %>%
  filter(hgnc_symbol != '') %>%
  full_join(data.frame(ensembl_gene_id = genesx)) %>%
  mutate(geneid = ifelse(is.na(hgnc_symbol),ensembl_gene_id,hgnc_symbol))
newgenesx = unname(setNames(idmap$geneid,idmap$ensembl_gene_id)[genesx])
xx = (genedat %>% filter(geneid %in% genesx) %>% select(geneid,all_diseases))

xx = strsplit(unname(setNames(xx$all_diseases,xx$geneid)[genesx]),', ')
names(xx) = newgenesx
xx = reshape2::melt(xx) %>%
  set_names(c('disease','L1'))
xx = xx %>%
  mutate(val = 1) %>%
  spread(disease,val,fill=0) 
rownames(xx) = xx$L1
xx$L1= NULL
geneannot = xx
geneannot$`<NA>`=NULL
xx= 1+gtex_sumsum[genesx,]
rownames(xx) = newgenesx
annotcolors = lapply(colnames(geneannot),function(i){setNames(c('gray95','gray25'),c('0','1'))})
names(annotcolors) = colnames(geneannot)
xx = xx[,colMeans(xx==1,na.rm=T)<1]
pheatmap::pheatmap(xx, cluster_rows = T, cluster_cols = T, color = c('gray95','darkred'), legend = F, border_color = 'gray80',cellwidth = 10,cellheight = 10, annotation_row = geneannot, annotation_colors = annotcolors, annotation_legend = F, filename = './results/temp/deneme.pdf')

xx = apply(gtexdat,1,function(x){
  x[which.max(abs(x))]
})
xx = xx[names(xx)%in%c(cl1genes_h1cat,cl2genes_h1cat,cl3genes_h1cat)]
clusterdat = setNames(rep(0,length(xx)),names(xx))
clusterdat[names(clusterdat)%in%cl1genes_h1cat] = 1
clusterdat[names(clusterdat)%in%cl2genes_h1cat] = 2
clusterdat[names(clusterdat)%in%cl3genes_h1cat] = 3
data.frame(expchange = unname(xx), gene = names(xx), clusterdat = unname(clusterdat)) %>%
  ggplot(aes(x = as.factor(clusterdat), y= expchange)) +
  geom_sina(size= 5, alpha= 0.7)