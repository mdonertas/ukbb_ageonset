source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)
mhcchr=6
mhcstart=28477797
mhcend=33448354

gwascat <-  read_tsv('../melike/projects/shared_data/GWASCatalog/20190607/data/raw/gwas_catalog_v1.0.2-associations_e96_r2019-05-03.tsv') %>%
  filter(`P-VALUE`<=5e-8)
gwasassocs <- lapply(unique(unlist(strsplit(gwascat$MAPPED_TRAIT,', '))),function(tr){
  rowx=grep(tr,gwascat$MAPPED_TRAIT,ignore.case = T)
  if(length(rowx)>0){
    unique(unlist(strsplit(unlist(strsplit(gwascat$MAPPED_GENE[rowx],', ')),' - ')))
  } else{
    NA
  }
})
names(gwasassocs) <- unique(unlist(strsplit(gwascat$MAPPED_TRAIT,', ')))

proxygenes <- sapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_proxyGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  setdiff(unique(x$proxy_hgnc),c('',NA))
})
names(proxygenes)=disIDs

eQTLgenes <- sapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_eQTLGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  setdiff(unique(x$eQTL_hgnc),c('',NA))
})
names(eQTLgenes)=disIDs

signifgenes <- sapply(disIDs,function(dis){
  union(proxygenes[[dis]],eQTLgenes[[dis]])
})
names(signifgenes)=disIDs

allingwascat <- unique(unlist(strsplit(unlist(strsplit(gwascat$MAPPED_GENE,', ')),' - ')))
allingwascat <- length(allingwascat)
resx <- lapply(signifgenes, function(disGene){
  t(sapply(gwasassocs,function(gcat){
    resx <- data.frame(x=length(unique(intersect(disGene,gcat))))
    resx$y <- length(unique(disGene)) - resx$x
    resx$z <- length(unique(gcat)) - resx$x
    resx$w <- allingwascat - resx$x - resx$y - resx$z
    resx
  }))
})

resx <- lapply(resx,function(resxx){
  resxxx <- as.data.frame(t(apply(resxx,1,function(myres){
    fi <- fisher.test(matrix(unlist(myres),ncol=2,byrow = T))
    c(unlist(myres),fi$p.value,fi$estimate)
  })))
  colnames(resxxx)[5:6]=c('p','odds')
  resxxx$p.adj = p.adjust(resxxx$p, method = 'fdr')
  resxxx$logodds = log2(resxxx$odds)
  resxxx$GWASCat = rownames(resxxx)
  resxxx
})

resx=reshape2::melt(resx,id.vars=colnames(resx[[1]]))%>%
  rename(disID=L1)%>%
  mutate(Disease = disCoding[disID])

gwascat <- resx%>%
  mutate(Percentage=x/(x+y))%>%
  rename(disease_and_gwascat=x,
         only_disease=y,
         only_gwascat=z,
         allothers=w)
write_tsv(gwascat,'results/genomicAnalysis/gwascat_noMHC.tsv')

myrank <- function(x){
  su=sort(unique(x))
  for (i in 1:length(su)) x[x==su[i]] = i
  return(x)
}
xx=gwascat %>%  
  filter(p.adj<=0.05 & abs(odds)>=2)%>%
  filter(disease_and_gwascat + only_disease >=5) %>%
  filter(disease_and_gwascat + only_gwascat >=5) %>%
  mutate(value=log2(odds)) %>%
  # group_by(Disease) %>%
  # top_n(5,odds) %>%
  # ungroup() %>%
  select(GWASCat,value,Disease)%>%
  spread(Disease,value,fill = 0)
# xx[xx==Inf] = 10
rownames(xx)=xx$GWASCat
xx$GWASCat=NULL
xx=as.matrix(xx)
xx = xx[which(rowMeans(xx==0)<1),which(colMeans(xx==0)<1)]
xx[xx==Inf]=unique(sort(xx,dec=T))[2]+1
pheatmap::pheatmap(xx,
                   color = c('white',colorRampPalette(brewer.pal(5,'Reds')[-1])(100)),
                   cutree_rows = 10,
                   # height =8.27, width = 11.69,
                   cutree_cols = 25,
                   # scale = 'column',
                   cellwidth = 10,cellheight = 10,
                   filename = 'results/genomicAnalysis/compare_with_GWASCat_allabove5_odds2_p005.pdf')
