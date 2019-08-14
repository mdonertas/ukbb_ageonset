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

gtexdat = readRDS('./data/processed/GTEx/rhomat.rds')
gtexp = readRDS('./data/processed/GTEx/padjmat.rds')

gtexsum = (abs(gtexdat)>=0.3 | gtexp<=0.1)
majtissues = sapply(strsplit(colnames( gtexsum),'-'),function(x)x[1])

gtex_sumsum = sapply(unique(majtissues),function(mt){
  if(sum(majtissues == mt)>1){
    xx = apply(gtexsum[,which(majtissues == mt)],1,function(x)max(x,na.rm=T))
  } else {
    xx = gtexsum[,which(majtissues == mt)]
  }
})

gtex_sumsum[is.na(gtex_sumsum)] = 0
# gtex_sumsum[(gtex_sumsum)==-Inf] = 0

martx = biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
genesx = intersect(rownames(gtex_sumsum), genedat$geneid)
xx = gtex_sumsum[genesx,]
# genesx = names(which(rowSums(abs(xx)>=0.3,na.rm=T)>=2))
genesx = names(which(rowSums(xx)>0))
xx = xx[genesx,]
onsetcl = setNames(genedat$ageonsetclusters,genedat$geneid)
numcat = setNames(genedat$numdiscat,genedat$geneid)
genannot = data.frame(cl1 = grepl('1',onsetcl[genesx])+1,cl2 = grepl('2',onsetcl[genesx])+1,cl3 = grepl('3',onsetcl[genesx])+1)
rownames(genannot) = genesx
genannot$numcat = 1+(numcat[rownames(genannot)]>1)
annotcolors = list(cl1 = setNames(c('gray90','gray25'),c(1:2)),
                   cl2 = setNames(c('gray90','gray25'),c(1:2)),
                   cl3 = setNames(c('gray90','gray25'),c(1:2)),
                   numcat = setNames(c('gray90','gray25'),c(1:2)))

idmap = biomaRt::getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),filters = 'ensembl_gene_id',values = genesx,mart = martx)
idmap = idmap %>%
  filter(hgnc_symbol != '') %>%
  full_join(data.frame(ensembl_gene_id = genesx)) %>%
  mutate(geneid = ifelse(is.na(hgnc_symbol),ensembl_gene_id,hgnc_symbol))
newgenesx = unname(setNames(idmap$geneid,idmap$ensembl_gene_id)[genesx])
rownames(xx)= newgenesx
rownames(genannot) = newgenesx
xx = xx + 1

genannot=genannot[which(genannot$numcat>1),]
genannot$numcat = NULL
rownames(genannot)[which(genannot$cl1 == 2 & genannot$cl2==1 & genannot$cl3==1)]
xx = xx[rownames(genannot),]
xx = xx[,colSums(xx==2)>0]
pheatmap::pheatmap(xx[rownames(genannot)[which(genannot$cl1 == 2 & genannot$cl2==1 & genannot$cl3==2)],],
                   cluster_rows = T, cluster_cols = T, 
                   # color = c('gray95','darkred'), 
                   legend = F, border_color = 'gray80',cellwidth = 10,cellheight = 10, 
                   # annotation_row = genannot, annotation_colors = annotcolors, annotation_legend = F,
                   filename = './results/temp/deneme2.pdf')

xx = apply(gtexdat,1,function(x){
  x[which.max(abs(x))]
})
xx= gtexdat[,1]
xx = xx[names(xx)%in%genedat$geneid]
clusterdat = onsetcl[names(xx)]
data.frame(expchange = unname(xx), gene = names(xx), clusterdat = unname(clusterdat)) %>%
  ggplot(aes(x = (clusterdat), y= expchange)) +
  geom_sina(size= 5, alpha= 0.7)

xxx = table(xx>0,clusterdat)
xxx[2,]/c(xxx[2,]+xxx[1,])

