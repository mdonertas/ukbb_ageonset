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

cl1genes= (genedat %>%
             filter(ageonsetclusters == '1') %>%
             filter(numdiscat > 0))$geneid

cl2genes = (genedat %>%
              filter(ageonsetclusters == '2') %>%
              filter(numdiscat > 0))$geneid

cl3genes= (genedat %>%
             filter(ageonsetclusters == '3') %>%
             filter(numdiscat > 0))$geneid

cl1genes_h1cat = (genedat %>%
                    filter(ageonsetclusters == '1') %>%
                    filter(numdiscat > 1))$geneid

cl2genes_h1cat = (genedat %>%
                    filter(ageonsetclusters == '2') %>%
                    filter(numdiscat > 1))$geneid

cl3genes_h1cat = (genedat %>%
                    filter(ageonsetclusters == '3') %>%
                    filter(numdiscat > 1))$geneid

####

senesfiles = list.files('./data/gkz555_supplemental_files/processed/')

senesdata = lapply(senesfiles,function(x){
  read_csv(paste('./data/gkz555_supplemental_files/processed/',x,sep='')) %>%
    select(1,2,4) %>%
    set_names(c('geneid','logFC','FDR'))
})
names(senesdata) = gsub('.csv','',senesfiles)

senesdata = lapply(senesdata,function(x)unique(filter(x,FDR<=0.1)))

genesx = unique(unname(unlist(sapply(senesdata,function(x)x$geneid))))
x=senesdata[[1]]
senesdata = sapply(senesdata,function(x){
  setNames(x$logFC,x$geneid)[genesx]
})
rowmeds = function(x){
  apply(x,1,function(x)median(x,na.rm=T))
}
martx = biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
genesx = intersect(rownames(senesdata), c(cl1genes_h1cat))
genesx = names(sort(abs(rowmeds(senesdata[genesx,])),dec=T))
# genesx = sample(rownames(gtex_sumsum),length(genesx))
idmap = biomaRt::getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),filters = 'ensembl_gene_id',values = genesx,mart = martx)
idmap = idmap %>%
  filter(hgnc_symbol != '') %>%
  full_join(data.frame(ensembl_gene_id = genesx)) %>%
  mutate(geneid = ifelse(is.na(hgnc_symbol),ensembl_gene_id,hgnc_symbol))
newgenesx = unname(setNames(idmap$geneid,idmap$ensembl_gene_id)[genesx])
geneannot = data.frame(cluster = as.character(ifelse(genesx%in%cl1genes_h1cat,'1', ifelse(genesx%in%cl2genes_h1cat,'2','3'))))
newgenesx[which(newgenesx == 'PINX1')[2]]='PINX1_2'
rownames(geneannot) = newgenesx
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
geneannot = cbind(geneannot,xx[newgenesx,])
geneannot$`<NA>`=NULL
annotcolors = lapply(setdiff(colnames(xx),NA),function(i){setNames(c('gray95','gray25'),c('0','1'))})
names(annotcolors) = setdiff(colnames(xx),NA)
annotcolors = c(annotcolors,cluster=list(setNames(ageonsetcolors[1:3],1:3)))
xx= senesdata[genesx, ]
rownames(xx) = newgenesx
keepx = which(rowSums(abs(xx)>=1,na.rm=T)!=0)
xx = xx[names(keepx),]
geneannot = geneannot[names(keepx),]
keepx = c('cluster',names(which(colSums(geneannot[,-1])!=0)))
geneannot = geneannot[,keepx]

geneannot = geneannot[,c('cluster',names(sort((readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$cluster)[colnames(geneannot)[-1]])))]

pheatmap::pheatmap(xx,breaks = seq(-10,10,length.out = 16), 
                   color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),
                                              'white',brewer.pal(8,'Reds')))(15), 
                   legend = T, border_color = 'gray80',
                   cellwidth = 10,cellheight = 10, 
                   annotation_row = geneannot,
                   annotation_colors = annotcolors, 
                   annotation_legend = F, 
                   # cluster_rows = T, cutree_rows = 6,
                   filename = './results/temp/cl1.pdf')

sc = senesdata[,grep('replicative',colnames(senesdata),ignore.case = T)]
sc = sc[complete.cases(sc),]
# seneschange = rowMeans( sc)

genesx = intersect(rownames(sc),genedat$geneid)
senesdat = sc[genesx,]
genesx = names(which(rowSums(abs(senesdat)>=1,na.rm=T)==2))
senesdat = senesdat[genesx,]
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
rownames(senesdat)= newgenesx
rownames(genannot) = newgenesx

pheatmap::pheatmap(t(senesdat),
                   breaks = seq(-8,8,by=1),
                   color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),
                                              'white',brewer.pal(8,'Reds')))(16),
                   legend = T, border_color = 'gray80',
                   cellwidth = 10,cellheight = 10,
                   annotation_col = genannot,
                   annotation_colors = annotcolors,
                   annotation_legend = F,
                  cutree_cols = 8,
                  # kmeans_k = 10,
                   filename = './results/temp/deneme.pdf')


xx =reshape2::melt(seneschange) %>%
  mutate(geneid = names(seneschange)) %>%
  rename(meanlogFC = value) %>%
  left_join(genedat) %>% 
  mutate(ageonsetclusters = ifelse(is.na(ageonsetclusters),'other',ageonsetclusters)) 

xx %>%
  mutate(ageonsetclusters = factor(ageonsetclusters,levels = c('1','2','3','1-2','1-3','2-3','1-2-3','other'))) %>%
  ggplot(aes(x = ageonsetclusters, y= abs(meanlogFC))) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))

