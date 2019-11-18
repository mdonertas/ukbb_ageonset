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

allgenesx = unique(c(cl1genes, cl2genes, cl3genes, cl1genes_h1cat, cl2genes_h1cat, cl3genes_h1cat))
martx = biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
genlen = biomaRt::getBM(attributes = c('ensembl_gene_id','start_position','end_position'), filters = 'ensembl_gene_id', values = allgenesx, mart = martx) %>%
  mutate(genlen = abs(end_position-start_position)) %>%
  rename(geneid = ensembl_gene_id) %>%
  select(geneid, genlen)
genlen$clgenes = 0
genlen$clgenes[genlen$geneid%in%cl1genes] = 1
genlen$clgenes[genlen$geneid%in%cl2genes] = 2
genlen$clgenes[genlen$geneid%in%cl3genes] = 3
genlen$clgenes_h1cat = 0
genlen$clgenes_h1cat[genlen$geneid%in%cl1genes_h1cat] = 1
genlen$clgenes_h1cat[genlen$geneid%in%cl2genes_h1cat] = 2
genlen$clgenes_h1cat[genlen$geneid%in%cl3genes_h1cat] = 3

genlen %>%
  ggplot(aes(x = as.factor(clgenes), y = genlen, fill = as.factor(clgenes))) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = 'white', outlier.shape = NA) +
  scale_y_log10() +
  stat_compare_means(label.x = 0.5, label.y = 6.5,hjust = 0) +
  scale_fill_manual(values = ageonsetcolors) +
  xlab ('Age of Onset Clusters') + ylab('Gene Length') +
  guides(fill = F)

ggsave('./results/functionalAnalysis/genelength.pdf', units = 'cm', width = 8, height = 8, useDingbats = F)
ggsave('./results/functionalAnalysis/genelength.png', units = 'cm', width = 8, height = 8)

genlen %>%
  ggplot(aes(x = as.factor(clgenes_h1cat), y = genlen)) +
  geom_boxplot() +
  scale_y_log10() +
  stat_compare_means()

####

source('../shared/functions/functions.R')

gnls=unique(signifGenes$geneid)
genex = setNames(rep(0,length(gnls)),gnls)
genex[names(genex) %in% cl1genes] = 1
genex[names(genex) %in% cl2genes] = 2
genex[names(genex) %in% cl3genes] = 3
bpres_cl1 = go_enrich.test(genelist = genex, selection = function(x)x==1) %>% 
  mutate(cluster = '1', ontology = 'BP')
bpres_cl2 = go_enrich.test(genelist = genex, selection = function(x)x==2) %>% 
  mutate(cluster = '2', ontology = 'BP')
bpres_cl3 = go_enrich.test(genelist = genex, selection = function(x)x==3) %>% 
  mutate(cluster = '3', ontology = 'BP')
mfres_cl1 = go_enrich.test(genelist = genex, selection = function(x)x==1, ontologyx = 'MF') %>% 
  mutate(cluster = '1', ontology = 'MF')
mfres_cl2 = go_enrich.test(genelist = genex, selection = function(x)x==2, ontologyx = 'MF') %>% 
  mutate(cluster = '2', ontology = 'MF')
mfres_cl3 = go_enrich.test(genelist = genex, selection = function(x)x==3, ontologyx = 'MF') %>% 
  mutate(cluster = '3', ontology = 'MF')
ccres_cl1 = go_enrich.test(genelist = genex, selection = function(x)x==1, ontologyx = 'CC') %>% 
  mutate(cluster = '1', ontology = 'CC')
ccres_cl2 = go_enrich.test(genelist = genex, selection = function(x)x==2, ontologyx = 'CC') %>% 
  mutate(cluster = '2', ontology = 'CC')
ccres_cl3 = go_enrich.test(genelist = genex, selection = function(x)x==3, ontologyx = 'CC') %>% 
  mutate(cluster = '3', ontology = 'CC')

gores = rbind(bpres_cl1,bpres_cl2,bpres_cl3,
              mfres_cl1,mfres_cl2,mfres_cl3,
              ccres_cl1,ccres_cl2,ccres_cl3) %>%
  filter(Annotated<=500 & Annotated>=5)
gores = gores %>% 
  mutate(p.adjusted = p.adjust(classicFisher, method='fdr'))
gores %>%
  dplyr::select(-genelist) %>%
  arrange(p.adjusted) %>%
  write_tsv('./results/functionalAnalysis/clusterXgenes.tsv')
p.val = 0.05

enrichment_table=dplyr::filter(gores, cluster=='1')
enrichment_table=enrichment_table[enrichment_table$p.adjusted<p.val,]
enrichment_table$Gr = c(rep('Detection of stimulus',6),
                        rep('Other',1),
                        rep('GPCR',1),
                        rep('Lipoprotein-related functions',3),
                        rep('Other',1),
                        rep('Lipoprotein-related functions',1),
                        rep('Morphogenesis Differentiation',1),
                        rep('Detection of stimulus',1),
                        rep('Other',1),
                        rep('Morphogenesis Differentiation',1),
                        rep('Other',1),
                        rep('Morphogenesis Differentiation',1),
                        rep('Detection of stimulus',1),
                        rep('GPCR',1),
                        rep('Morphogenesis Differentiation',1),
                        rep('Blood circulation',1),
                        rep('Other',2),
                        rep('Blood circulation',1),
                        rep('GPCR',1),
                        rep('Detection of stimulus',1),
                        rep('Other',1))
mymat=melt(setNames(enrichment_table$genelist,enrichment_table$GO.ID))
mygraph=graph_from_data_frame(mymat)
gr = setNames(enrichment_table$Gr,enrichment_table$GO.ID)[V(mygraph)$name]
verttype=gr
verttype[grep('ENSG',names(V(mygraph)))]='Gene'
verttype = setNames(verttype,names(V(mygraph)))
verttype=as.factor(verttype)
vertsize=setNames(rep(0.5,length(V(mygraph))),names(V(mygraph)))
vertsize[enrichment_table$GO.ID]=-log10(enrichment_table$p.adjusted)
vertsize=ceiling(vertsize)
labx=setNames(rep('',length(V(mygraph))),names(V(mygraph)))
labx[enrichment_table$GO.ID]=substr(enrichment_table$Term,1,41)
verttype = factor(verttype,levels = levels(verttype)[c(3,2,5,6,4,1,7)])
V(mygraph)$verttype = verttype
# vertsize[verttype=='Gene'] = 0.5
V(mygraph)$vertsize = vertsize
# gr = gr[complete.cases(gr)]
library(ggforce)
labx[verttype!='Other']=''
# xx = sapply(unique(setdiff(verttype,c('Gene','Other'))),function(x){
#   xx = verttype[!verttype%in%c('Gene','Other')]
#   sample(names(xx)[xx==x],1)
# })
# labx[xx]=names(xx)
labx[labx=='']=NA
cl1_go=ggnet2(mygraph, size = 0, edge.alpha = 0.8, edge.color = 'gray85') +
  geom_point(aes(fill = verttype), 
             size = vertsize,
             shape = c(23,21)[1+(verttype=='Gene')],
             alpha = c(0.8,0.8)[1+(verttype=='Gene')]) +
  geom_text(label=labx,size= 6/pntnorm,hjust=0,nudge_x = 0.01) +
  scale_fill_brewer(type = 'qual',palette = 6)+
  guides(fill=guide_legend('',override.aes = list(shape=21,size=2)))+
  theme_void() +
  theme(legend.position = 'bottom')
# system('mkdir ./results/functionalAnalysis')
ggsave('./results/functionalAnalysis/cl1_go.pdf',cl1_go, units = 'cm', width = 16.7,height = 10, useDingbats = F)
ggsave('./results/functionalAnalysis/cl1_go.png',cl1_go, units = 'cm', width = 16.7,height = 10)

enrichment_table=dplyr::filter(gores, cluster=='2')
enrichment_table=enrichment_table[enrichment_table$p.adjusted<p.val,]
# enrichment_table$Gr = c(rep('Detection of stimulus',6),
#                         rep('Other',1),
#                         rep('GPCR',1),
#                         rep('Lipoprotein-related functions',3),
#                         rep('Other',1),
#                         rep('Lipoprotein-related functions',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Detection of stimulus',1),
#                         rep('Other',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Other',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Detection of stimulus',1),
#                         rep('GPCR',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Blood circulation',1),
#                         rep('Other',2),
#                         rep('Blood circulation',1),
#                         rep('GPCR',1),
#                         rep('Detection of stimulus',1),
#                         rep('Other',1))
mymat=melt(setNames(enrichment_table$genelist,enrichment_table$GO.ID))
mygraph=graph_from_data_frame(mymat)
# gr = setNames(enrichment_table$Gr,enrichment_table$GO.ID)[V(mygraph)$name]
verttype=rep('GO.ID',length(V(mygraph)))
verttype[grep('ENSG',names(V(mygraph)))]='Gene'
verttype = setNames(verttype,names(V(mygraph)))
verttype=as.factor(verttype)
vertsize=setNames(rep(0.5,length(V(mygraph))),names(V(mygraph)))
vertsize[enrichment_table$GO.ID]=-log10(enrichment_table$p.adjusted)
vertsize=ceiling(vertsize)
labx=setNames(rep('',length(V(mygraph))),names(V(mygraph)))
labx[enrichment_table$GO.ID]=substr(enrichment_table$Term,1,41)
# verttype = factor(verttype,levels = levels(verttype)[c(3,2,5,6,4,1,7)])
V(mygraph)$verttype = verttype
# vertsize[verttype=='Gene'] = 0.5
V(mygraph)$vertsize = vertsize
# gr = gr[complete.cases(gr)]
library(ggforce)
# labx[verttype!='Other']=''
# xx = sapply(unique(setdiff(verttype,c('Gene','Other'))),function(x){
#   xx = verttype[!verttype%in%c('Gene','Other')]
#   sample(names(xx)[xx==x],1)
# })
# labx[xx]=names(xx)
# labx[labx=='']=NA
cl2_go = ggnet2(mygraph, size = 0, edge.alpha = 0.8, edge.color = 'gray85') +
  geom_point(aes(fill = verttype), 
             size = vertsize,
             shape = c(23,21)[1+(verttype=='Gene')],
             alpha = c(0.8,0.8)[1+(verttype=='Gene')]) +
  geom_text_repel(label=labx,size= 6/pntnorm,hjust=0,nudge_x = 0.01) +
  scale_fill_brewer(type = 'qual',palette = 6)+
  guides(fill=guide_legend('',override.aes = list(shape=21,size=2)))+
  theme_void() +
  theme(legend.position = 'bottom')
# system('mkdir ./results/functionalAnalysis')
ggsave('./results/functionalAnalysis/cl2_go.pdf',cl2_go, units = 'cm', width = 8,height = 8, useDingbats = F)
ggsave('./results/functionalAnalysis/cl2_go.png',cl2_go, units = 'cm', width = 8,height = 8)

enrichment_table=dplyr::filter(gores, cluster=='3')
enrichment_table=enrichment_table[enrichment_table$p.adjusted<p.val,]
# enrichment_table$Gr = c(rep('Detection of stimulus',6),
#                         rep('Other',1),
#                         rep('GPCR',1),
#                         rep('Lipoprotein-related functions',3),
#                         rep('Other',1),
#                         rep('Lipoprotein-related functions',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Detection of stimulus',1),
#                         rep('Other',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Other',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Detection of stimulus',1),
#                         rep('GPCR',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Blood circulation',1),
#                         rep('Other',2),
#                         rep('Blood circulation',1),
#                         rep('GPCR',1),
#                         rep('Detection of stimulus',1),
#                         rep('Other',1))
mymat=melt(setNames(enrichment_table$genelist,enrichment_table$GO.ID))
mygraph=graph_from_data_frame(mymat)
# gr = setNames(enrichment_table$Gr,enrichment_table$GO.ID)[V(mygraph)$name]
verttype=rep('GO.ID',length(V(mygraph)))
verttype[grep('ENSG',names(V(mygraph)))]='Gene'
verttype = setNames(verttype,names(V(mygraph)))
verttype=as.factor(verttype)
vertsize=setNames(rep(0.5,length(V(mygraph))),names(V(mygraph)))
vertsize[enrichment_table$GO.ID]=-log10(enrichment_table$p.adjusted)
vertsize=ceiling(vertsize)
labx=setNames(rep('',length(V(mygraph))),names(V(mygraph)))
labx[enrichment_table$GO.ID]=substr(enrichment_table$Term,1,41)
# verttype = factor(verttype,levels = levels(verttype)[c(3,2,5,6,4,1,7)])
V(mygraph)$verttype = verttype
# vertsize[verttype=='Gene'] = 0.5
V(mygraph)$vertsize = vertsize
# gr = gr[complete.cases(gr)]
library(ggforce)
# labx[verttype!='Other']=''
# xx = sapply(unique(setdiff(verttype,c('Gene','Other'))),function(x){
#   xx = verttype[!verttype%in%c('Gene','Other')]
#   sample(names(xx)[xx==x],1)
# })
# labx[xx]=names(xx)
# labx[labx=='']=NA
cl3_go = ggnet2(mygraph, size = 0, edge.alpha = 0.8, edge.color = 'gray85') +
  geom_point(aes(fill = verttype), 
             size = vertsize,
             shape = c(23,21)[1+(verttype=='Gene')],
             alpha = c(0.8,0.8)[1+(verttype=='Gene')]) +
  # geom_text(label=labx,size= 6/pntnorm,hjust=0,nudge_x = 0.01) +
  scale_fill_brewer(type = 'qual',palette = 6)+
  guides(fill=guide_legend('',override.aes = list(shape=21,size=2)))+
  theme_void() +
  theme(legend.position = 'bottom')
# system('mkdir ./results/functionalAnalysis')
ggsave('./results/functionalAnalysis/cl3_go.pdf',cl3_go, units = 'cm', width = 8,height = 8, useDingbats = F)
ggsave('./results/functionalAnalysis/cl3_go.png',cl3_go, units = 'cm', width = 8,height = 8)


gnls=unique(signifGenes$geneid)
genex = setNames(rep(0,length(gnls)),gnls)
genex[names(genex) %in% cl1genes_h1cat] = 1
genex[names(genex) %in% cl2genes_h1cat] = 2
genex[names(genex) %in% cl3genes_h1cat] = 3
bpres_cl1_h1 = go_enrich.test(genelist = genex, selection = function(x)x==1) %>% 
  mutate(cluster = '1', ontology = 'BP')
bpres_cl2_h1 = go_enrich.test(genelist = genex, selection = function(x)x==2) %>% 
  mutate(cluster = '2', ontology = 'BP')
bpres_cl3_h1 = go_enrich.test(genelist = genex, selection = function(x)x==3) %>% 
  mutate(cluster = '3', ontology = 'BP')
mfres_cl1_h1 = go_enrich.test(genelist = genex, selection = function(x)x==1, ontologyx = 'MF') %>% 
  mutate(cluster = '1', ontology = 'MF')
mfres_cl2_h1 = go_enrich.test(genelist = genex, selection = function(x)x==2, ontologyx = 'MF') %>% 
  mutate(cluster = '2', ontology = 'MF')
mfres_cl3_h1 = go_enrich.test(genelist = genex, selection = function(x)x==3, ontologyx = 'MF') %>% 
  mutate(cluster = '3', ontology = 'MF')
ccres_cl1_h1 = go_enrich.test(genelist = genex, selection = function(x)x==1, ontologyx = 'CC') %>% 
  mutate(cluster = '1', ontology = 'CC')
ccres_cl2_h1 = go_enrich.test(genelist = genex, selection = function(x)x==2, ontologyx = 'CC') %>% 
  mutate(cluster = '2', ontology = 'CC')
ccres_cl3_h1 = go_enrich.test(genelist = genex, selection = function(x)x==3, ontologyx = 'CC') %>% 
  mutate(cluster = '3', ontology = 'CC')

gores = rbind(bpres_cl1_h1,bpres_cl2_h1,bpres_cl3_h1,
              mfres_cl1_h1,mfres_cl2_h1,mfres_cl3_h1,
              ccres_cl1_h1,ccres_cl2_h1,ccres_cl3_h1) %>%
  filter(Annotated<=500 & Annotated>=5)
gores = gores %>% 
  mutate(p.adjusted = p.adjust(classicFisher, method='fdr'))
gores %>%
  dplyr::select(-genelist) %>%
  arrange(p.adjusted) %>%
  write_tsv('./results/functionalAnalysis/strict_clusterXgenes.tsv')

p.val = 0.05

enrichment_table=dplyr::filter(gores, cluster=='1')
enrichment_table=enrichment_table[enrichment_table$p.adjusted<p.val,]
enrichment_table$Gr = c(rep('Protein Binding',1),
                        rep('Cell Cycle',3),
                        rep('Other',2),
                        rep('Cell Cycle',5),
                        rep('Protein Binding',1),
                        rep('Cell Cycle',1),
                        rep('Protein Binding',1),
                        rep('Cell Cycle',1),
                        rep('Other',2),
                        rep('Cell Cycle',2))
mymat=melt(setNames(enrichment_table$genelist,enrichment_table$GO.ID))
mygraph=graph_from_data_frame(mymat)
gr = setNames(enrichment_table$Gr,enrichment_table$GO.ID)[V(mygraph)$name]
verttype=gr
verttype[grep('ENSG',names(V(mygraph)))]='Gene'
verttype = setNames(verttype,names(V(mygraph)))
verttype=as.factor(verttype)
vertsize=setNames(rep(0.5,length(V(mygraph))),names(V(mygraph)))
vertsize[enrichment_table$GO.ID]=-log10(enrichment_table$p.adjusted)
vertsize=ceiling(vertsize)
labx=setNames(rep('',length(V(mygraph))),names(V(mygraph)))
labx[enrichment_table$GO.ID]=substr(enrichment_table$Term,1,41)
verttype = factor(verttype,levels = levels(verttype)[c(2,1,4,3)])
V(mygraph)$verttype = verttype
# vertsize[verttype=='Gene'] = 0.5
V(mygraph)$vertsize = vertsize
# gr = gr[complete.cases(gr)]
library(ggforce)
labx[verttype!='Other']=''
# xx = sapply(unique(setdiff(verttype,c('Gene','Other'))),function(x){
#   xx = verttype[!verttype%in%c('Gene','Other')]
#   sample(names(xx)[xx==x],1)
# })
# labx[xx]=names(xx)
labx[labx=='']=NA
cl1_go=ggnet2(mygraph, size = 0, edge.alpha = 0.8, edge.color = 'gray85') +
  geom_point(aes(fill = verttype), 
             size = vertsize*2,
             shape = c(23,21)[1+(verttype=='Gene')],
             alpha = c(0.8,0.8)[1+(verttype=='Gene')]) +
  geom_text_repel(label=labx,size= 6/pntnorm,hjust=0,nudge_x = 0.01) +
  scale_fill_brewer(type = 'qual',palette = 6)+
  guides(fill=guide_legend('',override.aes = list(shape=21,size=2)))+
  theme_void() +
  theme(legend.position = 'bottom')
# system('mkdir ./results/functionalAnalysis')
ggsave('./results/functionalAnalysis/cl1_h1_go.pdf',cl1_go, units = 'cm', width = 8,height = 8, useDingbats = F)
ggsave('./results/functionalAnalysis/cl1_h1_go.png',cl1_go, units = 'cm', width = 8,height = 8)

enrichment_table=dplyr::filter(gores, cluster=='2')
enrichment_table=enrichment_table[enrichment_table$p.adjusted<p.val,]
# enrichment_table$Gr = c(rep('Detection of stimulus',6),
#                         rep('Other',1),
#                         rep('GPCR',1),
#                         rep('Lipoprotein-related functions',3),
#                         rep('Other',1),
#                         rep('Lipoprotein-related functions',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Detection of stimulus',1),
#                         rep('Other',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Other',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Detection of stimulus',1),
#                         rep('GPCR',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Blood circulation',1),
#                         rep('Other',2),
#                         rep('Blood circulation',1),
#                         rep('GPCR',1),
#                         rep('Detection of stimulus',1),
#                         rep('Other',1))
mymat=melt(setNames(enrichment_table$genelist,enrichment_table$GO.ID))
mygraph=graph_from_data_frame(mymat)
# gr = setNames(enrichment_table$Gr,enrichment_table$GO.ID)[V(mygraph)$name]
verttype=rep('GO.ID',length(V(mygraph)))
verttype[grep('ENSG',names(V(mygraph)))]='Gene'
verttype = setNames(verttype,names(V(mygraph)))
verttype=as.factor(verttype)
vertsize=setNames(rep(0.5,length(V(mygraph))),names(V(mygraph)))
vertsize[enrichment_table$GO.ID]=-log10(enrichment_table$p.adjusted)
vertsize=ceiling(vertsize)
labx=setNames(rep('',length(V(mygraph))),names(V(mygraph)))
labx[enrichment_table$GO.ID]=substr(enrichment_table$Term,1,41)
# verttype = factor(verttype,levels = levels(verttype)[c(3,2,5,6,4,1,7)])
V(mygraph)$verttype = verttype
# vertsize[verttype=='Gene'] = 0.5
V(mygraph)$vertsize = vertsize
# gr = gr[complete.cases(gr)]
library(ggforce)
# labx[verttype!='Other']=''
# xx = sapply(unique(setdiff(verttype,c('Gene','Other'))),function(x){
#   xx = verttype[!verttype%in%c('Gene','Other')]
#   sample(names(xx)[xx==x],1)
# })
# labx[xx]=names(xx)
# labx[labx=='']=NA
cl2_go = ggnet2(mygraph, size = 0, edge.alpha = 0.8, edge.color = 'gray85') +
  geom_point(aes(fill = verttype), 
             size = vertsize,
             shape = c(23,21)[1+(verttype=='Gene')],
             alpha = c(0.8,0.8)[1+(verttype=='Gene')]) +
  geom_text_repel(label=labx,size= 6/pntnorm,hjust=0,nudge_x = 0.01) +
  scale_fill_brewer(type = 'qual',palette = 6)+
  guides(fill=guide_legend('',override.aes = list(shape=21,size=2)))+
  theme_void() +
  theme(legend.position = 'bottom')
# system('mkdir ./results/functionalAnalysis')
ggsave('./results/functionalAnalysis/cl2_h1_go.pdf',cl2_go, units = 'cm', width = 8,height = 8, useDingbats = F)
ggsave('./results/functionalAnalysis/cl2_h1_go.png',cl2_go, units = 'cm', width = 8,height = 8)

enrichment_table=dplyr::filter(gores, cluster=='3')
enrichment_table=enrichment_table[enrichment_table$p.adjusted<p.val,]
# enrichment_table$Gr = c(rep('Detection of stimulus',6),
#                         rep('Other',1),
#                         rep('GPCR',1),
#                         rep('Lipoprotein-related functions',3),
#                         rep('Other',1),
#                         rep('Lipoprotein-related functions',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Detection of stimulus',1),
#                         rep('Other',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Other',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Detection of stimulus',1),
#                         rep('GPCR',1),
#                         rep('Morphogenesis Differentiation',1),
#                         rep('Blood circulation',1),
#                         rep('Other',2),
#                         rep('Blood circulation',1),
#                         rep('GPCR',1),
#                         rep('Detection of stimulus',1),
#                         rep('Other',1))
mymat=melt(setNames(enrichment_table$genelist,enrichment_table$GO.ID))
mygraph=graph_from_data_frame(mymat)
# gr = setNames(enrichment_table$Gr,enrichment_table$GO.ID)[V(mygraph)$name]
verttype=rep('GO.ID',length(V(mygraph)))
verttype[grep('ENSG',names(V(mygraph)))]='Gene'
verttype = setNames(verttype,names(V(mygraph)))
verttype=as.factor(verttype)
vertsize=setNames(rep(0.5,length(V(mygraph))),names(V(mygraph)))
vertsize[enrichment_table$GO.ID]=-log10(enrichment_table$p.adjusted)
vertsize=ceiling(vertsize)
labx=setNames(rep('',length(V(mygraph))),names(V(mygraph)))
labx[enrichment_table$GO.ID]=substr(enrichment_table$Term,1,41)
# verttype = factor(verttype,levels = levels(verttype)[c(3,2,5,6,4,1,7)])
V(mygraph)$verttype = verttype
# vertsize[verttype=='Gene'] = 0.5
V(mygraph)$vertsize = vertsize
# gr = gr[complete.cases(gr)]
library(ggforce)
# labx[verttype!='Other']=''
# xx = sapply(unique(setdiff(verttype,c('Gene','Other'))),function(x){
#   xx = verttype[!verttype%in%c('Gene','Other')]
#   sample(names(xx)[xx==x],1)
# })
# labx[xx]=names(xx)
# labx[labx=='']=NA
cl3_go = ggnet2(mygraph, size = 0, edge.alpha = 0.8, edge.color = 'gray85') +
  geom_point(aes(fill = verttype), 
             size = vertsize,
             shape = c(23,21)[1+(verttype=='Gene')],
             alpha = c(0.8,0.8)[1+(verttype=='Gene')]) +
  # geom_text(label=labx,size= 6/pntnorm,hjust=0,nudge_x = 0.01) +
  scale_fill_brewer(type = 'qual',palette = 6)+
  guides(fill=guide_legend('',override.aes = list(shape=21,size=2)))+
  theme_void() +
  theme(legend.position = 'bottom')
# system('mkdir ./results/functionalAnalysis')
ggsave('./results/functionalAnalysis/cl3_h1_go.pdf',cl3_go, units = 'cm', width = 8,height = 8, useDingbats = F)
ggsave('./results/functionalAnalysis/cl3_h1_go.png',cl3_go, units = 'cm', width = 8,height = 8)
