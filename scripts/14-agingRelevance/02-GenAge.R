source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)

proxyGenes <- sapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_proxyGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  setdiff(unique(x$proxy_hgnc),c('',NA))
})
names(proxyGenes)=disIDs

eqtlGenes <- sapply(paste('./data/processed/caseControl/a',disIDs,'/signif_gwasRes_eQTLGenes.rds',sep=''),function(x){
  x=readRDS(x)
  x=filter(x, !(CHR==mhcchr & BP>= mhcstart & BP<=mhcend))
  setdiff(unique(x$eQTL_hgnc),c('',NA))
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
  mutate(name = paste('cl',ageonsetclusters,'_',c('all','multicat')[1+monecat],sep=''))

# genegrs$name = c('cl1_o_o','cl1_m_o','cl1_m_m','cl1-2_m_o','cl1-2_m_m','cl1-2-3_m_m','cl1-3_m_m','cl2_o_o','cl2_m_o','cl2_m_m','cl2-3_m_o','cl2-3_m_m','cl3_o_o','cl3_m_o','cl3_m_m')
# genegrs = genegrs 

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

saveRDS(genegrs,file='./data/processed/clustergenes_hgnc.rds')
####

# genagehuman = read_csv('../melike/projects/shared_data/GenAge/20190813/data/raw/genage_human.csv')
# genagemodel = read_tsv('../melike/projects/shared_data/GenAge/20190813/data/raw/genage_models_orthologs_export.tsv')
# drugage = read_csv('../melike/projects/shared_data/GenAge/20190813/data/raw/drugage.csv')
# 
# humangenes = unique(genagehuman$symbol)
# modelgenes = unique(genagemodel$Symbol)
# wormgenes = unique((filter(genagemodel,`Model Organism`=="Caenorhabditis elegans"))$Symbol)
# micegenes = unique((filter(genagemodel,`Model Organism`=="Mus musculus"))$Symbol)
# flygenes = unique((filter(genagemodel,`Model Organism`=="Drosophila melanogaster"))$Symbol)
# yeastgenes = unique((filter(genagemodel,`Model Organism`=="Saccharomyces cerevisiae"))$Symbol)
# 
# DA_drugs_incSyn <- unique(drugage$compound_name)
# name2CID <- function(nm) {
#   library(RCurl)
#   library(jsonlite)
#   nm <- URLencode(nm, reserved = T)
#   name2cid <- getURL(paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", 
#                            nm, "/cids/JSON", sep = ""))
#   name2cid <- fromJSON(name2cid)
#   return(name2cid$IdentifierList$CID)
# }
# convertIDs <- function(drugList) {
#   cids <- sapply(drugList, name2CID)
#   cids <- cids[!sapply(cids, is.null)]
#   cids <- reshape2::melt(cids) %>% rename(CID = value, 
#                                           drugName = L1)
#   noCID <- setdiff(drugList, cids$drugName)
#   return(list(cids, noCID))
# }
# DA_drugs_incSyn_CIDs <- convertIDs(DA_drugs_incSyn)
# DA_drugs_incSyn_CIDlist <- DA_drugs_incSyn_CIDs[[1]]
# cid2chembl <- function(intab, out_fx) {
#   pubchem2chembl <- read_tsv("../melike/projects/shared_data/UniChem/20190813/data/raw/src1src22.txt") %>% 
#     setNames(., c("ChEMBLID", "CID")) %>% mutate(CID = as.character(CID))
#   fintable <- intab %>% mutate(CID = as.character(CID)) %>% 
#     left_join(pubchem2chembl) %>% unique()
#   return(fintable)
# }
# DA_drugs_incSyn_CHEMBLlist <- cid2chembl(DA_drugs_incSyn_CIDlist)
# interactions <- read_tsv("http://www.dgidb.org/data/interactions.tsv") %>% 
#   rename(ChEMBLID = drug_chembl_id) %>% dplyr::select(gene_name, 
#                                                       ChEMBLID) %>% unique() %>% na.omit()
# drugagegenes = unique((interactions %>% filter(ChEMBLID %in% setdiff(DA_drugs_incSyn_CHEMBLlist$ChEMBLID,NA)))$gene_name)
# 
# aginggenes = list(humangenes,modelgenes,micegenes,wormgenes,flygenes,yeastgenes,drugagegenes)
# 
# saveRDS(aginggenes,'./data/processed/agingrelevance/aginggenes.rds')
aginggenes = readRDS('./data/processed/agingrelevance/aginggenes.rds')
cellage = read_csv('./data/raw/cellage.csv')
aginggenes = c(aginggenes,list(cellage$`Gene Name` ))

eqtlgenes = readRDS('./data/processed/caseControl/a1071/gwasRes_eQTLGenes.rds')
eqtlgenes = setdiff(unique(eqtlgenes$eQTL_hgnc),NA)
proxygenes = readRDS('./data/processed/caseControl/a1071/gwasRes_proxyGenes.rds')
proxygenes = setdiff(unique(proxygenes$proxy_hgnc),NA)

allgenes = unique(c(eqtlgenes,proxygenes))
genegr = setNames(genegrs$genelist, genegrs$name)

# resx = lapply(genegr,function(x){
#   xx = t(sapply(aginggenes[c(1,2,7,8)],function(y){
#     a = length(unique(x%in%y))
#     b = length(unique(x)) - a
#     c = length(unique(y)) - a
#     d = length(unique(allgenes)) - a - b - c
#     matx = matrix(c(a,b,c,d), ncol=2, byrow = T)
#     fi = fisher.test(matx)
#     c(a,b,c,d,fi$p,fi$est)
#   }))
#   rownames(xx) = c('human','model','drug','senescence')
#   colnames(xx) = c('a','b','c','d','p','odds')
#   xx
# })
# names(resx) = paste(rep(1:3,2),rep(c('','_strict'),each=3),sep='')

aginggenes[[9]] = unique(unlist(aginggenes))

permres = lapply(genegr,function(x){
  xx = t(sapply(aginggenes[c(1,2,7,8,9)],function(y){
    xx = sapply(1:10000,function(i){100*mean(sample(allgenes,length(x))%in%y)})
    c(perm=mean(xx),val=100*mean(x%in%y),p=mean(xx>=(100*mean(x%in%y))),num = sum(x%in%y),genesize = length(x), agsize = length(intersect(y,allgenes)))
  }))
  rownames(xx) = c('human','model organism','drug','senescence','all-combined')
  xx
})
# names(permres) = paste(rep(1:3,2),rep(c('','_strict'),each=3),sep='')

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
library(ggthemes)

agingdb = setNames(c('GenAge-Human','GenAge-Model','DrugAge','CellAge','All Combined'),c('human','model organism','drug','senescence','all-combined'))
xx = xx %>% 
  mutate(cluster=gsub('cl','',cluster),
         type = c('Multidisease','Multicategory')[as.numeric(as.factor(type))]) %>%
  mutate(Aging = agingdb[Aging]) %>%
  mutate(cluster = factor(cluster, levels = c('1','2','3','1-2','1-3','2-3','1-2-3'))) %>%
  mutate(type = factor(type,levels = c('Multidisease','Multicategory'))) %>%
  mutate(prop = ifelse(prop==-Inf,NA,prop))
oddsplot = ggplot(xx, aes(x=cluster, y = prop, color = Aging, size = p<=0.05)) +
  geom_hline(yintercept = 0, color = 'darkred',linetype = 'dashed', size = 0.2, alpha = 0.5) +
  geom_text(data=dplyr::select(xx,cluster,numgenes,type),aes(label=paste(numgenes,'genes')), y = -1.6, color = 'gray25', size = 6/pntnorm,hjust=0,nudge_x = 0.4) +
  geom_point(alpha=0.7,aes(shape=p<0.05)) + 
  geom_text_repel(data = filter(xx,p<=0.05),aes(label = paste('(',num,')',sep='')), size = 6/pntnorm, box.padding = 0.25,min.segment.length = 0.25) +
  scale_size_manual(values=setNames(c(0.5,2),c(F,T))) +
  scale_shape_manual(values=setNames(c(4,19),c(F,T))) +
  scale_color_gdocs() +
  facet_wrap(~type, ncol = 2, nrow=1) +
  coord_flip() +
  xlab('') +
  ylab('Log2 Enrichment Score') +
  theme_pubr(base_size = 6) +
  theme(panel.grid.major.y = element_line(color = 'gray90',linetype='dotted'),
        legend.direction = 'horizontal', legend.position = 'bottom') +
  guides(size = F, shape=F,
         color = guide_legend('Database')) 
ggsave(filename = './results/agingRelevance/genage3.pdf',oddsplot,units = 'cm',width = 16,height = 6,useDingbats=F)
ggsave(filename = './results/agingRelevance/genage3.png',oddsplot,units = 'cm',width = 16,height = 6)

lapply(genegr,function(x){
  xx = lapply(aginggenes[c(1,2,7,8,9)],function(y){
    intersect(x,y)
  })
  names(xx) = c('human','model organism','drug','senescence','all-combined')
  xx
}) %>%
  reshape2::melt() %>%
  write_csv('./results/agingRelevance/genageoverlaps.csv')
beepr::beep()
