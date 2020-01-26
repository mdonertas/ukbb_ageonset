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
genegr = list(cl1genes,cl2genes,cl3genes,
              cl1genes_h1cat,cl2genes_h1cat,cl3genes_h1cat)

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
names(permres) = paste(rep(1:3,2),rep(c('','_strict'),each=3),sep='')

xx = reshape2::melt(permres) %>%
  spread(Var2,value) %>%
  separate(L1,into=c('cluster','type')) %>%
  mutate(type=ifelse(is.na(type),'all',type)) %>%
  rename(Aging = Var1) %>%
  mutate(prop = log2(val/perm))
oddsplot = xx %>% 
  filter(type =='strict') %>%
  ggplot() +
  geom_bar(aes(x = cluster, y= prop, fill=cluster), stat='identity', position = 'dodge') +
  facet_grid(~Aging) +
  scale_fill_manual(values = ageonsetcolors) +
  geom_label(aes(label = num, x= cluster),y=0.2)+
  # ylim(-3.5,3.5) +
  geom_text(data =filter(filter(xx,type=='strict'),p<=0.05),aes(x=cluster,y=-0.01,label=scales::pvalue(p,add_p = T,accuracy = 0.001)), angle=90,color='black',hjust=1) +
  # geom_text(data =filter(xx,type=='strict'),aes(x=cluster,y=-0.01,label=scales::pvalue(p,add_p = T,accuracy = 0.001)), angle=90,color='black',hjust=1) +
  xlab('Age of onset cluster') + ylab('Log2 Odds Ratio') +
  guides(fill = F)

odds_all = xx %>% 
  filter(type =='all') %>%
  ggplot() +
  geom_bar(aes(x = cluster, y= prop, fill=cluster), stat='identity', position = 'dodge') +
  facet_grid(~Aging) +
  scale_fill_manual(values = ageonsetcolors) +
  geom_label(aes(label = num, x= cluster),y=0.1)+
  # ylim(-1.5,1.5) +
  geom_text(data =filter(filter(xx,type=='all'),p<=0.05),aes(x=cluster,y=-0.01,label=scales::pvalue(p,add_p = T,accuracy = 0.001)), angle=90,color='black',hjust=1) +
  # geom_text(data =filter(xx,type=='all'),aes(x=cluster,y=-0.01,label=scales::pvalue(p,add_p = T,accuracy = 0.001)), angle=90,color='black',hjust=1) +
  xlab('Age of onset cluster') + ylab('Log2 Odds Ratio') +
  guides(fill = F)

oddsx = ggarrange(odds_all+ggtitle('Genes specific to one cluster'),oddsplot+ggtitle('Genes specific to one cluster and linked to multiple categories'),labels = 'auto',nrow=2)

system('mkdir ./results/agingRelevance/')
ggsave('./results/agingRelevance/genage2.pdf',oddsx,units ='cm',width = 15,height = 20,useDingbats=F)
ggsave('./results/agingRelevance/genage2.png',oddsx,units ='cm',width = 15,height = 20)

