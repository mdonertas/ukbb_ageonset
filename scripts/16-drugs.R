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

interactions <- read_tsv("http://www.dgidb.org/data/interactions.tsv") %>% 
  rename(ChEMBLID = drug_chembl_id) %>% dplyr::select(gene_name, 
                                                      ChEMBLID) %>% unique() %>% na.omit()

drug_target_list = lapply(unique(interactions$ChEMBLID),function(drug){
  unique(filter(interactions,ChEMBLID==drug)$gene_name)
})
names(drug_target_list)=unique(interactions$ChEMBLID)
genesx = intersect(genedat$geneid,interactions$gene_name)

cl1drugs = as.data.frame(t(sapply(drug_target_list,function(drugx){
  a = length(unique(intersect(drugx,cl1genes_h1cat)))
  b = length(unique(intersect(drugx,genesx))) - a
  c = length(unique(intersect(cl1genes_h1cat,genesx))) - a
  d = length(unique(genesx)) - a - b - c
  mydat=c(a,b,c,d)
  fi = fisher.test(matrix(mydat,nrow=2,byrow = F))
  c(mydat,fi$p.value,fi$estimate)
}))) %>%
  set_names(c('a','b','c','d','p','odds')) %>%
  mutate(padj = p.adjust(p,method ='fdr')) %>%
  mutate(ChEMBLID = names(drug_target_list))

drugage = read_csv('../melike/projects/shared_data/GenAge/20190813/data/raw/drugage.csv')
DA_drugs_incSyn <- unique(drugage$compound_name)
name2CID <- function(nm) {
  library(RCurl)
  library(jsonlite)
  nm <- URLencode(nm, reserved = T)
  name2cid <- getURL(paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                           nm, "/cids/JSON", sep = ""))
  name2cid <- fromJSON(name2cid)
  return(name2cid$IdentifierList$CID)
}
convertIDs <- function(drugList) {
  cids <- sapply(drugList, name2CID)
  cids <- cids[!sapply(cids, is.null)]
  cids <- reshape2::melt(cids) %>% rename(CID = value,
                                          drugName = L1)
  noCID <- setdiff(drugList, cids$drugName)
  return(list(cids, noCID))
}
DA_drugs_incSyn_CIDs <- convertIDs(DA_drugs_incSyn)
DA_drugs_incSyn_CIDlist <- DA_drugs_incSyn_CIDs[[1]]
cid2chembl <- function(intab, out_fx) {
  pubchem2chembl <- read_tsv("../melike/projects/shared_data/UniChem/20190813/data/raw/src1src22.txt") %>%
    setNames(., c("ChEMBLID", "CID")) %>% mutate(CID = as.character(CID))
  fintable <- intab %>% mutate(CID = as.character(CID)) %>%
    left_join(pubchem2chembl) %>% unique()
  return(fintable)
}
DA_drugs_incSyn_CHEMBLlist <- cid2chembl(DA_drugs_incSyn_CIDlist)
unique(filter(cl1drugs,p<=0.1)$ChEMBLID)
DA_drugs_incSyn_CHEMBLlist %>%
  filter(ChEMBLID %in% unique(filter(cl1drugs,p<=0.1)$ChEMBLID))
chembl2name <- function(nm){
  library(RCurl)
  library(jsonlite)
  nm <- URLencode(nm)
  dat <- getURL(paste("https://www.ebi.ac.uk/chembl/api/data/molecule/", 
                           nm, ".json", sep = ""))
  if(dat ==''){return(NA)}
  else{dat <- fromJSON(dat)
  return(unique(toupper(dat$pref_name)))}
}
chemid = unique(filter(cl1drugs,p<=0.1)$ChEMBLID)
chemname = sapply(chemid,chembl2name)
chemname = reshape2::melt(chemname) %>% set_names(c('Name','ChEMBLID')) 
chemtargets  = interactions %>%
  filter(ChEMBLID %in% chemid) %>%
  group_by(ChEMBLID) %>%
  summarise(genes = paste(sort(unique(gene_name)),collapse = ', '),
            cl1genes = paste(sort(unique(intersect(gene_name,cl1genes_h1cat))),collapse = ', '))

filter(cl1drugs,p<=0.1) %>%
  left_join(chemname) %>%
  left_join(chemtargets) %>%
  arrange(-odds,-a) %>%
  select(Name, cl1genes, odds, genes, ChEMBLID, odds, everything()) %>% View()

indat  = interactions %>% 
  filter(ChEMBLID%in%chemid & gene_name %in% genedat$geneid) %>%
  left_join(chemname) %>% select(-ChEMBLID) %>% unique()
myg = graph_from_data_frame(indat)
V(myg)$type = 'drug'
V(myg)$type[V(myg)$name %in% interactions$gene_name] = 'gene'
V(myg)$cl = '-'
V(myg)$cl[V(myg)$name %in% cl1genes_h1cat ] = 'cl1'
V(myg)$color=V(myg)$type
V(myg)$color[V(myg)$name %in% cl1genes_h1cat] = 'cl1'
labx= V(myg)$name
labx[!labx%in%genedat$geneid]=tolower(labx[!labx%in%genedat$geneid])
library(GGally)
drugnet = ggnet2(myg,size=0,edge.color='gray70')+
  geom_point(shape = c(18,19)[factor(V(myg)$type)], 
             size = c(3,2)[factor(V(myg)$type)],
             color = c('firebrick2','dodgerblue','azure3')[factor(V(myg)$color)]) +
  geom_text_repel(label=labx,box.padding = 0.01,size=6/pntnorm,
                  color = c('gray5','midnightblue','gray70')[factor(V(myg)$color)])

ggsave('./results/drugnet.pdf', drugnet,width = 16.7,height = 12,units = 'cm',useDingbats=F)
ggsave('./results/drugnet.png', drugnet,width = 16.7,height = 12,units = 'cm')

filter(cl1drugs,odds>0,a>0) %>%
  left_join(chemname) %>%
  left_join(chemtargets) %>%
  arrange(-odds,-a) %>%
  select(Name, cl1genes, odds, genes, ChEMBLID, odds, everything()) %>% select(cl1genes) %>% unique()
