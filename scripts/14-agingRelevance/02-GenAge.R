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

genagehuman = read_csv('../melike/projects/shared_data/GenAge/20190813/data/raw/genage_human.csv')
genagemodel = read_tsv('../melike/projects/shared_data/GenAge/20190813/data/raw/genage_models_orthologs_export.tsv')
drugage = read_csv('../melike/projects/shared_data/GenAge/20190813/data/raw/drugage.csv')

humangenes = unique(genagehuman$symbol)
modelgenes = unique(genagemodel$Symbol)
wormgenes = unique((filter(genagemodel,`Model Organism`=="Caenorhabditis elegans"))$Symbol)
micegenes = unique((filter(genagemodel,`Model Organism`=="Mus musculus"))$Symbol)
flygenes = unique((filter(genagemodel,`Model Organism`=="Drosophila melanogaster"))$Symbol)
yeastgenes = unique((filter(genagemodel,`Model Organism`=="Saccharomyces cerevisiae"))$Symbol)

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
interactions <- read_tsv("http://www.dgidb.org/data/interactions.tsv") %>% 
  rename(ChEMBLID = drug_chembl_id) %>% dplyr::select(gene_name, 
                                                      ChEMBLID) %>% unique() %>% na.omit()
drugagegenes = unique((interactions %>% filter(ChEMBLID %in% setdiff(DA_drugs_incSyn_CHEMBLlist$ChEMBLID,NA)))$gene_name)

calcprop = function(genelist){
  sapply(list(humangenes,modelgenes,micegenes,wormgenes,flygenes,yeastgenes,drugagegenes),function(y){
    mean(genelist%in%y)*100
  })
}

restable = t(sapply(list(cl1genes,cl1genes_h1cat,cl2genes,cl2genes_h1cat,cl3genes,cl3genes_h1cat),function(x){
  calcprop(x)
}))
restable = cbind(restable,sapply(list(cl1genes,cl1genes_h1cat,cl2genes,cl2genes_h1cat,cl3genes,cl3genes_h1cat),length))
colnames(restable)= c('human','anymodel','micegenes','wormgenes','flygenes','yeastgenes','drug-target', 'totnum')
rownames(restable) = c('cl1genes','cl1genes_h1cat','cl2genes','cl2genes_h1cat','cl3genes','cl3genes_h1cat')

restable[c(1,3,5),]
restable[c(2,4,6),]

permtables = lapply(list(cl1genes,cl1genes_h1cat,cl2genes,cl2genes_h1cat,cl3genes,cl3genes_h1cat),function(x){
  xx = t(sapply(1:1000,function(i){
    calcprop(sample(unique(genedat$geneid),length(x)))
  }))
  colnames(xx)= c('human','anymodel','micegenes','wormgenes','flygenes','yeastgenes','drug-target')
  xx
})

names(permtables) = c('cl1genes','cl1genes_h1cat','cl2genes','cl2genes_h1cat','cl3genes','cl3genes_h1cat')

prestable = t(sapply(names(permtables),function(nm){
  xx = permtables[[nm]]
  sapply(colnames(xx),function(setx){
    mean(xx[,setx]>=restable[nm,setx])
  })
}))


prestable[c(1,3,5),]
prestable[c(2,4,6),]
# human anymodel micegenes wormgenes flygenes yeastgenes drug-target
# cl1genes_h1cat 0.392    0.739     0.377     0.838        1          1       0.078
# cl2genes_h1cat 1.000    1.000     1.000     1.000        1          1       0.529
# cl3genes_h1cat 0.588    1.000     1.000     1.000        1          1       0.806