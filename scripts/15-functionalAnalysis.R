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

source('../shared/functions/functions.R')

go_enrich.test <- function(genelist, selection, description = "",ontologyx='BP',
                           ID = "Ensembl", nodesize = 1,corrmethod="fdr") {
  library(topGO)
  library("org.Hs.eg.db")
  genel = new("topGOdata", ontology = ontologyx, allGenes = genelist,
              geneSel = selection, description = description, annot = annFUN.org,
              mapping = "org.Hs.eg.db", ID = ID, nodeSize = nodesize)
  resultFisher = runTest(genel, algorithm = "classic", statistic = "fisher")
  sigterms = resultFisher@geneData["SigTerms"]
  genetable = GenTable(genel, classicFisher = resultFisher,
                       topNodes = sigterms)
  ps = c()
  for (i in genetable$classicFisher) {
    if (i == "< 1e-30")
      ps = c(ps, 1e-30) else {
        ps = c(ps, i)
      }
  }
  genetable_padjusted = cbind(genetable, p.adjust((ps),method=corrmethod))
  colnames(genetable_padjusted)[7]='p.adjusted'
  genetable_padjusted$genelist=lapply(genetable_padjusted$GO.ID,function(go)intersect(genesInTerm(genel,go)[[1]],names(which(sapply(genelist,selection)))))
  return(genetable_padjusted)
}
gnls=unique(signifGenes$geneid)
genex = setNames(rep(0,length(gnls)),gnls)
genex[names(genex) %in% cl1genes] = 1
genex[names(genex) %in% cl2genes] = 2
genex[names(genex) %in% cl3genes] = 3
gores_cl1 = go_enrich.test(genelist = genex, selection = function(x)x==1)
gores_cl2 = go_enrich.test(genelist = genex, selection = function(x)x==2)
gores_cl3 = go_enrich.test(genelist = genex, selection = function(x)x==3)

gnls=unique(signifGenes$geneid)
genex = setNames(rep(0,length(gnls)),gnls)
genex[names(genex) %in% cl1genes_h1cat] = 1
genex[names(genex) %in% cl2genes_h1cat] = 2
genex[names(genex) %in% cl3genes_h1cat] = 3
gores_cl1_h1 = go_enrich.test(genelist = genex, selection = function(x)x==1)
gores_cl2_h1 = go_enrich.test(genelist = genex, selection = function(x)x==2)
gores_cl3_h1 = go_enrich.test(genelist = genex, selection = function(x)x==3)


allgores = lapply(list(gores_cl1_h1,gores_cl2_h1,gores_cl3_h1),function(x){filter(x,p.adjusted<=0.1)%>%dplyr::select(-genelist)})
names(allgores) = c('cl1','cl2','cl3')

reshape2::melt(allgores, id.vars = colnames(allgores[[1]])) %>%
  dplyr::rename(cluster = `L1`) %>% 
  write_tsv('~/Desktop/clustergo.tsv')
