source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)
genegrs <- readRDS('./data/processed/clustergenes_hgnc.rds')
genegr <- setNames(genegrs$genelist,genegrs$name)
eqtlgenes = readRDS('./data/processed/caseControl/a1071/gwasRes_eQTLGenes.rds')
eqtlgenes = setdiff(unique(eqtlgenes$eQTL_hgnc),NA)
proxygenes = readRDS('./data/processed/caseControl/a1071/gwasRes_proxyGenes.rds')
proxygenes = setdiff(unique(proxygenes$proxy_hgnc),NA)

allgenes = unique(c(eqtlgenes,proxygenes))
rm(eqtlgenes,proxygenes)
####

interactions <- read_tsv("http://www.dgidb.org/data/interactions.tsv") %>% 
  rename(ChEMBLID = drug_chembl_id) %>% dplyr::select(gene_name, 
                                                      ChEMBLID) %>% unique() %>% na.omit()
drug_target_list = lapply(unique(interactions$ChEMBLID),function(drug){
  unique(filter(interactions,ChEMBLID==drug)$gene_name)
})
length(drug_target_list)
#6277
names(drug_target_list)=unique(interactions$ChEMBLID)
genesx = unique(intersect(allgenes,interactions$gene_name))
length(genesx)
# 2385
genegr$cl12combined_all=unique(c(genegr$cl1_all,genegr$`cl1-2_all`,genegr$cl2_all))
genegr$cl12combined_multicat =unique(c(genegr$cl1_multicat,genegr$`cl1-2_multicat`,genegr$cl2_multicat))
drugres = lapply(drug_target_list,function(drugx){
  xx = as.data.frame(t(sapply(genegr[13:14],function(geneset){
    a = unique(intersect(intersect(drugx,geneset),genesx))
    b = unique(intersect(setdiff(drugx,geneset),genesx))
    c = unique(intersect(setdiff(geneset,drugx),genesx))
    d = unique(setdiff(genesx,unique(c(a,b,c))))
    mydat=c(length(a),length(b),length(c),length(d))
    fi = fisher.test(matrix(mydat,nrow=2,byrow = T))
    c(mydat,fi$p.value,fi$estimate)
  }))) %>%
    set_names(c('a','b','c','d','p','odds')) %>%
    mutate(name = names(genegr)[13:14])
  xx
})
drugres = reshape2::melt(drugres,id.vars = colnames(drugres[[1]])) %>%
  rename(ChEMBLID =L1)
# drugres$ChEMBLID = rownames(drugres)
# drugage = read_csv('../melike/projects/shared_data/GenAge/20190813/data/raw/drugage.csv')
drugage = read_csv('./data/raw/drugage.csv')
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
  # pubchem2chembl <- read_tsv("../melike/projects/shared_data/UniChem/20190813/data/raw/src1src22.txt") %>%
  pubchem2chembl <- read_tsv("./data/raw/src1src22.txt") %>%
    setNames(., c("ChEMBLID", "CID")) %>% mutate(CID = as.character(CID))
  fintable <- intab %>% mutate(CID = as.character(CID)) %>%
    left_join(pubchem2chembl) %>% unique()
  return(fintable)
}
DA_drugs_incSyn_CHEMBLlist <- cid2chembl(DA_drugs_incSyn_CIDlist)
drugres = drugres %>%
separate(name, remove = F, sep='_',into = c('cluster','type')) %>%
# mutate(cluster = factor(gsub('cl','',cluster),levels=c('1','2','3','1-2','1-3','2-3','1-2-3'))) %>%
mutate(type = factor(c('Multidisease','Multicategory')[as.numeric(as.factor(type))],
levels = c('Multidisease','Multicategory')))

drugxx = drugres %>%
  filter(a+b <=10) %>%
  filter(type == 'Multicategory') %>%
  # mutate(padj = p.adjust(p, method='fdr')) %>%
  # filter(padj<=0.1)
  filter(p<=0.01 | odds == Inf)

DA_drugs_incSyn_CHEMBLlist %>%
  filter(ChEMBLID %in% drugxx$ChEMBLID)
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
chemid = drugxx$ChEMBLID
chemname = sapply(chemid,chembl2name)
chemname = reshape2::melt(chemname) %>% set_names(c('Name','ChEMBLID'))
# chemname = data.frame(Name = unname(chemname), ChEMBLID=names(chemname))
drugxx = drugxx %>%
  left_join(chemname) %>%
  mutate(Name = as.character(Name),
         ChEMBLID = as.character(ChEMBLID))
drugxx$Name[which(is.na(drugxx$Name))]=drugxx$ChEMBLID[which(is.na(drugxx$Name))]


chemtargets  = interactions %>%
  filter(ChEMBLID %in% chemid) %>%
  group_by(ChEMBLID) %>%
  summarise(genes = paste(sort(unique(gene_name)),collapse = ', '))

indat  = interactions %>% 
  filter(ChEMBLID%in%chemid & gene_name %in% allgenes) %>%
  left_join(chemname) %>% unique() %>%
  mutate(Name = as.character(Name),
         ChEMBLID = as.character(ChEMBLID))
indat$Name[which(is.na(indat$Name))]=indat$ChEMBLID[which(is.na(indat$Name))]


myg = indat %>%
  mutate(Name = ifelse(grepl('insulin',Name,ignore.case = T),'Insulin',Name))%>%
  select(-ChEMBLID) %>% unique( )%>%
  graph_from_data_frame()
# myg = graph_from_data_frame(indat)
V(myg)$type = 'drug'
V(myg)$type[V(myg)$name %in% interactions$gene_name] = 'gene'
V(myg)$color=V(myg)$type
V(myg)$color[V(myg)$name %in% genegr$cl1_multicat] = 'cl1'
V(myg)$color[V(myg)$name %in% genegr$`cl1-2_multicat`] = 'cl1-2'
labx= V(myg)$name
# labx[V(myg)$type=='drug']=''
labx[!labx%in%allgenes]=tolower(labx[!labx%in%allgenes])
library(GGally)
drugnet = ggnet2(myg,size=0,edge.color='gray70')+
  geom_point(shape = c(18,19)[factor(V(myg)$type)], 
             size = c(4,4,3,2)[factor(V(myg)$color)],
             color = setNames(c(ageonsetcolors[c('1','1-2')],'gray70','dodgerblue'),c('cl1','cl1-2','gene','drug'))[V(myg)$color]) +
  geom_text_repel(label=labx,
                  box.padding = 0.01,
                  size=6/pntnorm,
                  color = c('gray5','gray5','midnightblue',NA)[factor(V(myg)$color)]) +
  theme_void()

ggsave('./results/drug/drugnet.pdf', drugnet,width = 16,height = 12,units = 'cm',useDingbats=F)
ggsave('./results/drug/drugnet.png', drugnet,width = 16,height = 12,units = 'cm')

drugind = read_tsv('./data/raw/chembl_drug_indications.tsv') %>%
  rename(ChEMBLID = `Parent Molecule ChEMBL ID`)

xx = indat %>%
  group_by(Name,ChEMBLID) %>%
  summarise(genes_in_multicatCl1 = paste(sort(unique(intersect(gene_name,genegr$cl1_multicat))),collapse=', '),
            genes_in_multicatCl2 = paste(sort(unique(intersect(gene_name,genegr$cl2_multicat))),collapse=', '),
            genes_in_multicatCl12 = paste(sort(unique(intersect(gene_name,genegr$`cl1-2_multicat`))),collapse=', '),
            otherGenes = paste(sort(unique(setdiff(gene_name,unique(c(genegr$cl1_multicat,genegr$cl2_multicat,genegr$`cl1-2_multicat`))))),collapse=', ')) %>% setNames(c('Drug Name','ChEMBLID','Multicategory Cluster 1 Genes',                                             'Multicategory Cluster 2 Genes',                                    'Multicategory Cluster 1 & 2 Genes',                                                                      'Other Genes')) %>%
  left_join(drugind) %>%
  mutate_all(list(~na_if(.,"")))
write_tsv(xx,'results/drug/signifDrugs_multicat12combined_p001_oddsInf.tsv')
write_csv(xx,'results/drug/signifDrugs_multicat12combined_p001_oddsInf.csv')

xx = read_csv('results/drug/signifDrugs_multicat12combined_p001_oddsInf.csv')


unique((xx %>%
  ungroup() %>%
  # mutate(indication = ifelse(is.na(`EFO Terms`),'-',`EFO Terms`))%>%
  filter(`Max Phase for Indication`>=4) %>%
  select(`Drug Name`,`EFO Terms`) %>%
  unique())$`EFO Terms`)
  # mutate(drug = paste(`Drug Name`,' (',`ChEMBLID`,'; ',`Max Phase for Indication`,')',sep=''))%>%
  # select(`Drug Name`,`EFO Terms`,`ChEMBLID`) %>%
  # unique() %>%
  # na.omit() %>% head()
  # group_by(`EFO Terms`) %>%
  # summarise(n=length(unique(ChEMBLID)),
  #           Drugs = paste(sort(unique(`drug`)),collapse=', ')) %>%
  # View()

##### other drugs with the same indications
indicationlist = unique((xx %>%
                           ungroup() %>%
                           # mutate(indication = ifelse(is.na(`EFO Terms`),'-',`EFO Terms`))%>%
                           filter(`Max Phase for Indication`>=4) %>%
                           select(`Drug Name`,`EFO Terms`) %>%
                           unique())$`EFO Terms`)


library(ggthemes)
inddat = drugind %>%
  filter(`Max Phase for Indication`>=4) %>%
  filter(`EFO Terms` %in% indicationlist) %>%
  left_join(interactions) %>%
  select(ChEMBLID, `EFO Terms`, gene_name) %>%
  filter(gene_name %in% allgenes) %>%
  mutate(cl12 = gene_name %in% genegr$cl12combined_multicat) %>%
  mutate(hit_drugs = ChEMBLID %in% unique(drugxx$ChEMBLID)) 

ind_p1 = select(inddat, ChEMBLID, `EFO Terms`, hit_drugs) %>%
  unique() %>%
  mutate(`EFO Terms` = fct_reorder(`EFO Terms`, hit_drugs, function(x)mean(x,na.rm=T))) %>% 
  mutate(hit_drugs = c('other','significant hit')[1+hit_drugs]) %>%
  ggplot(aes(x = `EFO Terms`, fill = hit_drugs)) + 
  geom_bar(position = 'fill') +
  coord_flip() +
  xlab('') +
  ylab('% Approved Drugs') +
  scale_fill_manual(values = c('gray70','dodgerblue')) +
  guides(fill = guide_legend(title=NULL)) +
  theme(legend.position = 'bottom') 

ind_p2 = select(inddat, gene_name, `EFO Terms`, cl12) %>%
  unique() %>%
  mutate(`EFO Terms` = fct_reorder(`EFO Terms`, cl12, function(x)mean(x,na.rm=T))) %>% 
  mutate(cl12 = factor(c('other','multicat cl1 or cl2')[1+cl12],levels = c('other','multicat cl1 or cl2'))) %>%
  ggplot(aes(x = `EFO Terms`, fill = cl12)) + 
  geom_bar(position = 'fill') +
  coord_flip() +
  xlab('') +
  ylab('% Unique Targets') +
  scale_fill_manual(values = c('gray70','dodgerblue')) +
  guides(fill = guide_legend(title=NULL)) +
  theme(legend.position = 'bottom') 

indp=ggarrange(ind_p1,ind_p2, labels = 'auto', ncol = 1, nrow = 2)

ggsave('./results/drug/indication_perc.pdf',indp, units = 'cm', width = 16, height = 16, useDingbats = F)
ggsave('./results/drug/indication_perc.png', indp, units = 'cm', width = 16, height = 16)


druglist = xx$ChEMBLID


chembldat <- function(nm){
  library(RCurl)
  library(jsonlite)
  nm <- URLencode(nm)
  dat <- getURL(paste("https://www.ebi.ac.uk/chembl/api/data/molecule/", 
                      nm, ".json", sep = ""))
  if(dat ==''){return(NA)}
  else{dat <- fromJSON(dat)
  return(dat)}
}
chembldata = lapply(druglist,chembldat)

names(chembldata) = xx$`Drug Name`

atcnames = setNames(c('Alimentary tract and metabolism',
'Blood and blood forming organs',
'Cardiovascular system',
'Dermatologicals',
'Genito urinary system and sex hormones',
'Systemic hormonal preparations, excluding sex hormones and insulins',
'Antiinfective for systemic use',
'Antineoplastic and immunomodulating agents',
'Musculo-skeletal system',
'Nervous system',
'Antiparasitic products, insecticides and repellents',
'Respiratory system',
'Sensory organs',
'Various'), c('A', 'B', 'C', 'D', 'G', 'H', 'J', 'L', 'M', 'N', 'P', 'R', 'S', 'V'))

chembldata$CHEMBL56253$atc_classifications
atc = sapply(chembldata,function(x)unique(substr(x$atc_classifications,1,1))) %>%
  reshape2::melt() %>%
  unique() %>%
  mutate(atc = factor(atcnames[as.character(value)],levels = rev(atcnames))) %>%
  ggplot(aes(x = atc)) +
  geom_bar(aes(fill = value)) +
  coord_flip() +
  guides(fill = F) +
  xlab(NULL) + ylab('# of drugs') +
  scale_fill_gdocs() +
  ggtitle('ATC Classification Level I')
sapply(chembldata,function(x)unique(substr(x$atc_classifications,1,1))) %>%
  reshape2::melt() %>%
  unique() %>%
  mutate(atc = factor(atcnames[as.character(value)],levels = rev(atcnames))) %>%
  View()

xx %>%
select(`Drug Name`,`EFO Terms`, `Max Phase for Indication`) %>%
na.omit() %>%
set_names(c('drug','indications','phase')) %>%
  head()
filter(phase >=3) %>%
group_by(indications) %>%
summarise(n = length(unique(drug))) %>%
ggplot(aes(x = reorder(indications,n), y= n)) +
geom_bar(stat = 'identity') +
coord_flip()
