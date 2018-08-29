library(tidyverse)
library(data.table)

exclusions <- unique(c(read_delim('./data/processed/ukbb/gwas/remove/inplink_notin_bgen.fam',delim=' ',col_names = F)$X1,
                       read_delim('./data/processed/ukbb/gwas/remove/sampleQC_exc.fam',delim=' ',col_names = F)$X1,
                       read_delim('./data/processed/ukbb/gwas/remove/withdrawn.fam',delim=' ',col_names = F)$X1))

ukbfields <- read_tsv('./data/processed/ukbb/helperfiles/ukbfields.tsv')

cov_fields <- ukbfields %>%
  filter(field_id %in% c('eid',21000,22000,22009,54)) %>%
  mutate(fields = as.character(fields)) %>%
  filter(visit == 0)

cov_fields$names <- c('eid','centre','ethnicity','batch',paste('pc',1:40,sep=''))
cov_fieldDat <- setDF(fread('./data/raw/ukbb/ukb10904.csv',
                            select = cov_fields$fields,
                            col.names = cov_fields$names)) %>%
  filter(!eid %in% exclusions)

cov_fieldDat$batch <- gsub('-','UKBiLEVEAX',cov_fieldDat$batch)

cov_fieldDat[cov_fieldDat==-1 | cov_fieldDat ==-3] <- NA

traits <- readRDS('./data/processed/traits_clean/traitData_baseline.rds') %>%
  filter(!eid %in% exclusions)
traits <- traits %>%
  mutate(Sex = c('Female','Male')[Sex+1])%>%
  mutate(BMI = Weight/((`Standing height`/100)^2))
traits[traits==-1 | traits ==-3] <- NA

traits <- left_join(traits,cov_fieldDat)

pheno <- traits %>%
  select(eid, Sex, `Age when attended assessment centre`,BMI, centre, ethnicity, batch, paste('pc',1:20,sep=''))%>%
  rename(Age = `Age when attended assessment centre`)

colnames(pheno)[grep('pc',colnames(pheno))]=paste('ukbb_',grep('pc',colnames(pheno),v=T),sep='')

SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated.rds') %>%
  filter(!eid %in% exclusions)

diseaseCases <- SRdisease %>%
  group_by(node_id)%>%
  summarise(case=list(unique(as.character(eid))))

diseaseCols <- apply(diseaseCases,1,function(x){
  as.numeric(pheno$eid %in% x$case)
})

colnames(diseaseCols) <- paste('d',diseaseCases$node_id,sep='')

anydisease <- as.numeric(rowSums(diseaseCols)>0)

SRcancer <- readRDS('./data/processed/traits_clean/SRcancer_baseline.rds') %>%
  filter(!eid %in% exclusions)

cancerCases <- SRcancer %>%
  group_by(node_id)%>%
  summarise(case=list(unique(as.character(eid))))%>%
  mutate(node_id=paste('c',node_id,sep=''))

cancerCols <- apply(cancerCases,1,function(x){
  as.numeric(pheno$eid %in% x$case)
})

colnames(cancerCols) <- cancerCases$node_id

anycancer <- as.numeric(rowSums(cancerCols)>0)

anydisorcancer <- as.numeric(anycancer | anydisease)

disIDs <- unique((readRDS('./data/processed/traits_clean/SRdiseaseSet.rds') %>%
                    left_join(rename(read_tsv('./data/raw/ukbb/datacoding/coding6.tsv'),Disease=meaning)))$node_id)

ageCols <- lapply(disIDs,function(dis){
  filter(SRdisease,node_id==dis) %>%
    select(eid,Age,node_id)
})
ageCols <- reshape2::melt(ageCols,id.vars=c('eid','Age','node_id'))
ageCols$L1 <- NULL
ageCols <- ageCols%>%
  spread(node_id,Age)
rownames(ageCols) <- ageCols$eid
ageCols$eid <- NULL
ageCols <- as.matrix(ageCols)
colnames(ageCols) <- paste('aoo_',colnames(ageCols),sep='')

toadd <- setdiff(as.character(pheno$eid),rownames(ageCols))

ageCols <- rbind(ageCols,matrix(NA,nrow=length(toadd),ncol = ncol(ageCols),dimnames=list(toadd,colnames(ageCols))))
ageCols <- ageCols[as.character(pheno$eid),]

phenoData <- data.frame(pheno,anydisorcancer,anydisease,anycancer,diseaseCols[,paste('d',disIDs,sep='')],ageCols)

phenoData <- phenoData %>%
  mutate(FID=eid) %>%
  rename(IID=eid) %>%
  select(FID,IID,everything())

dim(phenoData)
# [1] 484666    263

phenoData %>%
  mutate(Sex = as.numeric(as.factor(Sex)),
         batch = as.numeric(as.factor(batch))) %>%
  write_delim('./data/processed/ukbb/gwas/pheno/phenoFile_noChar.pheno',delim=' ')
