library(tidyverse)

gwas2GRanges <- function(gwasRes, SNP = "RefSNP_id", start = "BP", chr = "CHR", cols2retain = c('ALLELE0', 'ALLELE1'), genome='hg19' ){
  library(tidyverse)
  library(GenomicRanges)
  gwasRes <- gwasRes %>%
    dplyr::rename(RefSNP_id = SNP,
                  start = start,
                  chr = chr) %>%
    dplyr::mutate(end=`start`)
  gwasRes <- gwasRes %>%
    dplyr::select( RefSNP_id, chr, start, end, cols2retain)
  snpinf <- makeGRangesFromDataFrame(gwasRes)
  if (substr(seqlevels(snpinf),1,3)!='chr') {
    seqlevels(snpinf)=paste('chr',seqlevels(snpinf),sep='')
  }
  genome(snpinf)=genome
  for( colnm in c('RefSNP_id', cols2retain)){
    values(snpinf)[[colnm]]=gwasRes[[colnm]]
  }
  return(snpinf)
}

martx=biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
genemap <- biomaRt::getBM(attributes = c('entrezgene_id','hgnc_symbol','ensembl_gene_id','description'),
                          mart = martx) %>%
  dplyr::rename(entrezgene = entrezgene_id)

snp2gene_eQTL <- function(gwasRes, eQTLfile, genemap){
  library(tidyverse)
  tissue <- strsplit(sapply(strsplit(eQTLfile,'/'),function(x)x[length(x)]),'[.]')[[1]][1]
  eqtl <- read_tsv(eQTLfile) %>%
    inner_join(gwasRes) %>%
    dplyr::select(SNP,CHR, BP, ALLELE1, ALLELE0, gene_id, slope, pval_beta, tissue)%>%
    unique()
  genemap <- genemap %>%
    unique()%>%
    dplyr::mutate(entrezgene=as.character(entrezgene))%>%
    dplyr::rename(gene_id = ensembl_gene_id)
  eqtl <- left_join(eqtl, genemap)
  eqtl$entrezgene[eqtl$entrezgene=='']=NA
  eqtl$hgnc_symbol[eqtl$hgnc_symbol=='']=NA
  eqtl$description[eqtl$description=='']=NA
  eqtl
}

gwasRes <- read_tsv('data/processed/ukbb/gwas/bolt/a1092.imp.stats')
gwas_as_GR <- gwas2GRanges(gwasRes, SNP = 'SNP',start = 'BP',chr = 'CHR',genome = 'hg19')

library(VariantAnnotation)
library(GenomicRanges)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
varInfo <- unique(granges(gwas_as_GR))
allvar <- locateVariants(varInfo, txdb,
                         AllVariants(promoter=PromoterVariants(upstream=1000,
                                                               downstream=1000)))
allvar <- allvar[which(!is.na(allvar$GENEID)),]
overs <- findOverlaps(gwas_as_GR,allvar)
for(colnm in colnames(mcols(gwas_as_GR))){
  mydat <-(mcols(gwas_as_GR)[[colnm]])[queryHits(overs)]
  mcols(allvar)[[colnm]] = '*'
  mcols(allvar)[[colnm]][subjectHits(overs)] <- mydat
}
geneids <- unique(allvar$GENEID)
geneids <- geneids[complete.cases(geneids)]
genemap <- genemap %>% unique() %>%
  mutate(entrezgene=as.character(entrezgene))
allvar <- as.tibble(allvar) %>%
  filter(LOCATION!='intergenic') %>%
  mutate(entrezgene=GENEID)%>%
  left_join(genemap)
rm(list=setdiff(ls(),c('allvar','gwas_as_GR','gwasRes','gwas2Granges','genemap')))
saveRDS(allvar,'./temp/allvar.rds')
proxyres <- allvar %>%
  dplyr::rename(SNP = RefSNP_id,
                CHR = seqnames,
                BP = start,
                proxy_entrez = entrezgene,
                proxy_hgnc = hgnc_symbol,
                proxy_ensembl = ensembl_gene_id,
                proxy_type = LOCATION,
                Ref = ALLELE1,
                Alt = ALLELE0) %>%
  dplyr::mutate(CHR = gsub('chr','',CHR)) %>%
  dplyr::select( SNP, CHR, BP, Ref, Alt, proxy_entrez, proxy_ensembl, proxy_hgnc, proxy_type) %>%
  unique()
print('proxy mapping is done')
# system('mkdir ./data/processed/genomicAnalysis')
saveRDS(proxyres, file='data/processed/genomicAnalysis/snp2gene_proxy.rds')
print('proxy file is saved')
rm(allvar)
rm(proxyres)
rm(gwas_as_GR)

filesx <- list.files('../melike/projects/shared_data/GTEx_v8/data/processed/signif_tissue_eQTLs_hg19_refaltamb',full.names=T)
eQTLres <- lapply(filesx,function(fx){
  snp2gene_eQTL(gwasRes,fx,genemap)
})
print('eqtl mapping is done')
saveRDS(eQTLres, file='data/processed/genomicAnalysis/snp2gene_eQTL_v1.rds')
print('eqtl_v1 data is saved')

eQTLres <- reduce(eQTLres,rbind) %>%
  dplyr::rename(Ref = ALLELE1, Alt=ALLELE0, eQTL_ensembl = gene_id, eQTL_entrez = entrezgene, eQTL_hgnc=hgnc_symbol,
                eQTL_slope=slope, eQTL_pval=pval_beta, eQTL_tissue=tissue) %>%
  dplyr::select(SNP, CHR, BP, Ref, Alt, eQTL_entrez, eQTL_hgnc, eQTL_ensembl, eQTL_slope, eQTL_pval, eQTL_tissue)%>%
  unique()
print('eqtl result is reformatted')
saveRDS(eQTLres, file='data/processed/genomicAnalysis/snp2gene_eQTL.rds')
print('eqtl data is saved')
eQTLGenes <- eQTLres %>%
  filter(eQTL_pval <= 5e-8) %>%
  mutate(CHR=as.character(CHR))

eQTLGenes2 <- eQTLGenes %>%
  group_by(SNP, eQTL_entrez) %>%
  summarise( numTissue = length(unique(eQTL_tissue)),
             direction = 2*(mean(eQTL_slope>0)-0.5))%>%
  unique()

eQTLGenes <- eQTLGenes %>%
  dplyr::select(SNP, CHR, BP, Ref, Alt, eQTL_entrez,eQTL_hgnc, eQTL_ensembl) %>%
  right_join(eQTLGenes2)
print('eqtl data is summarised')
saveRDS(eQTLGenes, file='data/processed/genomicAnalysis/snp2gene_eQTL_summary.rds')
print('eqttl summary data is saved')
print('finished')
