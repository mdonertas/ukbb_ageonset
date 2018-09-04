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

snp2gene_proxy <- function(variantInfo, genome='hg19', upstream=5000,downstream=5000){
  library(VariantAnnotation)
  library(GenomicRanges)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  varInfo <- unique(granges(variantInfo))
  allvar <- locateVariants(varInfo, txdb,
                           AllVariants(promoter=PromoterVariants(upstream=upstream,
                                                                 downstream=downstream)))
  overs <- findOverlaps(variantInfo,allvar)
  for(colnm in colnames(mcols(variantInfo))){
    mydat <-(mcols(variantInfo)[[colnm]])[queryHits(overs)]
    mcols(allvar)[[colnm]] = '*'
    mcols(allvar)[[colnm]][subjectHits(overs)] <- mydat
  }
  geneids <- unique(allvar$GENEID)
  geneids <- geneids[complete.cases(geneids)]
  martx=biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
  genemap <- biomaRt::getBM(attributes = c('entrezgene','hgnc_symbol','ensembl_gene_id','description'),
                            filters = c('entrezgene'),values = geneids,mart = martx) %>% unique()%>%
    mutate(entrezgene=as.character(entrezgene))
  allvar <- as.tibble(allvar) %>%
    mutate(entrezgene=GENEID)%>%
    left_join(genemap)
  return(allvar)
}

snp2gene_eQTL <- function(gwasRes, eQTLfile, tissue){
  eqtl <- read_tsv(eQTLfile) %>%
    dplyr::rename(BP = variant_pos,
                  CHR = chr,
                  ALLELE0 = alt,
                  ALLELE1 = ref) %>%
    dplyr::select(CHR, BP, ALLELE0, ALLELE1, gene_id, tss_distance, slope, pval_beta)%>%
    dplyr::mutate(gene_id=sapply(strsplit(gene_id,'[.]'),function(x)x[1]))%>%
    inner_join(gwasRes) %>%
    dplyr::select(SNP,CHR, BP, ALLELE1, ALLELE0, gene_id, tss_distance, slope, pval_beta)%>%
    mutate(tissue = tissue) %>%
    unique()
  martx=biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
  genemap <- biomaRt::getBM(attributes = c('entrezgene','hgnc_symbol','ensembl_gene_id','description'),
                            filters = c('ensembl_gene_id'),values = unique(eqtl$gene_id),mart = martx) %>%
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

proxyres <- snp2gene_proxy(gwas_as_GR, genome='hg19', upstream = 5000, downstream = 5000) %>%
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
  dplyr::select( SNP, CHR, BP, Ref, Alt, proxy_entrez, proxy_hgnc, proxy_ensembl, proxy_type) %>%
  unique() %>%
  filter(proxy_type!='intergenic')

saveRDS(proxyres, file='data/processed/pilot/snp2gene_proxy.rds')

filesx <- list.files('../melike/projects/shared_data/eQTL_GTEx_20180904/data/processed/signif_tissue_eQTLs/',full.names=T)
eQTLres <- lapply(filesx,function(fx){
  tis <- strsplit(sapply(strsplit(fx,'/'),function(x)x[length(x)]),'[.]')[[1]][1]
  snp2gene_eQTL(gwasRes,fx,tis)
})

eQTLres <- reshape2::melt(eQTLres,id.vars=colnames(eQTLres[[1]])) %>%
  dplyr::rename(Ref = ALLELE1, Alt=ALLELE0, eQTL_ensembl = gene_id, eQTL_entrez = entrezgene, eQTL_hgnc=hgnc_symbol,
         eQTL_slope=slope, eQTL_pval=pval_beta, eQTL_tissue=tissue) %>%
  dplyr::select(SNP, CHR, BP, Ref, Alt, eQTL_entrez, eQTL_hgnc, eQTL_ensembl, eQTL_slope, eQTL_pval, eQTL_tissue)%>%
  unique()

saveRDS(eQTLres, file='data/processed/pilot/snp2gene_eQTL.rds')