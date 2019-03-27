library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)
library(tidyverse)
# prepare coding - noncoding lists
snpinfo <- read_tsv('./data/processed/ukbb/gwas/bolt/a1071.imp.stats')
snpdat <- snpinfo %>%
  select(CHR,BP,everything()) %>%
  rename(chr = CHR, start = BP) %>%
  mutate(end = start) %>%
  select(chr,start,end,everything()) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

seqlevels(snpdat) = paste("chr", seqlevels(snpdat), sep = "")
genome(snpdat) = "hg19"
codingvar <- locateVariants(snpdat, txdb, CodingVariants())
ovx <- findOverlaps(codingvar,snpdat)
codingvarx <- snpdat[subjectHits(ovx),]$SNP
noncodingvarx <- setdiff(snpdat$SNP,codingvarx)

codingvar <- codingvar[queryHits(ovx),]
codingvar$SNP <- snpdat[subjectHits(ovx),]$SNP
system('mkdir -p data/processed/SNPinfo')
saveRDS(codingvarx,'./data/processed/SNPinfo/codingSNPIDs.rds')
saveRDS(noncodingvarx,'./data/processed/SNPinfo/noncodingSNPIDs.rds')
saveRDS(codingvar,'./data/processed/SNPinfo/codingSNPinfo.rds')