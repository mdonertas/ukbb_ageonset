library(tidyverse)
library(GenomicRanges)
library(Biostrings)
# disID = 'a1521'
# disdat = read_tsv(paste('./data/processed/ukbb/gwas/bolt/',disID,'.imp.stats',sep=''))
# disdat = disdat %>%
#   select(CHR,BP,SNP,ALLELE1,ALLELE0) %>%
#   mutate(CHR = paste('chr',CHR,sep='')) %>%
#   set_names(c('seqnames','start','ID','a1','a0')) %>%
#   mutate(end = start)
# disdat = makeGRangesFromDataFrame(disdat,keep.extra.columns = T)
# disdat = disdat[which(nchar(disdat$a1) == 1 & nchar(disdat$a0) == 1),]
# 
# library(rtracklayer)
# ch = import.chain("./data/raw/ucsc/liftover/hg19ToHg38.over.chain")
# disdat = unlist(liftOver(disdat, ch))
# print('hg19 to hg38 conversion is completed')
# saveRDS(disdat,'./data/processed/ukbb/snpCoord_hg38.rds')
# print('hg38 snp coordinates are saved')
unisx = grep('.rds',list.files('../melike/projects/shared_data/uniprot2genomicCoord/data/processed/uniprot_granges_hg38', full.names = T),v=T)
library(clustermq)
readunimap = function(unix){
  library(tidyverse)
  library(GenomicRanges)
  library(Biostrings)
  disdat = readRDS('./data/processed/ukbb/snpCoord_hg38.rds')
  uni = readRDS(unix)
  uni = uni[which(uni$ismatch),]
  uni = sort(uni)
  fo = findOverlaps(disdat,uni)
  xx = as.data.frame(disdat[queryHits(fo),]) %>%
    mutate(snpPos = start) %>%
    select(ID,snpPos,a0,a1)
  xx = data.frame(xx,as.data.frame(uni[subjectHits(fo),])) %>%
    select(-width)
  return(xx)
}

mapuni = function(i,unimap){
  x=unimap[i,]
  library(tidyverse)
  library(GenomicRanges)
  library(Biostrings)
  pos = ifelse(x$strand == '+', x$snpPos - x$start + 1, 3 - x$snpPos + x$start)
  if(substr(x$genSeq,pos,pos) == x$a0) {
    a1 = x$a1
    x$a1 = x$a0
    x$a0 = a1
    ismatch_seq = (substr(x$genSeq, pos, pos) == x$a1)
  } else if ( substr(x$genSeq, pos, pos) == x$a1){
    ismatch_seq = T
  } else if ( x$strand == '-'){
    x$a0 = as.character(reverseComplement(DNAStringSet(x$a0)))
    x$a1 = as.character(reverseComplement(DNAStringSet(x$a1)))
    if(substr(x$genSeq,pos,pos) == x$a0) {
      a1 = x$a1
      x$a1 = x$a0
      x$a0 = a1
      ismatch_seq = (substr(x$genSeq, pos, pos) == x$a1)
    } else if ( substr(x$genSeq, pos, pos) == x$a1){
      ismatch_seq = T
    } 
  } else if( substr(x$genSeq,pos,pos) %in% c(as.character(reverseComplement(DNAStringSet(x$a1))),as.character(reverseComplement(DNAStringSet(x$a0))))){
    x$a0 = as.character(reverseComplement(DNAStringSet(x$a0)))
    x$a1 = as.character(reverseComplement(DNAStringSet(x$a1)))
    if(substr(x$genSeq,pos,pos) == x$a0) {
      a1 = x$a1
      x$a1 = x$a0
      x$a0 = a1
      ismatch_seq = (substr(x$genSeq, pos, pos) == x$a1)
    } else if ( substr(x$genSeq, pos, pos) == x$a1){
      ismatch_seq = T
    } 
  } else { 
    ismatch_seq = F
  }
  genseq = tolower(x$genSeq)
  genseq = strsplit(genseq, '')[[1]]
  genseq2 = genseq
  genseq2[pos] = x$a0
  genseq[pos] = toupper(genseq[pos])
  genseq = paste(genseq,collapse ='')
  genseq2 = paste(genseq2,collapse ='')
  genChange = paste(genseq,'/',genseq2,sep='')
  newAA = unname(GENETIC_CODE[toupper(genseq2)])
  type = ifelse(x$AA == newAA,'syn','missense')
  xx = as.data.frame(x) %>%
    mutate(pos = pos, newGenSeq = toupper(genseq2), genChange =genChange, newAA = newAA,type = type, ismatch_seq =ismatch_seq) 
  return(xx)
}

# unimap = Q(readunimap,unix = unisx,n_jobs = 8000)
# unimap = reshape2::melt(unimap,id.vars = colnames(unimap[[1]])) %>%
#   select(-L1) %>%
#   as.tibble()
# saveRDS(unimap,'./temp/unimap_read.rds')
unimap = readRDS('./temp/unimap_read.rds')

# this step requires too much memory and I run it in multiple steps. 
unimap2 = Q(mapuni,i=1:nrow(unimap),n_jobs = 8000,const = list(unimap=unimap))
unimap2 = reshape2::melt(unimap2,id.vars = colnames(unimap2[[1]])) %>%
  select(-L1)
unimap3 = rbind(unimap3, unimap2)
saveRDS(unimap3,'./data/processed/codingVariants/unimap.rds')

unimap3 = unimap3[which(unimap3$ismatch_seq),]
# we filter out the following entry where the databases claim it is a variant from G to A/C/T but in UK BB it is from A-T variant.
# ID  snpPos a0 a1 seqnames   start     end strand
# 13498 rs2806234 7563750  A  T     chr6 7563748 7563750      +
#   protIndex            ensp AA genSeq accession     uni_id
# 13498       247 ENSP00000369129  A    GCG    P15924 DESP_HUMAN
# name ec ismatch pos newGenSeq genChange newAA type
# 13498 Desmoplakin  -    TRUE   3       GCA   gcG/gcA     A  syn
# ismatch_seq
# 13498       FALSE
saveRDS(unimap3,'./data/processed/codingVariants/unimap.rds')
