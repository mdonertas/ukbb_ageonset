disID=as.character(commandArgs(trailingOnly = T)[1])
library(tidyverse)
library(RFRlib)
library(clusterProfiler)
library(ggrepel)

gwasCat_report <- function(gwasCat_assocFile,genelist){
  gwascat <- read_tsv(gwasCat_assocFile)%>%
    filter(`P-VALUE`<=5e-8)
  genelist <- unique(genelist[complete.cases(genelist)])
  martx=biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
  geneInfo <- biomaRt::getBM(attributes = c('hgnc_symbol','description'),
                             filters = c('hgnc_symbol'),values = genelist,mart = martx) %>% unique()%>%
    mutate(description=sapply(strsplit(description,'[[]'),function(x)x[1])) 
  gwasCatRes <- geneInfo %>%
    group_by(hgnc_symbol,description) %>%
    summarise(MAPPED_TRAIT = list(gwascat$MAPPED_TRAIT[grep(hgnc_symbol,gwascat$MAPPED_GENE)])) %>%
    unique() %>%
    dplyr::rename(GWASCATALOG=MAPPED_TRAIT)
  gwasCatRes$GWASCATALOG <- sapply(gwasCatRes$GWASCATALOG,function(x){
    paste(unique(unlist(strsplit(x, ', '))),collapse = ', ')
  })
  gwasCatRes$GWASCATALOG[gwasCatRes$GWASCATALOG=='']=NA
  return(gwasCatRes)
}

system(paste('mkdir -p data/processed/caseControl/a',disID,sep=''))
system(paste('mkdir -p results/caseControl/a',disID,sep=''))
disCoding <- read_tsv('data/raw/ukbb/datacoding/coding6.tsv')
disname <- filter(disCoding,node_id==disID)$meaning
gwasRes <- read_tsv(paste('data/processed/ukbb/gwas/bolt/a',disID,'.imp.stats',sep='')) %>%
  mutate(CHR=as.character(CHR)) %>%
  rename(Ref = ALLELE1, Alt = ALLELE0)

## qqplot
qqplot <- gwasRes %>% 
  ggplot(aes(x=-log10((rank(P_BOLT_LMM_INF, ties.method="first")-.5)/(nrow(gwasRes)+1)),y=-log10(P_BOLT_LMM_INF)))+
  geom_vline(xintercept = -log10(5e-8), color='gray', linetype='dotted') +
  geom_hline(yintercept = -log10(5e-8), color='gray', linetype='dotted') +
  geom_hex(bins = 100, aes(fill= log10(..count..)))+
  geom_abline(slope=1, color='darkred',linetype='dashed')+
  theme_rfr(legend.pos = 'top') +
  xlab(expression(paste("Expected (",-log[10], " p-value)"))) +
  ylab(expression(paste("Observed (",-log[10], " p-value)"))) +
  scale_fill_gradient(low='gray60',high='gray25') +
  guides(fill = guide_colorbar(expression(paste(log[10],"Count"))))+
  ggtitle(disname)

ggsave(filename = paste('results/caseControl/a',disID,'/qqplot.pdf',sep=''),plot = qqplot,device ='pdf',width=7,height=6)
rm(qqplot)
## combine with gene info

### proxy
proxyGenes <- readRDS('./data/processed/genomicAnalysis/snp2gene_proxy.rds') %>%
  mutate(CHR=as.character(CHR)) %>%
  filter(SNP %in% gwasRes$SNP) %>%
  select(-proxy_type) %>%
  unique() %>%
  right_join(gwasRes)

write_tsv(proxyGenes, paste('data/processed/caseControl/a',disID,'/gwasRes_proxyGenes.tsv',sep=''))
saveRDS(proxyGenes,paste('data/processed/caseControl/a',disID,'/gwasRes_proxyGenes.rds',sep=''))

signifProxy <- proxyGenes %>%
  filter(P_BOLT_LMM_INF<=5e-8)

write_tsv(signifProxy,paste('data/processed/caseControl/a',disID,'/signif_gwasRes_proxyGenes.tsv',sep=''))
saveRDS(signifProxy,paste('data/processed/caseControl/a',disID,'/signif_gwasRes_proxyGenes.rds',sep=''))

### proxyGenes

proxyGenes_pvals <- proxyGenes %>%
  filter(!is.na(proxy_entrez)) %>%
  group_by ( proxy_entrez ) %>%
  summarise(pval = max(-log10(P_BOLT_LMM_INF)))%>%
  unique()
proxyGenes_pvals = sort(setNames(proxyGenes_pvals$pval,proxyGenes_pvals$proxy_entrez),decreasing = T)
proxyGenes_pvals [ proxyGenes_pvals == Inf ] = max( setdiff( proxyGenes_pvals, Inf ) )
proxyGenes_pvals = sort(proxyGenes_pvals, decreasing = T)

saveRDS(proxyGenes_pvals,paste('data/processed/caseControl/a',disID,'/proxy_gene_pvals.rds',sep=''))

kegg_proxy <- gseKEGG(geneList = proxyGenes_pvals, nPerm = 1000, minGSSize = 20, organism = 'hsa', keyType = 'kegg', maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = 'fdr', seed = T)

kegg_proxy@result %>%
  write_tsv(paste('results/caseControl/a',disID,'/proxyKEGG.tsv',sep=''))
saveRDS(kegg_proxy,paste('data/processed/caseControl/a',disID,'/proxyKEGG.rds',sep=''))

if (length(setdiff(unique(signifProxy$proxy_hgnc),c('',NA)))>=1){
  proxyGWASCatRep <- gwasCat_report(gwasCat_assocFile='../melike/projects/shared_data/GWASCatalog/20190607/data/raw/gwas_catalog_v1.0.2-associations_e96_r2019-05-03.tsv',genelist = setdiff(signifProxy$proxy_hgnc,c('',NA)))
  
  write_tsv(proxyGWASCatRep,paste('results/caseControl/a',disID,'/proxyGWASCatReport.tsv',sep=''))
  
  signifProxy2Label <- as.tibble(t(sapply(setdiff(unique(signifProxy$proxy_hgnc),c(NA,'')),function(gene){
    xx <- filter(signifProxy,proxy_hgnc==gene)
    c(P_BOLT_LMM_INF=xx$P_BOLT_LMM_INF[which.min(xx$P_BOLT_LMM_INF)],SNP=xx$SNP[which.min(xx$P_BOLT_LMM_INF)])
  }))) %>%
    mutate(hgnc = setdiff(unique(signifProxy$proxy_hgnc),c(NA,'')),
           P_BOLT_LMM_INF=as.numeric(P_BOLT_LMM_INF))
} else{
  signifProxy2Label <- data.frame(hgnc = NA, P_BOLT_LMM_INF=NA,SNP=NA)
}

rm(list=setdiff(ls(),c('disCoding','gwasRes','signifProxy2Label','disID','disname','gwasCat_report')))
### eQTL

eQTLGenes <- readRDS('./data/processed/genomicAnalysis/snp2gene_eQTL_summary.rds') %>%
  filter(SNP %in% gwasRes$SNP)%>%
  right_join(gwasRes)

write_tsv(eQTLGenes,paste('data/processed/caseControl/a',disID,'/gwasRes_eQTLGenes.tsv',sep=''))
saveRDS(eQTLGenes,paste('data/processed/caseControl/a',disID,'/gwasRes_eQTLGenes.rds',sep=''))

signifeQTL <- eQTLGenes %>%
  filter(P_BOLT_LMM_INF<=5e-8)

write_tsv(signifeQTL,paste('data/processed/caseControl/a',disID,'/signif_gwasRes_eQTLGenes.tsv',sep=''))
saveRDS(signifeQTL,paste('data/processed/caseControl/a',disID,'/signif_gwasRes_eQTLGenes.rds',sep=''))

eQTLGenes_pvals <- eQTLGenes %>%
  filter(!is.na(eQTL_entrez)) %>%
  group_by ( eQTL_entrez ) %>%
  summarise(pval = max(-log10(P_BOLT_LMM_INF)))%>%
  unique()
eQTLGenes_pvals = sort(setNames(eQTLGenes_pvals$pval,eQTLGenes_pvals$eQTL_entrez),decreasing = T)
eQTLGenes_pvals [ eQTLGenes_pvals == Inf ] = max( setdiff( eQTLGenes_pvals, Inf ) )
eQTLGenes_pvals = sort ( eQTLGenes_pvals, decreasing = T)

saveRDS(eQTLGenes_pvals,paste('data/processed/caseControl/a',disID,'/eqtl_gene_pvals.rds',sep=''))

kegg_eQTL <- gseKEGG(geneList = eQTLGenes_pvals, nPerm = 1000, minGSSize = 20, organism = 'hsa', keyType = 'kegg', maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = 'fdr', seed = T)

kegg_eQTL@result %>%
  write_tsv(paste('results/caseControl/a',disID,'/eQTLKEGG.tsv',sep=''))
saveRDS(kegg_eQTL,paste('data/processed/caseControl/a',disID,'/eQTLKEGG.rds',sep=''))

if (length(setdiff(unique(signifeQTL$eQTL_hgnc),c('',NA)))>=1){
  signifeQTL2Label <- as.tibble(t(sapply(setdiff(unique(signifeQTL$eQTL_hgnc),c(NA,'')),function(gene){
    xx <- filter(signifeQTL,eQTL_hgnc==gene)
    c(P_BOLT_LMM_INF=xx$P_BOLT_LMM_INF[which.min(xx$P_BOLT_LMM_INF)],SNP=xx$SNP[which.min(xx$P_BOLT_LMM_INF)])
  }))) %>%
    mutate(hgnc = setdiff(unique(signifeQTL$eQTL_hgnc),c(NA,'')),
           P_BOLT_LMM_INF=as.numeric(P_BOLT_LMM_INF))
  
  
  eQTLGWASCatRep <- gwasCat_report(gwasCat_assocFile='../melike/projects/shared_data/GWASCatalog/20190607/data/raw/gwas_catalog_v1.0.2-associations_e96_r2019-05-03.tsv',genelist = unique(setdiff(signifeQTL$eQTL_hgnc,c('',NA))))
  
  write_tsv(eQTLGWASCatRep,paste('results/caseControl/a',disID,'/eQTLGWASCatReport.tsv',sep=''))
} else{
  signifeQTL2Label <- data.frame(hgnc = NA, P_BOLT_LMM_INF=NA,SNP=NA)
}

rm(list=setdiff(ls(),c('disCoding','gwasRes','signifProxy2Label','disID','disname','gwasCat_report','signifeQTL2Label')))

# combine GWAScat results

if ((!all(is.na(signifProxy2Label))) & (!all(is.na(signifeQTL2Label)))){
  gwascat <- rbind(mutate(read_tsv(paste('results/caseControl/a',
                                         disID,'/eQTLGWASCatReport.tsv',sep='')),
                          type='eQTL'),
                   mutate(read_tsv(paste('results/caseControl/a',
                                         disID,'/proxyGWASCatReport.tsv',sep='')),
                          type='proxy')) %>%
    group_by(hgnc_symbol, description, GWASCATALOG) %>%
    summarise(type = paste(type,collapse=', '))
  
  write_tsv(gwascat,paste('results/caseControl/a',disID,'/GWASCatReport_combined.tsv',sep=''))
}

### merge all

don <- gwasRes %>% 
  mutate( CHR = as.numeric(CHR)) %>%
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(mutate(gwasRes, CHR= as.numeric(CHR)), ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  rename(P=P_BOLT_LMM_INF)%>%
  filter(-log10(P)>2)
don <- don %>%
  select(SNP,CHR,BP,P,tot,BPcum)%>%
  unique()%>%
  top_n(1,-P) %>%
  right_join(don)
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

if ( !all(is.na(signifeQTL2Label$SNP))){
  signifeQTL2Label <- signifeQTL2Label %>%
    left_join(don,by='SNP')
} else{
  signifeQTL2Label <- signifeQTL2Label %>%
    mutate(BPcum = NA,
           P = NA)
}
if ( !all(is.na(signifProxy2Label$SNP))){
  signifProxy2Label <- signifProxy2Label %>%
    left_join(don,by='SNP')
} else{
  signifProxy2Label <- signifProxy2Label %>%
    mutate(BPcum = NA,
           P = NA)
}
eqtlcol='steelblue4'
proxycol='sienna4'

p <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
  geom_point( aes(color=as.factor(CHR)), alpha=0.2, size=0.2) +
  scale_color_manual(values = rep(c("gray40", "gray70"), 22 )) +
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +
  geom_hline(yintercept = -log10(5e-8), color='darkred',linetype='dotted',size=0.5)+
  geom_point(data=filter(don,P<=5e-8),color='darkred',size=0.2,alpha=0.5)+
  theme_bw() +
  ggtitle(disname)+
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = 'bottom')+
  xlab('')+ylab(expression(paste(-log[10], " p-value)")))+
  guides(color=F)
ggsave(filename = paste('results/caseControl/a',disID,'/manhattan.jpeg',sep=''),plot = p,device = 'jpeg',dpi=300,width = 8,height = 3)

if (!all(is.na(signifeQTL2Label)) ){
  p <- p + geom_text_repel(data=signifeQTL2Label,
                           aes(label=hgnc,y=-log10(P),x=BPcum),size=1,box.padding = 0.01,
                           min.segment.length = 0.1,segment.alpha = 0.5,
                           segment.color = eqtlcol,segment.size = 0.3,color=eqtlcol)+
    geom_point(data=data.frame(color=c(eqtlcol,proxycol)),aes(fill=color),shape=21,size=0,stroke=0,x=axisdf$center[1],y=3)+
    scale_fill_manual(values=c(eqtlcol,proxycol),labels=c('eQTL','proximity'))+
    guides(color=F,fill=guide_legend('Type of Association', override.aes = list(stroke=0.5,size=2)))
}
if (!all(is.na(signifProxy2Label)) ){
  p <- p + geom_text_repel(data=signifProxy2Label,
                           aes(label=hgnc,y=-log10(P),x=BPcum),size=1,box.padding = 0.01,
                           min.segment.length = 0.1,segment.alpha = 0.5,
                           segment.color = proxycol,segment.size = 0.3,color=proxycol)+
    geom_point(data=data.frame(color=c(eqtlcol,proxycol)),aes(fill=color),shape=21,size=0,stroke=0,x=axisdf$center[1],y=3)+
    scale_fill_manual(values=c(eqtlcol,proxycol),labels=c('eQTL','proximity'))+
    guides(color=F,fill=guide_legend('Type of Association', override.aes = list(stroke=0.5,size=2)))}
ggsave(filename = paste('results/caseControl/a',disID,'/manhattan_annotated.jpeg',sep=''),plot = p,device = 'jpeg',dpi=300,width = 8,height = 4)
