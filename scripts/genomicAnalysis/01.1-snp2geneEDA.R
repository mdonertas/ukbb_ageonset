source('./scripts/00-setup.R')
proxy = readRDS('./data/processed/genomicAnalysis/snp2gene_proxy.rds')
x = unique(select(proxy,SNP,proxy_hgnc)) %>%
  filter(proxy_hgnc != '') %>%
  na.omit() %>%
  group_by(proxy_hgnc) %>%
  summarise( n = length(unique(SNP)))
summary(x$n)
length(unique(proxy$SNP))/length(unique(proxy$proxy_hgnc))
x$proxy_hgnc[which.max(x$n)]
proxy_snppergene_dist_p = ggplot(x,aes(x=n)) +
  geom_histogram(color='gray80',binwidth = 0.1) +
  scale_x_log10(labels = c(1,10,100,1000,10000),breaks = c(1,10,100,1000,10000)) +
  theme_bw() +
  xlab('Number of SNPs per gene (in log scale)')

eqtl = readRDS('./data/processed/genomicAnalysis/snp2gene_eQTL_summary.rds')
y = unique(select(eqtl,SNP,eQTL_hgnc)) %>%
  filter(eQTL_hgnc != '') %>%
  na.omit() %>%
  group_by(eQTL_hgnc) %>%
  summarise( n = length(unique(SNP)))
summary(y$n)
length(unique(eqtl$SNP))/length(unique(eqtl$eQTL_hgnc))
y$eQTL_hgnc[which.max(y$n)]
eqtl_snppergene_dist_p = ggplot(y,aes(x=n)) +
  geom_histogram(color='gray80',binwidth = 0.1) +
  scale_x_log10(labels = c(1,10,100,1000,10000),breaks = c(1,10,100,1000,10000)) +
  theme_bw() +
  xlab('Number of SNPs per gene (in log scale)')

allt = rename(x,gene=proxy_hgnc,proxy = n) %>%
  full_join(rename(y, gene = eQTL_hgnc, eQTL = n)) %>%
  filter(gene!='')

co = cor.test(allt$proxy,allt$eQTL,method = 's')

proxy_eqtl_scatter = ggplot(allt,aes(x = proxy, y = eQTL)) +
  geom_abline(slope = 1, intercept = 0, color = 'darkred', size = 0.5, linetype = 'dashed') +
  geom_point(size=0.01)+
  geom_smooth(method= 'lm', se = F, size = 0.7) +
  geom_rug(size = 0.05, alpha = 0.5) +
  scale_x_log10() + scale_y_log10() +
  xlab('# of SNPs per gene by proxy\n(in log10 scale)') + 
  ylab('# of SNPs per gene by eQTL data\n(in log10 scale)') +
  ggtitle(paste('Spearman\'s correlation rho =',round(co$est,2))) +
  theme(plot.title = element_text(size = 11))

ggsave('./results/genomicAnalysis/proxy_eqtl_scatter_allgenes.pdf',proxy_eqtl_scatter,units = 'cm',width = 8,height = 8, useDingbats = F)
ggsave('./results/genomicAnalysis/proxy_eqtl_scatter_allgenes.png',proxy_eqtl_scatter,units = 'cm',width = 8,height = 8)

densityPlot = allt %>%
  gather(type, number, -gene) %>%
  ggplot(aes(x = number, fill = type)) +
  geom_density(alpha = 0.7) +
  scale_x_log10() +
  scale_fill_manual(values = c('steelblue4','sienna4')) +
  xlab('# of SNPs per gene\n(in log10 scale)') +
  guides(fill = guide_legend('')) 

ggsave('./results/genomicAnalysis/proxy_eqtl_density_allgenes.pdf',densityPlot,units = 'cm',width = 8,height = 8, useDingbats = F)
ggsave('./results/genomicAnalysis/proxy_eqtl_density_allgenes.png',densityPlot,units = 'cm',width = 8,height = 8)

p1 = ggarrange(densityPlot,proxy_eqtl_scatter, labels = 'auto', align = 'hv')

ggsave('./results/genomicAnalysis/proxy_eqtl_allgenes.pdf',p1,units = 'cm',width = 16,height = 8, useDingbats = F)
ggsave('./results/genomicAnalysis/proxy_eqtl_allgenes.png',p1,units = 'cm',width = 16,height = 8)


# data:  allt$proxy and allt$eQTL
# S = 1.4365e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.1295103 

head(allt)

martx=biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
genesize <- biomaRt::getBM(attributes = c('hgnc_symbol','start_position',
                                          'end_position'), 
                           filter = 'hgnc_symbol', values = unique(allt$gene),
                           mart = martx)
genesize <- genesize %>%
  rename(gene = hgnc_symbol) %>%
  right_join(allt) %>%
  mutate(length = abs(end_position-start_position))

cor.test(genesize$proxy,genesize$length,method='s')
# Spearman's rank correlation rho
# 
# data:  genesize$proxy and genesize$length
# S = 1.7762e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8674424
cor.test(genesize$eQTL,genesize$length,method='s')
# Spearman's rank correlation rho
# 
# data:  genesize$eQTL and genesize$length
# S = 4.3328e+11, p-value = 0.0001073
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.03284316

promoter_proxy = unique(select(proxy,SNP,proxy_hgnc,proxy_type)) %>%
  filter(proxy_hgnc != '') %>%
  na.omit() %>%
  filter(proxy_type == 'promoter') %>%
  group_by(proxy_hgnc) %>%
  summarise( n = length(unique(SNP)))

promoter_proxy = promoter_proxy %>%
  rename(gene = proxy_hgnc, promoter_proxy = n) %>%
  right_join(genesize)

promoter_proxy %>%
  select(promoter_proxy,proxy,eQTL,length) %>%
  cor(method = 's', use = 'pairwise')
# promoter_proxy      proxy       eQTL     length
# promoter_proxy     1.00000000 0.42549223 0.07586191 0.20897216
# proxy              0.42549223 1.00000000 0.08706026 0.86744239
# eQTL               0.07586191 0.08706026 1.00000000 0.03284316
# length             0.20897216 0.86744239 0.03284316 1.00000000