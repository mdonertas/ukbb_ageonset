source('./scripts/00-setup.R')
traits <- readRDS('./data/processed/traits_clean/traitData_baseline_additions2.rds')
SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated_selected.rds')
prevDF = readRDS('./data/processed/traits_clean/SRdisease_prop_prevdf.rds')
disSet <- readRDS('./data/processed/traits_clean/SRdiseaseSet.rds')
disSet <- prevDF %>%
  filter(Disease %in% disSet$Disease)
disMat <- readRDS('./data/processed/traits_clean/disMat.rds')
disAgeMat <- readRDS('./data/processed/traits_clean/disAgeMat.rds')
RRtable <- readRDS('./data/processed/diseaseCooccur/RRtable.rds')
ORtable <- readRDS('./data/processed/diseaseCooccur/ORtable.rds')
cormat <- readRDS('./data/processed/diseaseCooccur/disCorrelations.rds')

rr2 <- RRtable %>%
  filter(uplevel == F & sublevel == F) %>%
  filter((RR > 1 & ll > 1) | (RR < 1 & ul < 1)) %>%
  mutate(RR = log2(RR)) %>%
  mutate(RR_rank = dense_rank(abs(RR)),
         RR_type = sign(RR)) %>%
  select(disA, disB, RR, RR_rank, RR_type)

or2 <- ORtable %>%
  filter(uplevel == F & sublevel == F) %>%
  filter((OR > 1 & ll > 1) | (OR < 1 & ul < 1)) %>%
  mutate(OR = log2(OR)) %>%
  mutate(OR_rank = dense_rank(abs(OR)),
         OR_type = sign(OR)) %>%
  select(disA, disB, OR, OR_rank, OR_type)

sumx <- full_join(rr2, or2)

or_rr_correlation = sumx %>%
  ggplot(aes(x = RR_rank, y = OR_rank)) + 
  geom_hex() +
  geom_abline(slope = 1, intercept = 0, color = 'darkred', size = 1) +
  coord_fixed() + 
  xlab('RR rank') + ylab('OR rank') +
  theme(legend.position = 'right') +
  ggtitle('Correlation between OR and RR') +
  scale_fill_viridis_c()

ggsave('./results/diseaseCooccur/or_rr_correlation.pdf', or_rr_correlation, 
       useDingbats = F, units = 'cm', width = 8, height = 6)
ggsave('./results/diseaseCooccur/or_rr_correlation.png', or_rr_correlation, 
       units = 'cm', width = 8, height = 6)

# decide to use RR 

cor2 <- reshape2::melt(cormat) %>%
  setNames(c('disA','disB','phi'))

sumx <- left_join(rr2, cor2)

## all 
xx <-  sumx %>%
  filter(sign(phi) == RR_type) %>%
  filter(abs(RR) >= 0)

xx <-  sumx %>%
  filter(sign(phi) == RR_type) %>%
  filter(disA %in% unique(c(xx$disA,xx$disB)) & disB %in% unique(c(xx$disA,xx$disB))) 

xxmat <- select(xx, disA, disB, RR) %>%
  spread(disB, RR, fill = 0) %>%
  as.data.frame()

rownames(xxmat) = xxmat$disA
xxmat$disA = NULL
xxmat = as.matrix(xxmat)
hcx <- hclust(dist(xxmat))
xx <- xx %>%
  mutate(disA = factor(disA, levels = hcx$labels[hcx$order]),
         disB = factor(disB, levels = hcx$labels[hcx$order]))

discolDF <- tibble(disA = names(disTreecl), disB = names(disTreecl), disCat = unname(disTreecl)) %>%
  mutate(disCol = discatcolors[disCat]) %>%
  filter(disA %in% xx$disA & disB %in% xx$disB) %>%
  mutate(disA = factor(disA, levels = levels(xx$disA)),
         disB = factor(disB, levels = levels(xx$disB)))

finres <- ggplot(xx, aes(x = disA, y = disB)) +
  geom_point(aes(color = RR, size = abs(phi)), shape = 15) + 
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(0.1,3), breaks = c(0.1,0.2,0.3)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 8),
        legend.position = 'top') +
  xlab('') + ylab('') +
  guides(size = guide_legend(expression(~phi))) +
  coord_fixed() +
  geom_tile(data = discolDF, aes(fill = disCol)) +
  scale_fill_identity() +
  theme(panel.grid.major = element_line(size = 0.1, linetype = 'solid'))
ggsave('./results/diseaseCooccur/disAssoc.pdf', finres, 
       useDingbats = F, units = 'cm', width = 30, height = 32)
ggsave('./results/diseaseCooccur/disAssoc.png', finres, 
       units = 'cm', width = 30, height = 32)

nolabel_all <- ggplot(xx, aes(x = disA, y = disB)) +
  geom_point(aes(color = RR, size = abs(phi)), shape = 15) + 
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(0.1,0.5), breaks = c(0.1,0.2,0.3)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'right') +
  xlab('') + ylab('') +
  guides(size = guide_legend(expression(~phi))) +
  coord_fixed() +
  geom_tile(data = discolDF, aes(fill = disCol)) +
  scale_fill_identity() +
  theme(panel.grid.major = element_line(size = 0.1, linetype = 'solid'),
        axis.line = element_blank(), 
        panel.border = element_rect(size = 0.5, color = 'gray35', fill = NA)) 
ggsave('./results/diseaseCooccur/disAssoc_nolabel.pdf', nolabel_all, 
       useDingbats = F, units = 'cm', width = 10, height = 8)
ggsave('./results/diseaseCooccur/disAssoc_nolabel.png', nolabel_all, 
       units = 'cm', width = 10, height = 8)
## RR>1
xx <-  sumx %>%
  filter(sign(phi) == RR_type) %>%
  filter(abs(RR) >= 1)

xx <-  sumx %>%
  filter(sign(phi) == RR_type) %>%
  filter(disA %in% unique(c(xx$disA,xx$disB)) & disB %in% unique(c(xx$disA,xx$disB))) 

xxmat <- select(xx, disA, disB, RR) %>%
  spread(disB, RR, fill = 0) %>%
  as.data.frame()

rownames(xxmat) = xxmat$disA
xxmat$disA = NULL
xxmat = as.matrix(xxmat)
hcx <- hclust(dist(xxmat))
xx <- xx %>%
  mutate(disA = factor(disA, levels = hcx$labels[hcx$order]),
         disB = factor(disB, levels = hcx$labels[hcx$order]))

discolDF <- tibble(disA = names(disTreecl), disB = names(disTreecl), disCat = unname(disTreecl)) %>%
  mutate(disCol = discatcolors[disCat]) %>%
  filter(disA %in% xx$disA & disB %in% xx$disB) %>%
  mutate(disA = factor(disA, levels = levels(xx$disA)),
         disB = factor(disB, levels = levels(xx$disB)))

finres <- ggplot(xx, aes(x = disA, y = disB)) +
  geom_point(aes(color = RR, size = abs(phi)), shape = 15) + 
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(0.1,3), breaks = c(0.1,0.2,0.3)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 8),
        legend.position = 'top') +
  xlab('') + ylab('') +
  guides(size = guide_legend(expression(~phi))) +
  coord_fixed() +
  geom_tile(data = discolDF, aes(fill = disCol)) +
  scale_fill_identity() +
  theme(panel.grid.major = element_line(size = 0.1, linetype = 'solid'))
ggsave('./results/diseaseCooccur/disAssoc_log2RR1.pdf', finres, 
       useDingbats = F, units = 'cm', width = 25, height = 27)
ggsave('./results/diseaseCooccur/disAssoc_log2RR1.png', finres, 
       units = 'cm', width = 25, height = 27)

## RR>2
xx <-  sumx %>%
  filter(sign(phi) == RR_type) %>%
  filter(abs(RR) >= 2)

xx <-  sumx %>%
  filter(sign(phi) == RR_type) %>%
  filter(disA %in% unique(c(xx$disA,xx$disB)) & disB %in% unique(c(xx$disA,xx$disB))) 

xxmat <- select(xx, disA, disB, RR) %>%
  spread(disB, RR, fill = 0) %>%
  as.data.frame()

rownames(xxmat) = xxmat$disA
xxmat$disA = NULL
xxmat = as.matrix(xxmat)
hcx <- hclust(dist(xxmat))
xx <- xx %>%
  mutate(disA = factor(disA, levels = hcx$labels[hcx$order]),
         disB = factor(disB, levels = hcx$labels[hcx$order]))

discolDF <- tibble(disA = names(disTreecl), disB = names(disTreecl), disCat = unname(disTreecl)) %>%
  mutate(disCol = discatcolors[disCat]) %>%
  filter(disA %in% xx$disA & disB %in% xx$disB) %>%
  mutate(disA = factor(disA, levels = levels(xx$disA)),
         disB = factor(disB, levels = levels(xx$disB)))

log2rr2 <- ggplot(xx, aes(x = disA, y = disB)) +
  geom_point(aes(color = RR, size = abs(phi)), shape = 15) + 
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(0.1,3), breaks = c(0.1,0.2,0.3)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 8),
        legend.position = 'top') +
  xlab('') + ylab('') +
  guides(size = guide_legend(expression(~phi))) +
  coord_fixed() +
  geom_tile(data = discolDF, aes(fill = disCol)) +
  scale_fill_identity() +
  theme(panel.grid.major = element_line(size = 0.1, linetype = 'solid'))
ggsave('./results/diseaseCooccur/disAssoc_log2RR2.pdf', log2rr2, 
       useDingbats = F, units = 'cm', width = 18, height = 20)
ggsave('./results/diseaseCooccur/disAssoc_log2RR2.png', log2rr2, 
       units = 'cm', width = 18, height = 20)

p = ggarrange(or_rr_correlation, nolabel_all, ncol = 2, nrow = 1, labels = c('a','b'), align = 'hv')
ggsave('./results/diseaseCooccur/figure1.pdf', p, useDingbats = F,
       units = 'cm', width = 18, height = 8)
ggsave('./results/diseaseCooccur/figure1.png', p,
       units = 'cm', width = 18, height = 8)
