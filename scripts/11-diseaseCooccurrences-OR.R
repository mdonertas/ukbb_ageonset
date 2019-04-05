source('./scripts/00-setup.R')
library(ggnet)
traits <- readRDS('./data/processed/traits_clean/traitData_baseline_additions2.rds')
SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated_selected.rds')
prevDF = readRDS('./data/processed/traits_clean/SRdisease_prop_prevdf.rds')
disSet <- readRDS('./data/processed/traits_clean/SRdiseaseSet.rds')
disSet <- prevDF %>%
  filter(Disease %in% disSet$Disease)
disMat <- readRDS('./data/processed/traits_clean/disMat.rds')
disAgeMat <- readRDS('./data/processed/traits_clean/disAgeMat.rds')

RRtable <- lapply(colnames(disMat), function(disAnm){
  disA <- disMat[,disAnm]
  xx <- as.tibble(t(sapply(colnames(disMat), function(disBnm){
    disB <- disMat[,disBnm]
    Nab <- sum((disA == 1) & (disB == 1))
    Nnab <- sum((disA == 0) & (disB == 1))
    Nanb <- sum((disA == 1) & (disB == 0))
    Nnanb <- sum((disA == 0) & (disB == 0))
    matx <- matrix(c(Nab, Nnab, Nanb, Nnanb), ncol=2, nrow=2)
    fi <- fisher.test(matx)
    c(Nab = Nab, Nnab = Nnab, Nanb = Nanb, Nnanb = Nnanb, OR = unname(fi$estimate), 
      ll = fi$conf.int[1], ul = fi$conf.int[2], p = fi$p.value)
  }))) %>%
    mutate(disB = colnames(disMat))
  return(xx)
})
names(RRtable) <- colnames(disMat)
RRtable <- reshape2::melt(RRtable, id.vars = colnames(RRtable[[1]])) %>%
  rename(disA=L1)
RRtable$sublevel <- sapply(1:nrow(RRtable),function(i)RRtable$disB[i]%in%subcomponent(disTree,RRtable$disA[i],'out')$name)
RRtable$uplevel <- sapply(1:nrow(RRtable),function(i)RRtable$disB[i]%in%subcomponent(disTree,RRtable$disA[i],'in')$name)
RRtable <- RRtable %>%
  select(disA, disB, everything())
saveRDS(RRtable, './data/processed/diseaseCooccur/ORtable.rds')
RRtable %>%
  arrange(-OR) %>%
  write_tsv('./data/processed/diseaseCooccur/ORtable.tsv')

RRtable <- readRDS('./data/processed/diseaseCooccur/ORtable.rds')

RRmat <- RRtable %>%
  filter(sublevel==F & uplevel == F) %>%
  arrange(-OR) %>%
  select(disA, disB, OR) %>%
  spread(disB,OR,fill = 1) %>%
  as.data.frame()

rownames(RRmat) = RRmat$disA
RRmat$disA = NULL
RRmat = as.matrix(RRmat)
# diag(RRmat)=NA
library(pheatmap)
# system('mkdir -p ./results/diseaseCooccur')
mymat <- log2(RRmat)
# mymat[mymat==Inf]=5
annotdata <- data.frame(categories = unname(disTreecl[rownames(mymat)]),
                        medianAge = sapply(rownames(mymat),function(dis)median(disAgeMat[,dis],na.rm=T)))
rownames(annotdata) <- rownames(mymat)
annotcols <- list(categories = discatcolors[intersect(names(discatcolors),unique(annotdata$categories))], 
                  medianAge = brewer.pal(8,'Purples'))
pheatmap(mymat, cellwidth = 10, cellheight = 10, 
         filename = './results/diseaseCooccur/ORmat_log2_all.pdf', 
         color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49), 
         breaks = seq(-5,5,length.out = 50),
         cutree_rows = 10,
         annotation_row = annotdata, 
         annotation_col = annotdata,
         annotation_colors = annotcols,
         cutree_cols = 10, cluster_rows = T, cluster_cols = T)
pheatmap(mymat, cellwidth = 10, cellheight = 10, 
         filename = './results/diseaseCooccur/ORmat_log2_all.png', 
         color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49), 
         breaks = seq(-5,5,length.out = 50),
         cutree_rows = 10,
         annotation_row = annotdata, 
         annotation_col = annotdata,
         annotation_colors = annotcols,
         cutree_cols = 10, cluster_rows = T, cluster_cols = T)

library(ggnet)

RRnet <- RRtable %>%
  filter(sublevel==F & uplevel == F & ((OR > 1 & ll > 1) | (OR < 1 & ul < 1))) %>%
  arrange(-OR) %>%
  select(disB, disA, OR) %>%
  mutate(OR = log2(OR)) %>%
  filter(OR >= 1  | OR <= -1)%>%
  mutate(sizex = abs(OR) / 5) %>%
  graph_from_data_frame(directed = T)
V(RRnet)$category = unname(disTreecl[V(RRnet)$name])
E(RRnet)$signx = c('darkred','midnightblue')[1+(E(RRnet)$OR < 0)]
RRnetp <- ggnet2(RRnet, arrow.size = 2, arrow.gap = 0.02, edge.size = 'sizex', size = 3, edge.color = 'signx', 
                 edge.alpha = 0.5, node.color = 'category', palette = discatcolors) + 
  geom_text_repel(label= V(RRnet)$name, size=6/pntnorm, box.padding = 0) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 6))
ggsave(filename = './results/diseaseCooccur/ORnet_log2_sig_log2OR1.pdf',useDingbats = F, units = 'cm', width = 18,height = 15)
ggsave(filename = './results/diseaseCooccur/ORnet_log2_sig_log2OR1.png', units = 'cm', width = 18,height = 15)

RRnetp <- ggnet2(RRnet, arrow.size = 2, arrow.gap = 0.02, edge.size = 'sizex', size = 3, edge.color = 'signx', 
                 edge.alpha = 0.5, node.color = 'category', palette = discatcolors) + 
  # geom_text_repel(label= V(RRnet)$name, size=6/pntnorm, box.padding = 0) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 6))
ggsave(filename = './results/diseaseCooccur/ORnet_log2_sig_log2OR1_nolabel.pdf',useDingbats = F, units = 'cm', width = 18,height = 15)
ggsave(filename = './results/diseaseCooccur/ORnet_log2_sig_log2OR1_nolabel.png', units = 'cm', width = 18,height = 15)

RRnet <- RRtable %>%
  filter(sublevel==F & uplevel == F & ((OR > 1 & ll > 1) | (OR < 1 & ul < 1))) %>%
  arrange(-OR) %>%
  select(disB, disA, OR) %>%
  mutate(OR = log2(OR)) %>%
  filter(OR >= 2  | OR <= -2)%>%
  mutate(sizex = abs(OR) / 5) %>%
  graph_from_data_frame(directed = T)

V(RRnet)$category = unname(disTreecl[V(RRnet)$name])
E(RRnet)$signx = c('darkred','midnightblue')[1+(E(RRnet)$OR < 0)]
RRnetp <- ggnet2(RRnet, arrow.size = 2, arrow.gap = 0.02, edge.size = 'sizex', size = 3, edge.color = 'signx', 
                 edge.alpha = 0.75, node.color = 'category', palette = discatcolors) + 
  geom_text_repel(label= V(RRnet)$name, size=6/pntnorm, box.padding = 0) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 8) )
ggsave(filename = './results/diseaseCooccur/ORnet_log2_sig_log2OR2.pdf',useDingbats = F, units = 'cm', width = 18,height = 15 )
ggsave(filename = './results/diseaseCooccur/ORnet_log2_sig_log2OR2.png', units = 'cm', width = 18,height = 15 )

RRnetp <- ggnet2(RRnet, arrow.size = 2, arrow.gap = 0.02, edge.size = 'sizex', size = 3, edge.color = 'signx', 
                 edge.alpha = 0.75, node.color = 'category', palette = discatcolors) + 
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 8) )
ggsave(filename = './results/diseaseCooccur/ORnet_log2_sig_log2OR2_nolabel.pdf',useDingbats = F, units = 'cm', width = 18,height = 15 )
ggsave(filename = './results/diseaseCooccur/ORnet_log2_sig_log2OR2_nolabel.png', units = 'cm', width = 18,height = 15 )

RRnet <- RRtable %>%
  filter(sublevel==F & uplevel == F & ((OR > 1 & ll > 1) | (OR < 1 & ul < 1))) %>%
  arrange(-OR) %>%
  select(disB, disA, OR) %>%
  mutate(OR = log2(OR)) %>%
  filter(OR >= 3  | OR <= -3)%>%
  mutate(sizex = abs(OR) / 5) %>%
  graph_from_data_frame(directed = T)
V(RRnet)$category = unname(disTreecl[V(RRnet)$name])
E(RRnet)$signx = c('darkred','midnightblue')[1+(E(RRnet)$OR < 0)]
RRnetp <- ggnet2(RRnet, arrow.size = 2, arrow.gap = 0.02, edge.size = 'sizex', size = 3, edge.color = 'signx', 
                 edge.alpha = 0.75, node.color = 'category', palette = discatcolors) + 
  geom_text_repel(label= V(RRnet)$name, size=6/pntnorm, box.padding = 0) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 8))
ggsave(filename = './results/diseaseCooccur/ORnet_log2_sig_log2OR3.pdf',useDingbats = F, units = 'cm', width = 18,height = 12)
ggsave(filename = './results/diseaseCooccur/ORnet_log2_sig_log2OR3.png', units = 'cm', width = 18,height = 12)

rm(list = ls())

