source('./scripts/00-setup.R')
library(ggnet)
traits <- readRDS('./data/processed/traits_clean/traitData_baseline_additions2.rds')
SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated_selected.rds')
prevDF = readRDS('./data/processed/traits_clean/SRdisease_prop_prevdf.rds')
disSet <- readRDS('./data/processed/traits_clean/SRdiseaseSet.rds')
disSet <- prevDF %>%
  filter(Disease %in% disSet$Disease)

# disMat <- SRdisease %>%
#   select(eid, Disease) %>%
#   mutate(value = 1) %>%
#   spread(Disease,value,fill=0) %>%
#   as.data.frame()
# 
# rownames(disMat) <- disMat$eid
# disMat$eid <- NULL
# disMat <- as.matrix(disMat)
# 
# disAgeMat <- SRdisease %>%
#   select(eid, Disease, Age) %>%
#   spread(Disease,Age,fill=NA) %>%
#   as.data.frame()
# 
# rownames(disAgeMat) <- disAgeMat$eid
# disAgeMat$eid <- NULL
# disAgeMat <- as.matrix(disAgeMat)
# 
# saveRDS(disMat,'./data/processed/traits_clean/disMat.rds')
# saveRDS(disAgeMat,'./data/processed/traits_clean/disAgeMat.rds')

disMat <- readRDS('./data/processed/traits_clean/disMat.rds')
disAgeMat <- readRDS('./data/processed/traits_clean/disAgeMat.rds')

# RRtable <- lapply(colnames(disMat), function(disAnm){
#   disA <- disMat[,disAnm]
#   xx <- as.tibble(t(sapply(colnames(disMat), function(disBnm){
#     disB <- disMat[,disBnm]
#     Nab <- sum((disA == 1) & (disB == 1))
#     Nnab <- sum((disA == 0) & (disB == 1))
#     Ta <- sum(disA == 1)
#     Tna <- sum(disA == 0)
#     Pexp <- Nab / Ta
#     Pnexp <- Nnab / Tna
#     RR <- Pexp / Pnexp
#     cif <- sqrt((((Ta - Nab) / Nab) / Ta) + (((Tna - Nnab) / Nnab) / Tna))
#     ll <- exp(log(RR) - (1.96 * cif))
#     ul <- exp(log(RR) + (1.96 * cif))
#     c(Nab=Nab,Nnab=Nnab,Ta=Ta,Tna=Tna,Pexp=Pexp,Pnexp=Pnexp,RR=RR,cif=cif,ll=ll,ul=ul)
#   }))) %>%
#     mutate(disB = colnames(disMat))
#   return(xx)
# })
# names(RRtable) <- colnames(disMat)
# RRtable <- reshape2::melt(RRtable, id.vars = colnames(RRtable[[1]])) %>%
#   rename(disA=L1)
# RRtable$sublevel <- sapply(1:nrow(RRtable),function(i)RRtable$disB[i]%in%subcomponent(disTree,RRtable$disA[i],'out')$name)
# RRtable$uplevel <- sapply(1:nrow(RRtable),function(i)RRtable$disB[i]%in%subcomponent(disTree,RRtable$disA[i],'in')$name)
# system('mkdir -p ./data/processed/diseaseCooccur')
# RRtable <- RRtable %>%
#   select(disA, disB, everything())
# saveRDS(RRtable, './data/processed/diseaseCooccur/RRtable.rds')
# RRtable %>%
#   arrange(-RR) %>%
#   write_tsv('./data/processed/diseaseCooccur/RRtable.tsv')

RRtable <- readRDS('./data/processed/diseaseCooccur/RRtable.rds')

RRmat <- RRtable %>%
  filter(sublevel==F & uplevel == F) %>%
  arrange(-RR) %>%
  select(disA, disB, RR) %>%
  spread(disB,RR,fill = 1) %>%
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
         filename = './results/diseaseCooccur/RRmat_log2_all.pdf', 
         color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49), 
         breaks = seq(-5,5,length.out = 50),
         cutree_rows = 10,
         annotation_row = annotdata, 
         annotation_col = annotdata,
         annotation_colors = annotcols,
         cutree_cols = 10, cluster_rows = T, cluster_cols = T)
pheatmap(mymat, cellwidth = 10, cellheight = 10, 
         filename = './results/diseaseCooccur/RRmat_log2_all.png', 
         color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49), 
         breaks = seq(-5,5,length.out = 50),
         cutree_rows = 10,
         annotation_row = annotdata, 
         annotation_col = annotdata,
         annotation_colors = annotcols,
         cutree_cols = 10, cluster_rows = T, cluster_cols = T)

library(ggnet)

RRnet <- RRtable %>%
  filter(sublevel==F & uplevel == F & ((RR > 1 & ll > 1) | (RR < 1 & ul < 1))) %>%
  arrange(-RR) %>%
  select(disB, disA, RR) %>%
  mutate(RR = log2(RR)) %>%
  filter(RR >= 1  | RR <= -1)%>%
  mutate(sizex = abs(RR) / 5) %>%
  graph_from_data_frame(directed = T)
V(RRnet)$category = unname(disTreecl[V(RRnet)$name])
E(RRnet)$signx = c('darkred','midnightblue')[1+(E(RRnet)$RR < 0)]
RRnetp <- ggnet2(RRnet, arrow.size = 2, arrow.gap = 0.02, edge.size = 'sizex', size = 3, edge.color = 'signx', 
                 edge.alpha = 0.5, node.color = 'category', palette = discatcolors) + 
  geom_text_repel(label= V(RRnet)$name, size=6/pntnorm, box.padding = 0) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 6))
ggsave(filename = './results/diseaseCooccur/RRnet_log2_sig_log2RR1.pdf',useDingbats = F, units = 'cm', width = 18,height = 15)
ggsave(filename = './results/diseaseCooccur/RRnet_log2_sig_log2RR1.png', units = 'cm', width = 18,height = 15)

RRnetp <- ggnet2(RRnet, arrow.size = 2, arrow.gap = 0.02, edge.size = 'sizex', size = 3, edge.color = 'signx', 
                 edge.alpha = 0.5, node.color = 'category', palette = discatcolors) + 
  # geom_text_repel(label= V(RRnet)$name, size=6/pntnorm, box.padding = 0) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 6))
ggsave(filename = './results/diseaseCooccur/RRnet_log2_sig_log2RR1_nolabel.pdf',useDingbats = F, units = 'cm', width = 18,height = 15)
ggsave(filename = './results/diseaseCooccur/RRnet_log2_sig_log2RR1_nolabel.png', units = 'cm', width = 18,height = 15)

RRnet <- RRtable %>%
  filter(sublevel==F & uplevel == F & ((RR > 1 & ll > 1) | (RR < 1 & ul < 1))) %>%
  arrange(-RR) %>%
  select(disB, disA, RR) %>%
  mutate(RR = log2(RR)) %>%
  filter(RR >= 2  | RR <= -2)%>%
  mutate(sizex = abs(RR) / 5) %>%
  graph_from_data_frame(directed = T)

V(RRnet)$category = unname(disTreecl[V(RRnet)$name])
E(RRnet)$signx = c('darkred','midnightblue')[1+(E(RRnet)$RR < 0)]
RRnetp <- ggnet2(RRnet, arrow.size = 2, arrow.gap = 0.02, edge.size = 'sizex', size = 3, edge.color = 'signx', 
                 edge.alpha = 0.75, node.color = 'category', palette = discatcolors) + 
  geom_text_repel(label= V(RRnet)$name, size=6/pntnorm, box.padding = 0) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 8) )
ggsave(filename = './results/diseaseCooccur/RRnet_log2_sig_log2RR2.pdf',useDingbats = F, units = 'cm', width = 18,height = 15 )
ggsave(filename = './results/diseaseCooccur/RRnet_log2_sig_log2RR2.png', units = 'cm', width = 18,height = 15 )

RRnetp <- ggnet2(RRnet, arrow.size = 2, arrow.gap = 0.02, edge.size = 'sizex', size = 3, edge.color = 'signx', 
                 edge.alpha = 0.75, node.color = 'category', palette = discatcolors) + 
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 8) )
ggsave(filename = './results/diseaseCooccur/RRnet_log2_sig_log2RR2_nolabel.pdf',useDingbats = F, units = 'cm', width = 18,height = 15 )
ggsave(filename = './results/diseaseCooccur/RRnet_log2_sig_log2RR2_nolabel.png', units = 'cm', width = 18,height = 15 )

RRnet <- RRtable %>%
  filter(sublevel==F & uplevel == F & ((RR > 1 & ll > 1) | (RR < 1 & ul < 1))) %>%
  arrange(-RR) %>%
  select(disB, disA, RR) %>%
  mutate(RR = log2(RR)) %>%
  filter(RR >= 3  | RR <= -3)%>%
  mutate(sizex = abs(RR) / 5) %>%
  graph_from_data_frame(directed = T)
V(RRnet)$category = unname(disTreecl[V(RRnet)$name])
E(RRnet)$signx = c('darkred','midnightblue')[1+(E(RRnet)$RR < 0)]
RRnetp <- ggnet2(RRnet, arrow.size = 2, arrow.gap = 0.02, edge.size = 'sizex', size = 3, edge.color = 'signx', 
                 edge.alpha = 0.75, node.color = 'category', palette = discatcolors) + 
  geom_text_repel(label= V(RRnet)$name, size=6/pntnorm, box.padding = 0) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 8))
ggsave(filename = './results/diseaseCooccur/RRnet_log2_sig_log2RR3.pdf',useDingbats = F, units = 'cm', width = 18,height = 12)
ggsave(filename = './results/diseaseCooccur/RRnet_log2_sig_log2RR3.png', units = 'cm', width = 18,height = 12)

rm(list = ls())

source('./scripts/00-setup.R')
library(ggnet)
traits <- readRDS('./data/processed/traits_clean/traitData_baseline_additions2.rds')
SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated_selected.rds')
prevDF = readRDS('./data/processed/traits_clean/SRdisease_prop_prevdf.rds')
disSet <- readRDS('./data/processed/traits_clean/SRdiseaseSet.rds')
disSet <- prevDF %>%
  filter(Disease %in% disSet$Disease)
disAgeMat <- readRDS('./data/processed/traits_clean/disAgeMat.rds')
disMatf <- matrix(cut((disAgeMat),breaks = seq(0,75,by = 10)),nrow = nrow(disAgeMat), ncol = ncol(disAgeMat),byrow = F,dimnames = list(rownames(disAgeMat),colnames(disAgeMat)))

# RRtable <- lapply(sort(na.omit(unique(c(disMatf)))),function(agerng){
#   disMat <- disMatf
#   disMat[disMatf == agerng] = 1
#   disMat[disMatf != agerng] = 0
#   disMat[is.na(disMatf)] = 0
#   disMat = disMat
#   RRtable <- lapply(colnames(disMat), function(disAnm){
#     disA <- as.numeric(disMat[,disAnm])
#     xx <- as.tibble(t(sapply(colnames(disMat), function(disBnm){
#       disB <- as.numeric(disMat[,disBnm])
#       Nab <- sum((disA == 1) & (disB == 1))
#       Nnab <- sum((disA == 0) & (disB == 1))
#       Ta <- sum(disA == 1)
#       Tna <- sum(disA == 0)
#       Pexp <- Nab / Ta
#       Pnexp <- Nnab / Tna
#       RR <- Pexp / Pnexp
#       cif <- sqrt((((Ta - Nab) / Nab) / Ta) + (((Tna - Nnab) / Nnab) / Tna))
#       ll <- exp(log(RR) - (1.96 * cif))
#       ul <- exp(log(RR) + (1.96 * cif))
#       c(Nab=Nab,Nnab=Nnab,Ta=Ta,Tna=Tna,Pexp=Pexp,Pnexp=Pnexp,RR=RR,cif=cif,ll=ll,ul=ul)
#     }))) %>%
#       mutate(disB = colnames(disMat))
#     return(xx)
#   })
#   names(RRtable) <- colnames(disMat)
#   RRtable <- reshape2::melt(RRtable, id.vars = colnames(RRtable[[1]])) %>%
#     rename(disA=L1)
#   RRtable$sublevel <- sapply(1:nrow(RRtable),function(i)RRtable$disB[i]%in%subcomponent(disTree,RRtable$disA[i],'out')$name)
#   RRtable$uplevel <- sapply(1:nrow(RRtable),function(i)RRtable$disB[i]%in%subcomponent(disTree,RRtable$disA[i],'in')$name)
#   RRtable$age <- agerng
#   return(RRtable)
# })
# RRtable <- reshape2::melt(RRtable,id.vars = colnames(RRtable[[1]]))
# RRtable <- RRtable %>%
#   select(-L1) %>%
#   select(disA, disB, everything())
# saveRDS(RRtable, './data/processed/diseaseCooccur/RRtable_agebins.rds')
# 
# RRtable %>%
#   arrange(-RR) %>%
#   write_tsv('./data/processed/diseaseCooccur/RRtable_agebins.tsv')
# 
# head(RRtable)

RRtable <- readRDS('./data/processed/diseaseCooccur/RRtable_agebins.rds')

for(agebn in unique(RRtable$age)){
  print(agebn)
  RRmat <- RRtable %>%
    filter(sublevel==F & uplevel == F & ((RR>1 & ll>1)|(RR<1 & ul<1))) %>%
    filter(age == agebn) %>%
    mutate(RR = log2(RR)) %>%
    filter(RR >= 2  | RR <= -2)%>%
    select(disA, disB, RR) %>%
    spread(disB,RR,fill = 0) %>%
    as.data.frame()
  
  rownames(RRmat) = RRmat$disA
  RRmat$disA = NULL
  RRmat = as.matrix(RRmat)
  # diag(RRmat)=NA
  library(pheatmap)
  # # system('mkdir -p ./results/diseaseCooccur')
  mymat <- RRmat
  if(min(mymat) == -Inf){mymat[mymat==-Inf]= unique(sort(mymat,dec=F))[2]-1}
  if(max(mymat) == Inf){mymat[mymat==Inf]= unique(sort(mymat,dec=T))[2]+1}
  annotdata <- data.frame(categories = unname(disTreecl[rownames(mymat)]),
                          medianAge = sapply(rownames(mymat),function(dis)median(disAgeMat[,dis],na.rm=T)))
  rownames(annotdata) <- rownames(mymat)
  annotcols <- list(categories = discatcolors[intersect(names(discatcolors),unique(annotdata$categories))],
                    medianAge = brewer.pal(8,'Purples'))
  pheatmap(mymat, cellwidth = 10, cellheight = 10,
           filename = paste('./results/diseaseCooccur/RRmat_log2_sig_log2RR2_',agebn,'.pdf',sep=''),
           color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49),
           breaks = seq(-5,5,length.out = 50),
           cutree_rows = 10,
           annotation_row = annotdata,
           annotation_col = annotdata,
           annotation_colors = annotcols,
           cutree_cols = 10, cluster_rows = T, cluster_cols = T)
  pheatmap(mymat, cellwidth = 10, cellheight = 10,
           filename = paste('./results/diseaseCooccur/RRmat_log2_sig_log2RR2_',agebn,'.png',sep=''),
           color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49),
           breaks = seq(-5,5,length.out = 50),
           cutree_rows = 10,
           annotation_row = annotdata,
           annotation_col = annotdata,
           annotation_colors = annotcols,
           cutree_cols = 10, cluster_rows = T, cluster_cols = T)
}
