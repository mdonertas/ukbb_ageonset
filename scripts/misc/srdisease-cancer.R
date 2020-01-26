source('./scripts/00-setup.R')
library(ggnet)
traits <- readRDS('./data/processed/traits_clean/traitData_baseline_additions2.rds')
SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated_selected.rds')
prevDF = readRDS('./data/processed/traits_clean/SRdisease_prop_prevdf.rds')
disSet <- readRDS('./data/processed/traits_clean/SRdiseaseSet.rds')
disSet <- prevDF %>%
  filter(Disease %in% disSet$Disease)
SRcancer <- readRDS('./data/processed/traits_clean/SRcancer_baseline.rds')

canMat <- SRcancer %>%
  select(eid, Cancer) %>%
  unique() %>%
  mutate(value = 1) %>%
  spread(Cancer,value,fill=0) %>%
  as.data.frame()
rownames(canMat) <- canMat$eid
canMat$eid <- NULL
canMat <- as.matrix(canMat)
disMat <- readRDS('./data/processed/traits_clean/disMat.rds')
eidx <- rownames(disMat)
canMat2 <- t(sapply(eidx,function(x){
  if (x %in% rownames(canMat)){canMat[x,]
  }else{
    rep(0,ncol(canMat))
  }}))
canMat2 <- canMat2[eidx,]
colnames(canMat2) <- colnames(canMat)
canMat <- canMat2
# disMat <- t(sapply(eidx,function(x){
#   if (x %in% rownames(disMat)){disMat[x,]
#   }else{
#     rep(0,ncol(disMat))
#   }}))
disMat <- disMat[eidx,]

RRtable <- lapply(colnames(disMat), function(disAnm){
  disA <- disMat[,disAnm]
  xx <- as.tibble(t(sapply(colnames(canMat), function(disBnm){
    disB <- canMat[,disBnm]
    Nab <- sum((disA == 1) & (disB == 1))
    Nnab <- sum((disA == 0) & (disB == 1))
    Ta <- sum(disA == 1)
    Tna <- sum(disA == 0)
    Pexp <- Nab / Ta
    Pnexp <- Nnab / Tna
    RR <- Pexp / Pnexp
    cif <- sqrt((((Ta - Nab) / Nab) / Ta) + (((Tna - Nnab) / Nnab) / Tna))
    ll <- exp(log(RR) - (1.96 * cif))
    ul <- exp(log(RR) + (1.96 * cif))
    c(Nab=Nab,Nnab=Nnab,Ta=Ta,Tna=Tna,Pexp=Pexp,Pnexp=Pnexp,RR=RR,cif=cif,ll=ll,ul=ul)
  }))) %>%
    mutate(disB = colnames(canMat))
  return(xx)
})
names(RRtable) <- colnames(disMat)
RRtable <- reshape2::melt(RRtable, id.vars = colnames(RRtable[[1]])) %>%
  rename(disA=L1)
RRtable$sublevel <- sapply(1:nrow(RRtable),function(i)RRtable$disB[i]%in%subcomponent(disTree,RRtable$disA[i],'out')$name)
RRtable$uplevel <- sapply(1:nrow(RRtable),function(i)RRtable$disB[i]%in%subcomponent(disTree,RRtable$disA[i],'in')$name)
RRtable <- RRtable %>%
  select(disA, disB, everything())
RRtable %>%
  filter(sublevel==F & uplevel == F & ((RR > 1 & ll > 1) | (RR < 1 & ul < 1))) %>%
  filter(Nab>=10) %>%
  arrange(-RR) %>%
  # select(disB, disA, RR) %>%
  mutate(RR = log2(RR)) %>%
  filter(RR >= 1  | RR <= -1)%>%
  filter(abs(RR)>=2)

cormat = apply(canMat,2,function(x){
  apply(disMat,2,function(y){
    cor(x,y)
  })
})

rr2 <- RRtable %>%
  filter(uplevel == F & sublevel == F) %>%
  filter((RR > 1 & ll > 1) | (RR < 1 & ul < 1)) %>%
  mutate(RR = log2(RR)) %>%
  mutate(RR_rank = dense_rank(abs(RR)),
         RR_type = sign(RR)) %>%
  select(disA, disB, RR, RR_rank, RR_type)
sumx = rr2
cor2 <- reshape2::melt(cormat) %>%
  setNames(c('disA','disB','phi'))
sumx <- left_join(rr2, cor2)
xx <-  sumx %>%
  filter(sign(phi) == RR_type)
xx = xx %>% filter(disB%in%unique(filter(xx,abs(RR) >= 2)$disB))

xxmat <- select(xx, disA, disB, RR) %>%
  spread(disB, RR, fill = 0) %>%
  as.data.frame()

rownames(xxmat) = xxmat$disA
xxmat$disA = NULL
xxmat = as.matrix(xxmat)
hcx <- hclust(dist(xxmat))
xx <- xx %>%
  mutate(disA = factor(disA, levels = hcx$labels[hcx$order]))

discolDF <- tibble(disA = names(disTreecl), disB = names(disTreecl), disCat = unname(disTreecl)) %>%
  mutate(disCol = discatcolors[disCat]) %>%
  filter(disA %in% xx$disA & disB %in% xx$disB) %>%
  mutate(disA = factor(disA, levels = levels(xx$disA)),
         disB = factor(disB, levels = levels(xx$disB)))

finres <- ggplot(xx, aes(x = disA, y = disB)) +
  geom_point(aes(color = RR, size = abs(phi)), shape = 15) + 
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  scale_size_continuous(range = c(0.01,3),limits = c(0,0.3),breaks = c(0.1,0.2,0.3)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 6),
        legend.position = 'top') +
  xlab('') + ylab('') +
  guides(size = guide_legend(expression(~phi))) +
  coord_fixed() +
  geom_tile(data = discolDF, aes(fill = disCol)) +
  scale_fill_identity() +
  theme(panel.grid.major = element_line(size = 0.1, linetype = 'solid'))
finres
ggsave('./results/diseaseCooccur/disAssoc_log2RR2_disCancer.pdf', finres, 
       useDingbats = F, units = 'cm', width = 25, height = 30)
ggsave('./results/diseaseCooccur/disAssoc_log2RR2_disCancer.png', finres, 
      units = 'cm', width = 25, height = 30)

disMat2 = sapply(unique(disTreecl),function(x){
  disnm = intersect(colnames(disMat),names(which(disTreecl==x)))
  as.numeric(rowSums(disMat[,disnm])>0)
})
apply(disMat2,2,function(x)cor(x,as.numeric(rowSums(canMat)>0)))


ageonset = readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$clustering

disMat2 = sapply(unique(ageonset),function(x){
  disnm = intersect(colnames(disMat),names(which(ageonset==x)))
  as.numeric(rowSums(disMat[,disnm])>0)
})
apply(disMat2,2,function(x)cor(x,as.numeric(rowSums(canMat)>0)))
disB = as.numeric(rowSums(canMat))
t(apply(disMat2,2,function(disA){
  Nab <- sum((disA == 1) & (disB == 1))
  Nnab <- sum((disA == 0) & (disB == 1))
  Ta <- sum(disA == 1)
  Tna <- sum(disA == 0)
  Pexp <- Nab / Ta
  Pnexp <- Nnab / Tna
  RR <- Pexp / Pnexp
  cif <- sqrt((((Ta - Nab) / Nab) / Ta) + (((Tna - Nnab) / Nnab) / Tna))
  ll <- exp(log(RR) - (1.96 * cif))
  ul <- exp(log(RR) + (1.96 * cif))
  c(Nab=Nab,Nnab=Nnab,Ta=Ta,Tna=Tna,Pexp=Pexp,Pnexp=Pnexp,RR=RR,cif=cif,ll=ll,ul=ul)
}))
# Nab  Nnab     Ta    Tna       Pexp      Pnexp        RR        cif        ll        ul
# [1,] 19968  8728 226300 128560 0.08823685 0.06789048 1.2996941 0.01234733 1.2686180 1.3315314
# [2,] 20336  8360 240886 113974 0.08442168 0.07335006 1.1509421 0.01248462 1.1231205 1.1794529
# [3,] 12844 15852 165910 188950 0.07741547 0.08389521 0.9227638 0.01138513 0.9024005 0.9435866
# [4,]  2066 26630  23775 331085 0.08689800 0.08043252 1.0803840 0.02182882 1.0351352 1.1276107
