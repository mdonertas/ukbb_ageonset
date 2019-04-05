source('./scripts/00-setup.R')
library(ggnet)
library(pheatmap)
traits <- readRDS('./data/processed/traits_clean/traitData_baseline_additions2.rds')
SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated_selected.rds')
prevDF = readRDS('./data/processed/traits_clean/SRdisease_prop_prevdf.rds')
disSet <- readRDS('./data/processed/traits_clean/SRdiseaseSet.rds')
disSet <- prevDF %>%
  filter(Disease %in% disSet$Disease)
disMat <- readRDS('./data/processed/traits_clean/disMat.rds')
sex = traits %>% 
  select(eid,Sex)

# disMatm <- disMat[rownames(disMat)%in%unique(filter(sex, Sex =='Male')$eid),]
# disMatf <- disMat[rownames(disMat)%in%unique(filter(sex, Sex =='Female')$eid),]
# 
# RRtable <- lapply(colnames(disMat), function(disAnm){
#   disAf <- disMatf[,disAnm]
#   disAm <- disMatm[,disAnm]
#   xx <- as.tibble(t(sapply(colnames(disMat), function(disBnm){
#     disBf <- disMatf[,disBnm]
#     disBm <- disMatm[,disBnm]
#     f_B <- sum(disBf)
#     f_AB <- sum(disBf == 1 & disAf ==1)
#     yk <- (f_B / f_AB) + 1
#     m_B <- sum(disBm)
#     m_AB <- sum(disBm == 1 & disAm ==1)
#     as <- (m_B / m_AB) + 1
#     SR <- log2(yk / as)
#     c(f_B = f_B, f_AB = f_AB, m_B = m_B, m_AB = m_AB, SR = SR)
#   }))) %>%
#     mutate(disB = colnames(disMat))
#   return(xx)
# })
# 
# names(RRtable) <- colnames(disMat)
# RRtable <- reshape2::melt(RRtable, id.vars = colnames(RRtable[[1]])) %>%
#   rename(disA=L1)
# RRtable$sublevel <- sapply(1:nrow(RRtable),function(i)RRtable$disB[i]%in%subcomponent(disTree,RRtable$disA[i],'out')$name)
# RRtable$uplevel <- sapply(1:nrow(RRtable),function(i)RRtable$disB[i]%in%subcomponent(disTree,RRtable$disA[i],'in')$name)
# RRtable <- RRtable %>%
#   select(disA, disB, everything())
# saveRDS(RRtable, './data/processed/diseaseCooccur/SRtable.rds')
# RRtable %>%
#   arrange(-SR) %>%
#   write_tsv('./data/processed/diseaseCooccur/SRtable.tsv')
RRtable=readRDS('./data/processed/diseaseCooccur/SRtable.rds')
SRmat = RRtable %>%
  filter( sublevel == F & uplevel == F) %>%
  filter( abs(SR)>=2) %>%
  select(disA, disB, SR) %>%
  spread(disB, SR, fill = 0) %>%
  as.data.frame()

rownames(SRmat) = SRmat$disA
SRmat$disA = NULL
SRmat = as.matrix(SRmat)

annotdata <- data.frame(categories = unname(disTreecl[rownames(SRmat)]),
                        medianAge = sapply(rownames(SRmat),function(dis)median(disMat[,dis],na.rm=T)))
rownames(annotdata) <- rownames(SRmat)
annotcols <- list(categories = discatcolors[intersect(names(discatcolors),unique(annotdata$categories))], 
                  medianAge = brewer.pal(8,'Purples'))
if(min(SRmat) == -Inf){SRmat[SRmat==-Inf]= unique(sort(SRmat,dec=F))[2]-1}
if(max(SRmat) == Inf){SRmat[SRmat==Inf]= unique(sort(SRmat,dec=T))[2]+1}
pheatmap(SRmat, cellwidth = 10, cellheight = 10, 
         filename = './results/diseaseCooccur/SRmat_log2_log2SR2.pdf',
         color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49),
         breaks = seq(-max(abs(SRmat)),max(abs(SRmat)),length.out = 50),
         cutree_rows = 10,
         # annotation_row = annotdata,
         # annotation_col = annotdata,
         # annotation_colors = annotcols,
         cutree_cols = 10,
         cluster_rows = T, cluster_cols = T)

pheatmap(SRmat, cellwidth = 10, cellheight = 10, 
         filename = './results/diseaseCooccur/SRmat_log2_log2SR2.png',
         color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49), 
         breaks = seq(-max(abs(SRmat)),max(abs(SRmat)),length.out = 50),
         cutree_rows = 10,
         # annotation_row = annotdata,
         # annotation_col = annotdata,
         # annotation_colors = annotcols,
         cutree_cols = 10,
         cluster_rows = T, cluster_cols = T)
