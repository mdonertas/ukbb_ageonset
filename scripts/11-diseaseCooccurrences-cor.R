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
disAgeMat <- readRDS('./data/processed/traits_clean/disAgeMat.rds')

co = cor(disMat)
saveRDS(co, './data/processed/diseaseCooccur/disCorrelations.rds')

mymat = co
diag(mymat) = 0

relmat = sapply(rownames(mymat),function(x){
  sapply(rownames(mymat),function(y){
    x%in%subcomponent(disTree, y,'out')$name | x%in%subcomponent(disTree, y,'in')$name
  })
})

mymat[relmat == T] = 0

annotdata <- data.frame(categories = unname(disTreecl[rownames(mymat)]),
                        medianAge = sapply(rownames(mymat),function(dis)median(disAgeMat[,dis],na.rm=T)))
rownames(annotdata) <- rownames(mymat)
annotcols <- list(categories = discatcolors[intersect(names(discatcolors),unique(annotdata$categories))], 
                  medianAge = brewer.pal(8,'Purples'))
pheatmap(mymat, cellwidth = 10, cellheight = 10, 
         filename = './results/diseaseCooccur/cormat_all.pdf', 
         color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49), 
         breaks = seq(-0.35,0.35,length.out = 50),
         cutree_rows = 10,
         annotation_row = annotdata, 
         annotation_col = annotdata,
         annotation_colors = annotcols,
         cutree_cols = 10, cluster_rows = T, cluster_cols = T)

pheatmap(mymat, cellwidth = 10, cellheight = 10, 
         filename = './results/diseaseCooccur/cormat_all.png', 
         color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49), 
         breaks = seq(-0.35,0.35,length.out = 50),
         cutree_rows = 10,
         annotation_row = annotdata, 
         annotation_col = annotdata,
         annotation_colors = annotcols,
         cutree_cols = 10, cluster_rows = T, cluster_cols = T)

mymat[abs(mymat)<0.1] = 0
mymat = mymat[rowMeans(mymat==0) !=1,rowMeans(mymat==0) !=1] 
annotdata <- data.frame(categories = unname(disTreecl[rownames(mymat)]),
                        medianAge = sapply(rownames(mymat),function(dis)median(disAgeMat[,dis],na.rm=T)))
rownames(annotdata) <- rownames(mymat)
annotcols <- list(categories = discatcolors[intersect(names(discatcolors),unique(annotdata$categories))], 
                  medianAge = brewer.pal(8,'Purples'))
pheatmap(mymat, cellwidth = 10, cellheight = 10, 
         filename = './results/diseaseCooccur/cormat_h01.pdf', 
         color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49), 
         breaks = seq(-0.35,0.35,length.out = 50),
         cutree_rows = 2,
         annotation_row = annotdata, 
         annotation_col = annotdata,
         annotation_colors = annotcols,
         cutree_cols = 2, cluster_rows = T, cluster_cols = T)

pheatmap(mymat, cellwidth = 10, cellheight = 10, 
         filename = './results/diseaseCooccur/cormat_h01.png', 
         color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49), 
         breaks = seq(-0.35,0.35,length.out = 50),
         cutree_rows = 2,
         annotation_row = annotdata, 
         annotation_col = annotdata,
         annotation_colors = annotcols,
         cutree_cols = 2, cluster_rows = T, cluster_cols = T)
