source('./scripts/00-setup.R')
# library(clustermq)
# SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated.rds')
# disAgeMat <- SRdisease %>%
#   select(eid, Disease, Age) %>%
#   spread(Disease,Age,fill=NA) %>%
#   as.data.frame()
# rownames(disAgeMat) <- disAgeMat$eid
# disAgeMat$eid <- NULL
# disAgeMat <- as.matrix(disAgeMat)
# disAgeMat <- cbind(disAgeMat,top = apply(disAgeMat,1,function(x)min(x,na.rm=T)))
# saveRDS(disAgeMat,'./data/processed/traits_clean/disAgeMat_all.rds')
# rm(list = ls())
# f = function(age,num,kx){
#   library(tidyverse)
#   disAgeMat <- readRDS('./data/processed/traits_clean/disAgeMat_all.rds')
#   permDat <- readRDS('./data/processed/traits_clean/traitData_baseline_additions2.rds') %>%
#     select(eid, `Age when attended assessment centre`) %>%
#     set_names(c('eid','Age')) %>%
#     mutate(Age = floor(Age)) 
#   resx = sapply(1:kx, function(k){
#     subx = as.character(intersect(rownames(disAgeMat),sample(filter(permDat,Age>=age)$eid,num,replace = F)))
#     colSums(floor(disAgeMat[subx,])==age,na.rm=T)
#   })
#   return(resx)
# }
# permRes = Q(f,age = 0:64,n_jobs = 65,const = list(num = 50000,kx=100))
# names(permRes) = 0:64
# permRes <- reshape2::melt(permRes) %>%
#   set_names(c('Disease','PermNum','value','Age')) %>%
#   mutate(Age = as.numeric(Age))
# # system('mkdir -p ./data/processed/ageonset')
# saveRDS(permRes,'./data/processed/ageonset/permRes_50000_all.rds')
# disSet = readRDS('./data/processed/traits_clean/SRdiseaseSet.rds')
# permRes <- readRDS('./data/processed/ageonset/permRes_50000_all.rds') %>%
#   mutate(category = factor(unname(disTreecl[as.character(Disease)]))) %>%
#   mutate(Disease = fct_reorder(Disease, as.numeric(category)))
# pdf('./results/ageonset/disRates_50000_95_all.pdf', width = 12, height = 8)
# for(i in 1:40){
#   print(permRes %>%
#           mutate(value = value/5) %>%
#           group_by(Disease, category, Age) %>%
#           summarise(mean = mean(value), sd = sd(value), 
#                     fq = quantile(value,0.025), tq = quantile(value,0.975),
#                     median = median(value)) %>%
#           ggplot(aes(x = Age)) +
#           geom_smooth(aes(y=median), color = 'gray60', se = F, method = 'loess') +
#           geom_point(aes(y=median, color = category)) +
#           geom_segment(aes(y = fq, yend = tq, xend = Age), color = 'gray70') + 
#           facet_wrap_paginate(~Disease, scales = 'free', nrow = 3, ncol = 4, page = i) +
#           ylab('Number of participants (in 10,000)') +
#           scale_color_manual(values = discatcolors))
# }
# dev.off()
# permRes <- permRes %>%
#   filter(Disease %in% c('Top',disSet$Disease))
# saveRDS(permRes,'./data/processed/ageonset/permRes_50000.rds')
permRes <- readRDS('./data/processed/ageonset/permRes_50000.rds') %>%
  mutate(category = factor(unname(disTreecl[as.character(Disease)]))) %>%
  mutate(Disease = fct_reorder(Disease, as.numeric(category)))

pdf('./results/ageonset/disRates_50000_sd.pdf', width = 12, height = 8)
for(i in 1:10){
  print(permRes %>%
          mutate(value = value/5) %>%
          group_by(Disease, category, Age) %>%
          summarise(mean = mean(value), sd = sd(value), fq = quantile(value,0.25), tq = quantile(value,0.75),
                    median = median(value)) %>%
          ggplot(aes(x = Age)) +
          geom_point(aes(y=mean, color = category)) +
          geom_segment(aes(y = mean-sd, yend = mean+sd, xend = Age), color = 'gray70') + 
          facet_wrap_paginate(~Disease, scales = 'free', nrow = 3, ncol = 4, page = i) +
          ylab('Number of participants (in 10,000)') +
          scale_color_manual(values = discatcolors))
}
dev.off()

pdf('./results/ageonset/disRates_50000_95.pdf', width = 12, height = 8)
for(i in 1:10){
  print(permRes %>%
          mutate(value = value/5) %>%
          group_by(Disease, category, Age) %>%
          summarise(mean = mean(value), sd = sd(value), 
                    fq = quantile(value,0.025), tq = quantile(value,0.975),
                    median = median(value)) %>%
          ggplot(aes(x = Age)) +
          geom_smooth(aes(y=median), color = 'gray60', se = F, method = 'loess') +
          geom_point(aes(y=median, color = category)) +
          geom_segment(aes(y = fq, yend = tq, xend = Age), color = 'gray70') + 
          facet_wrap_paginate(~Disease, scales = 'free', nrow = 3, ncol = 4, page = i) +
          ylab('Number of participants (in 10,000)') +
          scale_color_manual(values = discatcolors))
}
dev.off()

medRate <- permRes %>%
  group_by(Disease, Age, category) %>%
  summarise(medRate = median(value))

medRate <- medRate %>%
  select(-category) %>%
  spread(Age,medRate) %>%
  as.data.frame()

rownames(medRate) <- medRate$Disease
medRate$Disease <- NULL
medRate <- as.matrix(medRate)
medRate <- medRate/rowSums(medRate)

pcx <- prcomp(medRate)
pcres <- cbind(pcx$x) %>%
  as.data.frame() %>%
  mutate(Disease = rownames(pcx$x)) %>%
  left_join(group_by(permRes,Disease, Age, category) %>%
              summarise(medRate = median(value)) %>%
              ungroup() %>% group_by(Disease, category) %>%
              summarise(medAge = Age[which.max(medRate)]))

pc_catcolor <- pcres %>%
  ggplot(aes(x = PC1, y = PC2, color = category, size = medAge)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(0.5,2)) +
  scale_color_manual(values = discatcolors) + 
  guides(size = F) +
  xlab(paste('PC1 (',round(100*summary(pcx)$imp[2,1],2),'%)', sep='')) +
  ylab(paste('PC2 (',round(100*summary(pcx)$imp[2,2],2),'%)', sep='')) +
  coord_fixed() +
  theme(legend.position = 'right')

ggsave('./results/ageonset/pca_medRate_catcolored.pdf',plot = pc_catcolor, units = 'cm', width = 16, height = 8, useDingbats = F)
ggsave('./results/ageonset/pca_medRate_catcolored.png',plot = pc_catcolor, units = 'cm', width = 16, height = 8)

pc_agecolor <- pcres %>%
  ggplot(aes(x = PC1, y = PC2, color = medAge)) +
  geom_point() +
  scale_color_viridis_c(direction = -1) +
  xlab(paste('PC1 (',round(100*summary(pcx)$imp[2,1],2),'%)', sep='')) +
  ylab(paste('PC2 (',round(100*summary(pcx)$imp[2,2],2),'%)', sep='')) +
  coord_fixed() +
  guides(color = guide_colorbar('median Age'))

ggsave('./results/ageonset/pca_medRate_agecolored.pdf',plot = pc_agecolor, units = 'cm', width = 16, height = 12, useDingbats = F)
ggsave('./results/ageonset/pca_medRate_agecolored.png',plot = pc_agecolor, units = 'cm', width = 16, height = 12)


pheatmap::pheatmap(medRate, cluster_cols = F, 
                   color = c('white',colorRampPalette(brewer.pal(8,'Oranges')[-c(1:2)])(20)), 
                   filename = './results/ageonset/medAgeRate.pdf',
                   cellwidth = 1, cellheight = 10, show_colnames = F)

pheatmap::pheatmap(medRate, cluster_cols = F, 
                   color = c('white',colorRampPalette(brewer.pal(8,'Oranges')[-c(1:2)])(20)), 
                   filename = './results/ageonset/medAgeRate.png',
                   cellwidth = 1, cellheight = 10, show_colnames = F)

library(dendextend)
library(ggdendro)
library(factoextra)
library(jpeg)
library(RFRlib)
hc <- hclust(dist(medRate))
dend <- as.dendrogram(hc)
owd=getwd()
agedist = medRate
for(hcut in 2:25){
  system(paste('mkdir -p /nfs/research1/thornton/ukbb_ageonset/results/ageonset/',hcut,'cluster',sep=''))
  setwd(paste('/nfs/research1/thornton/ukbb_ageonset/results/ageonset/',hcut,'cluster',sep=''))
  print(hcut)
  if(hcut<=1){
    discl=cutree(hc,h=hcut)
  } else if(hcut>1){
    discl=cutree(hc,k=hcut)
  }
  
  par(mar=c(20, 4, 4, 2) + 0.1)
  colsx=c(brewer.pal(8,'Dark2'),brewer.pal(8,'Set2'),brewer.pal(8,'Set1'))[1:(luniq(discl))]
  
  dend <- dend %>% 
    set("branches_k_color", k = luniq(discl),
        value = colsx) %>% 
    set("labels_colors", k = luniq(discl),
        value = colsx) %>% 
    set("branches_lwd", 2)
  agedf <- lapply(unique(discl[hc$order]),function(i){
    dismat <- agedist[names(which(discl==i)),]
    data.frame(mnsx=if(class(dismat)=='matrix'){colMeans(dismat)}else{dismat},
               sdx=if(class(dismat)=='matrix'){apply(dismat,2,sd)}else{rep(0,length(seq(0,64,by=1)))},
               cl=i,
               numDiseases=luniq(names(which(discl==i))),
               age=seq(0,64,1))
  })
  agedf=reshape2::melt(agedf,id.vars=colnames(agedf[[1]]))%>%
    rename(clorder=L1)
  names(colsx)=unique(discl[hc$order])
  agedistplots <- lapply(unique(agedf$clorder),function(i){
    filter(agedf,clorder==i) %>%
      ggplot(aes(x=age))+
      geom_segment(aes(xend=age,y=mnsx-sdx,yend=mnsx+sdx),
                   color='gray60', size = 1)+
      geom_line(aes(y=mnsx),color='black',size=2)+
      theme_minimal()+
      xlab('')+ylab('')+
      theme(axis.text.y = element_blank(),
            panel.grid.major.x = element_line(color='gray70',size = 1),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color=colsx[i],linetype='dashed',size=2,fill=NA))+
      scale_x_continuous(breaks = seq(0,64,by=10),limits = c(0,64))
  })
  midpoints=sapply(unique(discl[hc$order]),function(i){
    mean(unname(which(discl[hc$order]==i)))
  })
  sapply(1:length(agedistplots),function(i){
    ggsave(paste(i,'.jpeg',sep=''),
           agedistplots[[i]],device = 'jpeg',width = 3,height = 2)
  })
  myimgdat=data.frame(midpoints=midpoints,
                      xstart=midpoints-3,
                      xend=midpoints+3,
                      ystart=rep(-0.65,luniq(discl)),
                      yend=rep(-0.85,luniq(discl)),
                      filename=paste(1:luniq(discl),'.jpeg',sep=''))
  mynewimgdat=myimgdat[1,]
  for(i in 2:nrow(myimgdat)){
    if((myimgdat$xstart[i]<mynewimgdat$xend[i-1]) & (myimgdat$ystart[i]==mynewimgdat$ystart[i-1])){
      newadd=myimgdat[i,]
      newadd$ystart=newadd$ystart-0.2
      newadd$yend=newadd$yend-0.2
      mynewimgdat=rbind(mynewimgdat,newadd)
    } else{
      mynewimgdat=rbind(mynewimgdat,myimgdat[i,])
    }
  }
  myimgdat=mynewimgdat
  p=fviz_dend(hc, k = luniq(discl), # Cut in four groups
              cex = 1.1, # label size
              k_colors = colsx,
              color_labels_by_k = T,
              rect = TRUE, # Add rectangle around groups
              rect_border = colsx,
              rect_fill = F)+
    ylim(-1,0.4)+
    theme(axis.text = element_blank(),
    axis.ticks = element_blank())+
    ylab('')+
    ggtitle('')
  for(i in 1:nrow(myimgdat)){
    p = p +
      annotation_raster(readJPEG(as.character(myimgdat$filename[i])),
                        xmin = myimgdat$xstart[i],
                        xmax = myimgdat$xend[i],
                        ymin = myimgdat$ystart[i],
                        ymax = myimgdat$yend[i],interpolate = T)
  }
  ggsave(paste(hcut,'.pdf',sep=''),p,width = 30,height = 12, useDingbats = F)
  ggsave(paste(hcut,'.png',sep=''),p,width = 30,height = 12)
}
setwd(owd)

topdat <- readRDS('./data/processed/ageonset/permRes_50000_all.rds') %>%
  filter(Disease == 'top')

p_topdat <- topdat %>%
  mutate(value = value/5) %>%
  group_by(Disease, Age) %>%
  summarise(mean = mean(value), sd = sd(value),
            fq = quantile(value,0.025), tq = quantile(value,0.975),
            median = median(value)) %>%
  ggplot(aes(x = Age)) +
  geom_smooth(aes(y=median), color = 'darkslategray', se = F, method = 'loess') +
  geom_segment(aes(y = fq, yend = tq, xend = Age), color = 'darkslategray') +
  geom_point(aes(y=median), color = 'midnightblue') +
  ylab('Number of participants (in 10,000)') +
  ggtitle('Age of onset for the first\nself-reported disease') 

ggsave('./results/ageonset/topdat.pdf', p_topdat, units = 'cm', width = 8, height = 8, useDingbats = F)
ggsave('./results/ageonset/topdat.png', p_topdat, units = 'cm', width = 8, height = 8)

sumx <- permRes %>%
  mutate(value = value) %>%
  group_by(Disease,category) %>%
  summarise(Ages = list(rep(Age,value))) %>%
  mutate(med = median(unlist(Ages)),
         fq = quantile(unlist(Ages), 0.25),
         tq = quantile(unlist(Ages), 0.75)) %>%
  select(-Ages) %>%
  ungroup()

allsummarised <- sumx %>%
  mutate(Disease = fct_reorder(fct_reorder(Disease,med) %>%
                                 fct_rev(), as.numeric(category))) %>%
  ggplot(aes(x = Disease, color = category)) +
  geom_pointrange(aes(y = med, ymin = fq, ymax = tq),fatten = 0.5) +
  scale_color_manual(values = discatcolors) +
  # facet_grid(.~category, scales = 'free_x') +
  # guides(color = F) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_blank(), 
        legend.position = 'right') +
  xlab('') + ylab('Age of Onset')

ggsave('./results/ageonset/allsummarised.pdf', allsummarised, units = 'cm', width = 18, height = 8, useDingbats = F)
ggsave('./results/ageonset/allsummarised.png', allsummarised, units = 'cm', width = 18, height = 8)

catsum <- sumx %>%
  mutate(category = fct_reorder(category,-med)) %>%
  ggplot(aes(x = category, y = med)) + 
  geom_boxplot(aes(fill = category), outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.5, width = 0.2) +
  scale_fill_manual(values = discatcolors)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'right') +
  xlab('') + ylab('Median Age of Onset')

ggsave('./results/ageonset/catsummarised.pdf', catsum, units = 'cm', width = 18, height = 8, useDingbats = F)
ggsave('./results/ageonset/catsummarised.png', catsum, units = 'cm', width = 18, height = 8)

allsummarised2 <- sumx %>%
  mutate(Disease = fct_reorder(Disease, med)) %>%
  ggplot(aes(x = Disease, color = category)) +
  geom_pointrange(aes(y = med, ymin = fq, ymax = tq),fatten = 0.5) +
  scale_color_manual(values = discatcolors) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_blank(), 
        legend.position = 'right') +
  xlab('') + ylab('Age of Onset')

ggsave('./results/ageonset/allsummarised2.pdf', allsummarised2, units = 'cm', width = 18, height = 8, useDingbats = F)
ggsave('./results/ageonset/allsummarised2.png', allsummarised2, units = 'cm', width = 18, height = 8)


px = fviz_dend(hc, # Cut in four groups
            cex = 1.1, # label size
            label_cols = unname(discatcolors[unname(disTreecl[hc$labels[hc$order]])]),
            color_labels_by_k = F,
            rect = F, # Add rectangle around groups
            rect_fill = F)+
  ylim(-1,0.4)+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  ylab('')+
  ggtitle('')
ggsave('./results/ageonset/clusterColoredCat.pdf',px,width = 50,height = 20, useDingbats = F, units = 'cm')
ggsave('./results/ageonset/clusterColoredCat.png',px,width = 50,height = 20, units = 'cm')
