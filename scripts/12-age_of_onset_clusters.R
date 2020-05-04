source('./scripts/00-setup.R')
library(cluster)
library(TSclust)

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
ageonsetcl = readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$cluster
permRes <- readRDS('./data/processed/ageonset/permRes_50000.rds') %>%
  mutate(cluster = as.character(ageonsetcl[as.character(Disease)])) %>%
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

i=1
as.character(unique(permRes$category)[i])
# [1] "gastrointestinal/abdominal"
permRes %>%
  filter(category == as.character(unique(permRes$category)[i])) %>%
  mutate(value = value/5) %>%
  group_by(Disease, category, Age,cluster) %>%
  summarise(mean = mean(value), sd = sd(value), 
            fq = quantile(value,0.025), tq = quantile(value,0.975),
            median = median(value)) %>% 
  ggplot(aes(x = Age)) +
  geom_smooth(aes(y=median), color = 'gray60', se = F, method = 'loess') +
  geom_point(aes(y=median, color = cluster)) +
  geom_segment(aes(y = fq, yend = tq, xend = Age), color = 'gray70', size = 0.3) + 
  facet_wrap(~Disease, scales = 'free', ncol = 4) +
  ylab('No of cases diagnosed / 10,000 people at a given age') +
  scale_color_manual(values = ageonsetcolors) +
  guides(color = guide_legend('Age-of-onset cluster')) +
  ggtitle(as.character(unique(permRes$category)[i]))
ggsave(paste('./results/ageonset/',gsub('/','_',as.character(unique(permRes$category)[i])),'.pdf',sep=''),units = 'cm', width = 24, height = 27, useDingbats = F)

i=2
as.character(unique(permRes$category)[i])
# [1] "immunological/systemic disorders"
permRes %>%
  filter(category == as.character(unique(permRes$category)[i])) %>%
  mutate(value = value/5) %>%
  group_by(Disease, category, Age,cluster) %>%
  summarise(mean = mean(value), sd = sd(value), 
            fq = quantile(value,0.025), tq = quantile(value,0.975),
            median = median(value)) %>% 
  ggplot(aes(x = Age)) +
  geom_smooth(aes(y=median), color = 'gray60', se = F, method = 'loess') +
  geom_point(aes(y=median, color = cluster)) +
  geom_segment(aes(y = fq, yend = tq, xend = Age), color = 'gray70', size = 0.3) + 
  facet_wrap(~Disease, scales = 'free', ncol = 4) +
  ylab('No of cases diagnosed / 10,000 people at a given age') +
  scale_color_manual(values = ageonsetcolors) +
  guides(color = guide_legend('Age-of-onset cluster')) +
  ggtitle(as.character(unique(permRes$category)[i]))
ggsave(paste('./results/ageonset/',gsub('/','_',as.character(unique(permRes$category)[i])),'.pdf',sep=''),units = 'cm', width = 24, height = 12, useDingbats = F)

i=3
as.character(unique(permRes$category)[i])
# [1] "haematology/dermatology"
permRes %>%
  filter(category == as.character(unique(permRes$category)[i])) %>%
  mutate(value = value/5) %>%
  group_by(Disease, category, Age,cluster) %>%
  summarise(mean = mean(value), sd = sd(value), 
            fq = quantile(value,0.025), tq = quantile(value,0.975),
            median = median(value)) %>% 
  ggplot(aes(x = Age)) +
  geom_smooth(aes(y=median), color = 'gray60', se = F, method = 'loess') +
  geom_point(aes(y=median, color = cluster)) +
  geom_segment(aes(y = fq, yend = tq, xend = Age), color = 'gray70', size = 0.3) + 
  facet_wrap(~Disease, scales = 'free', ncol = 4) +
  ylab('No of cases diagnosed / 10,000 people at a given age') +
  scale_color_manual(values = ageonsetcolors) +
  guides(color = guide_legend('Age-of-onset cluster')) +
  ggtitle(as.character(unique(permRes$category)[i]))
ggsave(paste('./results/ageonset/',gsub('/','_',as.character(unique(permRes$category)[i])),'.pdf',sep=''),units = 'cm', width = 24, height = 12, useDingbats = F)

i=4
as.character(unique(permRes$category)[i])
# "cardiovascular"
permRes %>%
  filter(category == as.character(unique(permRes$category)[i])) %>%
  mutate(value = value/5) %>%
  group_by(Disease, category, Age,cluster) %>%
  summarise(mean = mean(value), sd = sd(value), 
            fq = quantile(value,0.025), tq = quantile(value,0.975),
            median = median(value)) %>% 
  ggplot(aes(x = Age)) +
  geom_smooth(aes(y=median), color = 'gray60', se = F, method = 'loess') +
  geom_point(aes(y=median, color = cluster)) +
  geom_segment(aes(y = fq, yend = tq, xend = Age), color = 'gray70', size = 0.3) + 
  facet_wrap(~Disease, scales = 'free', ncol = 4) +
  ylab('No of cases diagnosed / 10,000 people at a given age') +
  scale_color_manual(values = ageonsetcolors) +
  guides(color = guide_legend('Age-of-onset cluster')) +
  ggtitle(as.character(unique(permRes$category)[i]))
ggsave(paste('./results/ageonset/',gsub('/','_',as.character(unique(permRes$category)[i])),'.pdf',sep=''),units = 'cm', width = 24, height = 22, useDingbats = F)

i=5
as.character(unique(permRes$category)[i])
# [1] "neurology/eye/psychiatry"
permRes %>%
  filter(category == as.character(unique(permRes$category)[i])) %>%
  mutate(value = value/5) %>%
  group_by(Disease, category, Age,cluster) %>%
  summarise(mean = mean(value), sd = sd(value), 
            fq = quantile(value,0.025), tq = quantile(value,0.975),
            median = median(value)) %>% 
  ggplot(aes(x = Age)) +
  geom_smooth(aes(y=median), color = 'gray60', se = F, method = 'loess') +
  geom_point(aes(y=median, color = cluster)) +
  geom_segment(aes(y = fq, yend = tq, xend = Age), color = 'gray70', size = 0.3) + 
  facet_wrap(~Disease, scales = 'free', ncol = 4) +
  ylab('No of cases diagnosed / 10,000 people at a given age') +
  scale_color_manual(values = ageonsetcolors) +
  guides(color = guide_legend('Age-of-onset cluster')) +
  ggtitle(as.character(unique(permRes$category)[i]))
ggsave(paste('./results/ageonset/',gsub('/','_',as.character(unique(permRes$category)[i])),'.pdf',sep=''),units = 'cm', width = 24, height = 27, useDingbats = F)

i=6
as.character(unique(permRes$category)[i])
# [1]  "musculoskeletal/trauma"
permRes %>%
  filter(category == as.character(unique(permRes$category)[i])) %>%
  mutate(value = value/5) %>%
  group_by(Disease, category, Age,cluster) %>%
  summarise(mean = mean(value), sd = sd(value), 
            fq = quantile(value,0.025), tq = quantile(value,0.975),
            median = median(value)) %>% 
  ggplot(aes(x = Age)) +
  geom_smooth(aes(y=median), color = 'gray60', se = F, method = 'loess') +
  geom_point(aes(y=median, color = cluster)) +
  geom_segment(aes(y = fq, yend = tq, xend = Age), color = 'gray70', size = 0.3) + 
  facet_wrap(~Disease, scales = 'free', ncol = 4) +
  ylab('No of cases diagnosed / 10,000 people at a given age') +
  scale_color_manual(values = ageonsetcolors) +
  guides(color = guide_legend('Age-of-onset cluster')) +
  ggtitle(as.character(unique(permRes$category)[i]))
ggsave(paste('./results/ageonset/',gsub('/','_',as.character(unique(permRes$category)[i])),'.pdf',sep=''),units = 'cm', width = 24, height = 32, useDingbats = F)

i=7
as.character(unique(permRes$category)[i])
# [1]  "respiratory/ent"
permRes %>%
  filter(category == as.character(unique(permRes$category)[i])) %>%
  mutate(value = value/5) %>%
  group_by(Disease, category, Age,cluster) %>%
  summarise(mean = mean(value), sd = sd(value), 
            fq = quantile(value,0.025), tq = quantile(value,0.975),
            median = median(value)) %>% 
  ggplot(aes(x = Age)) +
  geom_smooth(aes(y=median), color = 'gray60', se = F, method = 'loess') +
  geom_point(aes(y=median, color = cluster)) +
  geom_segment(aes(y = fq, yend = tq, xend = Age), color = 'gray70', size = 0.3) + 
  facet_wrap(~Disease, scales = 'free', ncol = 4) +
  ylab('No of cases diagnosed / 10,000 people at a given age') +
  scale_color_manual(values = ageonsetcolors) +
  guides(color = guide_legend('Age-of-onset cluster')) +
  ggtitle(as.character(unique(permRes$category)[i]))
ggsave(paste('./results/ageonset/',gsub('/','_',as.character(unique(permRes$category)[i])),'.pdf',sep=''),units = 'cm', width = 24, height = 22, useDingbats = F)

i=8
as.character(unique(permRes$category)[i])
# [1] "infections"
permRes %>%
  filter(category == as.character(unique(permRes$category)[i])) %>%
  mutate(value = value/5) %>%
  group_by(Disease, category, Age,cluster) %>%
  summarise(mean = mean(value), sd = sd(value), 
            fq = quantile(value,0.025), tq = quantile(value,0.975),
            median = median(value)) %>% 
  ggplot(aes(x = Age)) +
  geom_smooth(aes(y=median), color = 'gray60', se = F, method = 'loess') +
  geom_point(aes(y=median, color = cluster)) +
  geom_segment(aes(y = fq, yend = tq, xend = Age), color = 'gray70', size = 0.3) + 
  facet_wrap(~Disease, scales = 'free', ncol = 2) +
  ylab('No of cases diagnosed / 10,000 people at a given age') +
  scale_color_manual(values = ageonsetcolors) +
  guides(color = guide_legend('Age-of-onset cluster')) +
  ggtitle(as.character(unique(permRes$category)[i]))
ggsave(paste('./results/ageonset/',gsub('/','_',as.character(unique(permRes$category)[i])),'.pdf',sep=''),units = 'cm', width = 12, height = 12, useDingbats = F)

i=9
as.character(unique(permRes$category)[i])
# [1] "renal/urology"
permRes %>%
  filter(category == as.character(unique(permRes$category)[i])) %>%
  mutate(value = value/5) %>%
  group_by(Disease, category, Age,cluster) %>%
  summarise(mean = mean(value), sd = sd(value), 
            fq = quantile(value,0.025), tq = quantile(value,0.975),
            median = median(value)) %>% 
  ggplot(aes(x = Age)) +
  geom_smooth(aes(y=median), color = 'gray60', se = F, method = 'loess') +
  geom_point(aes(y=median, color = cluster)) +
  geom_segment(aes(y = fq, yend = tq, xend = Age), color = 'gray70', size = 0.3) + 
  facet_wrap(~Disease, scales = 'free', ncol = 4) +
  ylab('No of cases diagnosed / 10,000 people at a given age') +
  scale_color_manual(values = ageonsetcolors) +
  guides(color = guide_legend('Age-of-onset cluster')) +
  ggtitle(as.character(unique(permRes$category)[i]))
ggsave(paste('./results/ageonset/',gsub('/','_',as.character(unique(permRes$category)[i])),'.pdf',sep=''),units = 'cm', width = 24, height = 12, useDingbats = F)

i=10
as.character(unique(permRes$category)[i])
# [1] "endocrine/diabetes"
permRes %>%
  filter(category == as.character(unique(permRes$category)[i])) %>%
  mutate(value = value/5) %>%
  group_by(Disease, category, Age,cluster) %>%
  summarise(mean = mean(value), sd = sd(value), 
            fq = quantile(value,0.025), tq = quantile(value,0.975),
            median = median(value)) %>% 
  ggplot(aes(x = Age)) +
  geom_smooth(aes(y=median), color = 'gray60', se = F, method = 'loess') +
  geom_point(aes(y=median, color = cluster)) +
  geom_segment(aes(y = fq, yend = tq, xend = Age), color = 'gray70', size = 0.3) + 
  facet_wrap(~Disease, scales = 'free', ncol = 4) +
  ylab('No of cases diagnosed / 10,000 people at a given age') +
  scale_color_manual(values = ageonsetcolors) +
  guides(color = guide_legend('Age-of-onset cluster')) +
  ggtitle(as.character(unique(permRes$category)[i]))
ggsave(paste('./results/ageonset/',gsub('/','_',as.character(unique(permRes$category)[i])),'.pdf',sep=''),units = 'cm', width = 24, height = 12, useDingbats = F)

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


xx = sumx %>%
  mutate(category = fct_reorder(category,-med)) %>%
  group_by(category) %>%
  summarise(`Disease Category` = length(unique(Disease))) %>%
  mutate(`Disease Category` = paste(category,' (n=',`Disease Category`,')',sep='')) 
discatcolors2 = discatcolors
names(discatcolors2)=unname(setNames(xx$`Disease Category`,xx$category)[names(discatcolors2)])
catsum <- left_join(xx,mutate(sumx,category = fct_reorder(category,-med))) %>%
  mutate(`Disease Category` = fct_reorder(`Disease Category`,-med)) %>%
  ggplot(aes(x = `Disease Category`, y = med)) + 
  geom_boxplot(aes(fill = `Disease Category`), outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.5, width = 0.2) +
  scale_fill_manual(values = discatcolors2)+
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

# pam2 = function(x,k){
#   mydis = diss(x,'CORT')
#   px = cluster::pam(mydis,k,diss=T)
#   list(cluster = px$clustering)
# }
# gapstat = clusGap(medRate, FUNcluster = pam2, K.max = 20, B=1000)
# saveRDS(gapstat,'./data/processed/ageonset/gapstat.rds')
# k=sapply(c("firstSEmax", "Tibs2001SEmax", "globalSEmax",
#            "firstmax", "globalmax"),function(met){
#              maxSE(gapstat$Tab[,"gap"],gapstat$Tab[,"SE.sim"],method = met)
#            })
# k
# ggplot(as.data.frame(gapstat$Tab),aes(x=1:20))+
#   annotate(geom='segment',
#            x= k['Tibs2001SEmax'],
#            xend = k['Tibs2001SEmax'],
#            y = 0.85, yend = 0.81,
#            color = 'darkred',size=1,arrow=arrow(length = unit(5,'pt'))) +
#   geom_segment(aes(y=gap-SE.sim,yend=gap+SE.sim,xend=1:20)) +
#   geom_point(aes(y=gap))+
#   geom_line(aes(y=gap))+
#   xlab('k')+ylab('Gap') 
# ggsave('./results/ageonset/gapstat.pdf',units ='cm',width = 8,height = 6,useDingbats=F)
# ggsave('./results/ageonset/gapstat.png',units ='cm',width = 8,height = 6)
# mydis = diss(medRate,'CORT')
# pc = cluster::pam(mydis,k=k['Tibs2001SEmax'],diss=T)
# saveRDS(pc,'./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')
# saveRDS(mydis,'./data/processed/ageonset/distmat_cort.rds')

pc = readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')
mydis = readRDS('./data/processed/ageonset/distmat_cort.rds')

pam_cl = medRate %>%
  reshape2::melt() %>%
  set_names(c('Disease','Age','val')) %>%
  mutate(cl = pc$clustering[as.character(Disease)]) %>%
  arrange(Disease) %>%
  left_join(data.frame(Disease = names(disTreecl), Category = unname(disTreecl))) %>%
  left_join(reshape2::melt(table(pc$clustering)) %>%set_names(c('cl','disnum'))) %>% 
  mutate(cluster = paste('Cluster:',cl,', # of Diseases: ',disnum,sep='')) %>%
  ggplot(aes(x=Age,color = Category)) +
  geom_hline(yintercept = 0,color='black',size=0.1)+
  geom_line(aes(y=val,group=Disease),alpha=0.5,size=0.5)+
  facet_wrap(~cluster, ncol=2) +
  scale_color_manual(values=discatcolors)+
  theme(legend.position = 'right',
        panel.grid.major.x = element_line(color='gray60',size=0.2,linetype = 'dashed'),
        panel.grid.minor.x = element_line(color='gray70',size=0.1,linetype = 'dashed')) + 
  ylab('Disease Onset Rate') + labs(caption = 'clustered by PAM with CORT distance') 
ggsave('./results/ageonset/pam_cort_medRate.pdf',pam_cl, units = 'cm',width = 18,height = 12, useDingbats =F)  
ggsave('./results/ageonset/pam_cort_medRate.png',pam_cl, units = 'cm',width = 18,height = 12)  

overlapdat <- data.frame(Disease = names(pc$clustering), 
                         AgeOnsetCluster = as.factor(unname(pc$clustering))) %>%
  left_join(data.frame(Disease = names(disTreecl), Category = unname(disTreecl))) %>%
  group_by(AgeOnsetCluster, Category) %>%
  summarise(numDis = length(unique(Disease))) %>%
  spread(Category,numDis, fill = 0) %>%
  as.data.frame()
rownames(overlapdat) = overlapdat$AgeOnsetCluster
overlapdat$AgeOnsetCluster = NULL
overlapdat = as.matrix(overlapdat)
overlapdat = overlapdat[as.character(unique(pc$clustering)),]
overlapdat = 100 * overlapdat / rowSums(overlapdat)
pheatmap::pheatmap(overlapdat,
                   cluster_rows = F,
                   color = colorRampPalette(brewer.pal(8,'Oranges'))(max(overlapdat)),
                   number_color = 'black', cellwidth = 20, cellheight = 20, 
                   filename = paste('./results/ageonset/age_cat_overlap.pdf',sep=''),
                   display_numbers = T, number_format = '%.0f')
pheatmap::pheatmap(overlapdat,
                   cluster_rows = F,
                   color = colorRampPalette(brewer.pal(8,'Oranges'))(max(overlapdat)),
                   number_color = 'black', cellwidth = 20, cellheight = 20, 
                   filename = paste('./results/ageonset/age_cat_overlap.png',sep=''),
                   display_numbers = T, number_format = '%.0f')


data.frame(disease = names(pc$clustering), cluster = pc$clustering) %>%
  mutate(disease_category = disTreecl[as.character(disease)]) %>%
  arrange(cluster) %>% 
  write_tsv('./data/processed/ageonset/ageonset_clusters.tsv')
