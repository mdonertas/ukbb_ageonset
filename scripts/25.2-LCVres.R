source('./scripts/00-setup.R')
disIDs = gsub('a','',list.files('./results/caseControl/'))
disCoding <- setNames(disCoding$meaning,disCoding$node_id)
disCoding = disCoding[as.character(disIDs)]
aooclusters = readRDS('./data/processed/ageonset/clusters_pam_Tibs2001SEmax.rds')$cluster
signifSNPs <- lapply(paste('../ukbb_ageonset/data/processed/caseControl/a',
         disIDs, '/signif_gwasRes_proxyGenes.rds', sep = ''),
   function(x) {
     x = readRDS(x)
     x = filter(x, !(CHR == mhcchr & BP >= mhcstart & BP <= mhcend))
     select(x,SNP,CHR,BP,Ref,Alt) %>% unique()
   })
names(signifSNPs) = disIDs
signifSNPs = sapply(signifSNPs,function(x){
  filter(x, !SNP %in% unique(x$SNP[which(duplicated(x$SNP))])) %>%
    unique() %>%
    nrow()
})
pairs = list.files('./data/processed/LCV/',full.names = T)
lcvs = lapply(pairs,readRDS)
pairids = gsub('.rds','',sapply(strsplit(pairs,'/'),function(x)x[length(x)]))
names(lcvs) = pairids
lcvs2 = reshape2::melt(lcvs)
lcvs3 = lcvs2 %>%
  filter(L2 == 'h2.zscore') %>%
  mutate(L2 = paste(L2,1:2,sep='')) 
lcvs4 = lcvs2 %>%
  filter(L2 == 'pval.fullycausal') %>%
  mutate(L2 = paste(L2,1:2,sep='')) 
lcv = rbind(rbind(filter(lcvs2, !L2%in%c('h2.zscore','pval.fullycausal')),lcvs3) ,lcvs4)
rm(lcvs2,lcvs3,lcvs4)
lcv = lcv %>%
  unique() %>%
  spread(L2,value)

lcv = lcv %>%
  set_names(c('pair','gcp','gcp_se','heritability1','heritability2','fullcausal.p.1','fullcausal.p.2','pval_2tailed','rho_err','rho_est','z')) %>%
  separate(pair,into = c('dis1','dis2'), remove = F) %>%
  mutate(disease1 = disCoding[as.character(dis1)]) %>%
  mutate(disease2 = disCoding[as.character(dis2)]) 

lcv1 = lcv %>% filter(gcp < 0) %>% rename(disease2 = disease1, disease1 = disease2,
                                   heritability1 = heritability2,
                                   heritability2 = heritability1,
                                   fullcausal.p.1=fullcausal.p.2,
                                   fullcausal.p.2=fullcausal.p.1) %>%
  mutate(gcp = -gcp) 
lcv = rbind(filter(lcv, gcp>0),lcv1)
lcv$uplevel = sapply(1:nrow(lcv),function(i){
  (lcv$disease2[i] %in% subcomponent(disTree,as.character(lcv$disease1[i]),mode = 'out')$name)
})
lcv$downlevel = sapply(1:nrow(lcv),function(i){
  (lcv$disease1[i] %in% subcomponent(disTree,as.character(lcv$disease2[i]),mode = 'out')$name)
})

lcv2 = lcv %>%
  mutate(vertconnect = uplevel|downlevel) %>%
  filter(!vertconnect) %>%
  filter(heritability1>=7 & heritability2>=7) %>%
  mutate(padj = p.adjust(pval_2tailed, method = 'fdr'))

n = lcv2 %>%
  filter(padj<=0.01 & gcp>0.6) %>%
  select(disease1, disease2, everything() ) %>%
  graph_from_data_frame(directed = T)

V(n)$aoocluster = as.factor(aooclusters[V(n)$name])
V(n)$outdeg = degree(n,mode='all')[V(n)$name]
V(n)$namex = gsub('ent','ENT',V(n)$name)
library(ggnetwork)
ggplot(n) +
  geom_edges(aes(x = x, y=y, xend =xend, yend = yend),
             arrow = arrow(length = unit(3, "pt"), type = "closed"), color = 'gray80', alpha = 0.7) +
  # scale_color_gradient2(low = muted('blue'), mid = 'white', high = muted('red'), midpoint = 0) +
  geom_nodes(aes(x = x, y= y, color = as.factor(aoocluster), size = outdeg )) +
  scale_color_manual(values = muted(ageonsetcolors,l = 60)) +
  geom_nodetext_repel(aes(label = namex, x=x,y=y, color = as.factor(aoocluster)), size = 6/pntnorm,alpha=0.9) +
  theme_blank(base_size = 8)  +
  # guides(color = guide_colorbar('rho', barheight = unit(3,'mm'))) +
  guides(color = guide_legend('Age-of-onset Cluster'), size = guide_legend('Degree')) +
  theme(legend.position = 'bottom') +
  scale_size_continuous(range=c(1,4))
ggsave('./results/LCV/network.pdf',units = 'cm',width = 16,height = 14,useDingbats = F)
ggsave('./results/LCV/network.png',units = 'cm',width = 16,height = 14)

lcv2 %>% 
  filter(padj<=0.01 & gcp>0.6) %>%
  ggplot(aes(x = disease1, y = disease2)) +
  geom_point(aes(size = gcp, color = rho_est)) +
  theme_pubr(base_size = 6) +
  theme(axis.text.x= element_text(angle=90,hjust = 1,vjust=0.5),
        panel.grid.major = element_line(),legend.position = 'top') + 
  scale_y_discrete(position = 'right') +
  scale_color_gradient2(low=muted('blue'),mid = 'white',high = muted('red'),midpoint = 0) +
  scale_size_continuous(range = c(0.2,2.5)) +
  xlab(NULL) + ylab(NULL) + guides(color = guide_colorbar('rho',barheight = unit(2,'mm')), size = guide_legend('GCP'))
ggsave('./results/LCV/heatmap.pdf',units = 'cm',width = 12,height = 14, useDingbats = F)
ggsave('./results/LCV/heatmap.png',units = 'cm',width = 12,height = 14)

clstats = lcv2 %>% 
  filter(padj<=0.01 & gcp>0.6) %>%
  mutate( cluster1 = aooclusters[disease1],
          cluster2 = aooclusters[disease2]) %>%
  group_by(cluster1,cluster2) %>%
  summarise(n = n()) 

disset = unname(disCoding[unique(unlist(strsplit(names(lcvs),'_')))])
resx = sapply(1:3,function(i){
  sapply(1:3,function(j){
    gr1 = intersect(names(which(aooclusters == i)),disset)
    gr2 = intersect(names(which(aooclusters == j)),disset)
    lapply(gr1,function(g1){
      lapply(gr2,function(g2){
        exc = g2 %in% subcomponent(disTree, as.character(g1), 
                                   mode = 'out')$name | 
          g2 %in% subcomponent(disTree ,as.character(g1), mode = 'in')$name
        data.frame(g1,g2,exc)
      }) %>% reshape2::melt(id.vars = c('g1','g2','exc')) %>% select(-L1)
    }) %>% 
      reshape2::melt(id.vars = c('g1','g2','exc')) %>% 
      filter(!exc) %>% 
      select(-L1) %>%
      nrow()
  })
})

fi11=fisher.test(matrix(c(5,resx[1,1],17+9,sum(resx[1,2:3])),ncol=2))
fi12=fisher.test(matrix(c(17,resx[1,2],5+9,sum(resx[1,c(1,3)])),ncol=2))
fi13=fisher.test(matrix(c(9,resx[1,3],5+17,sum(resx[1,1:2])),ncol=2))
fi21=fisher.test(matrix(c(15,resx[2,1],38,sum(resx[1,2:3])),ncol=2))
fi22=fisher.test(matrix(c(30,resx[2,2],15+8,sum(resx[1,c(1,3)])),ncol=2))
fi23=fisher.test(matrix(c(8,resx[2,3],15+30,sum(resx[1,1:2])),ncol=2))
fi31=fisher.test(matrix(c(0,resx[3,1],7,sum(resx[1,c(2,3)])),ncol=2))
fi32=fisher.test(matrix(c(7,resx[3,2],0,sum(resx[1,c(1,3)])),ncol=2))
fi33=fisher.test(matrix(c(0,resx[3,3],7,sum(resx[1,c(1,2)])),ncol=2))

clstats$possible = c(resx[1,1],resx[1,2],resx[1,3],resx[2,1],resx[2,2],resx[2,3],resx[3,2])
clstats %>%
  mutate(perc = n/possible) %>%
  ggplot(aes(x = as.factor(cluster1), y = perc, fill = as.factor(cluster2))) +
  geom_bar(stat='identity',position=position_dodge2(preserve = "single",padding = 0.05)) +
  scale_fill_manual(values = ageonsetcolors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), limits = c(0,0.055)) +
  xlab('Disease 1 (Age-of-onset Cluster)') + 
  guides(fill = guide_legend('Disease 2\nAge-of-onset Cluster', title.hjust = 1, title.vjust = 0.5 )) +
  ylab(NULL) + ggtitle('Percent Causal Relationship\nDisease 1 -> Disease2') +
  theme_pubr(base_size = 6, legend = 'bottom') +
  theme(legend.key.size = unit(2,'mm')) +
  annotate('text',x=2.3,label='**', y= 0.025, size = 6/pntnorm) +
  annotate('text',x=3,label='**', y= 0.022, size = 6/pntnorm)
ggsave('./results/LCV/clustersummary.pdf',units = 'cm',width = 5,height = 5, useDingbats = F)
ggsave('./results/LCV/clustersummary.png',units = 'cm',width = 5,height = 5)

mynet2 = clstats %>%
  mutate(perc = n/possible) %>%
  mutate(pair = paste(cluster1,cluster2,sep='-')) %>%
  graph_from_data_frame(directed = T)
V(mynet2)$namex = V(mynet2)$name
ggplot(mynet2,aes(x = x, y= y)) +
  geom_edges(aes(x = x, y=y, xend =xend, yend = yend, size = perc),
             arrow = arrow(length = unit(5, "pt"), type = "closed"), color = 'gray80', self) +
  geom_nodes(size = 4, aes(color = namex)) +
  scale_color_manual(values = muted(ageonsetcolors,l = 60)) +
  scale_size_continuous(range = c(0.2,1)) +
  geom_nodetext(aes(label = namex),size = 6/pntnorm) +
  geom_edgetext(aes(xend=xend,yend=yend,label = round(100*perc,2)), size = 6/pntnorm) +
  theme_blank() +
  guides(size = F, color = F)
ggsave('./results/LCV/clustersummary.pdf',units = 'cm',width = 3,height = 3, useDingbats = F)

