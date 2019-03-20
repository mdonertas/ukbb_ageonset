source('./scripts/00-setup.R')
library(igraph)
traits <- readRDS('./data/processed/traits_clean/traitData_baseline_additions2.rds')
disCoding <- read_tsv('./data/raw/ukbb/datacoding/coding6.tsv')
disTree <- graph_from_data_frame(select(disCoding,parent_id,node_id),directed = T)
igraph_options(plot.layout=layout_as_tree)
SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated.rds')
prevDF = readRDS('./data/processed/traits_clean/SRdisease_prop_prevdf.rds')
disSet <- readRDS('./data/processed/traits_clean/SRdiseaseSet.rds')

# system('mkdir -p ./results/selectedDisease_EDA')
selectedNodes <- filter(disCoding,meaning%in%disSet$Disease)$node_id
pdf('./results/selectedDisease_EDA/disTree_selected.pdf',width=20,height = 5)
plot(disTree,
     vertex.size=c(0.75,0.75)[1+(V(disTree)$name%in%selectedNodes)],
     vertex.frame.color=NA,
     vertex.label=NA,
     asp=0.1,
     edge.arrow.size=0.2,
     edge.color='gray80',
     vertex.color=c('gray70','gray25')[1+(V(disTree)$name%in%selectedNodes)])
dev.off()

disSet <- prevDF %>%
  filter(Disease %in% disSet$Disease)

lmxdat <- disSet %>%
  mutate(mPrev = log10(1000 * mPrev), fPrev = log10(1000 * fPrev)) %>%
  mutate(mPrev = ifelse(is.na(mPrev),0,mPrev),fPrev = ifelse(is.na(fPrev),0,fPrev)) %>%
  na.omit() 

lmx <- lm(fPrev ~ mPrev, lmxdat)
lmxres <- lmx$residuals
names(lmxres) <- lmxdat$Disease

lmxdat <- lmxdat %>%
  mutate(resx = lmxres) %>%
  mutate(aresx = abs(resx)) %>%
  mutate(sexSpecific = aresx >= (3*sd(aresx)),
         sex = c('Male','Female')[1 + (fPrev > mPrev)]) %>%
  mutate(sexSpecName = ifelse(sexSpecific, Disease, NA))

disPrev <- lmxdat %>%
  ggplot(aes(x = mPrev, y = fPrev)) +
  geom_abline(slope = 1, intercept = 0, color = 'gray50', linetype = 'dashed') + 
  geom_point(aes(alpha = sexSpecific, color = sex)) +
  geom_smooth(method = 'lm') +
  scale_x_continuous(breaks = log10(c(1,10,30,100,300)), labels = c(1,10,30,50,300)) +
  scale_y_continuous(breaks = log10(c(1,10,30,100,300)), labels = c(1,10,30,50,300)) +
  xlab('Number of Cases in 1,000 Males') + 
  ylab('Number of Cass in 1,000 Females') + 
  scale_alpha_manual(values = c(0.4,1)) +
  geom_text_repel(aes(label = sexSpecName), size = 6 / pntnorm, box.padding = 0) +
  guides(alpha = F) +
  scale_color_manual("Sex", values=sexcolors[c('Female','Male')]) +
  coord_equal(xlim = c(log10(1),log10(400)),ylim = c(log10(1),log10(400))) 
ggsave('./results/selectedDisease_EDA/disPrevScatter.pdf', disPrev,  useDingbats=F, height=8, width = 8, units = 'cm')
ggsave('./results/selectedDisease_EDA/disPrevScatter.png', disPrev, height=8, width = 8, units = 'cm')

SRdisease_selected <- SRdisease %>%
  filter(Disease %in% disSet$Disease)
saveRDS(SRdisease_selected,'./data/processed/traits_clean/SRdisease_baseline_propagated_selected.rds')

xx <- traits %>% 
  select(eid, Sex) %>%
  right_join(SRdisease_selected) %>%
  na.omit() %>% 
  group_by(eid, Sex) %>%
  summarise(nSRdiseaseProp = length(unique(Disease))) %>%
  ungroup()

summary(xx$nSRdiseaseProp)

wilcox.test(xx$nSRdiseaseProp~as.factor(xx$Sex))

# Wilcoxon rank sum test with continuity correction
# 
# data:  xx$nSRdiseaseProp by as.factor(xx$Sex)
# W = 1.5835e+10, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

numDiseases_perInd <- traits %>% 
  select(eid, Sex) %>%
  right_join(SRdisease_selected) %>%
  na.omit() %>% 
  group_by(eid, Sex) %>%
  summarise(nSRdiseaseProp = length(unique(Disease))) %>%
  ungroup() %>%
  ggplot(aes(x = Sex, y = nSRdiseaseProp)) +
  geom_violin(aes(fill = Sex)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values=sexcolors[c('Female','Male')]) +
  ylab('Number of Self-Reported Diseases') + 
  xlab('')
ggsave('./results/selectedDisease_EDA/numDiseases_perInd.pdf', numDiseases_perInd,  useDingbats=F, height=8, width = 8, units = 'cm')
ggsave('./results/selectedDisease_EDA/numDiseases_perInd.png', numDiseases_perInd, height=8, width = 8, units = 'cm')

p1 <- ggarrange(numDiseases_perInd, disPrev, ncol = 2, nrow = 1, labels = 'auto', common.legend = T, legend = 'bottom', align = 'hv')

ggsave('./results/selectedDisease_EDA/figure1.pdf', p1,  useDingbats=F, height=9, width = 18, units = 'cm')
ggsave('./results/selectedDisease_EDA/figure1.png', p1, height=9, width = 18, units = 'cm')

disCoding <- read_tsv('./data/raw/ukbb/datacoding/coding6.tsv')
nms=setNames(disCoding$meaning,disCoding$node_id)
nms=c(nms,`0`='Top')
nmsdf=data.frame(node=unname(nms[as.character(disCoding$node_id)]),parent=unname(nms[as.character(disCoding$parent_id)]))
disTree <- graph_from_data_frame(select(nmsdf,parent,node),directed = T)
nodelabels=setNames(disCoding$meaning,disCoding$node_id)
disTreecl <- reshape::melt(sapply(neighbors(disTree,'Top')$name,function(nm){
  subcomponent(disTree,nm,'out')$name
}))%>%
  rename(node=value,
         cluster=L1)

disTreecl=setNames(disTreecl$cluster,disTreecl$node)

disset <- disSet$Disease
disTree <- induced_subgraph(disTree,c('Top',disset))
disTreecl=factor(unname(disTreecl[V(disTree)$name]))

level1=neighbors(disTree,'Top','out')$name
net=disTree
l <- layout_as_tree(net)
colsx <- c(brewer.pal(8,'Dark2'),brewer.pal(8,'Set1'))[as.numeric(disTreecl)]
colsx[is.na(colsx)]='black'
nms = V(disTree)$name
nms[!nms%in%level1]=NA
nms=sapply(strsplit(nms,'/'),function(x)x[1])
colsx2 <- setNames(colsx,nms)
colsx2 <- colsx2[!is.na(names(colsx2))]
pdf('./results/selectedDisease_EDA/disTree_selected2.pdf',width=18,height = 5,useDingbats = F)
plot(disTree2,
     layout=l,
     vertex.color=colsx,
     vertex.size=3,
     vertex.label = NA,
     vertex.frame.color=NA,
     asp=0.1,
     edge.arrow.size=0.2,
     edge.color='gray80',
     vertex.color=c('gray70'))
legend('bottom',col = colsx2,legend = names(colsx2),pch = 19, cex = 1, horiz = T, bty = 'n')
dev.off()

unique(colsx)
