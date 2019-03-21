source('./scripts/00-setup.R')
traits <- readRDS('./data/processed/traits_clean/traitData_baseline_additions2.rds')
igraph_options(plot.layout=layout_as_tree)
SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated.rds')
prevDF = readRDS('./data/processed/traits_clean/SRdisease_prop_prevdf.rds')
disSet <- readRDS('./data/processed/traits_clean/SRdiseaseSet.rds')

# system('mkdir -p ./results/selectedDisease_EDA')
selectedNodes <- filter(disCoding,meaning%in%disSet$Disease)$meaning
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

# SRdisease_selected <- SRdisease %>%
#   filter(Disease %in% disSet$Disease)
# saveRDS(SRdisease_selected,'./data/processed/traits_clean/SRdisease_baseline_propagated_selected.rds')

SRdisease_selected <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated_selected.rds')

xx <- traits %>% 
  select(eid, Sex) %>%
  right_join(SRdisease_selected) %>%
  na.omit() %>% 
  group_by(eid, Sex) %>%
  summarise(nSRdiseaseProp = length(unique(Disease))) %>%
  ungroup()

summary(xx$nSRdiseaseProp)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   3.000   5.000   5.407   7.000  46.000 
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

disset <- disSet$Disease
disTreesub <- induced_subgraph(disTree,c('Top',disset))
disTreeclsub=factor(unname(disTreecl[V(disTreesub)$name]))
l <- layout_as_tree(disTreesub)
colsx <- discatcolors[as.character(disTreeclsub)]
nms = V(disTreesub)$name
nms[!nms%in%disTree_level1]=NA
colsx2 <- discatcolors[na.omit(nms)]
names(colsx2)=sapply(strsplit(names(colsx2),'/'),function(x)x[1])
colsx[is.na(colsx)] = 'gray50'

pdf('./results/selectedDisease_EDA/disTree_selected2.pdf',width=18,height = 5,useDingbats = F)
plot(disTreesub,
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

agestrNumdis_5 <- traits %>% 
  select(eid, Sex) %>%
  right_join(SRdisease_selected) %>%
  mutate(AgeGr = cut(Age, breaks = seq(0,65,5))) %>%
  group_by(eid, Sex, AgeGr) %>%
  summarise(numDisease = length(unique(Disease))) %>%
  ungroup() %>%
  group_by(Sex, AgeGr) %>%
  summarise(Mean = mean(numDisease),
            Maximum = max(numDisease)) %>%
  na.omit()  %>%
  gather(key = 'Type', value = 'Number of Diseases', -Sex, -AgeGr)

agestrNumdis_5p <- agestrNumdis_5 %>%
  mutate(Type = factor(Type, levels = c('Mean','Maximum'))) %>%
  ggplot(aes(x = AgeGr, color = Sex, y = `Number of Diseases`, group = Sex)) +
  geom_point(size = 2) + geom_line(size = 2) +
  scale_color_manual("Sex", values=sexcolors[c('Female','Male')]) + 
  xlab('Age Group') +
  facet_wrap(~Type, scales = 'free_y', ncol = 1, nrow = 3)+
  theme(axis.text.x = element_text(angle=90))

ggsave('./results/selectedDisease_EDA/agestrNumdis_5.pdf', agestrNumdis_5p,  useDingbats=F, height=12, width = 8, units = 'cm')
ggsave('./results/selectedDisease_EDA/agestrNumdis_5.png', agestrNumdis_5p, height=12, width = 8, units = 'cm')

agestrNumdis_10 <- traits %>% 
  select(eid, Sex) %>%
  right_join(SRdisease_selected) %>%
  mutate(AgeGr = cut(Age, breaks = seq(0,70,10))) %>%
  group_by(eid, Sex, AgeGr) %>%
  summarise(numDisease = length(unique(Disease))) %>%
  ungroup() %>%
  group_by(Sex, AgeGr) %>%
  summarise(Mean = mean(numDisease),
            Maximum = max(numDisease)) %>%
  na.omit()  %>%
  gather(key = 'Type', value = 'Number of Diseases', -Sex, -AgeGr)

agestrNumdis_10p <- agestrNumdis_10 %>%
  mutate(Type = factor(Type, levels = c('Mean','Maximum'))) %>%
  ggplot(aes(x = AgeGr, color = Sex, y = `Number of Diseases`, group = Sex)) +
  geom_point(size = 2) + geom_line(size = 2) +
  scale_color_manual("Sex", values=sexcolors[c('Female','Male')]) + 
  xlab('Age Group') +
  facet_wrap(~Type, scales = 'free_y', ncol = 1, nrow = 3)+
  theme(axis.text.x = element_text(angle=90))

ggsave('./results/selectedDisease_EDA/agestrNumdis_10.pdf', agestrNumdis_10p,  useDingbats=F, height=12, width = 8, units = 'cm')
ggsave('./results/selectedDisease_EDA/agestrNumdis_10.png', agestrNumdis_10p, height=12, width = 8, units = 'cm')

cats <- unique(disTreecl)
catsp_age <- lapply(cats, function(cat){
  catdis <- names(disTreecl)[which(disTreecl %in% cat)]
  SRdisease_selected %>%
    filter(Disease%in%catdis) %>%
    left_join(select(traits,eid,Sex)) %>%
    mutate(AgeGr = cut(Age, breaks = seq(0,70,10))) %>%
    group_by(eid, Sex, AgeGr) %>%
    summarise(numDisease = length(unique(Disease))) %>%
    ungroup() %>%
    group_by(Sex, AgeGr) %>%
    summarise(Mean = mean(numDisease),
              Maximum = max(numDisease)) %>%
    na.omit()  %>%
    gather(key = 'Type', value = 'Number of Diseases', -Sex, -AgeGr) %>%
    ungroup() 
})
catsp_age <- reshape2::melt(catsp_age, id.vars = colnames(catsp_age[[1]])) 
catsp_age$Category <- cats[catsp_age$L1]

catsp_age_p <- catsp_age %>%
  mutate(Type = factor(Type, levels = c('Mean','Maximum'))) %>%
  filter(Type == 'Mean') %>%
  ggplot(aes(x = AgeGr, color = Sex, y = `Number of Diseases`, group = Sex)) +
  geom_point(size = 2) + geom_line(size = 2) +
  scale_color_manual("Sex", values=sexcolors[c('Female','Male')]) + 
  xlab('Age Group') + ylab('Averange Number of Diseases')+
  facet_wrap(~Category, scales = 'free_y', ncol = 3) +
  theme(axis.text.x = element_text(angle=90))

ggsave('./results/selectedDisease_EDA/age_cat_strdis.pdf', catsp_age_p,  useDingbats=F, height=15, width = 18, units = 'cm')
ggsave('./results/selectedDisease_EDA/age_cat_strdis.png', catsp_age_p, height=15, width = 18, units = 'cm')
