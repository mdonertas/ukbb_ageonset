source('./scripts/00-setup.R')
library(igraph)
traits <- readRDS('./data/processed/traits_clean/traitData_baseline_additions2.rds')
disCoding <- read_tsv('./data/raw/ukbb/datacoding/coding6.tsv')
disTree <- graph_from_data_frame(select(disCoding,parent_id,node_id),directed = T)
igraph_options(plot.layout=layout_as_tree)
SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated.rds')
prevDF = readRDS('./data/processed/traits_clean/SRdisease_prop_prevdf.rds')

# system('mkdir -p results/UKBB_disease_EDA')
pdf('./results/UKBB_disease_EDA/disTree.pdf',width=20,height = 5,useDingbats = F)
plot(disTree,
     vertex.size=0.75,
     vertex.frame.color=NA,
     vertex.label=NA,
     asp=0.1,
     edge.arrow.size=0.2,
     edge.color='gray60',
     vertex.color=c('gray30'))
dev.off()

select(traits, nSRdiseaseProp) %>%
  mutate(nSRdisease0inc = ifelse(is.na(nSRdiseaseProp),0,nSRdiseaseProp)) %>%
  summary()

# nSRdiseaseProp   nSRdisease0inc  
# Min.   : 1.00    Min.   : 0.000  
# 1st Qu.: 3.00    1st Qu.: 0.000  
# Median : 5.00    Median : 3.000  
# Mean   : 6.18    Mean   : 4.628  
# 3rd Qu.: 8.00    3rd Qu.: 7.000  
# Max.   :65.00    Max.   :65.000  
# NA's   :121826                   

select(traits, nSRdiseaseProp, Sex) %>%
  mutate(nSRdisease0inc = ifelse(is.na(nSRdiseaseProp),0,nSRdiseaseProp)) %>%
  wilcox.test(nSRdisease0inc ~ as.numeric(as.factor(Sex)), .)

# Wilcoxon rank sum test with continuity correction
# 
# data:  nSRdisease0inc by as.numeric(as.factor(Sex))
# W = 3.0138e+10, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

select(traits, nSRdiseaseProp, Sex) %>%
  mutate(nSRdisease0inc = ifelse(is.na(nSRdiseaseProp),0,nSRdiseaseProp)) %>%
  wilcox.test(nSRdiseaseProp ~ as.numeric(as.factor(Sex)), .)
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  nSRdiseaseProp by as.numeric(as.factor(Sex))
# W = 1.7347e+10, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

numDiseases_perInd <- select(traits, eid, nSRdiseaseProp, Sex) %>%
  na.omit() %>%
  ggplot(aes(x = Sex, y = nSRdiseaseProp)) +
  geom_violin(aes(fill = Sex)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values=sexcolors[c('Female','Male')]) +
  ylab('Number of Self-Reported Diseases') + 
  xlab('')
ggsave('./results/UKBB_disease_EDA/numDiseases_perInd.pdf', numDiseases_perInd,  useDingbats=F, height=8, width = 8, units = 'cm')
ggsave('./results/UKBB_disease_EDA/numDiseases_perInd.png', numDiseases_perInd, height=8, width = 8, units = 'cm')

lmxdat <- prevDF %>%
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
  scale_x_continuous(breaks = log10(c(1,10,100,300)), labels = c(1,10,50,300)) +
  scale_y_continuous(breaks = log10(c(1,10,100,300)), labels = c(1,10,50,300)) +
  xlab('Number of Cases in 1,000 Males') + 
  ylab('Number of Cass in 1,000 Females') + 
  scale_alpha_manual(values = c(0.4,1)) +
  geom_text_repel(aes(label = sexSpecName), size = 6 / pntnorm, box.padding = 0) +
  guides(alpha = F) +
  scale_color_manual("Sex", values=sexcolors[c('Female','Male')]) +
  coord_equal(xlim = c(log10(0.003),log10(400)),ylim = c(log10(0.003),log10(400))) 
ggsave('./results/UKBB_disease_EDA/disPrevScatter.pdf', disPrev,  useDingbats=F, height=8, width = 8, units = 'cm')
ggsave('./results/UKBB_disease_EDA/disPrevScatter.png', disPrev, height=8, width = 8, units = 'cm')

top10dis_sex <- prevDF %>%
  arrange(-fCase) %>%
  head(20) %>%
  mutate(fPrev = - fPrev) %>%
  mutate(fPrev = 1000 * fPrev, mPrev = 1000 * mPrev) %>%
  mutate(Disease = factor(Disease, levels = rev(reorder(Disease,nCases)))) %>%
  select(Disease, fPrev, mPrev) %>%
  gather( key = 'sex', value = 'Freq', -Disease) %>%
  mutate(sex = ifelse(sex == 'fPrev', 'Female','Male')) %>%
  ggplot(aes(y = Freq, x = Disease, fill = sex)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual("Sex", values=sexcolors[c('Female','Male')]) +
  coord_flip() +
  ylab('Number of Cases in 1,000') +
  xlab('')+
  scale_y_continuous(breaks = c(-400,-200,0,200,400),labels = c(400,200,0,200,400),limits = c(-450,450), position = 'right') +
  theme(axis.text = element_text(size = 8))

ggsave('./results/UKBB_disease_EDA/top10dis_bysex.pdf', top10dis_sex, useDingbats = F, height = 8, width = 18, units = 'cm')
ggsave('./results/UKBB_disease_EDA/top10dis_bysex.png', top10dis_sex, height = 8, width = 18, units = 'cm')

p1 <- ggarrange(ggarrange(numDiseases_perInd, disPrev, ncol = 2, nrow = 1, labels = 'auto', legend = 'none'), 
          top10dis_sex, ncol = 1, nrow = 2, labels = c('','c'), common.legend = T, legend = 'bottom', heights = c(1.25,1), align = 'h')

ggsave('./results/UKBB_disease_EDA/figure1.pdf', p1,  useDingbats=F, height=18, width = 18, units = 'cm')
ggsave('./results/UKBB_disease_EDA/figure1.png', p1, height=18, width = 18, units = 'cm')


