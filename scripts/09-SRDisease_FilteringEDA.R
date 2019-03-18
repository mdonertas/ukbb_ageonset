source('./scripts/00-setup.R')
prevDF = readRDS('./data/processed/traits_clean/SRdisease_prop_prevdf.rds')

mydat <- prevDF %>%
  mutate(mPrev = log10(1000 * mPrev), fPrev = log10(1000 * fPrev)) %>%
  mutate(mPrev = ifelse(is.na(mPrev),0,mPrev),fPrev = ifelse(is.na(fPrev),0,fPrev),
         sexSpec = ((mPrev < 0.001) | (fPrev < 0.001)) & (nCases > 2000)) %>%
  mutate(sexSpecName = sapply(strsplit(ifelse ( sexSpec == T, Disease, NA),'[(]'),function(x)x[1])) 

sum(mydat$nCases > 2000)

xx = prevDF %>%
  mutate(mPrev = ifelse(is.na(mPrev),0,mPrev),fPrev = ifelse(is.na(fPrev),0,fPrev))
sum(xx$fPrev>=0.001 & xx$mPrev>=0.001)

diseaseSel <- mydat %>%
  ggplot(aes(x = mPrev, y = fPrev)) +
  annotate('rect',xmin = log10(1), xmax = Inf, ymin = log10(1), ymax = Inf, fill = 'gray50', alpha = 0.5) +
  annotate('rect',xmin = -Inf, xmax = log10(1), ymin = log10(1), ymax = Inf, fill = sexcolors['Female'], alpha = 0.1) +
  annotate('rect',ymin = -Inf, ymax = log10(1), xmin = log10(1), xmax = Inf, fill = sexcolors['Male'], alpha = 0.1) +
  geom_vline(xintercept = log10(1), linetype = 'dashed', color = 'darkred') +
  geom_hline(yintercept = log10(1), linetype = 'dashed', color = 'darkred') +
  geom_point(aes(alpha = c('<2,000','>2,000')[1+(nCases > 2000)])) +
  scale_alpha_manual('# Cases', values = c(0.2,1)) +
  geom_smooth(method = 'lm') +
  scale_x_continuous(breaks = log10(c(1,10,100,300)), labels = c(1,10,50,300)) +
  scale_y_continuous(breaks = log10(c(1,10,100,300)), labels = c(1,10,50,300)) +
  xlab('Number of Cases in 1,000 Males') +
  ylab('Number of Cass in 1,000 Females') +
  annotate('text',x = log10(0.003), y = log10(300), label = 'Female-specific', color = sexcolors['Female'], hjust = 0, size = 8 / pntnorm) +
  annotate('text',y = log10(0.003), x = log10(1.5), label = 'Male-specific', color = sexcolors['Male'], hjust = 0, size = 8 / pntnorm) +
  annotate('text',x = log10(1.5), y = log10(300), label = 'Disease prevalence>=0.001\nin both males and females', color = 'black', hjust = 0, size = 8 / pntnorm) +
  annotate('text',y = log10(0.003), x = log10(0.003), label = 'Rare diseases', color = 'black', alpha = 0.5, hjust = 0, size = 8 / pntnorm) +
  coord_equal(xlim = c(log10(0.003),log10(400)),ylim = c(log10(0.003),log10(400))) +
  theme(legend.key.size = unit(0.1,'cm'))

ggsave('./results/UKBB_disease_EDA/disFiltering.pdf', diseaseSel,  useDingbats=F, height=12, width = 18, units = 'cm')
ggsave('./results/UKBB_disease_EDA/disFiltering.png', diseaseSel, height=12, width = 18, units = 'cm')

