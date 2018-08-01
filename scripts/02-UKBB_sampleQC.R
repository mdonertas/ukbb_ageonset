library(tidyverse)
library(data.table)
library(RFRlib)

fields <- read_tsv('./data/raw/ukbb/helperfiles/field.txt') %>%
  select(title, field_id) %>%
  rename(description = title) %>%
  mutate(field_id = as.character(field_id))

ukbfields <- colnames(read_csv('./data/raw/ukbb/ukb10904.csv', n_max = 0,
                               col_names = T))
ukbfields <- data.frame(fields = ukbfields)

ukbfields <- ukbfields %>%
  separate(fields, into = c('field_id', 'visit', 'num'), remove=F) %>%
  left_join(fields)

ukbfields[which(ukbfields$fields == 'eid'),] <- c('eid', 'eid', 0, 0, 'eid')
system('mkdir -p ./data/processed/ukbb/helperfiles')
  
write_tsv(ukbfields,'./data/processed/ukbb/helperfiles/ukbfields.tsv')  

## Exclusions

exc_fields <- ukbfields %>%
  filter(description %in% c('eid', 'Recommended genomic analysis exclusions',
                            'Genetic relatedness exclusions',
                            'Sex chromosome aneuploidy', 'Sex', 'Genetic sex',
                            'Missingness', 'Heterozygosity',
                            'Heterozygosity, PCA corrected',
                            'Outliers for heterozygosity or missing rate', 
                            'Ethnic background')) %>%
  mutate(fields = as.character(fields)) %>%
  filter(visit == 0)

ukb_exc <- setDF(fread('./data/raw/ukbb/ukb10904.csv',
                       select = exc_fields$fields,
                       col.names = exc_fields$description))
sum(is.na(ukb_exc$Heterozygosity))
# 14253
ukb_exc <- ukb_exc[!is.na(ukb_exc$Heterozygosity),]
nrow(ukb_exc)
# 488364
### Discordant sex
disc_sex <- ukb_exc$eid[which(ukb_exc$Sex != ukb_exc$`Genetic sex`)]
table(ukb_exc$Sex, ukb_exc$`Genetic sex`)
# 0      1
# 0 264623    143
# 1    235 223363
length(disc_sex)
# 378
mean(ukb_exc$eid %in% disc_sex) * 100
# [1] 0.07740128

### Sex chr aneuploidy
sexchr_aneup <- ukb_exc$eid[which(!is.na(ukb_exc$`Sex chromosome aneuploidy`))]
length(sexchr_aneup)
# 652
mean(ukb_exc$eid %in% sexchr_aneup) * 100
# [1] 0.133507
length(intersect(disc_sex,sexchr_aneup))
# 181
mean(disc_sex %in% sexchr_aneup) * 100
# [1] 47.8836
mean(sexchr_aneup %in% disc_sex) * 100
# [1] 27.76074

sex_exc <- union(disc_sex, sexchr_aneup)
length(sex_exc)
# 849

sex_exc_df <- data.frame(eid = sex_exc)
write_tsv(sex_exc_df, './data/processed/ukbb/sampleQC/discSex_exclusions.tsv')

ukb_exc <- ukb_exc[!ukb_exc$eid %in% sex_exc,]
dim(ukb_exc)
# [1] 487515     11

## Recommended genomic analysis exclusions:
# poor heterozygosity/missingness
poorheterozyg_or_missing <- ukb_exc$eid[which(!is.na(ukb_exc$`Recommended genomic analysis exclusions`))]
length(poorheterozyg_or_missing)
# 468

## Genetic relatedness exclusions
# Participant self-declared as having a mixed ancestral background
mixed_anc <- ukb_exc$eid[which(ukb_exc$`Genetic relatedness exclusions` == 1)]
length(mixed_anc)
# 689
# High heterozygosity rate (after correcting for ancestry) or high missing rate
highHet_or_missing <- ukb_exc$eid[which(ukb_exc$`Genetic relatedness exclusions` == 2)]
length(highHet_or_missing)
# 838

##Outliers for heterozygosity or missing rate

HetoMiss_outliers <- ukb_exc$eid[which(!is.na(ukb_exc$`Outliers for heterozygosity or missing rate`))]
length(HetoMiss_outliers)
# 963

all_exc <- list(`Recommended genomic analysis exclusions` = poorheterozyg_or_missing,
             `Participant self-declared as having a mixed ancestral background` = mixed_anc,
             `High heterozygosity rate (after correcting for ancestry) or high missing rate` = highHet_or_missing,
             `Outliers for heterozygosity or missing rate` = HetoMiss_outliers)
exc_overlap <- sapply(all_exc,function(x) {
  sapply(all_exc, function(y) {
    round(mean(x %in% y) * 100, 2)
  })
})
diag(exc_overlap) <- 0
system('mkdir -p results/SampleQC')
exc_overlapdf <- as.data.frame(exc_overlap)
exc_overlapdf$Name <- rownames(exc_overlapdf)
exc_overlapdf %>%
  gather(Type, value, -Name) %>%
  ggplot(aes(x = Name, y = Type, color = value))+
  geom_label(aes(label = value), size = 10)+
  theme_rfr(base_size = 20)+
  theme(axis.text.x = element_blank())+
  ggtitle('Percent Overlap Across Exclusion Categories')+
  ylab('')+ xlab('')+ guides(color=F)
ggsave('./results/SampleQC/Overlap_ExcCat.pdf', useDingbats=F, width = 17)
recomm_exc <- unique(unname(unlist(all_exc)))
length(recomm_exc)
# 2849

## Heterozygosity
eth_coding <- read_tsv('./data/raw/ukbb/datacoding/coding1001.tsv')
ethparent <- setNames(eth_coding$parent_id, eth_coding$coding)
ethparent[ethparent == 0] <- names(which(ethparent == 0))
ethname <- setNames(eth_coding$meaning, eth_coding$coding)
p1 <- ukb_exc %>%
  select(`eid`, `Heterozygosity, PCA corrected`,
         `Ethnic background`, `Missingness`) %>%
  mutate(logit_missingness = car::logit(Missingness)) %>%
  mutate(EthnicBackground_major = as.factor(ethname[as.character(unname(ethparent[as.character(`Ethnic background`)]))])) %>%
  mutate(EthnicBackground_minor = as.factor(ethname[as.character(`Ethnic background`)])) %>%
  mutate(RecommendedExc = factor(c('-', 'Exclude')[1 + (eid %in% recomm_exc)],
                                 levels=c('-', 'Exclude'))) %>%
  ggplot(aes(x = logit_missingness, y = `Heterozygosity, PCA corrected`,
             color = RecommendedExc))+
  geom_point(size = 0.4)+
  facet_wrap(~EthnicBackground_major)+
  theme_rfr()+
  geom_vline(xintercept = car::logit(0.05), linetype = 'dashed', 
             color = 'darkred')+
  geom_hline(yintercept = mean(ukb_exc$`Heterozygosity, PCA corrected`), 
             linetype = 'dashed', color = 'gray50')+
  scale_color_manual(values = c("#99999933",'gray20'))
ggsave('./results/SampleQC/Het_Miss_Plots.pdf', p1, width = 15, height = 10,
       useDingbats = F)  

## Exclude all of these

recomm_exc_df <- data.frame(eid = recomm_exc)
write_tsv(recomm_exc_df, 
          './data/processed/ukbb/sampleQC/recommended_exclusions.tsv')

ukb_exc <- ukb_exc[!ukb_exc$eid %in% recomm_exc,]
dim(ukb_exc)
# [1] 484666     13

# samples to include in the analysis:
samplePassQC <- data.frame(eid = ukb_exc$eid)
system('mkdir -p data/processed/ukbb/sampleQC')
write_tsv(samplePassQC, './data/processed/ukbb/sampleQC/sampleIDs_passedQC.tsv')
