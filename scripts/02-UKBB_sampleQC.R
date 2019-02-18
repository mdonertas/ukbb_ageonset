# libraries
library(tidyverse)
library(data.table)
library(RFRlib)

# definition of the fields
fields <- read_tsv('./data/raw/ukbb/helperfiles/field.txt') %>%
  select(title, field_id) %>%
  rename(description = title) %>%
  mutate(field_id = as.character(field_id))

# extract the field IDs from UKBB data
ukbfields <- colnames(read_csv('./data/raw/ukbb/ukb10904.csv', n_max = 0,
                               col_names = T))
ukbfields <- data.frame(fields = ukbfields)

# separate field IDs into an ID, visit and a number 
# to merge with field definitions
ukbfields <- ukbfields %>%
  separate(fields, into = c('field_id', 'visit', 'num'), remove=F) %>%
  left_join(fields)

# modify eid row so that it has the same structure as others
ukbfields[which(ukbfields$fields == 'eid'),] <- c('eid', 'eid', 0, 0, 'eid')

# create a folder to store helper files
system('mkdir -p ./data/processed/ukbb/helperfiles')

# write the field definitions into a file 
write_tsv(ukbfields,'./data/processed/ukbb/helperfiles/ukbfields.tsv')  

## Exclusions

# fields to decide exclusions
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
dim(ukb_exc)
# [1] 502617     11
withdrawn <- read_tsv('./data/raw/ukbb/withdrawn',col_names = 'eid')
ukb_exc <- ukb_exc[! ukb_exc$eid %in% withdrawn$eid,]
dim(ukb_exc)
# [1] 502543     11

sum(is.na(ukb_exc$Heterozygosity))
# 14248
ukb_exc <- ukb_exc[!is.na(ukb_exc$Heterozygosity),]
nrow(ukb_exc)
# 488295

### Discordant sex
disc_sex <- ukb_exc$eid[which(ukb_exc$Sex != ukb_exc$`Genetic sex`)]
table(ukb_exc$Sex, ukb_exc$`Genetic sex`)
# 0      1
# 0 264623    143
# 1    235 223363
length(disc_sex)
# 378
mean(ukb_exc$eid %in% disc_sex) * 100
# [1] 0.07741222

### Sex chr aneuploidy
sexchr_aneup <- ukb_exc$eid[which(!is.na(ukb_exc$`Sex chromosome aneuploidy`))]
length(sexchr_aneup)
# 651
mean(ukb_exc$eid %in% sexchr_aneup) * 100
# [1] 0.133321
length(intersect(disc_sex,sexchr_aneup))
# 181
mean(disc_sex %in% sexchr_aneup) * 100
# [1] 47.8836
mean(sexchr_aneup %in% disc_sex) * 100
# [1] 27.80338

sex_exc <- union(disc_sex, sexchr_aneup)
length(sex_exc)
# 848

system('mkdir -p ./data/processed/ukbb/sampleQC')

sex_exc_df <- data.frame(eid = sex_exc)
write_tsv(sex_exc_df, './data/processed/ukbb/sampleQC/discSex_exclusions.tsv')

## Recommended genomic analysis exclusions:
# poor heterozygosity/missingness
poorheterozyg_or_missing <- ukb_exc$eid[which(!is.na(ukb_exc$`Recommended genomic analysis exclusions`))]
length(poorheterozyg_or_missing)
# 469

## Genetic relatedness exclusions
# Participant self-declared as having a mixed ancestral background
mixed_anc <- ukb_exc$eid[which(ukb_exc$`Genetic relatedness exclusions` == 1)]
length(mixed_anc)
# 692
# High heterozygosity rate (after correcting for ancestry) or high missing rate
highHet_or_missing <- ukb_exc$eid[which(ukb_exc$`Genetic relatedness exclusions` == 2)]
length(highHet_or_missing)
# 840

##Outliers for heterozygosity or missing rate

HetoMiss_outliers <- ukb_exc$eid[which(!is.na(ukb_exc$`Outliers for heterozygosity or missing rate`))]
length(HetoMiss_outliers)
# 968

all_exc <- list(`Recommended genomic analysis exclusions` = poorheterozyg_or_missing,
                `Participant self-declared as having a mixed ancestral background` = mixed_anc,
                `High heterozygosity rate (after correcting for ancestry) or high missing rate` = highHet_or_missing,
                `Outliers for heterozygosity or missing rate` = HetoMiss_outliers,
                `Discordant Sex` = sex_exc)
saveRDS(all_exc,'./data/processed/ukbb/all_exclusions.rds')
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
ggsave('./results/SampleQC/Overlap_ExcCat.png', width = 17, device = 'png')
recomm_exc <- unique(unname(unlist(all_exc)))
length(recomm_exc)
# 3697

library(pheatmap)
newnames <- setNames(c('Rec. Exclusions','Mixed Ancestry','High hetero / missing','Hetero / missing outliers','Discordant Sex'),colnames(exc_overlap))
exc_overlap2 <- exc_overlap
colnames(exc_overlap2) <- unname(newnames[colnames(exc_overlap2)])
rownames(exc_overlap2) <- unname(newnames[rownames(exc_overlap2)])
pheatmap(exc_overlap2, cellwidth = 25, cellheight = 25, number_format = "%.2f",
         display_numbers = T, color = RColorBrewer::brewer.pal(8,'Oranges')[-8],
         file = './results/SampleQC/Overlap_ExcCat_Heatmap.pdf')
pheatmap(exc_overlap2, cellwidth = 25, cellheight = 25, number_format = "%.2f",
         display_numbers = T, color = RColorBrewer::brewer.pal(8,'Oranges')[-8],
         file = './results/SampleQC/Overlap_ExcCat_Heatmap.png')

mean(HetoMiss_outliers %in% poorheterozyg_or_missing)
mean(poorheterozyg_or_missing %in% HetoMiss_outliers)

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
ggsave('./results/SampleQC/Het_Miss_Plots.png', p1, width = 15, height = 10,
       device = 'png')  

## Exclude all of these

recomm_exc_df <- data.frame(eid = recomm_exc)
write_tsv(recomm_exc_df, 
          './data/processed/ukbb/sampleQC/recommended_exclusions.tsv')

ukb_exc <- ukb_exc[!ukb_exc$eid %in% recomm_exc,]
dim(ukb_exc)
# [1] 484598     11

# samples to include in the analysis:
samplePassQC <- data.frame(eid = ukb_exc$eid)
system('mkdir -p data/processed/ukbb/sampleQC')
write_tsv(samplePassQC, './data/processed/ukbb/sampleQC/sampleIDs_passedQC.tsv')
