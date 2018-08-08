# LOG

## 30.07.2018

1. Download data using ./scripts/01-downloadData.Rmd

## 01.08.2018

### Sample QC

processed using ./scripts/02-UKBB_sampleQC.R

1. Decode UKBB field names to field ID, visit, and number, and include the definition of field IDs. Mapping file is saved as './data/processed/ukbb/helperfiles/ukbfields.tsv'.
2. We first filtered out samples without genotypes (14253 eids, resulting in 488,364 eids).
3. Sample filtering is composed of 2 steps:
    1. Discordant sex (remaining: 487,515 samples)
    2. Genotype call rate & Heterozygosity (remaining: 484,666 samples).
    
Sample IDs that passed the QC process are saved as './data/processed/ukbb/sampleQC/sampleIDs_passedQC.tsv'.


#### Discordant sex

There are two columns for sex: 1) self reported sex and 2) genetic sex determined using the call intensities on sex chromosomes. There are multiple reasons why these two information do not correspond, such as sample mishandling, errors in data input, transgender individuals and sex-chromosome aneuploidies [1, 2]. Since we will probably need to stratify GWAS results by gender (because diseases show different characteristics in female and male patients), we preferred to be cautious about this issue and eliminated all cases where the genetic sex and self reported sex do not correspond and all cases where sex chromosome aneuploidy is detected. Specifically,

1. We used fields 31-0.0 (Sex) and 22001-0.0 (Genetic sex) to compile discordant information. 
    * Self reported sex: Male, Genetic sex: Female => 235
    * Self reported sex: Female, Genetic sex: Male => 143
    * In total: 378 cases, making 0.077% of all data
2. We used field 22019-0.0 (Sex chromosome aneuploidy) to exclude cases with sex chromosome aneuploidy. 
    * There were 652 cases of aneuploidy, making 0.134% of all data
    * 181 of these cases (27.76% of aneuploidy cases) are also detected as discordant information in the first step. This corresponds to 47.88% of discordant sex cases.
    
In total these two cases yielded 849 unique samples to be excluded. After the exclusion, the number of samples decreased to 487,515. The list of sample ids excluded from the analysis are saved as './data/processed/ukbb/sampleQC/discSex_exclusions.tsv'.

#### Genotype call rate & Heterozygosity

For the exclusions based on missingness and heterozygosity we only used the suggested exclusions by UK Biobank. Specifically,

1. We used the field 22010-0.0 (Recommended genomic analysis exclusions) and noted the cases with 'poor heterozygosity / missingness' => 468 samples
2. We used the field 22018-0.0 (Genetic relatedness exclusions)
    * 689 Cases with 'Participant self-declared as having a mixed ancestral background'
    * 838 Cases with 'High heterozygosity rate (after correcting for ancestry) or high missing rate'
3. We used the field 22027-0.0 (Outliers for heterozygosity or missing rate) - 963 cases

Before exclusion we checked the overlap across these categories since they are all related to missingness & heterozygosity, and this file is saved as './results/SampleQC/Overlap_ExcCat.pdf'. Overall, we saw the exclusion suggestions based on these criteria do not overlap. We then checked the scatter plots for logit(Missingness) vs. Heterozygosity for each Ethnic Background, in accordance with the identification of samples to exclude by the UK Biobank [2]. This plot is saved as
'./results/SampleQC/Het_Miss_Plots.pdf'. We confirmed these are in accordance with the article and we excluded the samples suggested by these fields, in total 2849 unique samples. Resulting in 484,666 samples. 2849 sample ids excluded by these criteria are saved as './data/processed/ukbb/sampleQC/recommended_exclusions.tsv'.

## 03.08.2018

### Prepare trait data for EDA

Processed using ./scripts/03-prepTraitData.R

Only values entered for the baseline visit are considered. Samples failing the initial QC are eliminated from the analysis.

#### Traits in UKBB

A file containing the following variables (only baseline measurements) are saved as './data/processed/traits_clean/traitData_baseline.tsv':
'eid', 'Sex', 'Age at recruitment', 'Age when attended assessment centre', 'Age at death', 'Standing height', 'Weight', 'Number of self-reported cancers', 'Number of self-reported non-cancer illnesses', 'Number of operations, self-reported', 'Number of treatments/medications taken', 'Sleep duration', 'Facial ageing', 'Maternal smoking around birth', 'Smoking status', 'Alcohol drinker status', "Father's age at death", "Mother's age at death", 'Overall health rating', 'Health satisfaction', 'Age when periods started (menarche)', 'Age at menopause (last menstrual period)', 'Non-accidental death in close genetic family'

#### Self Reported Non-Cancer Illnesses and Age of Diagnosis:

A file containing 'eid', 'diseaseID', 'Disease', 'node_id', 'parent_id', 'selectable', 'Age' is saved under './data/processed/traits_clean/SRdisease_baseline.tsv'.

#### Self Reported Cancers and Age of Diagnosis:

A file containing 'eid', 'cancerID', 'Cancer', 'node_id', 'parent_id', 'selectable', 'Age' is saved under './data/processed/traits_clean/SRcancer_baseline.tsv'.

## 06.08.2018

### Exploratory Data Analysis on UKBB Data

Done using ./scripts/04-UKBB_EDA.Rmd

## 07.08.2018

### Exploratory Data Analysis for the self-reported diseases

Done using ./scripts/05-UKBB_SRDiseases.Rmd

## 08.08.2018

### Exploratory Data Analysis for the self-reported diseases

Updated ./scripts/05-UKBB_SRDiseases.Rmd to account for the hierarchical relationship between diseases.

# References

[1] Anderson, C. A., Pettersson, F. H., Clarke, G. M., Cardon, L. R., Morris, A. P., & Zondervan, K. T. (2010). Data quality control in genetic case-control association studies. Nature Protocols, 5(9), 1564–1573. https://doi.org/10.1038/nprot.2010.116

[2] Bycroft, C., Freeman, C., Petkova, D., Band, G., Elliott, L. T., Sharp, K., … Marchini, J. (2017). Genome-wide genetic data on ~500,000 UK Biobank participants. BioRxiv, https://doi.org/10.1101/166298. https://doi.org/10.1101/166298
