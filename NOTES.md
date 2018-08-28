# BOLT-LMM

Notes are from [BOLT-LMM manual](https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/BOLT-LMM_v2.3.2_manual.pdf)

## Download

Download using 'https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/BOLT-LMM_v2.3.2.tar.gz'

## Install

```
tar -xvzf BOLT-LMM_v2.3.2.tar.gz
```
This package containes a standalone executable - thus, no installation step is required. This executable is now included in my path variable, thus I can simply use 'bolt' to run it.

## Files for BOLT-LMM

### Genotype file 

PLINK binary format (bed/bim/fam).
If all has the same prefix: **--bfile=prefix** otherwise use **--bed**, **--bim** and **--fam**. 

#### fam file

Need to fix sixth column.

### Reference genetic maps

Reference maps can be used to interpolate genetic map coordinates from SNP physical (base pair) positions in the event that PLINK bim file does not contain genetic coordinates (in units of Morgans). For UKBB data we need to use a reference genetic map. Depending on the genome build used in UKBB reference map can be selected from ../melike/soft/BOLT-LMM_v2.3.2/tables/genetic_map_hg#_withX.txt.gz 

```
--geneticMapFile=../melike/soft/BOLT-LMM_v2.3.2/tables/genetic_map_hg#_withX.txt.gz
```

### Imputed SNP dosages

When testing imputed SNPs, BOLT-LMM first performs its usual model-fitting on PLINK-format genotypes (supplied via --bfile or bed/bim/fam) and then applies the model to scan any provided imputed SNPs. The second step requires only a modest amount of additional computation and no additional RAM, as the it simply performs a genome scan of real-valued dosage SNPs against the residual phenotypes that BOLT-LMM computes during model-fitting. We currently recommend performing model-fitting on ~500K hard-called genotypes; this approach should sacrifice almost no statistical power while retaining computational efficiency. 

#### Imputed SNPs in BGEN format

To compute association statistics at SNPs in one or more BGEN data files, specify the .bgen file(s) with --bgenFile and the corresponding .sample file with --sampleFile. The --bgenMinMAF and --bgenMinINFO options allows limiting output to SNPs passing minimum allele frequency and INFO thresholds. (Note: the --bgenMinMAF filter is applied to the full BGEN file before any sample exclusions, whereas the MAF reported in BOLT-LMM’s output is computed in the subset of samples actually analyzed. Some SNPs may therefore pass the --bgenMinMAF filter but have lower reported MAF in the output file; if you wish to exclude such SNPs, you will need to post-process the results.)

Additionally, starting with BOLT-LMM v2.3.2, you may alternatively specify a list of whitespaceseparated .bgen / .sample file pairs using the --bgenSampleFileList option (instead of using the --bgenFile and --sampleFile options). This option enables analyses of data sets in which different BGEN files have different sample sets (e.g., the UK Biobank v3 imputation release; section 9.1).

```
--bgenSampleFileList
--bgenMinMAF
--bgenMinINFO
```

### Phenotypes

**--phenoFile** and **--phenoCol**: Whitespace-delimited file, specified with --phenoFile, with the first line containing column headers and subsequent lines containing records one per individual.

The first two columns bust be **FID** and **IID** (the PLINK identifiers of an individual). Any number of columns may follow; the column containing the phenotype to analyse is specified with **--phenoCol**. 

Values of -9 and NA are interpreted as missing data. All other values in the column should be numeric. The records in lines following the header line need not be in sorted order and need not match the fam file. BOLT-LMM will analyse only the individuals in the intersection of the genotype and phenotype files and will output a warning if these sets do not match. 

### Covariates 

Covariate data may be specified in a file **--covarFile** with the same format as the phenotype file. Same file for both phenotypes and covarietes can be used, but both --phenoFile and --covarFile should be specified. Each covariate must be specified using --covarCol (categorical covariates) and --qCovarCol (quantitative covariates) options. Categorical covariate values are allowed to be any text strings not containing whitespace; each unique text string in a column corresponds to a category. BOLT-LMM throws an error if a categorical covariate contains more than 10 distinct values, but the upper bound can be modified with --covarMaxLevels. Quantitative covariate values must be numeric with the exception of NA. In either case, both -9 and NA are interpreted as missing data. If groups of covariates of the same type are numbered sequentially, they may be specified using array shorthand (e.g. --qCovarCol=PC{1:10} for columns PC1, PC2, ... PC10).

### Missing data treatment

Individuals with missing phenotypes are ignored. By default, individuals with any missing covariates are also ignored - this approach is commonly used and referred to as 'complete case analysis'. As an alternative, 'missing indicator method' is possible using '--covarUseMissingIndic' option, which adds indicator variables demarcating missing status as additional covariates. 

Missing genotypes in plink data (--bfile or bed/bim/fam) are replaced with per-SNP averages. Imputed genotypes should not contain missing data; standard imputation software always produces genotype probability estimates even if uncertainty is high. 

### Genotype QC

BOLT-LMM automatically filter SNPs and individuals with missing rates exceeding thresholds of 0.1. These thresholds may be modified using --maxMissingPerSnp and --maxMissingPerIndiv. Note that filtering is not performed based on MAF or deviation from Hardy-Weinberg. Allele frequency and missingness of each SNP are included in the BOLT-LMM association test output, however, and we recommend checking these values and HW p-values when following up on significant associations. 

### User-specified filters

Individuals to remove from the analysis may be specified in one or more --remove files listing FID and IIDs (one individual per line). Similarly, SNPs to exclude from the analysis may be specified in one or more --exclude files listing SNP IDs (typically rs numbers).

Note that --exclude filters are not applied to imputed data; exclusions of specific imputed SNPs will need to be performed separately as a post-processing step. 

## Association analysis 

### Mixed model association tests 

#### BOLT-LMM
Association test on residuals from Bayesian modeling using a mixture-of-normals prior on SNP effect sizes. This approach can fit 'non-infinitesimal' traits with loci having moderate to large effects, allowing increased association power.

#### BOLT-LMM-inf
Standard (infinitesimal) mixed model association. This statistic approximates the standard approach used by eigendecomposition based software 

### Options for mixed model association test

1. *--lmm*: Performs default BOLT-LMM analysis, which consists of (1a) estimating heritability parameters, (1b) computing the BOLT-LMM-inf statistic, (2a) estimating Gaussian mixture parameters, and (2b) computing the BOLT-LMM statistic only if an increase in power is expected. If BOLT-LMM determines based on cross-validation that the non-infinitesimal model is likely to yield no increase in power, the BOLT-LMM (Bayesian) mixed model statistic is not computed.

2. *--lmmInfOnly*: Computes only infinitesimal mixed model association statistics (i.e., steps 1a and 1b).

3. *--lmmForceNonInf*: Computes both the BOLT-LMM-inf and BOLT-LMM statistics regardless of whether or not an increase in power is expected from the latter.


#### Reference LD score tables

A table of reference LD scores is needed to calibrate the BOLT-LMM statistic. Reference LD scores appropriate for analyses of Eurpoean-ancestry samples are provided in the tables subdirectory of the software. 

For analyses of non-European data, LD scores can be calculated using LDSC software on an ancestry-matched subset of the 1000 Genomes samples. 

By default, LD scores in the table are matched to SNPs in the PLINK data by rsIDs. The --LDscoresMatchBp option allows matching SNPs by base pair coordinate. 

#### Restricting SNPs used in the mixed model

If millions of SNPs are available from imputation, we suggest including at most 1 million SNPs at a time in the mixed model (using the --modelSnps option) when performing association analysis. Using an LD pruned set of at most 1 million SNPs should achieve near-optimal power and correction for confounding while reducing computational cost and improving convergence. Note that even when a file of --modelSnps is specified, all SNPs in the genotype data are still tested for association; only the random effects in the mixed model are restricted to the --modelSnps. Also note that BOLT-LMM automatically performs leave-one-chromosome-out (LOCO) analysis, leaving out SNPs from the chromosome containing the SNP being tested in order to avoid proximal contamination [4, 9].

### Standard linear regression

Setting the --verboseStats flag will output standard linear regression chi-square statistics and p-values in additional output columns CHISQ_LINREG and P_LINREG. Note that unlike mixed model association, linear regression is susceptible to population stratification, so you may wish to include principal components (computed using other software, e.g., PLINK2 or FastPCA [18] in EIGENSOFT v6.0+) as covariates when performing linear regression. Including PCs as covariates will also speed up convergence of BOLT-LMM’s mixed model computations.

## UK Biobank analysis

Example command line

```
./bolt \
--bed=ukb_cal_chr{1:22}_v2.bed \
--bim=ukb_snp_chr{1:22}_v2.bim \
--fam=ukb1404_cal_chr1_v2_CURRENT.fixCol6.fam \
--remove=bolt.in_plink_but_not_imputed.FID_IID.976.txt \
--remove=sampleQC/remove.nonWhite.FID_IID.txt \
--exclude=snpQC/autosome_maf_lt_0.001.txt \
--exclude=snpQC/autosome_missing_gt_0.1.txt \
--phenoFile=ukb4777.phenotypes.tab \
--phenoCol=height \
--covarFile=ukb4777.covars.tab.gz \
--covarCol=cov_ASSESS_CENTER \
--covarCol=cov_GENO_ARRAY \
--covarMaxLevels=30 \
--qCovarCol=cov_AGE \
--qCovarCol=cov_AGE_SQ \
--qCovarCol=PC{1:20} \
--LDscoresFile=tables/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=tables/genetic_map_hg19.txt.gz \
--lmmForceNonInf \
--numThreads=8 \
--statsFile=bolt_460K_selfRepWhite.height.stats.gz \
--bgenFile=ukb_imp_chr{1:22}_v2.bgen \
--bgenMinMAF=1e-3 \
--bgenMinINFO=0.3 \
--sampleFile=ukb1404_imp_chr1_v2_s487406.sample \
--statsFileBgenSnps=bolt_460K_selfRepWhite.height.bgen.stats.gz \
--verboseStats
```

