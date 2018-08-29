for i in {1..22}
do
cp /nfs/research1/ukbb/500K_release_v3/genotypes/EGAD00010001497/ukb_cal_chr$i\_v2.bed.gz /nfs/research1/thornton/ukbb_ageonset/data/raw/ukbb/genotypes/chr$i.bed.gz
done


for i in {1..22}
do
cp /nfs/research1/ukbb/500K_release_v3/genotypes/EGAD00010001497/ukb_snp_chr$i\_v2.bim.gz /nfs/research1/thornton/ukbb_ageonset/data/raw/ukbb/genotypes/chr$i.bim.gz
done


for i in {1..22}
do
cp /nfs/research1/thornton/ukbb_ageonset/data/processed/ukbb/gwas/fam4bolt/ukb30688_cal_chr$i\_v2_s488346.fam /nfs/research1/thornton/ukbb_ageonset/data/raw/ukbb/genotypes/chr$i.fam
done

gunzip /nfs/research1/thornton/ukbb_ageonset/data/raw/ukbb/genotypes/*.gz