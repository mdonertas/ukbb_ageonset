mkdir -p ./data/processed/ukbb/gwas/fam4bolt/
for i in {1..22}
do
echo ukb30688_cal_chr$i\_v2_s488346.fam
awk '{print $1, $2, $3, $4, $5, 1}' ./data/raw/ukbb/fam/ukb30688_cal_chr$i\_v2_s488346.fam > ./data/processed/ukbb/gwas/fam4bolt/ukb30688_cal_chr$i\_v2_s488346.fam
done
