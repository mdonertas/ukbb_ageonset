mkdir -p /nfs/research1/thornton/ukbb_ageonset/data/raw/ukbb/fam4bolt
for i in {1..22}
do
echo ukb30688_cal_chr$i\_v2_s488346.fam
awk '{print $1, $2, $3, $4, $5, 1}' /nfs/research1/thornton/ukbb_ageonset/data/raw/ukbb/fam/ukb30688_cal_chr$i\_v2_s488346.fam > /nfs/research1/thornton/ukbb_ageonset/data/raw/ukbb/fam4bolt/ukb30688_cal_chr$i\_v2_s488346.fam
done
