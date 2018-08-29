mkdir -p ./data/processed/ukbb/gwas/hwe

for i in {1..22}
do
echo $i
bsub -M 8000 -o /dev/null plink2 --bfile ./data/raw/ukbb/genotypes/chr$i --hardy --out ./data/processed/ukbb/gwas/hwe/chr$i
done

mkdir -p ./data/processed/ukbb/gwas/freq

for i in {1..22}
do
echo $i
bsub -M 8000 -o /dev/null plink2 --bfile ./data/raw/ukbb/genotypes/chr$i --freq --out ./data/processed/ukbb/gwas/freq/chr$i
done

