### Population analysis of SV and SNP
# vcf to PLINK
plink2 --vcf SNP.raw.vcf.gz --make-pgen --out SNP.raw &
plink2 --vcf SV.raw.vcf.gz --make-pgen --out SV.raw &
# MAF filter
/disk192/miaoj/bin/plink2 --pfile SNP.raw --make-pgen --maf 0.05 --out SNP &
/disk192/miaoj/bin/plink2 --pfile SV.raw --make-pgen --maf 0.01 --out SV &


# PCA for SNP
plink2 --pfile SNP --set-missing-var-ids @:#\$r,\$a --make-pgen --out SNP_tmp
plink2 --pfile SNP_tmp --indep-pairwise 50 10 0.2
plink2 --pfile SNP_tmp --extract plink2.prune.in --make-bed --out SNPprune
gcta64 --bfile SNPprune --make-grm --make-grm-alg 1 --out kinship
gcta64 --grm kinship --pca 10 --out SNP.pc

# PCA for SV
plink2 --pfile SV --set-missing-var-ids @:#\$r,\$a --make-pgen --out SV_tmp
plink2 --pfile SV_tmp --indep-pairwise 50 10 0.2
plink2 --pfile SV_tmp --extract plink2.prune.in --make-bed --out SVprune
gcta64 --bfile SVprune --make-grm --make-grm-alg 1 --out SV.kinship
gcta64 --grm SV.kinship --pca 10 --out SV.pc

# admixture
for K in {2..10}
do 
admixture --cv ../SNPprune.bed ${K} -j30 | tee snp.log${K}.out
admixture --cv ../SVprune.bed ${K} -j30 | tee sv.log${K}.out
done

# LD between SNP and SV
bcftools concat SNP.vcf.gz SV.vcf.gz -Oz -a -o merged.vcf.gz
plink2 --vcf merged.vcf.gz --make-pgen --out merged
plink --vcf merged.vcf.gz --double-id --ld-snp-list SV.id --r2 --ld-window 999999999 --ld-window-kb 50 --ld-window-r2 0.2 --out SVcor0.2 &
