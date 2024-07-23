### eQTL mapping

# RNAseq
fastq-dump ${SRAnumber} --gzip --split-3 -O ./data
gzip data/${SAMPLE}_1.fastq
gzip data/${SAMPLE}_2.fastq
# fastp
fastp -i ${SAMPLE}_1.fastq.gz -I ${SAMPLE}_2.fastq.gz -o ${SAMPLE}_1.fastqc.gz -O ${SAMPLE}_2.fastqc.gz -w 16
# STAR index
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir star \
--genomeFastaFiles Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
--sjdbGTFfile Sus_scrofa.Sscrofa11.1.109.gtf --sjdbOverhang 149 &
# STAR mapping
STAR --runThreadN 20 --genomeDir star --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
--readFilesCommand zcat --outFilterMismatchNmax 3 --readFilesIn ${SAMPLE}_1.fastqc.gz ${SAMPLE}_2.fastqc.gz --outFileNamePrefix ${SAMPLE}
# stringtie
stringtie -p 20 -e -B -G Sus_scrofa.Sscrofa11.1.109.gtf -o ${SAMPLE}.gtf -A ${SAMPLE}.tsv ${SAMPLE}Aligned.sortedByCoord.out.bam
# featurecounts
featureCounts -T 20 -p -t exon -g gene_id -a Sus_scrofa.Sscrofa11.1.109.gtf -o ${SAMPLE}.featureCounts.txt ${SAMPLE}Aligned.sortedByCoord.out.bam

# preprocess for eQTL mapping
# TSS position
python{
from io1 import gtf_to_tss_bed
annotation_gtf = "Sus_scrofa.Sscrofa11.1.109.gtf"
tss = gtf_to_tss_bed(annotation_gtf, feature='gene', exclude_chrs=[], phenotype_id='gene_id')
tss.to_csv("Sus11TSS.bed", index=False, sep="\t")
}
# preprocess for expression data
R{
library(edgeR)
library(RNOmni)
# remove samples with RNA, without WGS
snp.sample = scan("../../genotype/SNP.sample", what="character")
info = read.table("../../preprocess/sampleInfo.txt", header=T)
all.sample = unique(info$wgsID)
miss.sample = setdiff(all.sample, snp.sample) 
info1 = info[info$wgsID != miss.sample, ]

tissue = "M"
featuerCounts = read.table("../featureCounts.csv", header=T, sep=",")
tpm = read.table("../tpm_matirx.csv", header=T, sep=",")
if (!identical(featuerCounts[,1], tpm[,1])){print("######## ERROR ##########")}
# 2. remove samples without WGS
sample.info = read.table("../../preprocess/sampleInfo.txt", header=T)
rnaID.rm = sample.info[sample.info[,4] == "SAMEA111521213",1]
featuerCounts = featuerCounts[,!(colnames(featuerCounts) %in% rnaID.rm)]
tpm = tpm[,!(colnames(tpm) %in% rnaID.rm)]
sample.info = sample.info[sample.info[,4] != "SAMEA111521213",]
# 3. extract WGS according the samples in RNA
sample.tmp = sample.info[sample.info$tissue==tissue,]
rownames(sample.tmp) = sample.tmp[,4]
rnaID = sample.tmp[snp.sample,1]
featuerCounts.temp = featuerCounts[,c("GeneID", rnaID)]
tpm.temp = tpm[,c("GeneID",rnaID)]
# 4. replace sampleID in RNAseq
rnaID2wgsID = c(sample.tmp[,4], "GeneID")
names(rnaID2wgsID) = c(sample.tmp[,1], "GeneID")
colnames(featuerCounts.temp) = rnaID2wgsID[colnames(featuerCounts.temp)]
colnames(tpm.temp) = rnaID2wgsID[colnames(tpm.temp)]
rownames(featuerCounts.temp) = featuerCounts.temp[,1]
featuerCounts.temp = featuerCounts.temp[,2:ncol(featuerCounts.temp)]
rownames(tpm.temp) = tpm.temp[,1]
tpm.temp = tpm.temp[,2:ncol(tpm.temp)]
# 5. gene filter
count_threshold = 6
tpm_threshold = 0.1
sample_frac_threshold = 0.2
nsamples = ncol(tpm.temp)
tpm_th = rowSums(tpm.temp >= tpm_threshold)
count_th = rowSums(featuerCounts.temp >= count_threshold)
ctrl1 = tpm_th >= (sample_frac_threshold * nsamples)
ctrl2 = count_th >= (sample_frac_threshold * nsamples)
mask = ctrl1 & ctrl2
# 6. genrate TMM exp matrix
expr = DGEList(counts=featuerCounts.temp)
y = calcNormFactors(expr, method="TMM")
TMM = cpm(y,normalized.lib.sizes=T)
TMM_pass = TMM[mask,]
#expression values (TMM) were inverse normal transformed across samples.
inverse_normal_transform = function(x) {
    qnorm(rank(x) / (length(x)+1))
}
TMM_inv = t(apply(TMM_pass, MARGIN = 1, FUN = inverse_normal_transform))
# 7. genrate bed file for QTLtools
tss = read.table("../Sus11TSS.bed", header=T)
tss = tss[tss$gene_id %in% rownames(TMM_inv),]
colnames(tss) = c("#Chr", "start", "end", "pid")
tss$gid = "."
tss$strand = "+"
gene_expr = TMM_inv[tss$pid,]
bed.out = cbind(tss, gene_expr)
write.table(bed.out, paste0("tissue",tissue,".bed"), row.names=F, sep="\t", quote=F)
### COV file
cov.info = read.table("../preprocess/eQTLsampleInfo.txt", header=F)
sample2id = sample.tmp[,4]
names(sample2id) = sample.tmp[,2]
cov.info[,1] = sample2id[cov.info[,1]]
cov.out = rbind(c("id", "month", "sex", "breed"), cov.info)
colnames(cov.info) = c("id", "month", "sex", "breed")
cov.out = t(cov.out)
write.table(cov.out, "cov.txt", row.names=F, sep="\t", quote=F, col.names=F)
}
# add PCs from SNPs
plink2 --vcf SNP.FAANG.filter.vcf.gz --make-bed --out SNP
gcta64 --bfile SNP --autosome-num 18 --make-grm --make-grm-alg 1 --out kinship
gcta64 --grm kinship --pca 20 --out pc

# QTLtools for eQTL mapping
# permutation
cov=cov.pc.txt
vcf=merged.FAANG.vcf.gz
bed=tissueM.bed.gz
# permutation
for j in $(seq 1 10)
do
     echo "QTLtools cis  --vcf ${vcf} --bed ${bed} --cov ${cov} \
--permute 1000 --normal --seed 123456 --silent --chunk $j 10 --log M_permute_${j}.SV.log --out M_permute_${j}.SV.txt &" >> M_permutate.sh
done
nohup bash M_permutate.sh &> M_permutate.log &
# FDR correction
cat M_permute_*.SV.txt | gzip -c > M_permute.txt.gz
Rscript runFDR_cis.R M_permute.txt.gz 0.05 M_permute_all
# conditional pass
for j in $(seq 1 10)
do
     echo "QTLtools cis  --vcf ${vcf} --bed ${bed} --cov ${cov} \
--mapping M_permute_all.thresholds.txt --silent --chunk $j 10 --log M_conditional_${j}.log --out M_conditional_${j}.txt &" >> M_condition.sh
done
bash M_condition.sh &> M_condition.log &
cat M_conditional_*.txt > M_conditional_full.txt
# normal pass
QTLtools cis --vcf ${vcf} --bed ${bed} --cov ${cov}\
 --nominal 1 --normal --silent --log M_nominals.log --out M_nominals.txt &
# extract significant variants in normal pass
python SummyNominalPass.py -n M_nominals.txt -t M_permute_all.thresholds.txt -p M_nominal &
