### SV panel generation
# 01. only auto-chromosome.
bcftools view --min-ac 1 -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 -o filter.1.vcf.gz -Oz Sus-pan.vcf.gz
tabix filter.1.vcf.gz
# 02. remove alternative alleles that are not covered by any haplotype
bcftools view --trim-alt-alleles filter.1.vcf.gz > filter.2.vcf
# 03. Bubbles with more than 20% missing allele across haplotypes were discarded
bcftools filter --exclude 'F_MISSING > 0.2' filter.2.vcf -o filter.3.vcf
# 04. Only bubbles containing at least one SV with a length greater than 50 bp were kept
python /disk192/miaoj/PanGenome/PanCode/extract_sv_variants.py -i filter.3.vcf -o filtered.vcf
# 05. multi-allelic VCF to bi-allelic VCF
python3 annotate_vcf.py -vcf filtered.vcf -gfa Sus-pan.gfa.gz -o anno &> anno.log
cat anno.vcf | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' > SVanno.used.vcf
cat anno_biallelic.vcf | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' > SVanno_biallelic.used.vcf
rm anno.vcf anno_biallelic.vcf

### SV hotspot analysis
# extract SV
python3 extract_sv_variants.py -i SVanno_biallelic.used.vcf -o SV_biallelic.vcf
# vcf2bed
python3{
with open("SV_biallelic.vcf") as f, open("SV_biallelic.bed", "w") as fo:
	for line in f:
		if line.startswith("#"):
			continue
		else:
			cols = line.split("\t")
			CHR = cols[0]
			START = int(cols[1])
			reflen = len(cols[3])
			END = START + reflen - 1
			fo.write("\t".join([str(i) for i in [CHR, START, END]]) + "\n")
}
# hotspot
R{
library(primatR)
library(rtracklayer)
bed_file = "SV_biallelic.bed"
bed_granges <- import(bed_file, format="bed")
hotspotter(gr, bw=20000, pval = 1e-08, num.trial = 100)
}

### Accuracy of PanGenie
# PanGenie call SV
REF=Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
vcf=SVanno.used.vcf
bi_panel=SVanno_biallelic.used.vcf
PanGenie -i <(zcat ${SAMPLE}.qc_1.fastq.gz ${SAMPLE}.qc_2.fastq.gz) -r ${REF} -v ${vcf} -o ${SAMPLE} -j 30 -t 19 2> ${SAMPLE}.log &
cat ${SAMPLE}_genotyping.vcf | python3 convert-to-biallelic.py ${bi_panel} > ${SAMPLE}_biallelic.vcf
python extract_sv_variants.py -i ${SAMPLE}_biallelic.vcf -o ${SAMPLE}_biallelic.SV.vcf
bcftools view --min-ac 1 ${SAMPLE}_biallelic.SV.vcf > ${SAMPLE}_biallelic.SV.filter.vcf
cat ${SAMPLE}_biallelic.SV.filter.vcf | python prepare-for-truvari.py | bgzip -c > ${SAMPLE}.vcf.gz
tabix ${SAMPLE}.vcf.gz
# PBSV call SV
ngmlr -t 10 -r Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -q ${SAMPLE}_hifi.fastq.gz -o ${SAMPLE}.ngmlr.bam
samtools view -@ 1 -bS ${SAMPLE}.ngmlr.bam | samtools sort -@ 1 -o ${SAMPLE}.ngmlr.sort.bam -
samtools index ${SAMPLE}.ngmlr.sort.bam
pbsv discover ${SAMPLE}.ngmlr.sort.bam {SAMPLE}.svsig.gz
pbsv call Sus_scrofa.Sscrofa11.1.dna.toplevel.fa {SAMPLE}.svsig.gz ${SAMPLE}.pbsv.vcf
bcftools view -i 'SVTYPE!="DUP"' -o JH.pbsv.filter.noDup.vcf.gz -Oz JH.pbsv.filter.vcf.gz
# truvari
truvari bench -b bench/JH.pbsv.filter.noDup.vcf.gz -c ${SAMPLE}.vcf.gz -f Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -o Tru_pbsv/ -r 2000 -C 2000 --dup-to-ins -t --passonly &

