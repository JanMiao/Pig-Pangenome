### Mapping performance comparison of reference graph and pangenome
# reference graph
REF=Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
vg autoindex --workflow giraffe -r ${REF} -p ref

# fastp
fastp -i ${SAMPLE}_1.fq.gz -I ${SAMPLE}_2.fq.gz -o ${SAMPLE}_R1.fq.gz -O ${SAMPLE}_R2.fq.gz

# giraffe
vg giraffe -Z ${prefix}.giraffe.gbz -m ${prefix}.min -d ${prefix}.dist -f ${SAMPLE}_R1.fq.gz -f ${SAMPLE}_R2.fq.gzz > ${SAMPLE}.gam
vg stats -a ${SAMPLE}.gam > ${SAMPLE}.stats