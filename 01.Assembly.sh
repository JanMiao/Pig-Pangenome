### Assembly
# without HiC
hifiasm -o TC.asm -t16 $HIFI
# with HiC
hifiasm -o TC.asm -t16 --h1 $HIC1 --h2 $HIC2 $HIFI

### Quality Control
# remove contaminated contigs using NT database
python splitLongContig.py ${SAMPLE}-0.fa ${SAMPLE}-0.split.fa
blastn -db /disk191_2/miaoj/nt/nt -query ${SAMPLE}-0.split.fa -out ${SAMPLE}-0_NT.txt -num_threads 30 \
-outfmt '6 qseqid qlen sseqid slen pident qcovs evalue bitscore staxid' -max_target_seqs 1 -max_hsps 1 2> ${SAMPLE}-0_NT.log
python NT_add_topLevel.py ${SAMPLE}-0_NT.txt ${SAMPLE}-0_NT_TopLevel.txt
python RMcontaminationContigs.py ${SAMPLE}-0_NT_TopLevel.txt ${SAMPLE}-0.fa ${SAMPLE}-0.qc.fa ${SAMPLE}-0.contamination
# repeats mask
RepeatMasker -xsmall -species "Sus scrofa" -pa 40 -poly -html -gff -dir ${SAMPLE}-0 ${SAMPLE}-0.qc.fa
# Inspector
inspector.py -c ${SAMPLE}-0.qc.fa.masked -r ${SAMPLE}_hifi.fastq.gz -o ${SAMPLE}-0.asm/ --datatype hifi -t 64
inspector-correct.py -i ${SAMPLE}-0.asm/ --datatype pacbio-hifi -o ${SAMPLE}-0.asm.corrected/ --skip_structural -t 64
inspector.py -c ${SAMPLE}-0.asm.corrected/contig_corrected.fa -r ${SAMPLE}_hifi.fastq.gz -o ${SAMPLE}-0.cor.asm/ --datatype hifi -t 64

### Assembly stats
# stats
seqkit stats -a ${SAMPLE}-0.asm.corrected/contig_corrected.fa
# BUSCO
busco -i ${SAMPLE}-0.asm.corrected/contig_corrected.fa -c 30 -o ${SAMPLE}-0.new -m geno -l ${BUSCO_DIR}/mammalia_odb10



