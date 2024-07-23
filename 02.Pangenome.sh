### MC Pangenome
# 1. generate graph
nohup cactus-minigraph ./jobstore samples.txt Sus.gfa.gz --reference Sus11 --mapCores 20 --workDir ./temp &> Minigraph.map.log &
# 2. remap to the graph
nohup cactus-graphmap ./jobstore samples.txt Sus.gfa.gz Sus.paf --outputFasta Sus.gfa.fa  --reference Sus11 --mapCores 30 --workDir ./temp &> Minigraph.remap.log &
# 3. split chromosome
cactus-graphmap-split <jobStore> <seqFile> <inputGFA> <inputPAF> --reference --outDir
nohup cactus-graphmap-split ./jobstore samples.txt Sus.gfa.gz Sus.paf --outDir chroms --otherContig chrOther --refContigs $(for i in `seq 18`; do echo $i; done ; echo "X") --reference Sus11 --nodeStorage 10000 --betaInertia 0 --targetTime 1 --workDir ./temp &> Chroms.log &
# 4. cactus align
cp chroms/chromfile.txt .
nohup cactus-align-batch ./jobstore ./chromfile.txt chroms \
--alignCores 32 --maxCores 64 --alignOptions "--pangenome --maxLen 10000 --reference Sus11   --outVG" \
--defaultPreemptable --workDir ./temp &> Cactus.align.log &
# 5. build index
nohup cactus-graphmap-join ./jobstore --vg $(for j in $(for i in `seq 18`; do echo $i; done ; echo "X chrOther"); do echo chroms/${j}.vg; done) \
--hal $(for j in $(for i in `seq 18`; do echo $i; done ;echo "X chrOther"); do echo chroms/${j}.hal; done) \
--outDir Index/ --outName Sus-pan --reference Sus11 --giraffe clip filter --filter 2 --vcf --gbz --gfa --indexCores 32  --maxCores 64 &> Cactus.index1.log &

### Graph stats
gunzip -d -c Sus-pan.gfa.gz > Sus-pan.gfa
vg stats -zNElL Sus-pan.gfa > MC.stats
# non-reference nodes
graph_gfa=Sus.gfa
echo "nodeid Length Source Postion NonReference" > graph_len.tsv
awk '$1~/S/ { split($5,chr,":"); split($6,pos,":"); split($7,arr,":");print $2,length($3),chr[3],pos[3],arr[3] }' ${graph_gfa} >> graph_len.tsv
# variant stats
rtg vcfstats /disk191_4/miaoj/cactus/PPG3/Index/Sus-pan.raw.vcf.gz > Sus-pan.raw.txt
