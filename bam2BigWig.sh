#Genome="/home/imp/jingkui.wang/Genomes/C_elegans/ce10_index/ce10/ce10"
nb_cores=2
bams="${PWD}/alignments/BAMs_All"
OUT="${PWD}/bigWigNorm"
#DIR_rmdup="${PWD}/alignment/BAMs_All"
#echo $bams
#echo $OUT

mkdir -p $OUT
mkdir -p ${PWD}/logs

#chromSize="/groups/bell/jiwang/Genomes/C_elegans/ce10/ce10_sequence/ce10_chrom.sizes"
#chromSize="/groups/bell/jiwang/Genomes/Mouse/mm10_UCSC/Sequence/mm10_chrom_sizes.sizes"

for b in ${bams}/*.bam;
#for bam in ${DIR_input}/42964_filter_rmdup.bam
do
    echo $b
    wig="$(basename $b)"
    wig="${wig%.bam}"
    echo $wig
    qsub -q public.q -o ${PWD}/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N bamcoverage "if [ ! -e ${OUT}/${wig}.bw ]; then module unload deeptools; module load python/2.7.3; module load pysam/0.10.0; module load deeptools/2.2.3-python2.7.3;      bamCoverage -b ${b} -o $OUT/${wig}.bw --outFileFormat=bigwig --normalizeUsingRPKM --ignoreDuplicates --minMappingQuality 10; fi;"
    #qsub -q public.q -oqsub -q public.q -o ${PWD}/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N bam2bw "module load samtools; module load rseqc/2.6.1-python2.7.3; module load kent-ucsc/2.79; TMP=`mktemp`; bam2wig.py -i $bam -s $chromSize -o \$TMP -t 200000000; if [ ! -e \$TMP.bw ]; then wigToBigWig \$TMP.wig $chromSize \$TMP.bw; fi; mv \$TMP.bw $DIR_output/$wig.bw; rm \$TMP*; " 
    #samtools view -q 10 -b $file | bedtools intersect -a stdin -b /groups/bell/jiwang/Genomes/C_elegans/ce10/ce10_blacklist/ce10-blacklist.bed -v | samtools sort - ${DIR_output}/$newb; samtools index ${DIR_output}/${newb}.bam; samtools rmdup -s ${DIR_output}/${newb}.bam ${DIR_r     mdup}/${newb}_rmdup.bam;samtools index ${DIR_rmdup}/${newb}_rmdup.bam;" 
    #break
done