#####################################################################
#
# this script is to control the quality of chip-seq
# currently spp and deeptools were used to calculated metrics specifically 
# for ChIP-seq or ATAC-seq data  
#
#####################################################################
nb_cores=2
#cwd=`pwd`
DIR_cwd=$PWD
DIR_bam="${PWD}/alignments/BAMs_All"
echo $DIR_bam

DIR_SPP=$PWD/QCs/SPP_QCs;
DIR_deeptools=$PWD/QCs/deepTools;

mkdir -p $PWD/logs
mkdir -p ${DIR_SPP}
mkdir -p ${DIR_deeptools}

## spp for quality control
for file in ${DIR_bam}/*.bam;
do
    #echo $file
    output="$(basename $file)"
    output=${output%.bam}
    #echo $output
    qsub -q public.q -o $DIR_cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N sppQC "module unload R; module load R/2.15.3; export R_LIBS="/groups/vbcf-ngs/bin/R/R-2.15.3"; Rscript --vanilla /groups/vbcf-ngs/bin/spp/spp_package/run_spp.R -c=$file -odir=${DIR_SPP} -savp -out=${DIR_SPP}/${output}.txt" 
done

# deeptools for quality control
module load samtools;
BAMs=""
label=""
for b in ${DIR_bam}/*.bam; do
    if [ ! -e ${b}.bai ]; then
	samtools index ${b};
	echo $b;
    fi
        
    BAMs="$BAMs $b"
    bamname="$(basename $b)"
    label="$label ${bamname%.bam}"
done
#echo $BAMs
label=`echo $label | sed -e 's/^[ \t]*//'`
echo $label

mkdir -p ${DIR_bam}/logs

echo "module load python/2.7.3; module load deeptools/2.2.3-python2.7.3; multiBamSummary bins --bamfiles $BAMs -out $DIR_deeptools/multiBamSummary.npz -bs 1000 --minMappingQuality 10 --ignoreDuplicates --outRawCounts ${DIR_deeptools}/readCounts.tab; plotCorrelation -in ${DIR_deeptools}/multiBamSummary.npz -o ${DIR_deeptools}/correlation_pearson.pdf --corMethod pearson -p heatmap --removeOutliers --outFileCorMatrix ${DIR_deeptools}/correlation.tab --colorMap RdYlBu --plotNumbers; plotFingerprint -b $BAMs -l ${label} -plot ${DIR_deeptools}/fingerprint.pdf" > ${DIR_bam}/logs/run_deeptools_2cluster.sh

## sometimes have problem to run the job in cluster because the path or filename exceed 1024 character
qsub -q public.q -o $DIR_cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N deeptoolsQC "bash ${DIR_bam}/logs/run_deeptools_2cluster.sh "
