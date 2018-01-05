###########
## aim of script: map fastq (or fq) file to the selected gennome (mm9, mm10, ce11 or customized genomes already indexed) using bowtie2  
###########
while getopts ":hg:p" opts; do
    case "$opts" in
        "h") 
	    echo "script to map fastq or fq to the genome using bowite2"
	    echo "Usage: $0 -g ce11 (single_end)"
	    echo "$0 -g ce11 -p (paired_end) "
	    exit 0
	    ;;
	"g")
            genome="$OPTARG"
            ;;
	"p")
	    PAIRED="TRUE"
	    ;;
        "?")
            echo "Unknown option $opts"
            ;;
        ":")
            echo "No argument value for option $opts"
            ;;
        esac
done

## select the indexed genome
case "$genome" in 
    "ce10")
	echo "align to ce10 "
	Genome="/groups/bell/jiwang/Genomes/C_elegans/ce10/ce10_index_4Bowtie/ce10"
	;;
    "ce11")
	echo "alignment to ce11 "
	Genome="/groups/bell/jiwang/Genomes/C_elegans/ce11/ce11_index_4Bowtie/ce11"
	;;
    
    "mm9")
	echo "align to mm9"
	Genome="/groups/bell/jiwang/Genomes/Mouse/mm9_UCSC/mm9_index_4Bowtie/mm9"
	;;
    "mm10")
	echo "alignment to mm10"
	Genome="/groups/bell/jiwang/Genomes/Mouse/mm10_UCSC/mm10_index_4Bowtie/mm10"
	;;
    "mm9_G11D")
	echo "align to customized mm9_G11D"
	Genome="/groups/bell/jiwang/Genomes/Mouse/mm9_UCSC_G11D_reproter/mm9_G11D_index/mm9_G11D"
	;;
    
    "hg19")
	echo "align to hg19"
	Genome="/groups/bell/jiwang/Genomes/Human/hg19/hg19_index_4Bowtie/hg19"
	;;
    *)
	echo " No indexed GENOME Found !! "
	echo "Usage: $0 -g mm10 "
	exit 1;
	;;
esac

DIR_input="${PWD}/ngs_raw/FASTQs"
DIR_output="${PWD}/alignments/BAMs_All"

nb_cores=12

mkdir -p $DIR_output
mkdir -p $PWD/logs

if [ "$PAIRED" != "TRUE" ]; then
    echo "single-end fastq..."

    for file in ${DIR_input}/*.fastq;
    do
	fname="$(basename $file)"
        #echo $fname
	fname=${fname%.fastq}
	echo $fname
	qsub -q public.q -o ${PWD}/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N bowtie2Map "module load samtools/0.1.18; module load bowtie2/2.2.4;bowtie2 -q -p $nb_cores -x $Genome -U $file | samtools view -bSu - | samtools sort - ${DIR_output}/$fname; samtools index ${DIR_output}/$fname.bam;" 
    done
else
    echo "paired_end fastq ..."
    
    for seq1 in `ls ${DIR_input}/*.fastq |grep R1 `;
    do
	fname="$(basename $seq1)"
        SUFFIX=${fname#*_R1}
	fname=${fname%_R1*}
	echo $SUFFIX;
	echo $fname;
	seq2=${DIR_input}/${fname}_R2${SUFFIX};
	echo $seq1 
	echo $seq2
	qsub -q public.q -o ${PWD}/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N bowtie2Map "module load samtools/0.1.18; module load bowtie2/2.2.4; bowtie2 -q -p $nb_cores --no-mixed -X 2000 -x $Genome -1 $seq1 -2 $seq2  | samtools view -bSu - | samtools sort - ${DIR_output}/$fname; samtools index ${DIR_output}/$fname.bam;" 
    done
fi;