###########
## aim of script: map fastq (or fq) file to the selected gennome (mm9, mm10, ce11 or customized genomes already indexed) using bowtie2  
###########
while getopts ":hg:pc:" opts; do
    case "$opts" in
        "h") 
	    echo "script to map fastq or fq to the genome using bowite2"
	    echo "available genomes: ce10, ce11, mm10, mm9, hg19 and mm9_G11D (mm9 genome with G11 reporter) "
	    echo "Usage: "
	    echo "$0 -g ce11 -c 4 (single_end fastq aligned to ce11 and 2 cpus required)"
	    echo "$0 -g mm10 -c 8 (single_end fastq aligned to mm10 and 8 cpus required)"
	    echo "$0 -g mm10 -p (paired_end fastq aligned to mm10 and default number (6) cpus required)"
	    exit 0
	    ;;
	"g")
            genome="$OPTARG"
            ;;
	"p")
	    PAIRED="TRUE"
	    ;;
	"c")
	    nb_cores="$OPTARG"
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

if [ -z "$nb_cores" ]; then
    nb_cores=4
fi

DIR_input="${PWD}/ngs_raw/FASTQs"
DIR_output="${PWD}/alignments/BAMs_All"
DIR_logs=${PWD}/logs
jobName='bowtie2'

mkdir -p $DIR_output
mkdir -p $DIR_logs

if [ "$PAIRED" != "TRUE" ]; then
    echo "single-end fastq..."

    for file in ${DIR_input}/*.fastq;
    do
	fname="$(basename $file)"
        #echo $fname
	fname=${fname%.fastq}
	echo $fname
	if [ ! -e "${DIR_output}/$fname.bam" ]; then
	    # creat the script for each sample
	    script=$DIR_logs/${fname}_${jobName}.sh
	    cat <<EOF > $script
#!/usr/bin/bash	

#SBATCH --cpus-per-task=$nb_cores
#SBATCH --time=120
#SBATCH --mem=8000
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o $DIR_logs/${fname}.out
#SBATCH -e $DIR_logs/${fname}.err
#SBATCH --job-name $jobName
module load samtools/0.1.20-foss-2018b
module bowtie2/2.3.4.2-foss-2018b
bowtie2 -q -p $nb_cores -x $Genome -U $file | samtools view -bSu - | samtools sort - ${DIR_output}/${fname}; 
samtools index ${DIR_output}/${fname}.bam;

EOF
	    cat $script
	    #sbatch $script
	    
	fi
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
	
	if [ ! -e "${DIR_output}/$fname.bam" ]; then
	    # creat the script for each sample
	    script=$DIR_logs/${fname}_${jobName}.sh
	    cat <<EOF > $script
#!/usr/bin/bash	

#SBATCH --cpus-per-task=$nb_cores
#SBATCH --time=120
#SBATCH --mem=8000
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o $DIR_logs/${fname}.out
#SBATCH -e $DIR_logs/${fname}.err
#SBATCH --job-name $jobName
module load samtools/0.1.20-foss-2018b;
module load bowtie2/2.3.4.2-foss-2018b

bowtie2 -q -p $nb_cores --no-mixed -X 2000 -x $Genome -1 $seq1 -2 $seq2  | samtools view -bSu - | \
samtools sort - ${DIR_output}/$fname; 
samtools index ${DIR_output}/$fname.bam;

EOF
	    cat $script
	    sbatch $script

	fi
    done
fi;
