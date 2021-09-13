##########################################
# This script is to filter bam and find statistics for Bam files: 
# number of total, mapped, unique, after duplication removed reads
# Here the duplication removal is using picard that is taking really big memory 18-30G
# which is a little bit weird
# to adapt the script for axolotl genome, the duplication removal is done by samtools tools
# > samtools/1.10-foss-2018b in which the samtools index works
# 
#########################################
while getopts ":hD:" opts; do
    case "$opts" in
        "h")
	    echo "script to filter low quality mapped reads, remove duplicates for bams files for chip-seq or alike data (e.g. atac-seq)"
            echo "old version used samtools to remove duplicates; and the current version is using PICARD"
	    echo 'lastest version employed samtools again, becasue PICARD does not work for axolotl genome'
	    echo "Usage:"
	    echo "$0 (input in alignments/BAMs_All)"
            echo "$0 -D out/merged_bams"
	    exit 0
	    ;;
	"D")
	    DIR_input="$OPTARG"
	    ;;
	"?")
            echo "Unknown option $opts"
            ;;
        ":")
            echo "No argument value for option $opts"
            ;;
    esac
done

nb_cores=8
## mm10 or ce10
#blacklist="/groups/bell/jiwang/Genomes/C_elegans/ce10/ce10_blacklist/ce10-blacklist.bed" 
MAPQ_cutoff=30

if [ -z "$DIR_input" ]; then
    DIR_input="${PWD}/alignments/BAMs_All"
    echo "bam directory is $DIR_input"

else
    # check if provided directory ends with slash
    if [[ $DIR_input == */ ]]; then
        echo "parsing bam directory "
        DIR_input=${DIR_input%/}

    fi

    echo "bam directory is $DIR_input"
    # check if there are bam files in the directory or not
    if [ `ls -l $DIR_input/*.bam|wc -l` == "1" ]; then
        echo "no bam Found !!! "
        echo "wrong directory or no bams in there "
        exit 1;
    fi
fi


DIR_rmdup="${PWD}/alignments/BAMs_uniq"
DIR_rmdup_uniq="${PWD}/alignments/BAMs_uniq_rmdup"
DIR_stat="${PWD}/QCs/BamStat"
DIR_insertion=${PWD}/QCs/inserion_sizes
DIR_logs=${PWD}/logs

mkdir -p $DIR_logs
mkdir -p $DIR_rmdup
mkdir -p $DIR_rmdup_uniq
mkdir -p $DIR_stat
mkdir -p $DIR_insertion
mkdir -p $DIR_stat/logs

jobName='bam_filter_rmdup'
#picardPath=$EBROOTPICARD
#echo $picardPath

for file in ${DIR_input}/*.bam
#for file in ${DIR_input}/Embryo_Stage44_distal_93323_sorted.bam
do
    echo $file
    ff="$(basename $file)"
    ff="${ff%.bam}"
    fname=$ff;
    
    #echo $newb
    newb=$DIR_rmdup/${ff}_uniq
    newbb=$DIR_rmdup_uniq/${ff}_uniq_rmdup
    stat=${DIR_stat}/${ff}_stat.txt
    insertion_size=${DIR_insertion}/${ff}_insertion.size
    picardDup_QC=${DIR_stat}/${ff}_picardDup.qc.txt
    
    #echo $newb 
    #echo $newbb 
    #echo $stat
    
    # creat the script for each sample
    script=$DIR_logs/${fname}_${jobName}.sh
    cat <<EOF > $script
#!/usr/bin/bash

#SBATCH --cpus-per-task=$nb_cores
#SBATCH --time=300
#SBATCH --mem=24G

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o $DIR_logs/${fname}.out
#SBATCH -e $DIR_logs/${fname}.err
#SBATCH --job-name $jobName

module load samtools/1.10-foss-2018b
ml load picard/2.20.6--0-biocontainers
# filter low quality reads
#if [ $PAIRED != "TRUE" ]; then
#      samtools view -h -q $MAPQ_cutoff -b $file | samtools sort - $newb;
#   else
#      samtools view -h -F 1804 -q $MAPQ_cutoff -b $file | samtools sort - $newb;
      #samtools view -h -q 30 ${sample}.bam > ${sample}.rmMulti.bam
#fi
  
samtools view -h -q 30 $file > ${newb}.unsorted.bam
samtools sort -@ $nb_cores -o ${newb}.bam ${newb}.unsorted.bam
samtools index -c -m 14 ${newb}.bam 
rm ${newb}.unsorted.bam

# remove duplicates 
if [ ! -e $newbb.bam ]; then
   picard MarkDuplicates INPUT=$newb.bam OUTPUT=${newbb}.unsorted.bam METRICS_FILE=$picardDup_QC \
ASSUME_SORTED=true REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.2
fi;

if [ ! -e $newbb.bam.bai ]; then
   samtools sort -@ $nb_cores -o ${newbb}.bam ${newbb}.unsorted.bam
   samtools index -c -m 14 $newbb.bam;
   rm ${newbb}.unsorted.bam

fi;


echo 'done duplication removal...'

# save statistical number for each bam'
echo 'sample total mapped uniq rmdup_uniq' |tr ' ' '\t' > $stat 

total=\$(samtools view -@ $nb_cores -c $file); 
mapped=\$(samtools view -@ $nb_cores -c -F 4 $file); 
uniq=\$(samtools view -@ $nb_cores -c $newb.bam); 
uniq_rmdup=\$(samtools view -@ $nb_cores -c $newbb.bam); 
echo $ff \$total \$mapped \$uniq \$uniq_rmdup |tr ' ' '\t' >> $stat 

samtools stat -@ $nb_cores ${file} | grep ^IS | cut -f 2- > ${insertion_size}.txt
samtools stat -@ $nb_cores ${newbb}.bam | grep ^IS | cut -f 2- > ${insertion_size}_uniq.rmdup.txt

EOF
    
    #cat $script
    sbatch $script
    #break;
    
done
