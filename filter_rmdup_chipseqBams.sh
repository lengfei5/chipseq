##########################################
# This script is to filter bam and find statistics for Bam files: 
# number of total, mapped, unique, after duplication removed reads
# Here the duplication removal is using picard that is taking really big memory 18-30G
# which is a little bit weird
# to adapt the script for axolotl genome, the duplication removal is done by samtools tools
# > samtools/1.10-foss-2018b in which the samtools index works
# 
#########################################
while getopts ":hp" opts; do
    case "$opts" in
        "h")
	    echo "script to filter low quality mapped reads, remove duplicates for bams files for chip-seq or alike data (e.g. atac-seq)"
            echo "old version used samtools to remove duplicates; and the current version is using PICARD"
	    echo 'lastest version employed samtools again, becasue PICARD does not work for axolotl genome'
	    echo "Usage:"
	    echo "$0 (single_end bam)"
            echo "$0 -p (paired_end bam)"
	    exit 0
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

## mm10 or ce10
nb_cores=4
#blacklist="/groups/bell/jiwang/Genomes/C_elegans/ce10/ce10_blacklist/ce10-blacklist.bed" 
MAPQ_cutoff=30

DIR_input="${PWD}/alignments/BAMs_All"
DIR_uniq="${PWD}/alignments/BAMs_unique"
DIR_uniq_rmdup="${PWD}/alignments/BAMs_unique_rmdup"
DIR_stat="${PWD}/QCs/ALL/BamStat"
DIR_logs=${PWD}/logs

mkdir -p $DIR_logs
mkdir -p $DIR_uniq
mkdir -p $DIR_uniq_rmdup
mkdir -p $DIR_stat
mkdir -p $DIR_stat/logs

jobName='bamFiltering'
#picardPath=$EBROOTPICARD
echo $picardPath
for file in ${DIR_input}/*.bam;
do
    echo $file
    ff="$(basename $file)"
    ff="${ff%.bam}"
    fname=$ff;
    
    #newb="${newb%.bam}_filter"
    #echo $newb
    newb=$DIR_uniq/${ff}_uniq
    newbb=$DIR_uniq_rmdup/${ff}_uniq_rmdup
    stat=$DIR_stat/${ff}.txt
    #echo $newb 
    #echo $newbb 
    #echo $stat
    picardDup_QC=${DIR_stat}/${ff}_picardDup.qc.txt

    # creat the script for each sample
    script=$DIR_logs/${fname}_${jobName}.sh
    cat <<EOF > $script
#!/usr/bin/bash

#SBATCH --cpus-per-task=$nb_cores
#SBATCH --time=60
#SBATCH --mem=30000
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o $DIR_logs/${fname}.out
#SBATCH -e $DIR_logs/${fname}.err
#SBATCH --job-name $jobName

module load samtools/0.1.20-foss-2018b
module load bedtools/2.25.0-foss-2018b
module load picard/2.18.27-java-1.8

# filter low quality reads 
if [ ! -e $newb.bam ]; then
   if [ $PAIRED != "TRUE" ]; then
      samtools view -q $MAPQ_cutoff -b $file | samtools sort - $newb;
   else
      samtools view -F 1804 -q $MAPQ_cutoff -b $file | samtools sort - $newb;
   fi  
fi; 

if [ ! -e $newb.bam.bai ]; then 
   samtools index $newb.bam; 
fi;

# remove duplicates 
if [ ! -e $newbb.bam ]; then
   java -jar \$EBROOTPICARD/picard.jar MarkDuplicates INPUT=$newb.bam OUTPUT=$newbb.bam METRICS_FILE=$picardDup_QC \
ASSUME_SORTED=true REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.2 MAX_RECORDS_IN_RAM=250000 
fi;

if [ ! -e $newbb.bam.bai ]; then 
   samtools index $newbb.bam; 
fi;

echo 'done duplication removal...'

## save statistical number for each bam'
total=\$(samtools view -c $file); 
mapped=\$(samtools view -c -F 4 $file); 
unique=\$(samtools view -c $newb.bam); 
rmdup=\$(samtools view -c $newbb.bam); 
echo \"$ff \$total \$mapped \$unique \$rmdup\"|tr ' ' '\t' > $stat 

EOF
    cat $script
    sbatch $script
    #break; 
done
