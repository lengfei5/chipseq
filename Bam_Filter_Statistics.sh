###################
### This script is to filter bam and find statistics for Bam files: number of total, mapped, unique, after duplication removed reads
###################
## mm10 or ce10
Cores=2
blacklist="/groups/bell/jiwang/Genomes/C_elegans/ce10/ce10_blacklist/ce10-blacklist.bed" 
q_cutoff=10
DIR_input="${PWD}/alignments/BAMs_All"
DIR_uniq=${PWD}/alignments/BAMs_unique
DIR_uniq_rmdup=${PWD}/alignments/BAMs_unique_rmdup
DIR_stat="${PWD}/QCs/ALL/BamStat"
mkdir -p ${PWD}/logs
mkdir -p $DIR_uniq
mkdir -p $DIR_uniq_rmdup
mkdir -p $DIR_stat
mkdir -p $DIR_stat/logs
for file in ${DIR_input}/*.bam;
#for file in ${DIR_input}/515D10_H3K27me3_37157_8_C97P2ANXX.bam
do
    echo $file
    ff="$(basename $file)"
    ff="${ff%.bam}"
    #newb="${newb%.bam}_filter"
    #echo $newb
    newb=$DIR_uniq/${ff}_uniq
    newbb=$DIR_uniq_rmdup/${ff}_uniq_rmdup
    stat=$DIR_stat/${ff}.txt
    echo $newb 
    echo $newbb 
    echo $stat
    
    echo "module load samtools/0.1.18; 
module load bedtools/2.25.0; 
if [ ! -e $newb ]; then 
samtools view -q $q_cutoff -b $file | samtools sort - $newb; 
fi; 
if [ ! -e $newb.bam.bai ]; 
then 
samtools index $newb.bam; 
fi; 
if [ ! -e $newbb.bam ]; then 
samtools rmdup -s $newb.bam $newbb.bam; 
fi; 
if [ ! -e $newbb.bam.bai ]; then 
samtools index $newbb.bam; fi; 
total=\$(samtools view -c $file); 
mapped=\$(samtools view -c -F 4 $file); 
unique=\$(samtools view -c $newb.bam); 
rmdup=\$(samtools view -c $newbb.bam); 
echo \"$ff \$total \$mapped \$unique \$rmdup\"|tr ' ' '\t' > $stat" > ${DIR_stat}/logs/${ff}_bamstat_2cluster.sh
    
    qsub -q public.q -o ${PWD}/logs -j yes -pe smp $Cores -cwd -b y -shell y -N bamStat "bash ${DIR_stat}/logs/${ff}_bamstat_2cluster.sh "
done


