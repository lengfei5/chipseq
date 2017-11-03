#Genome="/groups/bell/jiwang/Genomes/C_elegans/ce10/ce10_index_4Bowtie/ce10"
Genome="/groups/bell/jiwang/Genomes/Mouse/mm10_UCSC/mm10_index_4Bowtie/mm10"
nb_cores=6
DIR_input="${PWD}/ngs_raw/FASTQs"
DIR_output="${PWD}/alignments/BAMs_All"
mkdir -p $DIR_output
mkdir -p $PWD/logs
#mkdir -p ${DIR_output}/logs

for file in ${DIR_input}/*.fastq;
do
    fname="$(basename $file)"
    #echo $fname
    fname=${fname%.fastq}
    echo $fname
    qsub -q public.q -o ${PWD}/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N mapping2genome "module load samtools/0.1.18;module load bowtie2/2.2.4;bowtie2 -q -p $nb_cores -x $Genome -U $file | samtools view -bSu - | samtools sort - ${DIR_output}/$fname; samtools index ${DIR_output}/$fname.bam " 
done