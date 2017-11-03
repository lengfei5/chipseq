###################
## bam2fastq script
###################
CORES=2
DIR_INPUT="${PWD}/ngs_raw/BAMs"
DIR_OUT="${PWD}/ngs_raw/FASTQs"
echo ${DIR_OUT};
mkdir -p ${DIR_OUT}
mkdir -p $PWD/logs

for file in $DIR_INPUT/*.bam;
do 
    echo "$file"
    #echo "$(dirname $file)"
    #echo "$(basename $file)"
    FILENAME="$(basename $file)";
    fname=${FILENAME%.bam};
    echo $fname
    qsub -q public.q -o $PWD/logs/ -j yes -N bam2fastq -cwd -b y -pe smp $CORES -shell y "module load bedtools; bamToFastq -i $file -fq ${DIR_OUT}/${fname}.fastq" 
done