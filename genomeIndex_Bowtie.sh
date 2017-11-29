nb_cores=10;
cwd=`pwd`

if [ ! -z "$1" ] && [ ! -z "$2" ]; then
    mkdir -p ${cwd}/logs
    genome="$1"
    outName="$2"
    OUT=${outName%/*}
    echo "directory for output $OUT"
    mkdir -p $OUT;
    qsub -q public.q -o ${cwd}/logs  -j yes -pe smp $nb_cores -cwd -b y -shell y -N genomeIndex "module load bowtie2/2.2.4; bowtie2-build $genome $outName; "
else
    echo "........"
    echo "$0 genome outName "
    echo "........"
fi