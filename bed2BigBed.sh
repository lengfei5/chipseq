###################
### This script is to convert bed file to BigBed file that can be visualized in track hubs in USCS genome brownser
###################
#if [ "$#" != "2" ]; then
#    echo -e "Usage of bed2BigBed : \t$0 bedfile(s) OutputDir \n "
#    exit 1;
#else
bedfiles="$@";
#DIR_output="$2";

#echo $bedfiles
#echo $DIR_output

#mkdir -p $DIR_output 
#mkdir -p ${PWD}/logs

chromSize="/groups/bell/jiwang/Genomes/C_elegans/ce11/ce11_sequence/ce11.chrom.sizes"

module load kent-ucsc/2.79;
for bed in $bedfiles;
do
    #echo $bed
    name="$(basename $bed)"
    name="${name%.bed}"
    echo $name
    cut -f 1-4 $bed | sort -k1,1 -k2,2n > ${name}.sorted.bed
    bedToBigBed ${name}.sorted.bed $chromSize ${name}.bb
    rm ${name}.sorted.bed
    
done
    
#fi