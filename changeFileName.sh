#######
## This script is to copy file and change file's name
## this should be done before mapping to the genome
#######
#DIR_Input=$PWD/ngs_raw/RAW_bams
DIR_Input="$PWD/alignments/BAMs_All"
DIR_OUT="$PWD/alignments/BAMs_All"
cwd=$PWD;
PARAM="${PWD}/Params_HagarCarin_R4889.txt"

mkdir -p $DIR_OUT
cd $DIR_Input
i=1;
while read -r line;
do 
    IFS=, read -r "ID" "CONDITION" "SAMPLE"  <<< "$line"
    if [[ $i -gt 1 ]]; then
	#echo $ID;
	old=`ls -l *.bai| grep $ID | awk '{print $9}'`
	if [ ! -z $old ]; then
	    #extension=`echo "$old" | cut -d'.' -f2`
	    #echo $extension;
	    extension="bam.bai"
	    new=${CONDITION}_${SAMPLE}_${ID}.${extension};
	    echo $old;
	    echo $new
	    #cp "$old" $DIR_OUT
	    mv "$old" "$new"
	    #mv "$new" $DIR_OUT
	fi
    fi

    i=$((i + 1));
    #IFS=$','
done < "$PARAM"
cd $cwd;
