#####################
# this script is to call peaks for ChIP-seq (or ATAC-seq) data 
# using MACS2 (sharp or broad peaks), Sicer (broad peaks) 
# it is quite tricky to specify Input files if bam files need different Input
#####################
while getopts ":hD:I:mbg:" opts; do
    case "$opts" in
        "h")
            echo "script to call chipseq peak using macs2 (sharp or broad optioins and  sicer for broad peaks"
	    echo "options: "
	    echo "-D  (specfiy the ABSOLUTE directory of bam files or by default alignments/BAMs_All)"
	    echo "-I (input, e.g. /groups/bell/jiwang/Projects/Jorge/Analysis_ChIP_seq/INPUT/merged_Input_49475_49911_49908.bam for mouse mm10)" 
	    echo "no input file needed by defaut"
            echo "-g (genome, i.e. mm10, ce11, hg19, am6)"
	    echo "-m (macs2 sharp peaks)"
	    echo "-b (macs2 broad peaks)"
	    
	    echo "Usage: "
            echo "$0 -m -g mm10 (if bam files in alignments/BAMs_All for chipseq)"
            echo "$0 -D /clustertmp/Jorge_R6118/alignments/BAMs_All_GCc_merged -m -g mm10 (bam files in XXX directory) "
            echo "$0 -b -g mm10 (call mm10 peaks with macs2 broad option) "
	    
            exit 0
            ;;
        "D")
            DIR_Bams="$OPTARG";
            ;;
        "I")
            INPUT="$OPTARG";
            ;;
	"m")
	    MACS2="TRUE"
	    ;;
	"f")
	    format="$OPTARG";
	    ;;
	"b")
	    MACS2_broad="TRUE"
	    ;;
	"g")
	    genome="$OPTARG";
	    ;;
        "?")
            echo "Unknown option $opts"
            ;;
        ":")
            echo "No argument value for option -D "
            exit 1;
            ;;
    esac
done

OUT=$PWD/calledPeaks

nb_cores=8;
cwd=`pwd`;

if [ -z "$DIR_Bams" ]; then
    DIR_Bams="$PWD/alignments/BAMs_All"
fi

if [ -z "$format" ]; then
    format=BAM
fi

if [ -z "$genome" ]; then
    echo "genome argument requried. ce11, mm10, hg19, am6"
    exit 1;
    
fi

if [ "$genome" == "mm10" ]; then
    species_macs="mm"
    chromSize="/groups/bell/jiwang/Genomes/Mouse/mm10_UCSC/Sequence/mm10_chrom_sizes.sizes"
fi

if [ "$genome" == "hg19" ]; then
    species_macs="hs";
    
fi

if [ "$genome" == "ce11" ]; then
    species_macs="ce";
    
fi

if [ "$genome" == "am6" ]; then
    species_macs=30000000000;
fi


pval=0.00001; # for macs sharp
fdr=0.05; # for macs broad and sicer

DIR_logs=$PWD/logs
jobName='macs2'
mkdir -p $DIR_logs
mkdir -p $OUT

echo $MACS2

if [ "$MACS2" == "TRUE" ]; then mkdir -p $OUT/macs2; fi;
if [ "$MACS2_broad" == "TRUE" ]; then mkdir -p $OUT/macs2_broad; fi

for sample in ${DIR_Bams}/*.bam; do
    samplename=`basename "$sample"`
    out=${samplename%.bam}
    fname=$out
    echo $sample $out
    
    #echo $outname
        
    # MACS2
    if [ "$MACS2" == "TRUE" ]; then
	cd $OUT/macs2
    fi;
    
    # MACS2 broad_peaks
    if [ "$MACS2_broad" == "TRUE" ]; then
	cd $OUT/macs2_broad
    fi
    
    # creat the script for each sample
    script=${fname}_${jobName}.sh
    cat <<EOF > $script
#!/usr/bin/bash

#SBATCH --cpus-per-task=$nb_cores
#SBATCH --time=240
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o $DIR_logs/${fname}.out
#SBATCH -e $DIR_logs/${fname}.err
#SBATCH --job-name $jobName

module load macs2/2.2.5-foss-2018b-python-3.6.6

EOF
    
    if [ "$MACS2" == "TRUE" ]; then
	#echo "peak calling with macs2"
	if [ -n "$INPUT" ]; then # with input
	    cat <<EOF >> $script
macs2 callpeak -t $sample -c $INPUT -n ${out}_macs2_pval_${pval} -f BAM \
-g $species_macs -p $pval --fix-bimodal -m 5 100 --call-summits

EOF

	# without input
	elif [ "$genome" == "am6" ]; then # axololt
	    cat <<EOF >> $script
	    macs2 callpeak -t $sample -n ${out}_macs2 -f BAMPE -g 20000000000 -p 0.001 	    
EOF

	else 
	    cat <<EOF >> $script
	    macs2 callpeak -t $sample -n ${out}_macs2 -f $format -g $species_macs --nomodel --shift -100 --extsize 200	    
EOF
	fi;
	
    fi;
    
    # MACS2 broad_peaks
    if [ "$MACS2_broad" == "TRUE" ]; then
	if [ -n "$INPUT" ]; then # with input
	    cat <<EOF >> $script
macs2 callpeak -t $sample -c $INPUT -n ${out}_macs2_broad_fdr_${fdr} -f BAM -g $species_macs --broad --broad-cutoff $fdr --fix-bimodal --extsize 200
EOF
       
	else # without input
	    cat <<EOF >> $script
macs2 callpeak -t $sample -n ${out}_macs2_broad_fdr_${fdr} -f BAM -g $species_macs --broad --broad-cutoff $fdr --fix-bimodal --extsize 200
EOF
    
	fi

    fi;
    

    cat $script
    sbatch $script
    
    cd $cwd #back to the main working directory
    
    #break;
    
done 


