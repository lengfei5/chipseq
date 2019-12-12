#####################
# this script is to call peaks for ChIP-seq (or ATAC-seq) data 
# using MACS2 (sharp or broad peaks), Sicer (broad peaks) 
# it is quite tricky to specify Input files if bam files need different Input
#####################
while getopts ":hD:I:mbSg:" opts; do
    case "$opts" in
        "h")
            echo "script to call chipseq peak using macs2 (sharp or broad optioins and  sicer for broad peaks"
	    echo "options: "
	    echo "-D  (specfiy the ABSOLUTE directory of bam files or by default alignments/BAMs_All)"
	    echo "-I (input, e.g. /groups/bell/jiwang/Projects/Jorge/Analysis_ChIP_seq/INPUT/merged_Input_49475_49911_49908.bam for mouse mm10)" 
	    echo "no input file needed by defaut"
            echo "-g (genome, i.e. mm10, ce11, hg19)"
	    echo "-m (macs2 sharp peaks)"
	    echo "-b (macs2 broad peaks)"
	    echo "-S (sicer)"
	    echo "Usage: "
            echo "$0 -m -g mm10 (if bam files in alignments/BAMs_All for chipseq)"
            echo "$0 -D /clustertmp/Jorge_R6118/alignments/BAMs_All_GCc_merged -m -g mm10 (bam files in XXX directory) "
            echo "$0 -b -g mm10 (call mm10 peaks with macs2 broad option) "
	    echo "$0 -S -g mm10 (call mm10 peaks with SICER) "
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
	"b")
	    MACS2_broad="TRUE"
	    ;;
	"S")
	    SICER="TRUE"
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

OUT=$PWD/Peaks
nb_cores=4;
cwd=`pwd`;

if [ -z "$DIR_Bams" ]; then
    DIR_Bams="$PWD/alignments/BAMs_All"
fi

if [ "$genome" == "mm10" ]; then
    species_macs="mm";
    species_sicer="mm10";
    chromSize="/groups/bell/jiwang/Genomes/Mouse/mm10_UCSC/Sequence/mm10_chrom_sizes.sizes"
fi

if [ "$genome" == "hg19" ]; then
    species_macs="hs";
    species_sicer="hg19";
fi

if [ "$genome" == "ce11" ]; then
    species_macs="ce";
    species_sicer="ce11";
fi

pval=0.00001; # for macs sharp
fdr=0.05; # for macs broad and sicer

# paras for sicer
window_sizes="200 500 1000";
redundancy_threshold=1;
fragment_size=140;
eff_genome_fraction=0.8;
gap_size=

DIR_logs=$PWD/logs
jobName='peakcalling'
mkdir -p $DIR_logs

if [ "$MACS2" == "TRUE" ]; then mkdir -p $OUT/macs2; fi;
if [ "$MACS2_broad" == "TRUE" ]; then mkdir -p $OUT/macs2_broad; fi
if [ "$SICER" == "TRUE" ]; then mkdir -p $OUT/sicer; fi

for sample in ${DIR_Bams}/*.bam; do
    samplename=`basename "$sample"`
    out=${samplename%.bam}
    fname=$out
    echo $sample $out
    #inputname=`basename "$input"`
    #outname=${CONDITION}_${ID}_${INPUT} 
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
#SBATCH --time=120
#SBATCH --mem=8000
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o $DIR_logs/${fname}.out
#SBATCH -e $DIR_logs/${fname}.err
#SBATCH --job-name $jobName

module load macs2/2.1.2.1-foss-2018b-python-2.7.15;

if [ "$MACS2" == "TRUE" ]; then
   #echo "peak calling with macs2"
   if [ -n "$INPUT" ]; then # with input
      macs2 callpeak -t $sample -c $INPUT -n ${out}_macs2_pval_${pval} -f BAM \
-g $species_macs -p $pval --fix-bimodal -m 5 100 --call-summits
   else # without input
      macs2 callpeak -t $sample -n ${out}_macs2_pval_${pval} -f BAM \
-g $species_macs -p $pval --fix-bimodal -m 5 100 --nomodel --extsize 200 --call-summits
      
   fi
	
fi;

# MACS2 broad_peaks
if [ "$MACS2_broad" == "TRUE" ]; then
    if [ -n "$INPUT" ]; then # with input
       macs2 callpeak -t $sample -c $INPUT -n ${out}_macs2_broad_fdr_${fdr} -f BAM -g $species_macs --broad --broad-cutoff $fdr --fix-bimodal --extsize 200
    else # without input
    macs2 callpeak -t $sample -n ${out}_macs2_broad_fdr_${fdr} -f BAM -g $species_macs --broad --broad-cutoff $fdr --fix-bimodal --extsize 200
    fi

fi;

EOF
    cat $script
    #sbatch $script
    cd $cwd #back to the main working directory
    break;
    
    # SICER (to correct and adapt for cbe) 
    if [ "$SICER" == "TRUE" ]; then
	if [ -z "$INPUT" ]; then
	    echo '--INPUT is required--'
	    exit;
	fi
	cd $OUT/sicer;
	bedsample=${out}
	bedinput=`basename $INPUT`
	inputExt=${INPUT##*.}
	bedinput=${bedinput%.*}
	echo $bedsample 
	echo $bedinput;
	echo $inputExt
	#echo "here !!"
	
	data=$OUT/sicer;
	
	#module load bedtools; 
	if [ ! -e ${bedsample}.bed ]; then bamToBed -i $sample > ${bedsample}.bed; fi; 
	if [ ! -e ${bedinput}.bed ]; then 
	    if [ "$inputExt" == "bed" ]; then
		cp $INPUT $data 
	    else
		bamToBed -i $input > ${bedinput}.bed; 
	    fi
	fi;
	#exit;

	for window_size in $window_sizes
	do
	    gap_size=$(expr $window_size \* 3);
	    odir=${data}/${bedsample}/${window_size};
	    
	    echo $gap_size 
	    echo $odir;
	    mkdir -p $odir; 
	    cd $odir;
	    #output of sicer (summary file): chrom, start, end, ChIP-island-read-count, control-island-read-count, p-value, fold-change, fdr-threshod
	    qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N sicer "module unload python; module load python/2.7.3; module load pythonlib; bash /groups/cochella/jiwang/local/SICER_V1.1/SICER/SICER.sh $OUT/sicer ${bedsample}.bed ${bedinput}.bed "." $species_sicer $redundancy_threshold $window_size $fragment_size $eff_genome_fraction $gap_size $fdr; "
	    
	done
	
	cd $cwd;
	
    fi
    
    # break;

done 


