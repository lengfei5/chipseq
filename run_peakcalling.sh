#####################
# this script is to call peaks for ChIP-seq (or ATAC-seq) data 
# using MACS2 (sharp or broad peaks), Sicer (broad peaks) 
# it is quite tricky to specify Input files if bam files need different Input
#####################
DIR_Bam="$PWD/alignments/BAMs_All"
OUT=$PWD/Peaks

# input file if exist
#INPUT="/groups/bell/jiwang/Projects/Jorge/INPUT/merged_Input_49475_49911_49908.bam"

nb_cores=6;
chromSize="/groups/bell/jiwang/Genomes/Mouse/mm10_UCSC/Sequence/mm10_chrom_sizes.sizes"
species_macs="hs";

species_sicer="mm10";
pval=0.00001;
fdr=0.01;
window_sizes="200 500 1000";
redundancy_threshold=1;
fragment_size=140;
eff_genome_fraction=0.8;
gap_size=

cwd=`pwd`;

mkdir -p $cwd/logs
mkdir -p $OUT/macs2

#mkdir -p $OUT/macs2_broad
#mkdir -p $OUT/peakranger
#mkdir -p $OUT/sicer

for sample in ${DIR_Bam}/*.bam; do
    # skip the header
    #if [ "$CONDITION" = "CONDITION" ]; then
    #echo "Header ..."
    #elif [ ! -z "$INPUT" ]; then 
    #sample=`ls -l ${DIR_Bam}/*.bam | grep $ID |awk '{print $9}'`
	#sample=${DIR_Bam}/${CONDITION}_${ID}.bam
	#echo $sample $INPUT
	
	#IFS=; read -a inputs <<< "$INPUT"
	#inputs=$(echo $INPUT | tr ";" " ")
	#inputs=(${INPUT//;/})
	#echo $sample $INPUT
	#input=`ls -l ${DIR_Bam}/*.bam | grep $INPUT |awk '{print $9}'`
    #echo $sample $INPUT
    
    samplename=`basename "$sample"`
    out=${samplename%.bam}
    echo $out
    #inputname=`basename "$input"`
    #outname=${CONDITION}_${ID}_${INPUT} 
    #echo $outname
        
    ## MACS2
    cd $OUT/macs2
    #echo "peak calling with macs2"
    if [ -n "$INPUT" ]; then  
    qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N macs2 "module load macs/2.1.0; macs2 callpeak -t $sample -c $INPUT -n ${out}_macs2_pval_${pval} -f BAM -g $species_macs -p $pval --fix-bimodal -m 5 100 -B --call-summits" 
    else
    qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N macs2 "module load macs/2.1.0; macs2 callpeak -t $sample -n ${out}_macs2_pval_${pval} -f BAM -g $species_macs -p $pval --fix-bimodal -m 5 100 -B --call-summits" 
    fi
    cd $cwd
    #break;
	
	## MACS2 broad_peaks
	#cd $OUT/macs2_broad
	#qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N macs2_broad "module load macs/2.1.0; macs2 callpeak -t $sample -c $input -n ${outname}_macs2_broad_fdr_${fdr} -f BAM -g $species_macs -q $fdr --broad --fix-bimodal -m 5 100 -B " 
	#cd $cwd
	
        ## peakranger
	#cd $OUT/peakranger
        #qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N peakranger "module load peakranger; peakranger ccat --format bam $sample $input ccat_$outname "
	#cd $cwd
	
        ## sicer
	#cd $OUT/sicer;
	#bedsample=${CONDITION}_${ID};
	#bedinput="INPUT"_${INPUT};
	#echo $bedsample $bedinput;
	#data=$OUT/sicer;
	
	#module load bedtools; 
	#if [ ! -e ${bedsample}.bed ]; then bamToBed -i $sample > ${bedsample}.bed; fi; 
	#if [ ! -e ${bedinput}.bed ]; then bamToBed -i $input > ${bedinput}.bed; fi;
	
	#for window_size in $window_sizes
	#do
	 #   gap_size=$(expr $window_size \* 3);
	  #  odir=${data}/${bedsample}/${window_size};
	    
	   # echo $gap_size 
	    #echo $odir;
	    #mkdir -p $odir; 
	    #cd $odir;
	    #output of sicer (summary file): chrom, start, end, ChIP-island-read-count, control-island-read-count, p-value, fold-change, fdr-threshod
	    #qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N sicer "module unload python; module load python/2.7.3; module load pythonlib; bash /groups/cochella/jiwang/local/SICER_V1.1/SICER/SICER.sh $OUT/sicer ${bedsample}.bed ${bedinput}.bed "." $species_sicer $redundancy_threshold $window_size $fragment_size $eff_genome_fraction $gap_size $fdr; "
	    
	#done
	#cd $cwd;
	
	#break;
    
    #fi
	
done #< $PARAMETER


