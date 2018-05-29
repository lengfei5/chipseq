###############################################
# This script is to correct GG bias sample by sample for ChIP-seq Bams 
# currently deeptools is used to do so 
# Input is bam file and associated peaks, because GC correction is done by excluding those peak regions
# otherwise it would be easily over-corrected
# 
###############################################
cwd=`pwd`;
nb_cores=2;

#INPUT="/groups/bell/jiwang/Projects/Jorge/INPUT/merged_Input_49475_49911_49908.bam"

OUT="$PWD/alignments/GC_correction_4Bams_peaksExcluded"
OUT_Bam="$PWD/alignments/BAMs_All_GCc"
OUT_BW="$PWD/alignments/BigWigs_All"

DIR_Bam="$PWD/alignments/BAMs_All"
DIR_peak="$PWD/Peaks/macs2_broad"
genome2b="/groups/bell/jiwang/Genomes/Mouse/mm10.2bit"
genome_effect=2150570000

mkdir -p ${PWD}/logs
mkdir -p $OUT;
mkdir -p $OUT_Bam;
mkdir -p $OUT_BW;
mkdir -p ${OUT}/logs

bams=( $(echo ${DIR_Bam}/*.bam) );
#echo ${bams[@]};
if [ -n "$INPUT" ]; then 
    bams+=($INPUT);
fi

#echo ${bams[0]};
#echo ${bams[1]}
echo ${#bams[@]};
#echo ${bams[@]}

### convert bam file to bigWig (optioal: corrected GC before)
for b in "${bams[@]}" 
do
    #echo $b;
    ff=`basename $b`;
    ff=${ff%.bam};
    cc=${ff}
    cc2=${ff}_GGc
    #prot=`echo $ff | cut -d"_" -f2`;
    echo $ff
    #peak=`ls ${DIR_peak}/*_macs2_pval_10_5_CENTER.bed | grep AN312_ | grep $prot `;
    peak=`ls ${DIR_peak}/*macs2_broad_fdr_0.1_peaks.broadPeak | grep $ff `;
    #echo $prot
    echo $peak
    echo $b
        
    sscript="${OUT}/logs/${ff}_Bam_GCc.sh";
    if [ ! -e ${OUT_Bam}/${cc2}.bam ]; then
	## GC correction with bam output
	if [ -n $peak ]; then
	    echo -e "module unload deeptools; module load libcurl;  module load deeptools/2.2.3-python2.7.3;
            
            computeGCBias -b $b --effectiveGenomeSize $genome_effect -g $genome2b -l 200 -bl $peak --GCbiasFrequenciesFile ${OUT}/${ff}_freq_peakExcluded.txt --biasPlot ${OUT}/${ff}_freq_before_GCc.png;                                                              correctGCBias -b $b --effectiveGenomeSize $genome_effect -g $genome2b -freq ${OUT}/${ff}_freq_peakExcluded.txt -o ${OUT_Bam}/${cc2}.bam; 
            computeGCBias -b ${OUT_Bam}/${cc2}.bam --effectiveGenomeSize $genome_effect -g $genome2b -l 200 --GCbiasFrequenciesFile ${OUT}/${cc2}_freq_after_GCc.txt --biasPlot ${OUT}/${cc2}_freq_after_GCc.png;" > $sscript;
	else
	    echo -e "module unload deeptools; module load libcurl;  module load deeptools/2.2.3-python2.7.3;
            
            computeGCBias -b $b --effectiveGenomeSize $genome_effect -g $genome2b -l 200 --GCbiasFrequenciesFile ${OUT}/${ff}_freq_peakExcluded.txt --biasPlot ${OUT}/${ff}_freq_before_GCc.png;                                                                        correctGCBias -b $b --effectiveGenomeSize $genome_effect -g $genome2b -freq ${OUT}/${ff}_freq_peakExcluded.txt -o ${OUT_Bam}/${cc2}.bam;                                                                                                                    computeGCBias -b ${OUT_Bam}/${cc2}.bam --effectiveGenomeSize $genome_effect -g $genome2b -l 200 --GCbiasFrequenciesFile ${OUT}/${cc2}_freq_after_GCc.txt --biasPlot ${OUT}/${cc2}_freq_after_GCc.png;" > $sscript        
	fi
	
        ## convert bam to bigWig
	if [ ! -e ${OUT_BW}/${cc}.bw ]; then
	    echo -e "module unload deeptools; module load python/2.7.3; module load pysam/0.8.3; module load deeptools/2.2.3-python2.7.3;
             bamCoverage -b $b -o ${OUT_BW}/${cc}.bw --outFileFormat=bigwig --normalizeUsingRPKM --ignoreDuplicates --minMappingQuality 10 --extendReads 200;" >> $sscript;
	fi;
	
	if [ ! -e ${OUT_BW}/${cc2}.bw ]; then 
	    echo -e "bamCoverage -b ${OUT_Bam}/${cc2}.bam -o ${OUT_BW}/${cc2}.bw --outFileFormat=bigwig --normalizeUsingRPKM --ignoreDuplicates --minMappingQuality 10 --extendReads 200; " >> $sscript
	fi;
	
	qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N BamGCc "bash $sscript" 
	
    fi
    
    #break;

done





