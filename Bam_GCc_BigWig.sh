##########
### this script is to correct GG bias sample by sample for ChIP-seq Bams
##########
cwd=`pwd`;
nb_cores=12;

OUT="$PWD/heatmaps_PartB/GC_correction"
DIR_Bam_merged="/groups/bell/jiwang/Projects/Jorge/PART_B/BAMs_merged"
DIR_peak="/groups/bell/jiwang/Projects/Jorge/PART_B/Peaks_Union_BEDs_tmp"

INPUT="/groups/bell/jiwang/Projects/Jorge/INPUT/merged_Input_49475_49911_49908.bam"
#DIR_Bam_rep="/groups/bell/jiwang/Projects/Jorge/PART_A/BAMs_replicates"
genome2b="/groups/bell/jiwang/Genomes/Mouse/mm10.2bit"
genome_effect=2150570000

mkdir -p ${PWD}/logs
mkdir -p $OUT
mkdir -p ${OUT}/logs

bams=( $(echo ${DIR_Bam_merged}/*.bam) );
#echo ${bams[@]};
bams+=($INPUT);
#echo ${bams[0]};
#echo ${bams[1]}

### convert bam file to bigWig (optioal: corrected GC before)
for b in "${bams[@]}" 
do
    #echo $b;
    ff=`basename $b`;
    ff=${ff%.bam};
    cc=${ff}
    cc2=${ff}_GGc_peakExcluded
    prot=`echo $ff | cut -d"_" -f2`;
    cond=`echo $ff | cut -d"_" -f1`;
    cond_prot=${cond}_${prot}
    echo $ff
    #peak=`ls ${DIR_peak}/*_macs2_pval_10_5_CENTER.bed | grep AN312_ | grep $prot `;
    peak=`ls ${DIR_peak}/*_pval_10_5_PEAK.bed | grep "$cond_prot" `;
    echo $prot
    echo $peak
    echo $b
    #if [ "$prot" == "Ring1B" ]; then
    #echo $prot 
    #echo $peak
    #computeGCBias -b $b --effectiveGenomeSize $genome_effect -g $genome2b -l 200 --GCbiasFrequenciesFile ${OUT}/${ff}_freq_before_GCcorrection.txt --biasPlot ${OUT}/${ff}_freq_before_GCcorrection.png;
    #computeGCBias -b ${OUT}/${cc}.bam --effectiveGenomeSize $genome_effect -g $genome2b -l 200 --GCbiasFrequenciesFile ${OUT}/${ff}_freq_after_GCcorrection.txt --biasPlot ${OUT}/${ff}_freq_after_GCcorrection.png; 
    #if [ ! -e ${OUT}/${cc}.bam ]; then
    #computeGCBias -b $b --effectiveGenomeSize $genome_effect -g $genome2b -l 200 --GCbiasFrequenciesFile ${OUT}/${ff}_freq.txt --biasPlot ${OUT}/${ff}_freq_before_GCc.png;
    #correctGCBias -b $b --effectiveGenomeSize $genome_effect -g $genome2b -freq ${OUT}/${ff}_freq.txt -o ${OUT}/${cc}.bam;
    #computeGCBias -b ${OUT}/${cc}.bam --effectiveGenomeSize $genome_effect -g $genome2b -l 200 --GCbiasFrequenciesFile ${OUT}/${cc}_freq.txt --biasPlot ${OUT}/${cc}_freq.png;
    #fi;

   echo -e "module unload deeptools; module load libcurl;  module load deeptools/2.2.3-python2.7.3;
if [ ! -e ${OUT}/${cc2}.bam -a -n "$peak" ]; then
computeGCBias -b $b --effectiveGenomeSize $genome_effect -g $genome2b -l 200 -bl $peak --GCbiasFrequenciesFile ${OUT}/${ff}_freq_peakExcluded.txt --biasPlot ${OUT}/${ff}_freq_before_GCc.png;                                                         
correctGCBias -b $b --effectiveGenomeSize $genome_effect -g $genome2b -freq ${OUT}/${ff}_freq_peakExcluded.txt -o ${OUT}/${cc2}.bam; 
computeGCBias -b ${OUT}/${cc2}.bam --effectiveGenomeSize $genome_effect -g $genome2b -l 200 --GCbiasFrequenciesFile ${OUT}/${cc2}_freq.txt --biasPlot ${OUT}/${cc2}_freq.png;                                                                            
fi;

module unload deeptools; module load python/2.7.3; module load pysam/0.8.3; module load deeptools/2.2.3-python2.7.3;
if [ ! -e ${OUT}/${cc}.bw ]; then bamCoverage -b $b -o $OUT/${cc}.bw --outFileFormat=bigwig --normalizeUsingRPKM; fi;
if [ ! -e ${OUT}/${cc2}.bw -a -n "$peak" ]; then bamCoverage -b ${OUT}/${cc2}.bam -o $OUT/${cc2}.bw --outFileFormat=bigwig --normalizeUsingRPKM; fi; " > ${OUT}/logs/make_matrix4heatmap_${ff}.sh

    qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N GCcorr "bash ${OUT}/logs/make_matrix4heatmap_${ff}.sh" 
    
    #break;
    #fi;
done





 



