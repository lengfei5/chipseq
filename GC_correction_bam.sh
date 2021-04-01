###############################################
# This script is to correct GG bias sample by sample for ChIP-seq aligned bams 
# currently deeptools is used to do so 
# Input is bam file and associated peaks, because GC correction is done by excluding those peak regions
# otherwise it would be easily over-corrected
# the initial script was modified for cbe
#
###############################################
cwd=`pwd`;
nb_cores=2;

DIR_Bam="$PWD/alignments/BAMs_All" # bam folder 
DIR_peak="$PWD/Peaks/macs2_broad" # peak folder
#INPUT="/groups/bell/jiwang/Projects/Jorge/INPUT/merged_Input_49475_49911_49908.bam"

OUT="$PWD/alignments/GC_correction_4Bams_peaksExcluded"
OUT_Bam="$PWD/alignments/BAMs_All_GCc"
OUT_BW="$PWD/alignments/BigWigs_All"

# genome parameters
genome2b="/groups/bell/jiwang/Genomes/Mouse/mm10.2bit"
genome_effect=2150570000

mkdir -p ${PWD}/logs
mkdir -p $OUT;
mkdir -p $OUT_Bam;
mkdir -p $OUT_BW;
mkdir -p ${OUT}/logs

#bams=( $(echo ${DIR_Bam}/*.bam) );
#echo ${bams[@]};
#if [ -n "$INPUT" ]; then 
#    bams+=($INPUT);
#fi
#echo ${bams[0]};
#echo ${#bams[@]};

## convert bam file to bigWig (optioal: corrected GC before)
for b in ${DIR_Bam}/*.bam 
do
    #echo $b;
    ff=`basename $b`;
    ff=${ff%.bam};
    cc=${ff}
    cc2=${ff}_GGc
    #prot=`echo $ff | cut -d"_" -f2`;
    echo $ff

    # select the peak associated with bam file (be careful the name !!!) 
    peak=`ls ${DIR_peak}/*macs2_broad_fdr_0.1_peaks.broadPeak | grep $ff `;

    echo $peak
    echo $b

    # make job submission script
    script="${OUT}/logs/${ff}_Bam_GCc.sh"
cat <<EOF > $script
#!/usr/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=180
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o $OUT/logs/$fname.out
#SBATCH -e $OUT/logs/$fname.err
#SBATCH --job-name gc_cor

ml load deeptools/3.3.1-foss-2018b-python-3.6.6

# compute GC bias
computeGCBias -b $b --effectiveGenomeSize $genome_effect -g $genome2b \
-l 200 -bl $peak --GCbiasFrequenciesFile ${OUT}/${ff}_freq_peakExcluded.txt \
--biasPlot ${OUT}/${ff}_freq_before_GCc.png; \

# correct GC bias
correctGCBias -b $b --effectiveGenomeSize $genome_effect -g $genome2b \
-freq ${OUT}/${ff}_freq_peakExcluded.txt -o ${OUT_Bam}/${cc2}.bam; \
 
# compute again GC bias after correction
computeGCBias -b ${OUT_Bam}/${cc2}.bam --effectiveGenomeSize $genome_effect \
-g $genome2b -l 200 --GCbiasFrequenciesFile ${OUT}/${cc2}_freq_after_GCc.txt \
--biasPlot ${OUT}/${cc2}_freq_after_GCc.png;

# convert bam to bigWig
bamCoverage -b $b -o ${OUT_BW}/${cc}.bw --outFileFormat=bigwig \
--normalizeUsingRPKM --ignoreDuplicates --minMappingQuality 10 --extendReads 200
	
bamCoverage -b ${OUT_Bam}/${cc2}.bam -o ${OUT_BW}/${cc2}.bw --outFileFormat=bigwig \
--normalizeUsingRPKM --ignoreDuplicates --minMappingQuality 10 --extendReads 200;

EOF

    cat $script;
    sbatch $script
    #break;

done





