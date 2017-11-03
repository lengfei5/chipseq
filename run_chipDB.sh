########
### this script is to do differential analysis for ChIP-seq under multiple conditions at different levels: peak overlapping 
### after downsampling; differential analysis (DB) with merged bams within peak regions; DB with replicates withinin regions; 
### window-based method (csaw) to finely detect DB regions.  
########
hostname;
date;
## Inputs
cwd=`pwd`;
DIR_Bam_merged="/groups/bell/jiwang/Projects/Jorge/PART_A/BAMs_merged"
DIR_Bam_rep="/groups/bell/jiwang/Projects/Jorge/PART_A/BAMs_replicates"
DIR_Bed_peak="/groups/bell/jiwang/Projects/Jorge/PART_A/BEDs_Peaks_Summits"

INPUT="/groups/bell/jiwang/Projects/Jorge/INPUT/merged_Input_49475_49911_49908.bam"
q_cutoff=10;
Conditions="AN312_ 515D10_ 515D10H3_"
Samples="Eed Ezh2 Suz12 Ring1B Cbx7 H3K27me3 H2A119Ub"
union="pval_10_5_PEAK_Union.bed"

DB_PEAK="TRUE";
USE_MergedBam="FALSE"
PEAK_COMP="FALSE";
DB_WINDOW="FALSE";

nb_cores=6;
chromSize="/groups/bell/jiwang/Genomes/Mouse/mm10_UCSC/Sequence/mm10_chrom_sizes.sizes"
species_macs="mm"; 
species_sicer="mm10";
pval=0.00001; fdr=0.01;

## Output
OUT=$PWD/chipDB_PartA; 
mkdir -p $OUT
#mkdir -p $OUT/peakcomp/macs2
mkdir -p $OUT/peakcomp;
mkdir -p $OUT/DB_peak;
mkdir -p $OUT/DB_window; 
mkdir -p $cwd/logs

module load samtools/1.3.1;

for ss in $Samples; do
    ## DB with merged bams and peaks
    if [ $DB_PEAK == "TRUE" ]; then
	echo "start DB for merged peaks with merged bam files"
	mkdir -p $OUT/DB_peak/$ss
	cd $OUT/DB_peak/$ss;
	peaks=$(ls $DIR_Bed_peak/*${union} |grep ${ss});
	echo $peaks;
	
	bams_merged=$(ls $DIR_Bam_merged/*.bam |grep ${ss});
	bams_rep=$(ls $DIR_Bam_rep/*.bam |grep ${ss});
	echo ${bams_rep[@]}
			
	bbs=$(echo $bams_rep $INPUT);
	norms=${ss}_cntStat_chipseq_used4normalization.txt;
        #counts=${ss}_readcounts_chipeaks.txt;
	echo "sampleName total mapped unique unique.rmdup " |tr ' ' '\t' > $norms;
        	      
	echo " module load bedtools/2.25.0; 
               module unload samtools
               #module load samtools/1.3.1;
               module load samtools/0.1.18; 
               for b in $bbs; do
                 if [[ ! -e \${b}.bai ]]; then 
                    samtools index \$b; 
                 fi; 
                 file=\$(basename \$b); 
                 ff=\${file%.bam}; 
                 ## count reads within peak regions
                 bedtools multicov -bams \$b -bed $peaks -q $q_cutoff  > \${ff}_readcounts_chipeaks.txt
                 newb=\${ff}_uniq
                 newbb=\${ff}_uniq_rmdup
                 
                 ## save the unique and non-duplicated read counts as normalization factors.
                 #ff=\"$ffs $file\"; 
                 #nb_mapped=$(samtools view -c \$b); mapped=\"$mapped $nb_mapped\";
                 if [ ! -e \${newb}.bam ]; then 
                   samtools view -q $q_cutoff -b \$b | samtools sort - \$newb; 
                 fi; 
                 if [ ! -e \${newb}.bam.bai ]; then 
                 samtools index \$newb.bam; 
                 fi; 
                 if [ ! -e \${newbb}.bam ]; then 
                 samtools rmdup -s \${newb}.bam \${newbb}.bam; 
                 fi; 
                 if [ ! -e \${newbb}.bam.bai ]; then 
                 samtools index \${newbb}.bam; fi; 
                 total=\$(samtools view -c \$b); 
                 mapped=\$(samtools view -c -F 4 \$b); 
                 unique=\$(samtools view -c \$newb.bam); 
                 rmdup=\$(samtools view -c \$newbb.bam);
                  
                 echo \"\$ff \$total \$mapped \$unique \$rmdup\" |tr ' ' '\t' >> $norms;
	      done   
              #mkdir -p Olds; mv *.bam* Olds; " > ${OUT}/DB_peak/${ss}/${ss}_chipeak_readcounts.sh
	
	qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N chipseq_DB "bash $OUT/DB_peak/${ss}/${ss}_chipeak_readcounts.sh" 
	
	#break; 
	
	if [ "$USE_MergedBam" == "TRUE" ]; then
	    for b in ${bams_merged[@]}; do
		if [[ ! -e ${b}.bai ]]; then
		    samtools index $b;
		fi
		#bams_both="$bams_both $b"
		file=`basename $b`;file=${file%.bam}; 
		#echo qsub -q public.q -o $cwd/logs -j yes -pe smp 1 -cwd -b y -shell y -N readcounts "module load bedtools/2.25.0; bedtools multicov -bams $b -bed $common_peaks > ${file}_readcounts_withinpeaks.bed" 
		ffs="$ffs $file";
		nb_mapped=$(samtools view -c $b); echo $nb_mapped;
		mapped="$mapped $nb_mapped";
	    done
	fi
	
	#echo $ffs |tr " " "\t" > ${ss}_Stat_mappedreads.txt
	#echo $mapped |tr " " "\t" >> ${ss}_Stat_mappedreads.txt
        #echo ${bams_both[@]};
	#echo $files
	#echo "bedtools multicov -bams ${bams_both[@]} -bed $common_peaks > ${ss}_read_counts_within_peaks.bed";
	cd $cwd;
    fi;
   
    ## peak comparison part
    if [[ $PEAK_COMP == "TRUE" ]]; then
	mkdir -p $OUT/peakcomp/$ss
	cd $OUT/peakcomp/$ss
        ## filter the unique reads and remove duplicates; the read number for downsamples
	reads4down=$(echo "10^10"|bc);echo $reads4down
	for cond in $Conditions; do
            bam=$(ls $DIR_Bam_merged/*.bam | grep ${ss} | grep ${cond});
	    echo $bam
	    newb=$(basename $bam);
	    newb=${newb%.bam}_uniq_rmdup.bam;
	    echo $newb
       	    if [ ! -e $newb ]; then 
		samtools rmdup -s $bam $newb;
            fi
	    nb_reads=$(samtools view -c $newb);
	    echo $nb_reads;
	    if [[ $reads4down -gt $nb_reads ]]; then
		reads4down=$nb_reads
	    fi
      	    #echo $reads4down;	
	done
        ## downsample and peak calling
	echo $reads4down;
	mkdir -p macs2
	for bb in *_uniq_rmdup.bam; do
	    #samplename=`basename "$bb"`
	    pp=${bb%.bam}
	    newbed=${pp}_downsampled.bed
	    echo $pp
	    #qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N macs2 "module load macs/2.1.0; macs2 randsample -t $bb -n $reads4down -o $newbed; macs2 callpeak -t $newbed -c $INPUT -n ${pp}_macs2_pval_${pval} --outdir macs2/ -g $species_macs -p $pval --fix-bimodal -m 5 100 -B"
    	done
	cd $cwd;
    fi
    
     
    #break;
    done

echo "Code ends"