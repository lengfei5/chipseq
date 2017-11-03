#####
## this script is to demultiplex the bam files from vbcf
####
nb_cores=8
DIR_cwd=`pwd`
DIR_download=$PWD
DIR_FC=$PWD/ngs_raw/RAWs
DIR_OUT=$PWD/ngs_raw/BAMs

mkdir -p $DIR_OUT;
mkdir -p $DIR_cwd/logs;
#cd ngs_raw
#mkdir -p RAWS

## download the raw data 
#module load R/3.3.0
#FC="CAEBDANXX"
#url="http://ngs.vbcf.ac.at/data/47414_ATCACG_CA88VANXX_3_20161216B_20161216.bam"
#RAW="../CAFVMANXX_8_20170120B_20170123.bam"
#wget -c --no-check-certificate --auth-no-challenge --user 'Jingkui.Wang' --password '8A1P2EWa1K' http://gecko:9100/data/CAKNHANXX_2_20170221B_20170223.bam
#wget -c --no-check-certificate --auth-no-challenge --user 'Jingkui.Wang' --password '8A1P2EWa1K' $url
#bash /clustertmp/jiwang/Julien_R4224/getCSFNGS.sh -p $NGS_USER -f $FC
#cd ..

#cd $DIR_FC;
## demultiplex raw data
#for RAW in ${DIR_FC}/*.bam; do
#    echo $RAW
#    qsub -q public.q -o $DIR_cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N demultiplex  "/groups/vbcf-ngs/bin/funcGen/jnomicss.sh illumina2BamSplit --inputFile $RAW"
#done
cd $DIR_cwd