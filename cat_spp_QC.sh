###########
## this script is to generate spp quality control summary
##########
txtfiles="$@";
echo "filename numReads estFragLen corr_estFragLen phantomPeak corr_phantomPeak argmin_corr min_corr NSC(1.05) RSC(0.8) QualityTag(-2:verylow;-1:low;0:medium;1:high;2:veryhigh)" | tr " " "\t" 
#cat QCs/ALL/SPP_QCs/*.txt
#sed -i '1ifilename\tnumReads\testFragLen\tcorr_estFragLen\tphantomPeak\tcorr_phantomPeak\targmin_corr\tmin_corr\tNSC(1.05)\tRSC(0.8)\tQualityTag(-2:verylow;-1:low;0:medium;1:high;2:veryhigh)' summary_spp_QCs.txt
cat $txtfiles;