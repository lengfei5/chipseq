##################################################
## Project: general purpused quality control specific for ChIP-seq data
## Script purpose: 
## Usage example: Rscript ~/scripts/Chipseq/chipseqQC.R XXX(current path)
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Wed Jan 10 11:17:03 2018
## tested in cluster
##################################################
source('/home/imp/jingkui.wang/scripts/Chipseq/functions_chipSeq.R')

simplify.filename = function(x){
  test = unlist(strsplit(as.character(x), '[.]'))
  return(paste0(test[1:(length(test))-1], sep='', collapse = '.'))
}

#DIR.cwd = "../R5618_atac/"
DIR.cwd = getwd();
DIR.spp = paste0(DIR.cwd, "/QCs/SPP_QCs")
DIR.stat = paste0(DIR.cwd, "/QCs/ALL/BamStat")
DIR.out = paste0(DIR.cwd, "/QCs")

cat("current directory --", DIR.cwd, "\n");
cat("spp directory --", DIR.spp, "\n");
cat("stat directory --", DIR.spp, "\n")
cat("QC summary directory --", DIR.out, "\n");

if(dir.exists(DIR.spp))
{
  spp.files = list.files(path = DIR.spp, pattern = "*.txt", full.names = TRUE)
  if(length(spp.files)>0)
  {
    spp = c()
    for(n in 1:length(spp.files))
    {
      spp = rbind(spp, read.table(spp.files[n], sep="\t", header = FALSE))
    }
    
  }else{
    cat("no SPP txt files FOUND !!! \n")
  }
  
  colnames(spp) = c("filename",  "numReads",  "estFragLen", 
                    "corr_estFragLen",  "phantomPeak",  "corr_phantomPeak", 
                    "argmin_corr",  "min_corr", "NSC",  "RSC",  "Quality")
  
}else{
  cat("NO Folder found for SPP !!! \n")
}

if(dir.exists(DIR.stat))
{
  stat.files = list.files(path = DIR.stat, pattern = "*.txt", full.names = TRUE)
  jj = grep("picardDup", stat.files)
  if(length(jj)>0) stat.files = stat.files[-jj]
  if(length(stat.files)>0)
  {
    stat = c()
    for(n in 1:length(stat.files)) stat = rbind(stat, read.table(stat.files[n], sep = "\t", header = FALSE))
    
  }else{
    cat("No mapping stat files FOUND !!! \n")
  }
}else{
  cat("NO Folder found for mapping stat !!! \n")
}

#stat = read.table('BamStat_summary.txt', sep='\t', header = FALSE)
colnames(stat) = c('filename','nb.RawRead', 'nb.MappedRead', 'nb.UniqueRead', 'nb.UniqueRmdupRead')
#spp = read.table('../R4646/summary_spp_QCs.txt', sep = '\t', header = TRUE)

xx = sapply(spp[,1], simplify.filename)
spp[,1] = xx

mm = match(spp[,1], stat[,1])
stat = data.frame(stat, spp[, -c(1,2)], stringsAsFactors = FALSE)

max.nb.sample = 12;

#source("functions_chipSeq.R")
pdf(paste0(DIR.out, "/QCs_Summary_chipseq.pdf"), width = 20, height = 12)
par(cex = 2.0, las = 1, mgp = c(1.6,0.5,0), mar = c(3,22,2,0.8)+0.1, tcl = -0.3)

nn = nrow(stat)%/%max.nb.sample;

# stat = matrix(NA, ncol=12, nrow=49)
for(n in 0:nn)
{
  index = intersect(c(1:nrow(stat)), seq((max.nb.sample*n+1), (max.nb.sample*(n+1)), by = 1))
  if(length(index)>0){
    #print(index)
    PLOT.Quality.Controls.Summary(stat[index, ])
  }
}


dev.off()


