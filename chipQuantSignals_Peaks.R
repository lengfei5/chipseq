#############################
##################################################
## Project: ZFP445 binding 
## Script purpose: to correlat ZFP445 binding peaks with H3K27me3 and H3K9me3
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Tue Jan 23 11:56:33 2018
##################################################
library(Rsubread)

##################################################
## Section: specify input, output and parameters
##################################################
DIR.bed = getwd();
DIR.bam = "alignments/BAMs_All";
DIR.res = "compare_H3K27me3_H3K9me3"
if(!dir.exists(DIR.res)) system(paste0('mkdir -p ', DIR.res))

Merge.bed.file = TRUE

# import chipseq peqks and bams 
bam.list = list.files(path=DIR.bam, pattern = "*.bam$", full.names = TRUE, ignore.case = TRUE)
bed.list = list.files(path=DIR.bed, pattern = "*.bed$", full.names = TRUE, ignore.case = TRUE)

if(Merge.bed.file){
  peaks = c();
  for(n in 1:length(bed.list))
  {
    bed = read.table(bed.list[n], sep = "\t", header = FALSE)
    
    if(length(unique(bed[,4])) != nrow(bed)) stop("peak names are Redundant !");
    if(length(grep('random', basename(bed.list[n])))>0) bed[, 4] = paste0("random_", seq(1:nrow(bed)));
    peaks = rbind(peaks, bed); 
  }
  peaks = data.frame(peaks, stringsAsFactors = FALSE)
  colnames(peaks) = c("chr.peak", "start.peak", "end.peak", "peak.name", "score.peak", "strand.peak")[c(1:ncol(peaks))]
}

######################################################
## Quantify chipseq signals for all peaks regions
######################################################
Quantify.signals.within.peaks = TRUE
Count.Reads.for.Peak.Regions = TRUE
Change.Sample.Names = FALSE
if(Quantify.signals.within.peaks)
{
  if(Count.Reads.for.Peak.Regions)
  {
    ## peak regions (configure your peak regions into a data.frame)
    jj = match(unique(peaks$peak.name), peaks$peak.name)
    df = peaks[jj, ];
    SAF = data.frame(GeneID=df$peak.name, Chr=df$chr.peak, Start=df$start.peak, End=df$end.peak, Strand="+", stringsAsFactors = FALSE)
    
    ## count reads for those peak regions using 
    fc <- featureCounts(files=bam.list, annot.ext = SAF, countMultiMappingReads = FALSE, minMQS = 10, 
                        ignoreDup = TRUE, strandSpecific = 0, juncCounts = FALSE, nthreads = 20)
    stat = fc$stat;
    counts = fc$counts;
    counts.annot = fc$annotation
    
    rpkm = matrix(NA, ncol = ncol(counts), nrow = nrow(counts))
    colnames(rpkm) = colnames(counts)
    row.names(rpkm) = rownames(counts)
    kk = which(stat$Status=="Assigned" | stat$Status== "Unassigned_NoFeatures")
    
    for(n in 1:ncol(counts))
    {
      #n = 1 
      jj = which(colnames(stat) == colnames(counts)[n])
      rpkm[, n] = (counts[, n]+0.25)/counts.annot$Length/sum(stat[kk, jj])*10^9
    }
    
    #rpkm = log2(rpkm)
    res = data.frame(SAF, Length=counts.annot$Length[match(SAF$GeneID, counts.annot$GeneID)], 
                     rpkm[match(SAF$GeneID, rownames(rpkm)), ], stringsAsFactors = FALSE)
    
    ## change colnames
    if(Change.Sample.Names)
    {
      names = colnames(res)[-c(1:7)]
      for(n in 1:length(names))
      {
        #n = 1
        test = unlist(strsplit(as.character(names[n]), "[.]"))
        test = test[3]
        if(length(grep("Input", test))>0){
          ttest = unlist(strsplit(as.character(test), "_"))
          names[n] = paste0(ttest[c(1,3, 2)], collapse = "_")
        }
        if(length(grep("_H3K", test))>0){
          ttest = unlist(strsplit(as.character(test), "_"))
          names[n] = paste0(ttest[c(1, 4, 2)], collapse = "_")
        }
        if(length(grep("tbx3", test))>0) names[n] = test
      }
      colnames(res)[-c(1:7)] = names
      
      kk = grep("6291|6344|Input", colnames(res)) ## all inputs for histone modifications are the same one.
      if(length(kk)>0) res = res[, -kk]
    }
    
    ## save the peak information and quantified different rpkm signals within peaks
    write.table(res, file=paste0(DIR.res, "/rpkm_in_chipeaks.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    #save(peaks, res, file='Rdata/Peaks_merged_macs2_p_5_filtered_N2_gene_assignment_TSS_WBcel235_analysis_hisModif_tbx_signals_v1.Rdata')
  }
}


######################################################
### test the correlations
######################################################
Search4correlations = FALSE
if(Search4correlations)
{
  dd = as.matrix(res[, -c(1:6)])
  rownames(dd) = res$GeneID
  colnames(dd) = c("H3K27me3", "H3K9me3_rep1", "H3K9me3_rep2", "Input", "ZFP445")
  dd = log2(dd)
  groups = rep(NA, nrow(dd))
  groups[grep('MACS', rownames(dd))] = 1
  groups[grep('MACS', rownames(dd), invert = TRUE)] = 2
  
  pdfname = paste0(DIR.res, "/comparison_plots.pdf")
  pdf(pdfname, width = 12, height = 8)
  ## boxplot for peaks and randome locations
  par(cex = 1.8, las = 1, mgp = c(1.6,0.5,0), mar = c(12,3,2,0.8)+0.1, tcl = -0.3)
  par(mfrow=c(1,1))
  plot(c(0, 12), c(-8, 6), type = 'n', axes = FALSE, ylab = "log2 (rpkm)", xlab=NA)
  for(n in 1:ncol(dd))
  {
    add = TRUE
    if(n>1) add = TRUE
    boxplot(dd[,n] ~ groups, col= c("blue", "darkgray"), at=c((2*n-1), 2*n), las=2,
            names=paste0(colnames(dd)[n], c("_peak", "_random")), add =add)
  }
  
  cex = 0.25
  ## scatterplot for all 
  par(pty="s")
  pairs(dd, cex=cex)
  ## scatterplot for peaks
  kk = grep("MACS", rownames(dd))
  pairs(dd[kk,], cex=cex)
  ## scatterplot for randome sequence
  kk = grep("MACS", rownames(dd), invert = TRUE)
  pairs(dd[kk,], cex=cex)
  
  dev.off()
  
}








