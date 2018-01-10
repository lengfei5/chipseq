#setwd("~/clustertmp/Jorge_Arturo.Zepeda_Martinez/")
library(ChIPseeker)
library(rtracklayer)
#library(UpSetR)
library("ChIPpeakAnno")
library("Vennerable") #install.packages("Vennerable", repos="http://R-Forge.R-project.org", type="source")
library("ggplot2")
library("GenomicFeatures")

venn_cnt2venn <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts") - 1
  SetNames=colnames(venn_cnt)[1:n]
  Weight=venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  
  Venn(SetNames=SetNames, Weight=Weight)
}

venn_cnt2barplot <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts")

  cnt <- venn_cnt[,(n+1):length(colnames(venn_cnt))]
  rownames(cnt) <- cbind(apply(venn_cnt[,1:(n-1)], 1, paste, collapse="|"))
  cnt <- cnt[rownames(cnt) != "0|0|0", ]
  cnt <- melt(cnt)
  colnames(cnt) <- c("intersections", "condition", "peaks")
  legend.title <- paste(colnames(venn_cnt[,1:(n-1)]), collapse="|")

  p <- ggplot(cnt, aes(x=condition, y=peaks, fill=intersections)) +  geom_bar(stat="identity", colour="grey50", alpha = 0.9) + scale_fill_brewer(palette="Set3") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank()) + labs(fill=legend.title)

  plot(p)  
}


filterPeaks <- function(files, pval = 10, qval=0, remove.blacklist=FALSE, barplot = FALSE, fileType="macs2", input1B=TRUE, addChr=FALSE)
{
  # files = macs2.files;remove.blacklist=FALSE; barplot = FALSE; fileType="macs2"; input1B=TRUE; addChr=FALSE; pval = 100; qval=0;
  bar <- NULL
  peaks <- NULL
  for (f in files)
  {
    p <- readPeakFile(f, as = "GRanges")
    if (addChr)
    {
      seqlevels(p) <- paste0("chr", seqlevels(p))
    }
    if (input1B)
    {
      start(p) <- start(p) - 1
    }
    all <- length(p)
    if (fileType == "macs14")
    {
      p <- p[mcols(p)[,"X.10.log10.pvalue."] > pval]
    } else if (fileType == "macs2") {
      if (qval == 0)
      {
        p <- p[mcols(p)[,"X.log10.pvalue."] > pval]
      } else {
        p <- p[mcols(p)[,"X.log10.qvalue."] > qval]
      }
    }
    p100 <- length(p)
    if(remove.blacklist)
    {
      p <- p[!overlapsAny(p, blacklist)]
      p.bl <- p100 - length(p)
    }
    all <- all - p100
    p100 <- length(p)
    
    if(remove.blacklist){
      bar <- cbind(bar, c(all, p100, p.bl))
    }else{
      bar <- cbind(bar, c(all, p100))
    }
    
    peaks <- c(peaks, p)
  }
  
  names(peaks) <- names(files)
  
  if (barplot)
  {
    if (qval !=0)
    {
      rownames(bar) <- c("all", paste0("q", pval), "blacklisted")
    } else {
      rownames(bar) <- c("all", paste0("p", pval), "blacklisted")
    }
    colnames(bar) <- names(peaks)
    p <- ggplot(melt(bar), aes(x=X2, y=value, fill=X1)) + geom_bar(stat="identity", colour="black") + scale_fill_brewer() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.title=element_blank())
    plot(p)
  }
  
  return(peaks)  
  
}


plotV <- function(peaks, main = "")
{
    if (length(peaks) > 5)
    {
        peaks <- peaks[1:5] #up to five!!!
    }
    ol.peaks <- makeVennDiagram(peaks, NameOfPeaks=names(peaks), maxgap=0, minoverlap =1, main=main, connectedPeaks="keepAll")
 
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    #plot(v, doWeights=FALSE)
    #plot(v, doEuler = TRUE)
    venn_cnt2barplot(ol.peaks$vennCounts)
    return(ol.peaks)
}


doAllPeak <- function(set)
{
  if (set == "macs14")
  {
    allPeaks <- filterPeaks(macs14.files, barplot=TRUE)
    allPeaks.50 <- filterPeaks(macs14.files, pval=50, blacklist=blacklist)
  } else if (set == "macs2")
  {
    allPeaks <- filterPeaks(macs2.files, pval = 10, blacklist=blacklist, barplot=TRUE, fileType="macs2")
    allPeaks.50 <- filterPeaks(macs2.files, pval = 5, blacklist=blacklist, barplot=TRUE, fileType="macs2")
  } else if (set == "macs2qval")
  {
    allPeaks <- filterPeaks(macs2.files, qval = 10, blacklist=blacklist, barplot=TRUE, fileType="macs2")
    allPeaks.50 <- filterPeaks(macs2.files, qval = 5, blacklist=blacklist, barplot=TRUE, fileType="macs2")
  } else if (set == "rep")
  {
    allPeaks <- filterPeaks(rep.files, blacklist=blacklist, barplot=TRUE)
    allPeaks.50 <- filterPeaks(rep.files, pval=50, blacklist=blacklist)
  }
  
  promoter <- getPromoters(TxDb=txdb, upstream=5000, downstream=5000) #ChIPseeker method
  tagMatrixList <- lapply(allPeaks, getTagMatrix, windows=promoter)
  
  tagMatrixList <- tagMatrixList[!sapply(lapply(tagMatrixList, nrow), is.null)]
  tagMatrixList[sapply(tagMatrixList, nrow) == 0] <- NULL
  
  #promoter <- promoters(genes(txdb), upstream=5000, downstream=5000) #GRanges
  #tagMatrixList <- lapply(allPeaks, getTagMatrix, windows=promoter)
  #plotAvgProf(tagMatrixList, xlim=c(-5000, 5000))
  #tagHeatmap(tagMatrixList, xlim=c(-5000, 5000), color=NULL)
  
  peakAnnoList <- lapply(allPeaks, annotatePeak, TxDb=txdb, tssRegion=c(-2000, 200), verbose=FALSE)
  
  return(list(allPeaks, allPeaks.50, peakAnnoList, tagMatrixList))
  
}


doAllPlot <- function(set, allPeaksList)
{
    allPeaks <- allPeaksList[[1]]
    allPeaks.50 <- allPeaksList[[2]]
    peakAnnoList <- allPeaksList[[3]]
    tagMatrixList <- allPeaksList[[4]]
    
    if (set == "rep")
      {
        plotV(allPeaks, main = "Ints (< 1e-100)") #different order due to error otherwise ...
        plotV(allPeaks.50, main = "Ints (< 1e-50)")
      } else {
    
        plotV(allPeaks, main = "Ints (< 1e-100)") #different order due to error otherwise ...
        plotV(allPeaks.50, main = "Ints (< 1e-50)")
      }

  #covplot(allPeaks[[9]], weightCol="V5", chrs = c(paste0("chr", c(1:19, "X", "Y"))))


  #peakAnnoList <- list("Normal"=annotatePeak(allPeaks[[1]],  TxDb=txdb,tssRegion=c(-2000, 200), verbose=FALSE), "TSS"=annotatePeak(allPeaks[[1]],  TxDb=txdb,tssRegion=c(0, 0), verbose=FALSE), "5kb"=annotatePeak(allPeaks[[1]],  TxDb=txdb,tssRegion=c(-5000,5000), verbose=FALSE), "3kb"=annotatePeak(allPeaks[[1]],  TxDb=txdb, verbose=FALSE))

    print(plotAnnoBar(peakAnnoList))
    print(plotDistToTSS(peakAnnoList))
    
    for (n in names(peakAnnoList))
      {
        par(mfrow=c(1,1))
        vennpie(peakAnnoList[[n]])
        text(x=0, y=-1, n)
        par(mfrow=c(1,1))
        plotAnnoPie(peakAnnoList[[n]])
#        par(mfrow=c(1,1))
#        upsetplot(peakAnnoList[[n]]) #only in TB-3.2.1-dev #  vennpie=TRUE,
      }
}

dowload.file.for.making.TxDb = function()
{
  #txdb = 
    ####  Followings are also from Thomas' codes
    ##blacklist: R slightly different the UCSC command line ...
    #download.file("http://www.broadinstitute.org/~anshul/projects/mouse/blacklist/mm9-blacklist.bed.gz", "mm9-blacklist.bed.gz")
    #download.file("http://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz", "mm9ToMm10.over.chain.gz")
    #system("gunzip mm9ToMm10.over.chain.gz")
    #ch <- import.chain("mm9ToMm10.over.chain")
    #original <- import("mm9-blacklist.bed.gz")
    #blacklist <- liftOver(x=original, chain=ch)
    #blacklist <- unlist(blacklist)
    #saveRDS(file="mm10-blacklist.rds", blacklist)
    #blacklist <- import("mm10-blacklist.bed")
  
  #TxDB
  txdb <- makeTxDbFromUCSC(genome="ce10", tablename="refGene")
  saveDb(txdb, file="refGene_ce10.sqlite")
  #download.file("http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz", "refGene20150925.txt.gz")
  #download.file("http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refFlat.txt.gz", "refFlat20150925.txt.gz")
  #download.file("http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refLink.txt.gz", "refLink20150925.txt.gz")
  #saveDb(txdb, file="refGene20150925.sqlite")
  #export(genes(txdb), "refGene20150925.gene.bed", format='BED')
  #export(promoters(genes(txdb), upstream=0, downstream=0), "refGene20150925.tss.bed", format='BED')
  #txdb <- loadDb("refGene20150925.sqlite")
  #txdb <- makeTxDbFromGFF("/groups/bioinfo/thomas/ucsc/mm10/mm10_refSeq.gene.noRandomChr.gtf", format="gtf") #exclude random chromosomes
  #txdb <- makeTxDbFromGFF(pipe('grep -P -v "chr\\S+_random" /groups/bioinfo/thomas/ucsc/mm10/mm10_refGene'), format="gtf") #exclude random chromosomes
  #txdb
  
}

#for (set in c("macs14", "macs2", "macs2qval"))
#  {
#    pdf(paste0("peaksAnalysis.", set, ".pdf"), width = 10, height=10)
#    doAllPeak(set=set)
#    dev.off()
#  }




#ChipQC
#library(ChIPQC)
#samples <- read.delim("Jorge_Arturo.Zepeda_Martinez.chipqc.txt")
#experiment = ChIPQC(samples)
#experiment
#ChIPQCreport(experiment)
find.samples.conditions = function(x, ID='samples')
{
  if(ID=='samples') 
  { 
    j = 2
    temp = unlist(strsplit(as.character(x), '_'));
    return(temp[j])
  }
  if(ID=='conditions') 
  {
    j = 1;
    temp = unlist(strsplit(as.character(x), '_'));
    temp = temp[j]
    temp = unlist(strsplit(as.character(temp), '[.]'));
    return(temp[length(temp)])
    #kk = which(temp=='macs2')
    #return(paste(temp[c(n1:(kk-n2))], sep='', collapse = '_'))
  }
}

PLOT.Quality.Controls.Summary = function(stat, selected.samples=c(1:nrow(stat)), group.samples.conditions=FALSE)
{
  par(mfrow=c(2,2))
  
  ss = stat[selected.samples, ]
  if(group.samples.conditions)
  {
    samples = sapply(ss$filename, find.samples.conditions, ID='samples')
    conditions =  sapply(ss$filename, find.samples.conditions, ID='conditions')
    ss$conditions = conditions
    ss = data.frame(ss, stringsAsFactors = FALSE)
    ss = ss[with(ss, order(conditions, samples)), ]
    cols.conditions = data.frame(unique(conditions), c(1:length(unique(conditions))))
    cols = cols.conditions[match(ss$conditions, cols.conditions[,1]), 2]
  }else{
    cols = c(1:nrow(ss))
  }

  barplot(as.numeric(ss$nb.UniqueRead)/10^6, col=cols, main="nb of uniquely-mapped Reads", horiz=TRUE, names.arg = ss$filename, 
          las=1, xlim = c(0, 50))
  abline(v=c(10, 20), col='red', lty=1, lwd=2.0);abline(v=c(20, 45), col='red', lty=1, lwd=2.0)
  
  pbc = ss$nb.UniqueRmdupRead/ss$nb.UniqueRead
  barplot(pbc, col=cols, main="NRF", horiz=TRUE, names.arg = ss$filename, 
          las=1, xlim = c(0, 1))
  abline(v=c(0.9), col='red', lty=1, lwd=2.0);
  
  
  plot(ss$NSC, ss$RSC, type='n', xlab='NSC', ylab='RSC', ylim=range(ss$RSC), xlim=range(ss$NSC))
  colfunc <- colorRampPalette(c("red", "green"))
  cols = colfunc(5)
  #jj = which(stat$samples !='Input')
  #stat = stat[jj,]
  for(n in 1:nrow(ss))
  {
    points(ss$NSC[n], ss$RSC[n], type='p', col=cols[(ss$Quality[n]+3)], pch=16, cex=2.0)
    #else  points(stat$NSC.1.05.[n], stat$RSC.0.8.[n], type='p', col=cols[(stat$QualityTag..2.verylow..1.low.0.medium.1.high.2.veryhigh.[n]+3)], pch=17, cex=2.0)
    text(ss$NSC[n], ss$RSC[n], n, cex=1.6, offset=0.5, pos = 1)
  }
  abline(h=c(0.8), col='red', lty=3, lwd=2.0);abline(v=c(1.05), col='red', lty=3, lwd=2.0)
  legend('topleft', legend=rev(c('very low', 'low', 'medium', 'high', 'very high')), col=rev(cols), pch=16, bty='n', cex=1.5)
  
  plot(ss$NSC, ss$RSC, type='n',ylim=c(0, 2.0), xlim=c(1, 1.5), xlab=NA, ylab=NA, axes=FALSE)
  
  if(nrow(ss)<=40){
    legend('topleft', legend=paste(c(1:nrow(ss)), ss$filename, sep='-'), col=cols[ss$Quality+3], pch=16, bty='n', cex=1.2)
  }else{
    legend('topleft', legend=paste(c(1:40), ss$filename[1:40], sep='-'), col=cols[ss$Quality+3][1:40], pch=16, bty='n', cex=1.)
    legend('topright', legend=paste(c(41:nrow(ss)), ss$filename[41:nrow(ss)], sep='-'), col=cols[ss$Quality+3][41:nrow(ss)], pch=16, bty='n', cex=1.2)
  }
  
}

PeakAnnotation.customized = function(pp, window.size=2000, annotation='wormbase')
{
  peaks = pp;
  ### import Refseq annotation
  #annot = read.table('/Volumes/groups/cochella/jiwang/Projects/Julien/R4224/Refseq_annotation.txt', header = TRUE, sep='\t')
  #save(annot, file='Refseq_annotation.Rdata')
  chrom.size = read.delim('ce10_chrom.sizes', sep='\t', header = FALSE)
  if(annotation=='wormbase')
  {
    #annot = read.delim('/Volumes/groups/cochella/jiwang/annotations/BioMart_WBcel235.txt', sep='\t', header = TRUE)
    #colnames(annot)[c(3,)]
    #save(annot, file='/Volumes/groups/cochella/jiwang/annotations/BioMart_WBcel235.Rdata')
    load(file='/Volumes/groups/cochella/jiwang/annotations/BioMart_WBcel235.Rdata')
    ## filter the annotation
    sels = which(annot$Chromosome.scaffold.name != "MtDNA" & annot$Status..gene.=="KNOWN");
    annot = annot[sels, ]
    mm = match(unique(annot$Gene.stable.ID), annot$Gene.stable.ID)
    annot = annot[mm, ]
    aa = data.frame(annot$Chromosome.scaffold.name,
                    annot$Gene.Start..bp., annot$Gene.End..bp., 
                    annot$Strand, annot$Gene.stable.ID, annot$Gene.name, stringsAsFactors = FALSE);
    colnames(aa) = c('chr', 'start', 'end', 'strand', 'wormbase.ID', 'gene')
    aa$strand[which(aa$strand>0)] = '+';
    aa$strand[which(aa$strand<0)] = '-';
    aa$chr = paste0('chr', aa$chr)
  }
  if(annotation=='refseq')
  {
    load(file='Refseq_annotation.Rdata');
    aa = data.frame(annot$ce10.refGene.chrom, 
                    annot$ce10.refGene.txStart, annot$ce10.refGene.cdsEnd, 
                    annot$ce10.refGene.strand, annot$ce10.refGene.name, annot$ce10.refGene.name2, stringsAsFactors = FALSE);
    colnames(aa) = c('chr', 'start', 'end', 'strand', 'RefID', 'gene')
  }
  
  aa = makeGRangesFromDataFrame(aa, keep.extra.columns=TRUE)
  
  targets = NULL;
  for(n in 1:length(peaks))
  #for(n in 1:100)
  {
    #n = 2
    cat('peak index ....', n, '\n')
    px = data.frame(peaks[n], stringsAsFactors = FALSE); 
    
    ll = chrom.size[which(chrom.size[,1]==px[,1]),2]
    wsize = window.size;
    find = 1
    while(find>0)
    {
      start = px[,2]-wsize;end=px[,3]+wsize;
      
      if(start<0) start=0;
      if(end>ll) end=ll;
      dd = data.frame(px[,1], start, end, stringsAsFactors = FALSE); colnames(dd) = c('chr', 'start', 'end')
      newp = makeGRangesFromDataFrame(dd)
      
      if(overlapsAny(newp, aa))
      {
        tt = aa[overlapsAny(aa, newp)]
        for(m in 1:length(tt))
        {
          targets = rbind(targets, data.frame(data.frame(px), data.frame(tt[m]), window.size = wsize, stringsAsFactors = FALSE))
        }
       find = -1; 
      }else{
        wsize = wsize + window.size;
      }
    }
    #resize(px, fix = '', width = 10)
    #reduce(px, drop.empty.ranges=FALSE, min.gapwidth=1L, with.revmap=FALSE, with.inframe.attrib=FALSE, ignore.strand=FALSE)
  }
  
  return(targets)
}

Peak.GPS = function(pp)
{
  #pp = res[407, ]
  summit = floor(mean(as.numeric(pp[match(c('start.peak', 'end.peak'), names(pp))])));
  tss = as.numeric(pp[which(names(pp)=='TSS')]);
  strand = pp[which(names(pp)=='strand.gene')];
  start = as.numeric(pp[which(names(pp)=='start.gene')]);
  end = as.numeric(pp[which(names(pp)=='end.gene')]);
  if(is.na(tss)){
    if(strand=='+'){tss = start;
    }else{tss = end;}
  }
  dist.tss = summit - tss; if(strand =='-') dist.tss = - dist.tss;
  promoter.1kb = (summit< (tss+1000) & summit> (tss-1000))
  promoter.2kb = ((summit< (tss+2000) & summit > (tss+1000)) | (summit> (tss-2000) & summit<(tss-1000)))
  gene.body = (summit< end & summit> start)
  if(strand=='+')
  {
    dowstream.3kb = (summit > end & summit < (end+3000))
  }else{
    dowstream.3kb = (summit > (start -3000) & summit < (start))
  }
  test = c(promoter.1kb, promoter.2kb, gene.body, dowstream.3kb)
  if(sum(test)==0) {test = c(test, TRUE) 
  }else {test = c(test, FALSE)}
  
  return(c(dist.tss, test))
}


