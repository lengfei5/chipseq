library("ChIPseeker");
library("rtracklayer")
#library(UpSetR);library("ChIPpeakAnno")
library("Vennerable") #install.packages("Vennerable", repos="http://R-Forge.R-project.org", type="source")
library("ggplot2");
library("GenomicFeatures");
library("GenomicRanges")
#library("DiffBind");library("IRanges");
source("functions_chipSeq.R") ## these functions are from Thomas Burkard

### Import peaks and design matrix
version.analysis = 'peaks_GCc_merged_R6118'
peakDir = "../../NGS_requests/R6118_chipseq/Peaks/macs2_broad"

#design.file = ""
outDir = "../Results/Peaks_R6118"
if(!dir.exists(outDir)) dir.create(outDir)

xlist<-list.files(path=peakDir, pattern = "*macs2_broad_fdr_0.1_peaks.xls", full.names = TRUE)
#xlist = xlist[grep('_H', xlist)]

### select samples to check using design matrix
Import.Design.Matrix = FALSE
if(Import.Design.Matrix){
    ## improt design matrix
    design = read.delim(paste0("ChIPseq_design_matrix_all.txt"), header = TRUE, as.is = 2)
      index = c()
      for(n in 1:nrow(design))
        {
              jj = grep(design$SampleID[n], xlist)
                  if(length(jj)==1){index = c(index, jj)
                                  }else{cat("NOT FOUND sample for ", design$sampleID[n], "\n")}
            }
      macs2.files = xlist[index]
      bname = basename(macs2.files)
      design.matrix = data.frame(macs2.files, bname, design, stringsAsFactors = FALSE)
      colnames(design.matrix) = c('file.path', 'file.name', 'sampleID', 'condition', 'factor')
      #ff = data.frame(ff, sapply(bname, find.samples.conditions, ID='samples'))
      o1 = with(design.matrix, order(factor, condition))
      design.matrix = design.matrix[with(design.matrix, order(factor, condition)), ];
  }else{
      ## make design matrix based on the peak files
      macs2.files = xlist;
        bname = basename(macs2.files)
        xx = c()
        yy = c()
        for(n in 1:length(bname))
          {
                #n = 1;
                test = bname[n];
                    test = unlist(strsplit(as.character(test), "_"))
                    xx = c(xx, test[1])
                    #test = unlist(strsplit(as.character(test), "[.]"))[-1]
                    #test = paste0(test, collapse = ".")
                    yy = c(yy, test[2])
                    #bb = sapply(bname, function(x) gsub("..", ".", x))
              }
        design = data.frame(yy, xx)

        #design = data.frame(sapply(bname, find.samples.conditions, ID='conditions'), sapply(bname, find.samples.conditions, ID='samples'),
        #                    stringsAsFactors = FALSE)
        design.matrix = data.frame(macs2.files, bname, design, stringsAsFactors = FALSE)
        colnames(design.matrix) = c('file.path', 'file.name', 'condition', 'factor')
        design.matrix = design.matrix[with(design.matrix, order(factor, condition)), ];
        #condition= sapply(bname, find.sample.names)
        #macs2.files <- paste("/Volumes/groups/bell/jiwang/Projects/Jorge/ChIP_seq_Jorge/peakcalling/macs2/", xlist, sep='')
    }

factor.condition = paste0(design.matrix$factor, '_', design.matrix$condition)
design.matrix = data.frame(design.matrix, factor.condition, stringsAsFactors = FALSE)

##################################################
##################################################
## Section: Quality control after peak calling by macs2
## check the peak overlapping
##################################################
##################################################
## import macs peaks as Genomic Range object
peaks.list = c()
for(n in 1:nrow(design.matrix)){
    cat(n, '\n')
      xx = readPeakFile(design.matrix$file.path[n], as = "GRanges");
      peaks.list = c(peaks.list, xx)
      #eval(parse(text = paste0("pp.", n, " = xx")))
  }

pdf(paste0(outDir, "/Comparison_Peaks_All_Replicates_overlapping_macs2_", version.analysis, ".pdf"), width = 12, height = 8)

source("functions_chipSeq.R") #
#par(mfrow=c(1,2))
Comparison.overlapping.peaks(design.matrix, peaks.list, toCompare="factor.condition")

dev.off()

##################################################
##################################################
## Section: Peak extension, Peak-to-gene assignment, save tables and annotation plots
##################################################
##################################################
ChIPseq.peak.gene.assignment = TRUE
Save.Peak.Annotation = TRUE;

Chr2save = paste0('chr', c(1:19, 'X', 'Y'));
chroms = Chr2save;
Merge.close.peaks = TRUE
merge.dist = 2000;

if(ChIPseq.peak.gene.assignment)
  {
      library(csaw)
        library(AnnotationDbi)
        library(AnnotationHub)

        ## import gene annotation
        #txdb <- makeTxDbFromUCSC(genome="mm10", tablename="refGene")
        #saveDb(txdb, file="refGene_mm10.sqlite")
        txdb <- loadDb("../refGene_mm10.sqlite")

        ah = AnnotationHub()
        #ce = ah[["AH52238"]]
        #mm <- ah[ah$species=="Mus musculus" & ah$rdataclass=='OrgDb']
        mm = ah[["AH57974"]]
        entrez2symbol <- unique(select(mm, keys=keys(mm, keytype="ENTREZID"), columns=c("SYMBOL"), keytype="ENTREZID"))

        #blacklist = readPeakFile(, as = "GRanges")
        pdfname = paste0(outDir, "/Peaks_Annotation_",version.analysis, ".pdf")
        pdf(pdfname, width = 10, height = 6)

        for(n in 1:nrow(design.matrix))
          {
                # n = 1
                cat(design.matrix$file.name[n], "\n")

                    #grep('H3K27me3', macs2.files)
                    #n = 6;
                    p = readPeakFile(design.matrix$file.path[n], as = "GRanges");
                    #p10 <- p[mcols(p)[,"X.log10.pvalue."] > pval];
                    cat(length(p), '\n')
                    if(design.matrix$factor[n] == "H2AK119" | design.matrix$factor[n] == "H2AUb"){
                            merged <- mergeWindows(p, tol=5000L, ignore.strand = TRUE)
                          }else{
                                  merged <- mergeWindows(p, tol=2000L, ignore.strand = TRUE)
                                }
                    cat(length(merged$region))
                    dmeta = as.data.frame(p)[, -c(1:7)]
                    dmeta[, 1] = 10^(-dmeta$X.log10.pvalue.); dmeta[, 2] = log2(dmeta[,2])
                    colnames(dmeta) = c("PValue", "log.fc", "log10.qvalue", "name")
                    tabcom <- combineTests(merged$id,  dmeta, pval.col = 1, fc.col = 2 )

                    pp = merged$region
                    elementMetadata(pp) = tabcom
                    #p = reduce(p); p10 = reduce(p10);
                    #cat(length(p), length(p10), '\n')
                    #pp = reduce(pp, drop.empty.ranges=FALSE, min.gapwidth=1000L, with.revmap=TRUE, with.inframe.attrib=FALSE)

                    peakAnnots = annotatePeak(pp, TxDb=txdb, tssRegion = c(-3000, 3000))
                    #xx = data.frame(peakAnnots)
                    print(plotAnnoBar(peakAnnots, title = design.matrix$file.name[n]))
                    print(plotDistToTSS(peakAnnots, title = design.matrix$file.name[n]))

                    annotatedPeak = as.data.frame(peakAnnots)
                    ### Use the ChIPpeakAnno (we did NOT use here)
                    #annotatedPeak <- annotatePeakInBatch(allPeaksList[[1]][[1]], AnnotationData = genes(txdb), output="both")
                    #annotatedPeak <- annotatePeakInBatch(allPeaks[[n]], AnnotationData = promoters(genes(txdb), upstream=0, downstream=0))
                    ii = grep('geneId', colnames(annotatedPeak))
                    mm = match(annotatedPeak[,ii], entrez2symbol$ENTREZID)
                    df = data.frame(annotatedPeak, entrez2symbol[mm, ], stringsAsFactors = FALSE)
                    kk = which(!is.na(match(df[,1], chroms)));
                    df = df[kk, ]

                    name.annotation = paste0(outDir, "/", sub('\\.xls$', '', design.matrix$file.name[n]), '_combineTests_geneAssignment.txt')
                    write.table(df, file = name.annotation, row.name = FALSE, col.name = TRUE, quote=FALSE, sep='\t')

              }

        dev.off()

    }
