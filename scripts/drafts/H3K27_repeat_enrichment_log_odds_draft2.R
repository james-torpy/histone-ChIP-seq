### H3K27me3_repeat_enrichment.R ###

# This script takes a bam file from a ChIP-seq BWA alignment and overlaps
# the bam GRanges with the Repbase repeats annotation and a fake repeats
# annotation to check for H3K27me3 enrichment in repeat regions


##########################################################################
### 0. Define variables/paths ###
##########################################################################

# load packages needed:
library(tibble)
library(dplyr)
library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)
library("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome)
library(plyr)
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "histone-ChIP-seq"
sampleName <- "HGSOC_SAMN01761041/SRR600559-H3K27me3"
descrip <- paste0(sampleName, "_enrichment")
Type <- "all"


# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/", expName, "/")
resultsDir <- paste0(projectDir, "/results")
refDir <- paste0(projectDir, "/refs/")
RobjectDir <- paste0(projectDir, "/Robjects/",
                     expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/plots/", descrip, "/")

inDir <- paste0(resultsDir, "/bwa/", sampleName, "/")

DE_dir <- paste0(homeDir, 
  "projects/hgsoc_repeats/RNA-seq/Robjects/exp9//DEverify/htseq_EdgeR_primary_HGSOC_vs_FT/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", RobjectDir))


##########################################################################
### 1. Load in bam, repeat and gencode annotations ###
##########################################################################

inFile <- paste0(inDir, "/SRR600559-H3K27me3.sorted.rmdup.bam")

# assign column names for dataframe:
what = c("qname","rname","strand","pos","qwidth")
# flag unmapped sequences to leave out:
flag = scanBamFlag(isUnmappedQuery=FALSE)
# define parameters of bam scan:
param = ScanBamParam(what=what,flag=flag)

#define the inFile:
bam=scanBam(inFile,param=param)

#assign the imported bam file to the variable name 'results':
results=bam[[1]]

# remove all entries with '_' or '-' in the name to get rid of alternate 
# chromosome entries and/or ERCC spike-ins:
correct_chromosomes_underscore=grepl("_",bam[[1]]$rname)
correct_chromosomes_dash=grepl("-",bam[[1]]$rname)
both=correct_chromosomes_underscore|correct_chromosomes_dash
correct_chromosomes=!both

#convert the inFile to GRanges object:
bam_gr=GRanges(
  seqnames = results$rname[correct_chromosomes],
  ranges = IRanges(start=results$pos[correct_chromosomes], 
  width=results$qwidth[correct_chromosomes]),
  strand = results$strand[correct_chromosomes],
  read_id  = results$qname[correct_chromosomes]
)

saveRDS(bam_gr, paste0(RobjectDir, "/SRR600559-H3K27me3_bam_gr.rds"))

# load repeats annotation:
rp_annot <- import(paste0(refDir, "custom3rep.final.gff"))

# load gencode annotation:
gc <- import(paste0(refDir, "gencode_v24_hg38_annotation.gtf"))

# remove unwanted chromosome constructs:
rp_annot <- rp_annot[grep("[0-9].[0-9]|M", seqnames(rp_annot), 
  invert = T)]

# assign length of all chromosomes to the variable name 'seq_lengths':
seq_lengths <- seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]

rp_annot$seq_lengths <- rep(NA, length(rp_annot))

for ( v in names(seq_lengths) ) {
  print(v)
  rp_annot$seq_lengths[as.character(seqnames(rp_annot)) == v] <- seq_lengths[v]
}

# add 500 to the start of each ranges if the start is 500 bp or more from
# the start of the chromosome:
start(ranges(rp_annot))[start(ranges(rp_annot)) >= 500
  ] <- start(ranges(rp_annot))[start(ranges(rp_annot)) >= 500] - 500

# add 500 to the end of each ranges if the end is 500 bp or more from
# the end of the chromosome:
end(ranges(rp_annot))[end(ranges(rp_annot)) <= 
  (rp_annot$seq_lengths - 500)] <- 
    end(ranges(rp_annot))[end(ranges(rp_annot)) <= 
      (rp_annot$seq_lengths - 500)] + 500

saveRDS(rp_annot, paste0(RobjectDir, "/rp_exp.rds"))

save.image(file = paste0(RobjectDir, "/rp_annot_expanded.rds"))


##########################################################################
### 2. Calculate log odds ratios of H3K27me3 enrichment in repeat and
# control regions ###
##########################################################################

olaps <- findOverlaps(bam_gr, rp_annot)

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
chrGR<-GRanges(seqnames=names(chrs),IRanges(1,chrs))

oddsRatio <- function(region,snpDatabase){
  
  #for a given nucleotide in region, what are the odds it has a SNP
  totalNuc<-sum(as.numeric(width(region)))
  reducedsnpDatabase=reduce(snpDatabase)
  
  #mat<-findOverlaps(region,reducedsnpDatabase)
  #total hits is the length of intersect
  strand(region)="*"
  strand(reducedsnpDatabase)="*"
  
  totalHits<-sum(width(intersect(region,reducedsnpDatabase)))
  
  firstRatio<-totalHits/(totalNuc-totalHits)
  
  nonRegion<-sum(as.numeric(width(chrGR)))-totalNuc	
  
  genomeHits<-sum(as.numeric(width(reducedsnpDatabase)))-totalHits
  secondRatio<-genomeHits/(nonRegion-genomeHits)
  
  odds<-firstRatio/secondRatio
  
  se<-sqrt(1/totalNuc+1/totalHits+1/nonRegion+1/genomeHits)
  cat(odds)
  cat("\n")
  cat(se)
  cat("\n")
  result<-list()
  result[[1]]=odds
  result[[2]]=se
  
  #p-value
  SNPenrichment <-
    matrix(c(as.numeric(totalHits), as.numeric(totalNuc), as.numeric(genomeHits), as.numeric(nonRegion)),
           nrow = 2,
           dimnames = list(Regions = c("WithSNP", "Total"),
                           Genome = c("WithSNP", "Total")))
  result[[3]]=chisq.test(SNPenrichment)$p.value
  
  return(result)
  
}

# split into GRangesList by repeat ID:
rp_annot <- split(rp_annot, rp_annot$ID)

# apply odds ratio function to each element:
rp_odds <- lapply(rp_annot, oddsRatio, bam_gr)

saveRDS(rp_odds, file = paste0(RobjectDir, "/rp_odds.rds"))


##########################################################################
### 3. Plot log odds ratios of H3K27me3 enrichment in repeat regions ###
##########################################################################

# create function to filter rp_odds by selected p-value and plot them on
# barplot:
plot_odds <- function(odds, pval) {
  # filter by p-value:
  rp_odds_temp <- lapply(odds, function(x) {
    if ( x[[3]] < pval & x[[3]] != 0 ) {
      return(data.frame(x[[1]], x[[2]], x[[3]]))
    } else {
      return(NULL)
    }
  })
  rp_odds_sig <- rp_odds_temp[!sapply(rp_odds_temp, is.null)]
  
  # collapse list to data frame:
  odds_df <- do.call("rbind", rp_odds_sig)
  colnames(odds_df) <- c("log_odds", "se", "p_value")
  odds_names <- unlist(
    lapply(rp_odds_sig, function(x) {
      return(x$ID[1])
    })
  )
  odds_df$ID <- names(rp_odds_sig)
  
  # plot all significant log odds:
  pdf(paste0(plotDir, "/H3K27me3_enrichment_log_odds_pval", pval, ".pdf"))
  p <- ggplot(odds_df, aes(x = ID, y = log_odds))
  p <- p + geom_bar(stat = "identity")
  p <- p + geom_errorbar(aes(ymin=log_odds+se, 
                             ymax=log_odds-se))
  p <- p + coord_cartesian(ylim = c(0.1, 10)) 
  p <- p + ylab("log odds")
  p <- p + xlab("repeat loci")
  p <- p + scale_y_log10()
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                            vjust = 0.6))
  print(p)
  dev.off()
}

# plot with p<0.05:
plot_odds(rp_odds, 0.05)
# plot with p<0.1:
plot_odds(rp_odds, 0.1)

# keep only DE genes and plot any with p<0.1:
DE <- readRDS(file = paste0("/share/ScratchGeneral/jamtor/",
        "projects/hgsoc_repeats/RNA-seq/Robjects/exp9/",
        "DEverify/htseq_EdgeR_primary_HGSOC_vs_FT/", 
        "custom3_DEsigReps.rds"))
DE <- rownames(DE[[1]])

DE_odds <- rp_odds[DE]
plot_odds(DE_odds, 0.1)

