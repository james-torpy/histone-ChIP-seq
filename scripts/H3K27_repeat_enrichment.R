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
### 1. Load in bam and repeat annotations ###
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

rp_annot <- import(paste0(refDir, "custom3rep.final.gff"))

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
### 2. Count overlaps between bam and repeats annotation ###
##########################################################################

olaps <- findOverlaps(bam_gr, rp_annot)

# add hits to bam_gr:
qHits <- queryHits(olaps)
sHits <- subjectHits(olaps)
    
# create repeatHits entry and add hits to olaps elements:
bam_gr$repeatHits <- rep(NA, length(bam_gr))
bam_gr$repeatHits[qHits] <- as.character(rp_annot[sHits]$ID)
    
# create repeat_peaks gr by removing NAs, dust and trf hits:
repeat_peaks <- bam_gr[!is.na(bam_gr$repeatHits)]
repeat_peaks <- repeat_peaks[repeat_peaks$repeatHits!="dust" &
  repeat_peaks$repeatHits!="trf"]
    
# count occurance of peaks in each repeat:
peak_nos <- table(repeat_peaks$repeatHits)

saveRDS(peak_nos, file = paste0(RobjectDir, "/repeat_peaks.rds"))


##########################################################################
### 3. Create fake repeats annotation for repeats with hits ###
##########################################################################

# redistribute chromosome names randomly to create fake annotations:
for ( i in 1:10 ) {
  fakie <- rp_annot
  seqnames(fakie) <- sample(seqnames(fakie))
  if (i==1) {
    fake_annot <- list(fakie)
  } else {
    fake_annot[[i]] <- fakie
  }
}

# check if any ranges fall outside of new chromosome lengths:
writeLines("\n")
print("Do any ranges fall outside of new chromosome lengths?")
lapply(fake_annot, any(end(ranges(fake_annot)) > fake_annot$seq_lengths))

saveRDS(fake_annot, paste0(RobjectDir, "/fake_rp_exp.rds"))


##########################################################################
### 4. Count overlaps between bam and fake annotations ###
##########################################################################

fake_bam_gr <- bam_gr

fake_peak_nos <- lapply(fake_annot, function(x) {

  fake_olaps <- findOverlaps(fake_bam_gr, x)

  fake_qHits <- queryHits(fake_olaps)
  fake_sHits <- subjectHits(fake_olaps)
    
  # create repeatHits entry and add hits to bam:

  fake_bam_gr$repeatHits <- rep(NA, length(fake_bam_gr))
  fake_bam_gr$repeatHits[fake_qHits] <- as.character(rp_annot[fake_sHits
    ]$ID)
    
  # create repeat_peaks gr by removing NAs, dust and trf hits:
  fake_repeat_peaks <- fake_bam_gr[!is.na(fake_bam_gr$repeatHits)]
  fake_repeat_peaks <- fake_repeat_peaks[fake_repeat_peaks$repeatHits!=
    "dust" & fake_repeat_peaks$repeatHits!="trf"]
    
  # count occurance of peaks in each repeat:
  return(table(fake_repeat_peaks$repeatHits))
})

# return the element-wise means of all fake repeat peaks dfs
fake_peak_means <- colSums(laply(fake_peak_nos, as.matrix))/length(fake_peak_nos)

saveRDS(fake_peak_means, file = paste0(RobjectDir, "/fake_repeat_peaks.rds"))


##########################################################################
### 4. Plot real and fake H3K27me3 peaks of DE repeats  ###
##########################################################################

DE <- rownames(readRDS(paste0(DE_dir, "/custom3_DEsigReps.rds"))[[1]])

DE_peak_nos <- data.frame(peak_nos[DE], fake_peak_means[DE])
colnames(DE_peak_nos) <- c("repeat_id", "repeat_nos", "random_nos")

melted_peaks <- melt(DE_peak_nos, id.vars="repeat_id")

pdf(paste0(plotDir, "/DE_repeat_peak_counts.pdf"))
p <- ggplot(melted_logs, aes(x = melted_peaks$repeat_id, y = melted_peaks$value))
p <- p + geom_bar(aes(fill = variable), position = "dodge", stat = "identity")
p <- p + ylab("Peak numbers")
p <- p + xlab("Repeat loci")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
p
dev.off()

pdf(paste0(plotDir, "/DE_repeat_peak_log10_counts.pdf"))
p <- ggplot(melted_logs, aes(x = melted_peaks$repeat_id, y = melted_peaks$value))
p <- p + geom_bar(aes(fill = variable), position = "dodge", stat = "identity")
p <- p + ylab("log10 peak numbers")
p <- p + xlab("Repeat loci")
p <- p + scale_y_log10()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
p
dev.off()
