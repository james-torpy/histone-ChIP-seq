### histone_ChIP-seq_repeat_enrichment.R ###

# This script takes a macs2 output file from macs2 output and overlaps 
# the results with the Repbase repeats annotation  to check for histone 
# methylation enrichment in repeat regions


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
sampleName <- "HGSOC_SAMN01761041/SRR600956-H3K4me3"
descrip <- paste0(sampleName, "_enrichment")
Type <- "all"

pos_ctl <- as.character(read.table(paste0(refDir, "/HGSOC_H3K27me3_H3K4me3_pos_ctls_curry2017.txt"),
                      header=F)[,1][1:50])
neg_ctl <- c("OVGP1", "CLEC4M", "ANXA8")

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/", expName, "/")
resultsDir <- paste0(projectDir, "/results")
refDir <- paste0(projectDir, "/refs/")
RobjectDir <- paste0(projectDir, "/Robjects/",
                     expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/plots/", descrip, "/")

inDir <- paste0(resultsDir, "/macs2/HGSOC_SAMN01761041/SRR600956-H3K4me3/")

DE_dir <- paste0(homeDir, 
                 "projects/hgsoc_repeats/RNA-seq/Robjects/exp9//DEverify/htseq_EdgeR_primary_HGSOC_vs_FT/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", RobjectDir))


##########################################################################
### 1. Load in bam, repeat and gencode annotations ###
##########################################################################

inFile <- paste0(inDir, 
                 'HGSOC_SAMN01761041_SRR600956-H3K4me3_peaks.narrowPeak')

results <- read.table(inFile, header = F, sep = "\t")

#convert the inFile to GRanges object:
data_gr=GRanges(
  seqnames = results$V1,
  ranges = IRanges(start=results$V2, 
                   end=results$V3),
  strand = rep("*", nrow(results)),
  read_id  = results$V4,
  fold_change = results$V7,
  pval = results$V8
)

# define lengths of chromosomes:
seq_lengths <- seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]

# create function to remove unwanted references and add to either side of 
# ranges of annotation:
exp_annot <- function(annot, leng) {
  # remove unwanted chromosome constructs:
  annot <- annot[grep("[0-9].[0-9]|M", seqnames(annot), 
                      invert = T)]
  
  # assign length of all chromosomes to the variable name 'seq_lengths':
  annot$seq_lengths <- rep(NA, length(annot))
  
  for ( v in names(seq_lengths) ) {
    print(v)
    annot$seq_lengths[as.character(seqnames(annot)) == v] <- seq_lengths[v]
  }
  
  # add length to the start of each ranges if the start is that length or more from
  # the start of the chromosome:
  start(ranges(annot))[start(ranges(annot)) >= leng
                       ] <- start(ranges(annot))[start(ranges(annot)) >= leng] - leng
  
  # add length to the end of each ranges if the end is length bp or more from
  # the end of the chromosome:
  end(ranges(annot))[end(ranges(annot)) <= 
                       (annot$seq_lengths - leng)] <- 
    end(ranges(annot))[end(ranges(annot)) <= 
                         (annot$seq_lengths - leng)] + leng
  return(annot)
}


if ( !file.exists(paste0(RobjectDir, "/rp_exp.rds")) ) {
  # load repeats annotation:
  rp_annot <- import(paste0(refDir, "custom3rep.final.gff"))
  
  # expand ranges by 500 bp either side:
  rp_annot <- exp_annot(rp_annot, 500)
  
  # split annot into GRangesLists by IDs/names:
  rp_annot <- split(rp_annot, rp_annot$ID)
  
  saveRDS(rp_annot, paste0(RobjectDir, "/rp_exp.rds"))
} else {
  rp_annot <- readRDS(paste0(RobjectDir, "/rp_exp.rds"))
}

if ( !file.exists(paste0(RobjectDir, "/ctls_exp.rds")) ) {
  # load gencode annotation:
  gc <- import(paste0(refDir, "gencode_v24_hg38_annotation.gtf"))
  
  # expand ranges by 500 bp either side:
  gc <- exp_annot(gc, 500)
  saveRDS(gc, paste0(RobjectDir, "/gc_exp.rds"))
  
  # isolate control ranges from gc:
  ctls <- gc[gc$gene_name %in% pos_ctl|gc$gene_name %in% neg_ctl]
  #values(ctls)[ctls$gene_name %in% pos_ctl]$ctl_type <- "positive"
  #values(ctls)[ctls$gene_name %in% neg_ctl]$ctl_type <- "negative"
  values(ctls) <- subset(values(ctls), select=gene_name)
  colnames(values(ctls)) <- "ID"
  
  saveRDS(ctls, paste0(RobjectDir, "/ctls_exp.rds"))
  
} else {
  ctls <- readRDS(paste0(RobjectDir, "/ctls_exp.rds"))
}

# fetch DE genes:
DE <- readRDS(file = paste0("/share/ScratchGeneral/jamtor/",
                            "projects/hgsoc_repeats/RNA-seq/Robjects/exp9/",
                            "DEverify/htseq_EdgeR_primary_HGSOC_vs_FT/", 
                            "custom3_DEsigReps_names.rds"))


##########################################################################
### 2. Calculate fold change from random regions of peaks in repeats and
# controls ###
##########################################################################

mean_peaks <- function(Annot, Data, Pval=0.05) {
  olaps <- findOverlaps(Annot, Data)
  
  peaks <- Data[subjectHits(olaps)]
  peaks$ID <- Annot$ID[queryHits(olaps)]
  peaks <- unique(peaks)
  sig_peaks <- peaks[10^-(peaks$pval) < Pval]
  all_peaks <- c(sig_peaks, peaks[peaks$ID %in% c(as.character(pos_ctl),
                  neg_ctl)])
  peak_split <- split(all_peaks, all_peaks$ID)
  fc <- lapply(peak_split, function(x) {
    return(data.frame(mean(x$fold_change), mean(x$pval)))
  })
  
  df <- do.call("rbind", fc)
  df$ID <- rownames(df)
  
  return(df)
}

rp_df <- mean_peaks(rp_annot, data_gr)
ctl_df <- mean_peaks(ctls, data_gr)

for ( i in 1:length(ctls)) {
  if ( !(ctls[i]$ID %in% ctl_df$ID) ) {
    ctl_df[nrow(ctl_df)+1,] <- c(0, NA, ctls[i]$ID)
  }
}

pos_df <- ctl_df[ctl_df$ID %in% pos_ctl,][1:3,]
pos_df$type <- "positive"
neg_df <- ctl_df[ctl_df$ID %in% neg_ctl,]
neg_df$type <- "negative"
rp_df$type <- "repeat"

# add identifier column to each data frame and bind them
# together:
all_df <- rbind(pos_df, neg_df, rp_df)


##########################################################################
### 3. Plot histone methylation enrichment in repeat regions ###
##########################################################################








rp_olaps <- findOverlaps(rp_annot, data_gr)


rp_peaks <- data_gr[subjectHits(rp_olaps)]
rp_peaks$ID <- rp_annot$ID[queryHits(rp_olaps)]
rp_peaks <- unique(rp_peaks)
rp_peaks <- rp_peaks[rp_peaks$pval > 1.3]

rp_peak_split <- split(rp_peaks, rp_peaks$ID)
rp_fc <- lapply(rp_peak_split, function(x) {
  return(data.frame(mean(x$fold_change), mean(x$pval)))
})
rp_df <- do.call("rbind", rp_fc)
rp_df$ID <- rownames(rp_df)


ctl_olaps <- findOverlaps(ctls, data_gr)

ctl_peaks <- data_gr[subjectHits(ctl_olaps)]
ctl_peaks$ID <- ctls$ID[queryHits(ctl_olaps)]
ctl_peaks <- unique(ctl_peaks)











# create function to filter peaks by selected p-value and plot them on
# barplot:
plot_odds <- function(odds, pval = 0.05, con = F, diff_ex = F) {
  # filter by p-value:
  odds_temp <- lapply(odds, function(x) {
    if ( x[[3]] < pval & x[[3]] != 0 ) {
      return(data.frame(x[[1]], x[[2]], x[[3]]))
    } else {
      return(NULL)
    }
  })
  odds_sig <- odds_temp[!sapply(odds_temp, is.null)]
  
  # collapse list to data frame:
  odds_df <- do.call("rbind", odds_sig)
  colnames(odds_df) <- c("log_odds", "se", "p_value")
  odds_names <- unlist(
    lapply(odds_sig, function(x) {
      return(x$ID[1])
    })
  )
  odds_df$ID <- names(odds_sig)
  
  # if diff_ex == TRUE, keep only DE genes:
  if (diff_ex) {
    odds_df <- odds_df[rownames(odds_df) %in% DE,]
  }
  
  # if con == TRUE, collect odds for plot:
  if (con) {
    
    # plot all significant log odds:
    if (DE) {
      pdf(paste0(plotDir, "/H3K27me3_enrichment_positive_ctl_log_odds_pval", pval, "_DE.pdf"))
    } else {
      pdf(paste0(plotDir, "/H3K27me3_enrichment_positive_ctl_log_odds_pval", pval, ".pdf"))
    }
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
    
  } else {
    
    ctl_temp <- c(ctl_odds[names(ctl_odds) %in% pos_ctl][1:3], 
                  ctl_odds[names(ctl_odds) %in% neg_ctl])
    ctl_df <- do.call("rbind", ctl_temp)
    colnames(ctl_df) <- c("log_odds", "se", "p_value")
    ctl_names <- unlist(
      lapply(ctl_odds, function(x) {
        return(x$ID[1])
      })
    )
    ctl_df$ID <- ctl_names
    
    # add identifier column to each data frame and bind them
    # together:
    odds_df$type <- "repeat"
    ctl_df$type[1:3] <- "positive_control"
    ctl_df$type[4:nrow(ctl_df)] <- "negative_control"
    all_df <- rbind(ctl_df, odds_df)
    
    # plot all significant log odds:
    if (DE) {
      pdf(paste0(plotDir, "/H3K27me3_enrichment_log_odds_pval", pval, "_DE.pdf"))
    } else {
      pdf(paste0(plotDir, "/H3K27me3_enrichment_log_odds_pval", pval, ".pdf"))
    }
    pdf(paste0(plotDir, "/H3K27me3_enrichment_log_odds_pval", pval, ".pdf"))
    p <- ggplot(all_df, aes(x = ID, y = log_odds))
    p <- p + geom_bar(aes(fill = type), stat = "identity")
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
}

# plot with p<0.05:
plot_odds(rp_odds, 0.05)
# plot with p<0.1:
plot_odds(rp_odds, 0.1)

# plot all controls:
plot_odds(ctl_odds, 0.1, con = TRUE)

# plot DE genes only:
plot_odds(rp_odds, pval = 0.1, diff_ex = T)
