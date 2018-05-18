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
expName <- "FT_SAMN03496012"
sampleTypes <- c("H3K4me1", "H3K27me3")
Type <- "all"

pos_ctl <- ""
neg_ctl <- ""
pVal <- 0.1
peak_type <- "narrow"
diff_exp <- TRUE

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/histone-ChIP-seq/")
resultsDir <- paste0(projectDir, "/results")
refDir <- paste0(projectDir, "/refs/")


inDir <- paste0(resultsDir, "/macs2/", expName, "/")

sampleNames <- list.files(inDir, pattern = "Peak")
sampleNames <- gsub(
  "_peaks.*$", "", sampleNames
)

for ( s_name in sampleNames) {
  
  descrip <- paste0(s_name, "_enrichment")
  RobjectDir <- paste0(projectDir, "/Robjects/histone-ChIP-seq/", expName, "/", descrip, "/")
  plotDir <- paste0(resultsDir, "/R/plots/", expName, "/", descrip, "/")
  
  system(paste0("mkdir -p ", plotDir))
  system(paste0("mkdir -p ", RobjectDir))
  
  ##########################################################################
  ### 1. Load in bam, repeat and gencode annotations ###
  ##########################################################################
  
  if ( peak_type == "narrow" ) {
    inFile <- paste0(inDir, s_name, "_peaks.narrowPeak")
  } else if ( peak_type == "broad" ) {
    inFile <- paste0(inDir, s_name, "_peaks.broadPeak")
  }
  
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
  
  if ( !exists("rp_annot") ) {
    if ( !file.exists(paste0(RobjectDir, "/rp_exp.rds")) ) {
      # load repeats annotation:
      rp_annot <- import(paste0(refDir, "custom3rep.final.gff"))
      
      # expand ranges by 500 bp either side:
      rp_annot <- exp_annot(rp_annot, 500)
      
      # split annot into GRangesLists by IDs/names:
      #rp_annot <- split(rp_annot, rp_annot$ID)
      
      saveRDS(rp_annot, paste0(RobjectDir, "/rp_exp.rds"))
    } else {
      rp_annot <- readRDS(paste0(RobjectDir, "/rp_exp.rds"))
    }
  }
  
  if ( !exists("ctls") ) {
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
  }
  
  # fetch DE genes:
  if ( !exists("DE") ) {
    DE <- readRDS(file = paste0(homeDir,
                                "projects/hgsoc_repeats/RNA-seq/Robjects/exp9/",
                                "DEverify/htseq_EdgeR_primary_HGSOC_vs_FT/", 
                                "custom3_DEsigReps_names.rds"))
  }
  
  
  
  ##########################################################################
  ### 2. Calculate fold change from random regions of peaks in repeats and
  # controls ###
  ##########################################################################
  
  #define standard error function:
  std_e <- function(x) sd(x)/sqrt(length(x))
  
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
      return(data.frame(mean(x$fold_change), std_e(x$fold_change), mean(x$pval)))
    })
    
    df <- do.call("rbind", fc)
    df$ID <- rownames(df)
    
    return(df)
  }
  
  rp_df <- mean_peaks(rp_annot, data_gr, Pval = pVal)
  ctl_df <- mean_peaks(ctls, data_gr, Pval = pVal)
  
  for ( i in 1:length(ctls)) {
    if ( !(ctls[i]$ID %in% ctl_df$ID) ) {
      ctl_df[nrow(ctl_df)+1,] <- c(0, NA, NA, ctls[i]$ID)
    }
  }
  
  pos_df <- ctl_df[ctl_df$ID %in% pos_ctl,][1:3,]
  pos_df$type <- "positive"
  neg_df <- ctl_df[ctl_df$ID %in% neg_ctl,]
  neg_df$type <- "negative"
  rp_df$type <- "repeat"
  
  # add identifier column to each data frame and bind them
  # together:
  
  if ( diff_exp ) {
    rp_df <- rp_df[rp_df$ID %in% DE,]
    
    all_df <- rbind(pos_df, neg_df, rp_df)
    colnames(all_df) <- c("mean_fold_change", "fc_standard_error", "p-value", "ID", "type")
    all_df$mean_fold_change <- round(as.numeric(all_df$mean_fold_change), 3)
    all_df$fc_standard_error <- round(as.numeric(all_df$fc_standard_error), 3)
    all_df$fc_standard_error[is.na(all_df$fc_standard_error)] <- 0
    
    all_df$ID <- factor(all_df$ID, levels = c(neg_ctl, pos_df$ID, rp_df$ID))
    
    ##########################################################################
    ### 3. Plot histone methylation enrichment in repeat regions ###
    ##########################################################################
    
    # plot fold changes:
    pdf(paste0(plotDir, "fold_change_from_random_regions_p", pVal, "_DE.pdf"))
    p <- ggplot(all_df, aes(x = ID, y = mean_fold_change))
    p <- p + geom_bar(aes(fill = type), stat = "identity")
    p <- p + geom_errorbar(aes(ymin=mean_fold_change+fc_standard_error, 
                               ymax=mean_fold_change-fc_standard_error))
    #p <- p + coord_cartesian(ylim = c(0.1, 10)) 
    p <- p + ylab("fold_change_from_random_regions")
    p <- p + xlab("repeat loci")
    #p <- p + scale_y_log10()
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                              vjust = 0.6))
    print(p)
    dev.off()
    
  } else {
    all_df <- rbind(pos_df, neg_df, rp_df)
    colnames(all_df) <- c("mean_fold_change", "fc_standard_error", "p-value", "ID", "type")
    all_df$mean_fold_change <- round(as.numeric(all_df$mean_fold_change), 3)
    all_df$fc_standard_error <- round(as.numeric(all_df$fc_standard_error), 3)
    all_df$fc_standard_error[is.na(all_df$fc_standard_error)] <- 0
    
    
    ##########################################################################
    ### 3. Plot histone methylation enrichment in repeat regions ###
    ##########################################################################
    
    # plot fold changes:
    pdf(paste0(plotDir, "fold_change_from_random_regions_p", pVal, ".pdf"))
    p <- ggplot(all_df, aes(x = ID, y = mean_fold_change))
    p <- p + geom_bar(aes(fill = type), stat = "identity")
    p <- p + geom_errorbar(aes(ymin=mean_fold_change+fc_standard_error, 
                               ymax=mean_fold_change-fc_standard_error))
    #p <- p + coord_cartesian(ylim = c(0.1, 10)) 
    p <- p + ylab("fold_change_from_random_regions")
    p <- p + xlab("repeat loci")
    #p <- p + scale_y_log10()
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                              vjust = 0.6))
    print(p)
    dev.off()
  }
}
