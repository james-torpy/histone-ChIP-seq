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
sampleName <- "chapman-rothe_2013"
descrip <- paste0(sampleName, "_enrichment")

pos_ctl <- read.table(paste0(refDir, "/HGSOC_H3K27me3_H3K4me3_pos_ctls_curry2017.txt"),
                      header=F)[,1]
neg_ctl <- c("GAPDH", "CD47", "CCNE1")

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/", expName, "/", sampleName, "/")
resultsDir <- paste0(projectDir, "/results")
refDir <- paste0(projectDir, "/refs/")
RobjectDir <- paste0(projectDir, "/Robjects/",
                     expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/plots/", descrip, "/")

inDir <- paste0(resultsDir, "/macs2/")

DE_dir <- paste0(homeDir,
                 "projects/hgsoc_repeats/RNA-seq/Robjects/exp9//DEverify/htseq_EdgeR_primary_HGSOC_vs_FT/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", RobjectDir))


##########################################################################
### 1. Load in peak files ###
##########################################################################

if ( !file.exists(paste0(RobjectDir, "peak_gr.rds")) ) {
  
  print("Creating peaks GRanges object...")
  inFiles <- grep(
    "gapped", list.files(inDir, pattern = "Peak", full.names = T), invert=T, value = T
  )
  
  for ( i in 1:length(inFiles) ) {
    print(i)
    id <- gsub(
      "_peaks.*$", "", basename(inFiles[i])
    )
    if (i==1) {
      peaks <- list(read.table(inFiles[i], header = F, sep = "\t"))
      names(peaks)[i] <- id
    } else {
      peaks[[i]] <- read.table(inFiles[i], header = F, sep = "\t")
      names(peaks)[i] <- id
    }
  }
  
  #convert the inFile to GRanges object:
  peak_gr <- lapply(peaks, function(x) {
    return(
      GRanges(
        seqnames = x$V1,
        ranges = IRanges(start=x$V2, 
                         end=x$V3),
        strand = rep("*", nrow(x)),
        read_id  = x$V4
      )
    )
  })
  saveRDS(peak_gr, paste0(RobjectDir, "/peak_gr.rds"))
  
} else {
  print("Loading peaks GRanges object...")
  peak_gr <- readRDS(paste0(RobjectDir, "/peak_gr.rds"))
}



##########################################################################
### 2. Load in repeat and gencode annotations ###
##########################################################################

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
  
  print("Creating expanded repeat annotation...")
  # load repeats annotation:
  rp_annot <- import(paste0(refDir, "custom3rep.final.gff"))
  
  # expand ranges by 1000 bp either side:
  rp_annot <- exp_annot(rp_annot, 1000)
  
  # split annot into GRangesLists by IDs/names:
  rp_annot <- split(rp_annot, rp_annot$ID)
  
  saveRDS(rp_annot, paste0(RobjectDir, "/rp_exp.rds"))
} else {
  print("Loading expanded repeat annotation...")
  rp_annot <- readRDS(paste0(RobjectDir, "/rp_exp.rds"))
}

if ( !file.exists(paste0(RobjectDir, "/ctls_exp.rds")) ) {
  
  print("Creating expanded control annotation...")
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
  
  # split annot into GRangesLists by IDs/names:
  ctls <- split(ctls, ctls$ID)
  
  saveRDS(ctls, paste0(RobjectDir, "/ctls_exp.rds"))
  
} else {
  print("Loading expanded control annotation...")
  ctls <- readRDS(paste0(RobjectDir, "/ctls_exp.rds"))
}

save.image(file = paste0(RobjectDir, "/annot_expanded.rds"))


##########################################################################
### 3. Calculate log odds ratios of H3K27me3 enrichment in repeat and
# control regions ###
##########################################################################

# fetch DE genes:
DE <- readRDS(file = paste0(DE_dir, "custom3_DEsigReps_names.rds"))

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
chr_gr<-GRanges(seqnames=names(chrs),IRanges(1,chrs))

oddsRatio <- function(database, region){
  
  # for a given nucleotide in region, what are the odds it overlaps with 
  # database ranges:
  region_length <- sum(as.numeric(width(region)))
  reduced_database=reduce(database)
  
  strand(region)="*"
  strand(reduced_database)="*"
  
  # calculate the length of all region hits with database:
  region_hits <- sum(width(intersect(region,reduced_database)))
  
  # calculate odds of a region hit with database:
  first_ratio <- region_hits/(region_length - region_hits)
  
  # fetch regions other than region hits:
  non_region <- sum(as.numeric(width(chr_gr))) - region_length
  
  # fetch non-region hits with database:
  non_region_hits <- 
  	sum(as.numeric(width(reduced_database))) - region_hits

  # calculate odds of a non-region hit with database:
  second_ratio <- non_region_hits/(non_region - non_region_hits)
  
  odds <- first_ratio/second_ratio
  
  # calculate standard error for odds:
  se <- sqrt(1/region_length + 1/region_hits + 1/non_region + 
  	1/non_region_hits)
  cat(odds)
  cat("\n")
  cat(se)
  cat("\n")
  result<-list()
  result[[1]]=odds
  result[[2]]=se
  
  # calculate p-value for odds:
  database_enrichment <-
    matrix(c(as.numeric(region_hits), as.numeric(region_length), 
    	as.numeric(non_region_hits), as.numeric(non_region)),
           nrow = 2,
           dimnames = list(Regions = c("Intersecting_repeats", 
           	"Total"),
                           Genome = c("Intersecting_repeats", 
                           	"Total")))
  result[[3]]=chisq.test(database_enrichment)$p.value
  
  return(result)
}

# apply odds ratio function to each element:
rp_odds_peaks1 <- lapply(rp_annot, oddsRatio, peak_gr[[1]])
saveRDS(rp_odds_peaks1, file = paste0(RobjectDir, "/repeats_", 
	names(peak_gr)[1], "_overlap_odds.rds"))

rp_odds_peaks2 <- lapply(rp_annot, oddsRatio, peak_gr[[2]])
saveRDS(rp_odds_peaks2, file = paste0(RobjectDir, "/repeats_", 
	names(peak_gr)[2], "_overlap_odds.rds"))

ctl_odds_peaks1 <- lapply(ctls, oddsRatio, peak_gr[[1]])
saveRDS(ctl_odds_peaks1, file = paste0(RobjectDir, "/repeats_", 
	names(peak_gr)[1], "_ctl_odds.rds"))

ctl_odds_peaks2 <- lapply(ctls, oddsRatio, peak_gr[[2]])
saveRDS(ctl_odds_peaks2, file = paste0(RobjectDir, "/repeats_", 
	names(peak_gr)[2], "_ctl_odds.rds"))


########################################################################
### 3. Plot log odds ratios of H3K27me3 enrichment in repeat regions #
########################################################################

# create function to filter rp_odds by selected p-value and plot them 
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
      pdf(paste0(plotDir, "H3K27me3_enrichment_positive_ctl_log_odds_pval", pval, "DE.pdf"))
    } else {
      pdf(paste0(plotDir, "H3K27me3_enrichment_positive_ctl_log_odds_pval", pval, ".pdf"))
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
                 pdf(paste0(plotDir, "/H3K27me3_enrichment_log_odds_pval", pval, DE.pdf"))
                 } else {
                   pdf(paste0(plotDir, "/H3K27me3_enrichment_log_odds_pval", pval, pdf"))
                 }
                              pdf(paste0(plotDir, "/H3K27me3_enrichment_log_odds_pval", pval, pdf"))
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
  
} else if ( bam_or_peak == "peak" ) {
  rp_olaps <- findOverlaps(rp_annot, data_gr)
  ctl_olaps <- findOverlaps(ctls, data_gr)
  
}

