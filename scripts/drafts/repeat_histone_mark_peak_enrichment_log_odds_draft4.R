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

p_thresh <- 0.1

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/", expName, "/", sampleName, "/")
resultsDir <- paste0(projectDir, "/results")
refDir <- paste0(projectDir, "/refs/")
RobjectDir <- paste0(projectDir, "/Robjects/",
                     expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/plots/")
tableDir <- paste0(resultsDir, "/R/tables/")

inDir <- paste0(resultsDir, "/macs2/")

DE_dir <- paste0(homeDir,
                 "projects/hgsoc_repeats/RNA-seq/results/R/exp9/plots/DEplots/htseq_EdgeR_primary_HGSOC_vs_FT/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", tableDir))
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

if ( !file.exists(paste0(RobjectDir, "/repeats_peak_overlap_odds.rds")) ) {
  
  # fetch summed length of all chromosomes
  chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
  chr_gr<-GRanges(seqnames=names(chrs),IRanges(1,chrs))
  
  oddsRatio <- function(database, region){
    
    # for a given nucleotide in region, what are the odds it overlaps with 
    # database ranges:
    region_length <- sum(as.numeric(width(region)))
    reduced_database <- reduce(database)
    
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
    print("Odds are:")
    cat(odds)
    cat("\n")
    print("Std error is:")
    cat(se)
    cat("\n")
    
    # calculate p-value for odds:
    database_enrichment <-
      matrix(c(as.numeric(region_hits), as.numeric(region_length), 
               as.numeric(non_region_hits), as.numeric(non_region)),
             nrow = 2,
             dimnames = list(Regions = c("Intersecting_repeats", 
                                         "Total"),
                             Genome = c("Intersecting_repeats", 
                                        "Total")))
    pval=chisq.test(database_enrichment)$p.value
    print("P-value is:")
    cat(pval)
    cat("\n")
    
    result <- data.frame(odds, se, pval)
    colnames(result) <- c("odds", "std_error", "pval")
    
    return(result)
  }
  
  # apply odds ratio function to each element:
  for ( j in 1:length(peak_gr) ) {
    if (j==1) {
      
      rp_peak_odds <- list(lapply(rp_annot, oddsRatio, peak_gr[[j]]))
      names(rp_peak_odds)[j] <- names(peak_gr)[j]
      ctl_peak_odds <- list(lapply(ctls, oddsRatio, peak_gr[[j]]))
      names(ctl_peak_odds)[j] <- names(peak_gr)[j]
      
    } else {
      
      rp_peak_odds[[j]] <- lapply(rp_annot, oddsRatio, peak_gr[[j]])
      names(rp_peak_odds)[j] <- names(peak_gr)[j]
      ctl_peak_odds[[j]] <- lapply(ctls, oddsRatio, peak_gr[[j]])
      names(ctl_peak_odds)[j] <- names(peak_gr)[j]
      
    }
  }
  
  # bind odds results for repeats and controls into one data frame per 
  # sample/annotation combination:
  rp_test <- lapply(rp_peak_odds, function(x) {
    return(do.call("rbind", x))
  })
  
  saveRDS(rp_peak_odds, file = paste0(RobjectDir, "/repeats_peak_overlap_odds.rds"))
  saveRDS(ctl_peak_odds, file = paste0(RobjectDir, "/ctl_peak_overlap_odds.rds"))
  
} else {
  rp_peak_odds <- readRDS(file = paste0(RobjectDir, "/repeats_peak_overlap_odds.rds"))
  ctl_peak_odds <- readRDS(file = paste0(RobjectDir, "/ctl_peak_overlap_odds.rds"))
}


########################################################################
### 3. Plot log odds ratios of H3K27me3 enrichment in repeat regions #
########################################################################

# create function to filter odds data frames by selected p-value and plot them 
# with barplot:
plot_odds <- function(odds, pval = p_thresh, sig_only = F, cont = NA) {
  
  # filter by odds:
  odds_DE <- odds[odds$odds != 0,]
  # filter by p-value:
  odds_sig <- odds_DE[odds_DE$pval < p_thresh,]
  
  # if ctl exists include controls for plot:
  if ( is.list(cont) ) {
    
    # fetch bind positive and negative control data frames and annotate
    # before binding with odds_sig:
    odds_sig$type <- "repeat"
    
    pos_temp <- cont[[k]][rownames(cont[[k]]) %in% pos_ctl,]
    pos_temp$type <- "positive_control"
    neg_temp <- cont[[k]][rownames(cont[[k]]) %in% neg_ctl,]
    neg_temp$type <- "negative_control"

    ctl_df <- rbind(pos_temp, neg_temp)
    
    # filter by odds:
    ctl_DE <- ctl_df[ctl_df$odds != 0,]
    # filter by p-value:
    ctl_sig <- ctl_DE[ctl_DE$pval < p_thresh,]
    
    
    # plot controls and log odds:
    if (sig_only) {
      
      all_df <- rbind(ctl_DE, odds_DE)
      all_df$ID <- rownames(all_df)
      
      p <- ggplot(all_df, aes(x = ID, y = odds))
      p <- p + geom_bar(aes(fill = type), stat = "identity")
      p <- p + geom_errorbar(aes(ymin=odds+std_error, 
                                 ymax=odds-std_error))
      p <- p + scale_y_log10()
      p <- p + ylab("log10 odds")
      p <- p + xlab("repeat loci")
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                vjust = 0.6))
      
    } else {
      
      odds_DE$type[odds_DE$pval <= pval] <- 
        paste0(odds_DE$type[odds_DE$pval <= pval], "_significant")
      odds_DE$type[odds_DE$pval > pval] <- 
        paste0(odds_DE$type[odds_DE$pval > pval], "_non-significant")
      
      all_df <- rbind(ctl_DE, odds_DE)
      all_df$ID <- rownames(all_df)

      p <- ggplot(all_df, aes(x = ID, y = odds))
      p <- p + geom_bar(aes(fill = type), stat = "identity")
      p <- p + geom_errorbar(aes(ymin=odds+std_error, 
                                 ymax=odds-std_error))
      p <- p + scale_y_log10() 
      p <- p + ylab("log10 odds")
      p <- p + xlab("repeat loci")
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                vjust = 0.6))
      
    }
  result <- all_df
  
  } else {

    # plot log odds:
    if (sig_only) {
      
      odds_sig$ID <- rownames(odds_sig) 
      
      p <- ggplot(odds_sig, aes(x = ID, y = odds))
      p <- p + geom_bar(stat = "identity")
      p <- p + geom_errorbar(aes(ymin=odds+std_error, 
                                 ymax=odds-std_error))
      p <- p + scale_y_log10()
      p <- p + ylab("log10 odds")
      p <- p + xlab("repeat loci")
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                vjust = 0.6))
      result <- odds_sig
      
    } else {
      
      odds_DE$ID <- rownames(odds_DE) 
      odds_DE$thresh <- NA
      odds_DE$thresh[odds_DE$pval <= pval] <- "significant"
      odds_DE$thresh[odds_DE$pval > pval] <- "non-significant"
      
      p <- ggplot(odds_DE, aes(x = ID, y = odds))
      p <- p + geom_bar(aes(fill = thresh), stat = "identity")
      p <- p + geom_errorbar(aes(ymin=odds+std_error, 
                                 ymax=odds-std_error))
      p <- p + scale_y_log10() 
      p <- p + ylab("log10 odds")
      p <- p + xlab("repeat loci")
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                vjust = 0.6))
      result <- odds_DE
    }
  }
k <<- k+1
return(list(p, result))
}
  
# plot DE odds with p < 0.05:
k=1
odds_p0.05 <- lapply(rp_peak_odds, plot_odds, pval = 0.05)

# plot DE odds with p < 0.1:
k=1
odds_p0.1 <- lapply(rp_peak_odds, plot_odds, pval = 0.1)

# plot DE odds and controls with p < 0.05:
k=1
odds_ctls_p0.05 <- lapply(rp_peak_odds, plot_odds, pval = 0.05, cont = ctl_peak_odds)

# plot DE odds and controls with p < 0.1:
k=1
odds_ctls_p0.1 <- lapply(rp_peak_odds, plot_odds, pval = 0.05, cont = ctl_peak_odds)

# plot DE odds of repeats with DE RNA expression:
RNA_DE <- read.table(file = paste0(DE_dir, "/sig_rep_list.txt"))[,1]

RNA_DE_odds <- lapply(rp_peak_odds, function(x) {
  return(x[rownames(x) %in% RNA_DE,])
})

k=1
RNA_DE_odds_p0.1 <- lapply(RNA_DE_odds, plot_odds, pval = 0.1)

k=1
RNA_DE_odds_p0.1_conts <- lapply(RNA_DE_odds, plot_odds, pval = 0.1, cont = ctl_peak_odds)

# save selected plots and tables:
pdf(paste0(plotDir, "/HGSOC_H3K27me3_SRR600559_repeat_peak_log_odds_p0.05.pdf"))
odds_p0.05$HGSOC_H3K27me3_SRR600559[[1]]
dev.off()

pdf(paste0(plotDir, "/HGSOC_H3K4me3_SRR600956_repeat_peak_log_odds_p0.05.pdf"))
print(odds_p0.05$HGSOC_H3K4me3_SRR600956[[1]])
dev.off()




pdf(paste0(plotDir, "/HGSOC_H3K27me3_SRR600559_repeat_peak_log_odds_p0.05.pdf"))
odds_p0.05$HGSOC_H3K27me3_SRR600559[[1]]
dev.off()

pdf(paste0(plotDir, "/HGSOC_H3K4me3_SRR600956_repeat_peak_log_odds_p0.05.pdf"))
print(odds_p0.05$HGSOC_H3K4me3_SRR600956[[1]])
dev.off()


