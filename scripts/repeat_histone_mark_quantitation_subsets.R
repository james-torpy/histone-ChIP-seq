### repeat_histone_mark_read_enrichment_log_odds.R ###

# This script takes a bam file from a histone mark ChIP-seq alignment and 
# calculates the log odds of mark enrichment in repeat regions


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
library(plyr)
library(ggplot2)
library(scales)
library(reshape2)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "histone-ChIP-seq"
sampleName <- "chapman-rothe_2013_hg19"
descrip <- paste0(sampleName, "SICER_peak_enrichment")

# specify regions to include for marks (body, upstream, up_and_downstream)
incl_marks <- "body_no_promoter"
# specify how many bp up/downstream to expand repeats annotation by:
exp_no <- 2000

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/", expName, "/", sampleName, "/")
resultsDir <- paste0(projectDir, "/results")
refDir <- paste0(projectDir, "/refs/")
RobjectDir <- paste0(projectDir, "/Robjects/")
newRobjectDir <- paste0(RobjectDir, "/", incl_marks, "/")
plotDir <- paste0(resultsDir, "/R/plots/")
tableDir <- paste0(resultsDir, "/R/tables/")

inDir <- paste0(resultsDir, "/epic/")
DE_dir <- paste0(homeDir,
                 "projects/hgsoc_repeats/RNA-seq/results/R/exp9/plots/DEplots/htseq_EdgeR_primary_HGSOC_vs_FT/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", tableDir))
system(paste0("mkdir -p ", RobjectDir))
system(paste0("mkdir -p ", newRobjectDir))

pos_ctl <- as.character(read.table(file=paste0(DE_dir, "/top_CPMR.txt"))[,1])
neg_ctl <- as.character(read.table(file=paste0(DE_dir, "/bottom_CPMR.txt"))[,1])

# set up parallel workers:
no_cores <- 3
print(paste0("No. cores are: ", no_cores))
cl <- makeCluster(no_cores)             


##########################################################################
### 1. Load in ChIP reads files ###
##########################################################################

if ( !file.exists(paste0(RobjectDir, "peak_gr.rds")) ) {
  
  in_files <- list.files(inDir, pattern = ".bed", full.names = T, 
    recursive = T)
  
  s_ids <- gsub("\\.bed", "", basename(in_files))
  
  for ( i in 1:length(in_files) ) {
    if (i==1) {
      files_list <- list(in_files[i])
    } else {
      files_list[[i]] <- in_files[i]
    }
  }
  
  load_bed <- function(x) {

    library(rtracklayer)
    library(GenomicRanges)
    
    # read in bed as table:
    tab <- read.table(x, sep = "\t")
    
    # convert to GRanges object:
    result <- GRanges(
      seqnames = tab$V1,
      ranges = IRanges(start=tab$V2, end=tab$V3),
      strand = "*",
      FDR  = tab$V4
    ) 
    result <- result[grep("K|G|MT", seqnames(result), invert=T)]
    return(result)
    
  }

  # load in bams in parallel:
  peak_gr <- parLapply(cl, in_files, load_bed)
  names(peak_gr) <- s_ids

  saveRDS(peak_gr, paste0(RobjectDir, "/peak_gr.rds"))
  
} else {
  print("Loading reads GRanges object...")
  peak_gr <- readRDS(paste0(RobjectDir, "/peak_gr.rds"))
}


##########################################################################
### 2. Load in repeat and gencode annotations ###
##########################################################################

# create function to remove unwanted references and add to either side of 
# ranges of annotation:
exp_annot <- function(annot, Length, pos, is_ctl=F) {
  
  library("BSgenome.Hsapiens.UCSC.hg38")
  # define lengths of chromosomes:
  seq_lengths <- seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
  
  # remove unwanted chromosome constructs:
  annot <- annot[grep("[0-9].[0-9]|MT|G|K", seqnames(annot), invert = T)]
  
  # reduce ranges:
  annot <- reduce(annot)
  
  if ( is_ctl ) {
    ranges(annot) <- IRanges(start=min(start(annot)), end=max(end(annot)))
    annot <- unique(annot)
  }

  # assign length of all chromosomes to the variable name 'seq_lengths':
  annot$seq_lengths <- rep(NA, length(annot))
  
  for ( v in names(seq_lengths) ) {
    print(v)
    annot$seq_lengths[as.character(seqnames(annot)) == v] <- seq_lengths[v]
  }
  
  if ( pos == "body" ) {

    # add length to the start of each ranges if the start is that length or more from
    # the start of the chromosome:
    start(ranges(annot))[start(ranges(annot)) >= Length] <- 
      start(ranges(annot))[start(ranges(annot)) >= Length] - Length
    
    # add length to the end of each ranges if the end is length bp or more from
    # the end of the chromosome:
    end(ranges(annot))[end(ranges(annot)) <= (annot$seq_lengths - Length)] <- 
      end(ranges(annot))[end(ranges(annot)) <= (annot$seq_lengths - Length)] + Length

    return(annot)

  } else if ( pos == "body_no_promoter" ) {

    return(annot) 

  } else {
    
    # add length upstream of each range to another gr:
    start_annot <- annot
    # make end position equal to original start position:
    end(ranges(start_annot)) <- start(ranges(start_annot))
    # make start position the original start position - Length if the original start is at least Length
    # away from start of chromosome:
    start(start_annot)[start(ranges(start_annot)) >= Length] <- 
      start(ranges(start_annot))[start(ranges(start_annot)) >= Length] - Length
    
    if ( pos == "upstream" ) {

      return(start_annot)
      
    } else if ( pos == "up_and_downstream" ) {
      
      # add length downstream of each range to another gr:
      end_annot <- annot
      
      # make start position equal to original end position:
      start(ranges(end_annot)) <- end(ranges(end_annot))
      # make end position the original end position + Length if the end is at least Length
      # away from end of chromosome:
      far_enuff <- end(ranges(end_annot)) <= (end_annot$seq_lengths - Length)
      end(ranges(end_annot))[far_enuff] <- end(ranges(end_annot))[far_enuff] + Length
      
      # make end position the end of the chromosome if the original end is not at least Length
      # away from end of chromosome:
      not_far_enuff <- end(ranges(end_annot)) > (end_annot$seq_lengths - Length)
      end(ranges(end_annot))[not_far_enuff] <- seq_lengths[not_far_enuff]
      
      # return combined start and end annotations:
      return(c(start_annot, end_annot))
      
    }
  }
}

if ( !file.exists(paste0(newRobjectDir, "/rp_", incl_marks, ".rds")) ) {
  
  print("Creating expanded repeat annotation...")
  # load repeats annotation:
  rp_annot <- import(paste0(refDir, "custom3rep.hg19.gtf"))
  
  # split annot into GRangesLists by IDs/names:
  rp_annot <- split(rp_annot, rp_annot$ID)
  
  # expand ranges by 500 bp either side:
  n=1
  rp_annot <- lapply(rp_annot, exp_annot, Length = exp_no,
                       pos = incl_marks)
  
  saveRDS(rp_annot, paste0(newRobjectDir, "/rp_", incl_marks, ".rds"))
  
} else {
  print("Loading expanded repeat annotation...")
  rp_annot <- readRDS(paste0(newRobjectDir, "/rp_", incl_marks, ".rds"))
}

if ( !file.exists(paste0(newRobjectDir, "/CPMR_ctls_", incl_marks, ".rds")) ) {
  
  print("Creating expanded control annotation...")
  # load gencode annotation:
  gc <- import(paste0(refDir, "gencode_ercc.v19.annotation.gtf"))
  
  # isolate control ranges from gc:
  ctls <- gc[gc$gene_name %in% pos_ctl|gc$gene_name %in% neg_ctl]
  values(ctls) <- subset(values(ctls), select=gene_name)
  colnames(values(ctls)) <- "ID"
  
  # split annot into GRangesLists by IDs/names:
  ctls <- split(ctls, ctls$ID)
  
  # extract regions of interest from ranges:
  ctls <- lapply(ctls, exp_annot, Length = exp_no, 
                    pos = incl_marks, is_ctl = T)
  
  saveRDS(ctls, paste0(newRobjectDir, "/CPMR_ctls_", incl_marks, ".rds"))

} else {
  print("Loading expanded control annotation...")
  ctls <- readRDS(paste0(newRobjectDir, "/CPMR_ctls_", incl_marks, ".rds"))
}

save.image(file = paste0(newRobjectDir, "/annot_expanded.rds"))

# check rp_annot and ctls generated properly against IGV:
# if ( incl_marks == "upstream" ) {
#   AluJo_eg <- rp_annot$AluJo[seqnames(rp_annot$AluJo) == "chr8"]
#   AluJo_eg <- AluJo_eg[!is.na(end(ranges(AluJo_eg)))]
#   AluJo_eg <- AluJo_eg[start(ranges(AluJo_eg)) > 27559000]
#   AluJo_eg <- AluJo_eg[end(ranges(AluJo_eg)) < 27562000]
#   
#   ACTB_eg <- ctls$ACTB[seqnames(ctls$ACTB) == "chr7"]
#   ACTB_eg <-ACTB_eg[!is.na(end(ranges(ACTB_eg)))]
#   ACTB_eg <-ACTB_eg[start(ranges(ACTB_eg)) > 5564000]
#   ACTB_eg <-ACTB_eg[end(ranges(ACTB_eg)) < 5567000]
# }


##########################################################################
### 3. Quantitate histone mark enrichment in repeat and control 
# regions ###
##########################################################################

if ( !file.exists(paste0(newRobjectDir, "/repeats_peak_overlap_data.rds")) ) {
  
  count_repeat_peaks <- function(database, region, is_ctl = F){
    library(rtracklayer)
    library(GenomicRanges)
    library("BSgenome.Hsapiens.UCSC.hg38")
    
    r_region <- reduce(region)
    r_database <- reduce(database)
    
    # count number of ranges in each gene annotation to normalise results:
    seq_no <- length(r_database)
    
    # count region hits with database:
    hits <- findOverlaps(r_region, r_database)
    
    return(list(length(hits), region[queryHits(hits)], database[subjectHits(hits)], seq_no))
  }
  
  # apply odds ratio function to each element with parLapply):
  for ( j in 1:length(peak_gr) ) {
    if (j==1) {
  
      rp_peak_data <- list(lapply(rp_annot, count_repeat_peaks, peak_gr[[j]]))
      names(rp_peak_data)[j] <- names(peak_gr)[j]
      ctl_peak_data <- list(lapply(ctls, count_repeat_peaks, peak_gr[[j]], is_ctl = T))
      names(ctl_peak_data)[j] <- names(peak_gr)[j]
        
    } else {
        
      rp_peak_data[[j]] <- lapply(rp_annot, count_repeat_peaks, peak_gr[[j]])
      names(rp_peak_data)[j] <- names(peak_gr)[j]
      ctl_peak_data[[j]] <- lapply(ctls, count_repeat_peaks, peak_gr[[j]], is_ctl = T)
      names(ctl_peak_data)[j] <- names(peak_gr)[j]
      
    }
  }

  saveRDS(rp_peak_data, file = paste0(newRobjectDir, "/repeats_peak_overlap_data.rds"))
  saveRDS(ctl_peak_data, file = paste0(newRobjectDir, "/ctl_peak_overlap_data.rds"))
  
} else {
  rp_peak_data <- readRDS(file = paste0(newRobjectDir, "/repeats_peak_overlap_data.rds"))
  ctl_peak_data <- readRDS(file = paste0(newRobjectDir, "/ctl_peak_overlap_data.rds"))
}

# # check counts generated properly against IGV:
# if ( incl_marks == "upstream" ) {
#   ######
#   eg <- region[seqnames(region) == "chr8"]
#   eg <- eg[!is.na(end(ranges(eg)))]
#   eg <- eg[start(ranges(eg)) > 27559000]
#   eg <- eg[end(ranges(eg)) < 27565000]
# 
#   eg2 <- database[seqnames(database) == "chr8"]
#   eg2 <- eg2[!is.na(end(ranges(eg2)))]
#   eg2 <- eg2[start(ranges(eg2)) > 27559000]
#   eg2 <- eg2[end(ranges(eg2)) < 27565000]
#   ######
#   test <- rp_peak_data[[2]]$AluJo
#   AluJo_eg2 <- test[[2]][seqnames(test[[2]]) == "chr8"]
#   AluJo_eg2 <- AluJo_eg2[!is.na(end(ranges(AluJo_eg2)))]
#   AluJo_eg2 <- AluJo_eg2[start(ranges(AluJo_eg2)) > 27559000]
#   AluJo_eg2 <- AluJo_eg2[end(ranges(AluJo_eg2)) < 27565000]
# 
#   AluJo_eg3 <- test[[3]][seqnames(test[[3]]) == "chr8"]
#   AluJo_eg3 <- AluJo_eg3[!is.na(end(ranges(AluJo_eg3)))]
#   AluJo_eg3 <- AluJo_eg3[start(ranges(AluJo_eg3)) > 27559000]
#   AluJo_eg3 <- AluJo_eg3[end(ranges(AluJo_eg3)) < 27565000]
#   ######
# }


########################################################################
### 3. Plot counts of peak enrichment in repeat regions #
########################################################################

# create function for formatting and converting counts to percentage
# per sequence:
format_counts <- function(peak_data) {
  
  # fetch counts only, normalise to number of ranges per gene annotation
  # and collapse into data frame:
  Percent <- lapply(peak_data, function(x) {
    peak_percent <- lapply(x, function(y) {
      return( (y[[1]]/y[[4]])*100 )
    })
    vec <- unlist(peak_percent)
    result <- data.frame(names(vec), vec)
    colnames(result) <- c("id", "percent")
    return(result)
  })
  
  # collapse data frames for all peak types into one df:
  df <- do.call("cbind", Percent)
  
  # get rid of all but one id column:
  ind <- grep("id", colnames(df))
  ind <- ind[2:length(ind)]
  return(df[,-ind])

}

# format and melt data for barplot:
rp_df <- format_counts(rp_peak_data)
rp_df$group <- "repeat"
ctl_df <- format_counts(ctl_peak_data)
ctl_df$group <- "control"
plot_df <- rbind(ctl_df, rp_df)
ids <- gsub("\\_SRR.*$", "", s_ids)
colnames(plot_df) <- c("id", ids[1], ids[2], "group")

m_df <- melt(plot_df)

# fetch DE repeats from RNA-seq with FDR < 0.1 and 0.3
RNA_DE_0.1 <- read.table(file = paste0(DE_dir, "/sig_reps_FDR_0.1.txt"))
RNA_DE_0.3 <- read.table(file = paste0(DE_dir, "/sig_reps_FDR_0.3.txt"))

DE_list <- list(rownames(RNA_DE_0.1)[RNA_DE_0.1$logFC > 0])
DE_list[[2]] <- rownames(RNA_DE_0.3)[RNA_DE_0.3$logFC > 0]
DE_list[[3]] <- rownames(RNA_DE_0.1)[RNA_DE_0.1$logFC < 0]
DE_list[[4]] <- rownames(RNA_DE_0.3)[RNA_DE_0.3$logFC < 0]
names(DE_list) <- c("up_0.1", "up_0.3", "down_0.1", "down_0.3")

# adjust groups column to add control and repeat type information:
m_df$group[m_df$id %in% pos_ctl] <- "positive_control"
m_df$group[m_df$id %in% neg_ctl] <- "negative_control"
m_df$group[m_df$id %in% DE_list$up_0.3] <- "up_repeat_FDR<0.3"
m_df$group[m_df$id %in% DE_list$up_0.1] <- "up_repeat_FDR<0.1"
m_df$group[m_df$id %in% DE_list$down_0.3] <- "down_repeat_FDR<0.3"
m_df$group[m_df$id %in% DE_list$down_0.1] <- "down_repeat_FDR<0.1"
m_df$group[m_df$group == "repeat"] <- "non_DE_repeat"
colnames(m_df) <- c("gene", "group", "peak", "percent")

# plot on barplot:
p <- ggplot(m_df, aes(x=group, y=percent, fill=peak))
p <- p + geom_bar(stat = "identity", position=position_dodge(0.9))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
p <- p + ylab("Percentage of sequence copy number in genome")
p

