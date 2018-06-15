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
library(pheatmap)
library(parallel)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "histone-ChIP-seq"
sampleName <- "chapman-rothe_2013_hg19"
descrip <- "absolute_values_chapmanDE_pos_ctl"
exp_no <- 500

Posit <- "upstream"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/", expName, "/", sampleName, "/")
resultsDir <- paste0(projectDir, "/results")
refDir <- paste0(projectDir, "/refs/")
RobjectDir <- paste0(projectDir, "/Robjects/")
plotDir <- paste0(resultsDir, "/R/plots/")

inDir <- paste0(resultsDir, "/epic/")
ref_dir <- paste0(projectDir, "/refs/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", tableDir))
system(paste0("mkdir -p ", RobjectDir))

pos_chapman_ctl <- as.character(read.table(file=paste0(ref_dir, 
  "/chapman-roethe_top_DE_RNA_symbols_ids.txt"))[,1])

##########################################################################
### 1. Load in ChIP reads files ###
##########################################################################

in_files <- list.files(inDir, pattern = ".bed", full.names = T, 
                       recursive = T)


s_ids <- gsub("\\.bed", "", basename(in_files))


if ( !file.exists(paste0(RobjectDir, "peak_gr.rds")) ) {
  
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
  peak_gr <- lapply(in_files, load_bed)
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
exp_annot <- function(annot, Length, Posit, is_ctl=F) {
  
  library("BSgenome.Hsapiens.UCSC.hg38")
  # define lengths of chromosomes:
  seq_lengths <- seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
  
  # remove unwanted chromosome constructs:
  annot <- annot[grep("[0-9].[0-9]|MT|G|K", seqnames(annot), invert = T)]
  
  if ( length(annot) > 0 ) {
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
    
    if ( Posit == "body" ) {
      
      # add length to the start of each ranges if the start is that length or more from
      # the start of the chromosome:
      start(ranges(annot))[start(ranges(annot)) >= Length] <- 
        start(ranges(annot))[start(ranges(annot)) >= Length] - Length
      
      # add length to the end of each ranges if the end is length bp or more from
      # the end of the chromosome:
      end(ranges(annot))[end(ranges(annot)) <= (annot$seq_lengths - Length)] <- 
        end(ranges(annot))[end(ranges(annot)) <= (annot$seq_lengths - Length)] + Length
      
      return(annot)
      
    } else if ( Posit == "body_no_promoter" ) {
      
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
      
      if ( Posit == "upstream" ) {
        
        return(start_annot)
        
      } else if ( Posit == "up_and_downstream" ) {
        
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
}

if ( !file.exists(paste0(RobjectDir, "/ctls_", Posit, "_", 
                           exp_no, "_bp_", descrip, ".rds")) ) {
      
      print("Creating expanded control annotation...")
      # load gencode annotation:
      gc <- import(paste0(refDir, "gencode_ercc.v19.annotation.gtf"))
      
      # isolate control ranges from gc:
      ctls <- gc[gc$gene_name %in% pos_chapman_ctl]
      values(ctls) <- subset(values(ctls), select=gene_name)
      colnames(values(ctls)) <- "ID"
      
      # split annot into GRangesLists by IDs/names:
      ctls <- split(ctls, ctls$ID)
      
      # extract regions of interest from ranges:
      ctls <- lapply(ctls, exp_annot, Length = exp_no, 
                        Posit = Posit, is_ctl = T)
      
      # remove NULL values:
      if ( any(unlist(lapply(ctls, is.null))) ) {
        ctls <- ctls[-which(unlist(lapply(ctls, is.null)))]
      }
      
      saveRDS(ctls, paste0(RobjectDir, "/ctls_", Posit, "_", 
                           Subset, "_", exp_nos[o], "_bp_", descrip, ".rds"))
      
    } else {
      print("Loading expanded control annotation...")
      ctls <- readRDS(paste0(RobjectDir, "/ctls_", Posit, "_", 
                             Subset, "_", exp_nos[o], "_bp_", descrip, ".rds"))
    }