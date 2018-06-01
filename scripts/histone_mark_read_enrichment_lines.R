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
library(scales)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "histone-ChIP-seq"
sampleName <- "gayther_2015_kalender_2017"
descrip <- paste0(sampleName, "_enrichment")

p_thresh <- 0.1

# specify regions to include for marks (body, upstream, up_and_downstream)
incl_marks <- "upstream"
# specify how many bp up/downstream to expand repeats annotation by:
exp_no <- 500

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/", expName, "/", sampleName, "/")
resultsDir <- paste0(projectDir, "/results")
refDir <- paste0(projectDir, "/refs/")
RobjectDir <- paste0(projectDir, "/Robjects/", sampleName, "/")
plotDir <- paste0(resultsDir, "/R/plots/")
tableDir <- paste0(resultsDir, "/R/tables/")

inDir <- paste0(resultsDir, "/bwa/")
DE_dir <- paste0(homeDir,
                 "projects/hgsoc_repeats/RNA-seq/results/R/exp9/plots/DEplots/htseq_EdgeR_primary_HGSOC_vs_FT/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", tableDir))
system(paste0("mkdir -p ", RobjectDir))

pos_ctl <- read.table(file=paste0(DE_dir, "/top_50_up_genes_allsymbols_DE.txt"))$symbol
neg_ctl <- read.table(file=paste0(DE_dir, "/top_50_down_genes_allsymbols_DE.txt"))$symbol
                      

##########################################################################
### 1. Load in ChIP reads files ###
##########################################################################

if ( !file.exists(paste0(RobjectDir, "read_gr.rds")) ) {
  
  in_files <- grep(
    "input", list.files(
      inDir, pattern = "sorted.bam$", full.names = T
    ), invert=T, value = T
  )
  
  s_ids <- gsub("\\.sorted.bam", "", basename(in_files))
  
  # specify scambam parameters:
  what <- c("qname", "rname", "strand", "pos", "qwidth")
  param <- ScanBamParam(what=what)
  
  # load in bams:
  for ( i in 1:length(in_files) ) {
    if (i==1) {
      bams <- list(scanBam(in_files[i], param=param))
    } else {
      bams[[i]] <- scanBam(in_files[i], param=param)
    }
  }
  
  # convert the bams to GRanges object:
  read_gr <- lapply(bams, function(x) {
    x <- x[[1]]
    x$qname <- as.character(x$qname)
    result <- GRanges(
          seqnames = x$rname,
          ranges = IRanges(start=x$pos, width=x$qwidth),
          strand = x$strand,
          read_id  = x$qname
      ) 
    result <- result[grep("K|G|MT", seqnames(result), invert=T)]
    return(result)
  })

  names(read_gr) <- s_ids

  saveRDS(read_gr, paste0(RobjectDir, "/read_gr.rds"))
  
} else {
  print("Loading reads GRanges object...")
  read_gr <- readRDS(paste0(RobjectDir, "/read_gr.rds"))
}


##########################################################################
### 2. Load in repeat and gencode annotations ###
##########################################################################

# define lengths of chromosomes:
seq_lengths <- seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]

# create function to remove unwanted references and add to either side of 
# ranges of annotation:
exp_annot <- function(annot, Length, Position) {
  # remove unwanted chromosome constructs:
  annot <- annot[grep("[0-9].[0-9]|M", seqnames(annot), 
                      invert = T)]
  
  # assign length of all chromosomes to the variable name 'seq_lengths':
  annot$seq_lengths <- rep(NA, length(annot))
  
  for ( v in names(seq_lengths) ) {
    print(v)
    annot$seq_lengths[as.character(seqnames(annot)) == v] <- seq_lengths[v]
  }
  
  
  if ( Position == "body" ) {
    # add length to the start of each ranges if the start is that length or more from
    # the start of the chromosome:
    start(ranges(annot))[start(ranges(annot)) >= Length] <- 
      start(ranges(annot))[start(ranges(annot)) >= Length] - Length
    
    # add length to the end of each ranges if the end is length bp or more from
    # the end of the chromosome:
    end(ranges(annot))[end(ranges(annot)) <= (annot$seq_lengths - Length)] <- 
      end(ranges(annot))[end(ranges(annot)) <= (annot$seq_lengths - Length)] + Length
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
    
    if ( Position == "upstream" ) {
      
      return(start_annot)
      
    } else if ( Position == "up_and_downstream" ) {
      
      # add length downstream of each range to another gr:
      end_annot <- annot
      
      # make start position equal to original end position:
      start(ranges(end_annot)) <- end(ranges(end_annot))
      # make end position the original end position + Length if the send is at least Length
      # away from end of chromosome:
      end(ranges(end_annot))[end(ranges(end_annot)) <= (end_annot$seq_lengths - Length)] <- 
        end(ranges(end_annot))[end(ranges(end_annot)) <= (end_annot$seq_lengths - Length)] + Length
      
      # return combined start and end annotations:
      return(c(start_annot, end_annot))
      
    }
  }
}

if ( !file.exists(paste0(RobjectDir, "/rp_", incl_marks, ".rds")) ) {
  
  print("Creating expanded repeat annotation...")
  # load repeats annotation:
  rp_annot <- import(paste0(refDir, "custom3rep.final.gff"))
  
  # expand ranges by 500 bp either side:
  rp_annot <- exp_annot(rp_annot, exp_no, incl_marks)
  
  # split annot into GRangesLists by IDs/names:
  rp_annot <- split(rp_annot, rp_annot$ID)
  
  saveRDS(rp_annot, paste0(RobjectDir, "/rp_", incl_marks, ".rds"))
} else {
  print("Loading expanded repeat annotation...")
  rp_annot <- readRDS(paste0(RobjectDir, "/rp_", incl_marks, ".rds"))
}

if ( !file.exists(paste0(RobjectDir, "/ctls_", incl_marks, ".rds")) ) {
  
  print("Creating expanded control annotation...")
  # load gencode annotation:
  gc <- import(paste0(refDir, "gencode_v24_hg38_annotation.gtf"))
  
  # expand ranges by 500 bp either side:
  gc <- exp_annot(gc, exp_no, incl_marks)
  saveRDS(gc, paste0(RobjectDir, "/gc_", incl_marks, ".rds"))
  
  # isolate control ranges from gc:
  ctls <- gc[gc$gene_name %in% pos_ctl|gc$gene_name %in% neg_ctl]
  values(ctls) <- subset(values(ctls), select=gene_name)
  colnames(values(ctls)) <- "ID"
  
  # split annot into GRangesLists by IDs/names:
  ctls <- split(ctls, ctls$ID)
  
  saveRDS(ctls, paste0(RobjectDir, "/ctls_", incl_marks, ".rds"))
  
} else {
  print("Loading expanded control annotation...")
  ctls <- readRDS(paste0(RobjectDir, "/ctls_", incl_marks, ".rds"))
}

save.image(file = paste0(RobjectDir, "/annot_expanded.rds"))


##########################################################################
### 3. Calculate log odds ratios of H3K27me3 enrichment in repeat and
# control regions ###
##########################################################################

if ( !file.exists(paste0(RobjectDir, "/repeats_read_overlap_odds.rds")) ) {
  
  # fetch summed length of all chromosomes
  chrs <- seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
  chr_gr <- GRanges(seqnames=names(chrs),IRanges(1,chrs))
  
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
      matrix( c(as.numeric(region_hits), as.numeric(non_region_hits), 
      as.numeric(region_length), as.numeric(non_region)), nrow = 2,
      ncol = 2 )
    colnames(database_enrichment) = c("hits", "non-hits")
    rownames(database_enrichment) = c("region", "non-region")

    pval <- chisq.test(database_enrichment)$p.value
    print("P-value is:")
    cat(pval)
    cat("\n")
    
    result <- data.frame(odds, se, pval)
    colnames(result) <- c("odds", "std_error", "pval")
    
    return(result)
  }
  
  # apply odds ratio function to each element:
  for ( j in 1:length(read_gr) ) {
    if (j==1) {
      
      rp_read_odds <- list(lapply(rp_annot, oddsRatio, read_gr[[j]]))
      names(rp_read_odds)[j] <- names(read_gr)[j]
      ctl_read_odds <- list(lapply(ctls, oddsRatio, read_gr[[j]]))
      names(ctl_read_odds)[j] <- names(read_gr)[j]
      
    } else {
      
      rp_read_odds[[j]] <- lapply(rp_annot, oddsRatio, read_gr[[j]])
      names(rp_read_odds)[j] <- names(read_gr)[j]
      ctl_read_odds[[j]] <- lapply(ctls, oddsRatio, read_gr[[j]])
      names(ctl_read_odds)[j] <- names(read_gr)[j]
      
    }
  }
  
  # bind odds results for repeats and controls into one data frame per 
  # sample/annotation combination:
  rp_test <- lapply(rp_read_odds, function(x) {
    return(do.call("rbind", x))
  })
  
  saveRDS(rp_read_odds, file = paste0(RobjectDir, "/repeats_read_overlap_odds.rds"))
  saveRDS(ctl_read_odds, file = paste0(RobjectDir, "/ctl_read_overlap_odds.rds"))
  
} else {
  rp_read_odds <- readRDS(file = paste0(RobjectDir, "/repeats_read_overlap_odds.rds"))
  ctl_read_odds <- readRDS(file = paste0(RobjectDir, "/ctl_read_overlap_odds.rds"))
}

