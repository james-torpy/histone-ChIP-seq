### repeat_histone_mark_odds_100_chapmanDE_ctls.R ###

# This script takes a bam file from a histone mark ChIP-seq alignment 
# and calculates the log odds of H3K4me3 enrichment at promoters of
# control genes


##########################################################################
### 0. Define variables/paths ###
##########################################################################

# load packages needed:
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library("BSgenome.Hsapiens.UCSC.hg38")

# define starting variables:
project <- "hgsoc_repeats"
exp_name <- "histone-ChIP-seq"
sample_name <- "chapman-rothe_2013_hg19"
descrip <- paste0("histone_mark_odds_100_chapmanDE_ctls")

bp_range <- 500
posit <- "up_and_downstream"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "/projects/", project, "/", 
	exp_name, "/", sample_name, "/")
results_dir <- paste0(project_dir, "/results/")
ref_dir <- paste0(project_dir, "/refs/")
Robject_dir <- paste0(project_dir, "/Robjects/")
plot_dir <- paste0(results_dir, "/R/plots/")
table_dir <- paste0(results_dir, "/R/tables/")

in_dir <- paste0(results_dir, "/epic/")
DE_dir <- paste0(home_dir,
	"/projects/hgsoc_repeats/RNA-seq/results/R/exp9/plots/DEplots/", 
	"htseq_EdgeR_primary_HGSOC_vs_FT/")

system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", table_dir))
system(paste0("mkdir -p ", Robject_dir))

pos_chapman_ctl <- as.character(read.table(file=paste0(ref_dir, 
	"/chapman-roethe_top_DE_RNA_symbols_ids.txt"))[,1])
neg_chapman_ctl <- as.character(read.table(file=paste0(ref_dir, 
    "/chapman-roethe_bottom_DE_RNA_symbols_ids.txt"))[,1])


##########################################################################
### 1. Load in ChIP reads files ###
##########################################################################

in_files <- grep(
	"H3K4me3", list.files(
		in_dir, pattern = ".bed", full.names = T, 
    	recursive = T
    ), value = T
)

s_ids <- gsub("\\.bed", "", basename(in_files))


##########################################################################
### 2. Create peak_gr ###
##########################################################################

load_bed <- function(x) {
        
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

if ( !file.exists(paste0(Robject_dir, "/", descrip, 
    	"_peak_gr.rds")) ) {

	for ( i in 1:length(in_files) ) {
        if (i==1) {
          files_list <- list(in_files[i])
        } else {
          files_list[[i]] <- in_files[i]
        }
      }

	# load in bams in parallel:
    peak_gr <- lapply(in_files, load_bed)
    names(peak_gr) <- s_ids
  
    saveRDS(peak_gr, paste0(Robject_dir, "/", descrip, 
    	"_peak_gr.rds"))
  
} else {

    print("Loading reads GRanges object...")
    peak_gr <- readRDS(paste0(Robject_dir, "/", descrip, 
    	"_peak_gr.rds"))

}


##########################################################################
### 3. Create control annotation and expand ###
##########################################################################
    
# create function to remove unwanted references and add to either 
# side of ranges of annotation:
exp_annot <- function(annot, Length) {

	print(n)

	# define lengths of chromosomes:
	seq_lengths <- seqlengths(Hsapiens)[!grepl("_",
		names(seqlengths(Hsapiens)))]
	# remove unwanted chromosome constructs:
	annot <- annot[grep("[0-9].[0-9]|MT|G|K", seqnames(annot), 
		invert = T)]
	  
    if ( length(annot) > 0 ) {
      # reduce ranges:
      annot <- reduce(annot)
    
    ranges(annot) <- IRanges(start=min(start(annot)), 
    	end=max(end(annot)))
    annot <- unique(annot)
        
    # assign length of all chromosomes to the variable name 'seq_lengths':
    annot$seq_lengths <- rep(NA, length(annot))
    
    for ( v in names(seq_lengths) ) {
      annot$seq_lengths[as.character(seqnames(annot)) == v] <- seq_lengths[v]
    }
      
    # add length upstream of each range to another gr:
    start_annot <- annot
    # make end position equal to original start position:
    end(ranges(start_annot)) <- start(ranges(start_annot))
    # make start position the original start position - Length if the original start is at least Length
    # away from start of chromosome:
    start(start_annot)[start(ranges(start_annot)) >= Length] <- 
      start(ranges(start_annot))[start(ranges(start_annot)) >= Length] - Length

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

    if ( not_far_enuff ) {
    	#n <<- n+1
    	return(NULL)
    } else {
    	# return combined start and end annotations:
    	#n <<- n+1
    	return(c(start_annot, end_annot))
    }
  }
}

if ( !file.exists(paste0(Robject_dir, "/ctl_annot_", bp_range, "_bp_", 
	posit, "_", descrip, ".rds")) ) {
      
    print("Creating expanded control annotation...")
    # load gencode annotation:
    gc <- import(paste0(ref_dir, "/gencode_ercc.v19.annotation.gtf"))
    
    # isolate control ranges from gc:
    ctls <- gc[gc$gene_name %in% pos_chapman_ctl|gc$gene_name %in% 
    	neg_chapman_ctl]
    values(ctls) <- subset(values(ctls), select=gene_name)
    colnames(values(ctls)) <- "ID"
    
    # split annot into GRangesLists by IDs/names:
    ctls <- split(ctls, ctls$ID)
    
    # extract regions of interest from ranges:
    ctls <- lapply(ctls, exp_annot, Length = bp_range)
    
    # remove NULL values:
    if ( any(unlist(lapply(ctls, is.null))) ) {
      ctls <- ctls[-which(unlist(lapply(ctls, is.null)))]
    }
    
    saveRDS(ctls, paste0(RobjectDir, "/ctls_", Posit,
                           "_", exp_nos[o], "_bp_", descrip, ".rds"))
      
} else {

    print("Loading expanded control annotation...")
    ctls <- readRDS(paste0(Robject_dir, "/ctl_annot_", bp_range, "_bp_", 
		posit, "_", descrip, ".rds"))

}







