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
descrip <- paste0("absolute_values_chapmanDE_ctls")
exp_nos <- c(500)
posits <- "up_and_downstream"
Subset=""

# specify regions to include for marks (body, upstream, up_and_downstream)
# Posit <- "upstream"
# # specify how many bp up/downstream to expand repeats annotation by:
# exp_no <- 2000
# # specify whether subsets of annotations will be used (for debugging) - 
# # either "body_no_promoter", "upstream" or "":

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/", expName, "/", sampleName, "/")
resultsDir <- paste0(projectDir, "/results")
refDir <- paste0(projectDir, "/refs/")
RobjectDir <- paste0(projectDir, "/Robjects/")
plotDir <- paste0(resultsDir, "/R/plots/")
tableDir <- paste0(resultsDir, "/R/tables/")

inDir <- paste0(resultsDir, "/epic/")
ref_dir <- paste0(projectDir, "/refs/")
DE_dir <- paste0(homeDir,
                 "projects/hgsoc_repeats/RNA-seq/results/R/exp9/plots/DEplots/htseq_EdgeR_primary_HGSOC_vs_FT/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", tableDir))
system(paste0("mkdir -p ", RobjectDir))

pos_chapman_ctl <- as.character(read.table(file=paste0(ref_dir, 
                                                       "/chapman-roethe_top_DE_RNA_symbols_ids.txt"))[,1])
neg_chapman_ctl <- as.character(read.table(file=paste0(ref_dir, 
                                                       "/chapman-roethe_non-DE_RNA_symbols_ids.txt"))[,1])
neg_chapman_ctl2 <- as.character(read.table(file=paste0(ref_dir, 
                                                        "/chapman-roethe_bottom_DE_RNA_symbols_ids.txt"))[,1])


# set up parallel workers:
no_cores <- 3
print(paste0("No. cores are: ", no_cores))
cl <- makeCluster(no_cores)             

for (o in 1:length(exp_nos)) {
  
  clusterExport(cl, "o")
  res <- lapply(posits, function(Posit) {  
    
    exp_nos <- c(500)
    
    
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
    
    # define starting variables:
    project <- "hgsoc_repeats"
    expName <- "histone-ChIP-seq"
    sampleName <- "chapman-rothe_2013_hg19"
    descrip <- paste0(sampleName, "SICER_peak_enrichment")
    Subset=""
    
    
    # specify regions to include for marks (body, upstream, up_and_downstream)
    # Posit <- "upstream"
    # # specify how many bp up/downstream to expand repeats annotation by:
    # exp_no <- 2000
    # # specify whether subsets of annotations will be used (for debugging) - 
    # # either "body_no_promoter", "upstream" or "":
    
    # define directories:
    #homeDir <- "/Users/jamestorpy/clusterHome/"
    homeDir <- "/share/ScratchGeneral/jamtor/"
    projectDir <- paste0(homeDir, "/projects/", project, "/", expName, "/", sampleName, "/")
    resultsDir <- paste0(projectDir, "/results")
    refDir <- paste0(projectDir, "/refs/")
    RobjectDir <- paste0(projectDir, "/Robjects/")
    plotDir <- paste0(resultsDir, "/R/plots/")
    tableDir <- paste0(resultsDir, "/R/tables/")
    
    inDir <- paste0(resultsDir, "/epic/")
    ref_dir <- paste0(projectDir, "/refs/")
    DE_dir <- paste0(homeDir,
                     "projects/hgsoc_repeats/RNA-seq/results/R/exp9/plots/DEplots/htseq_EdgeR_primary_HGSOC_vs_FT/")
    
    
    system(paste0("mkdir -p ", plotDir))
    system(paste0("mkdir -p ", tableDir))
    system(paste0("mkdir -p ", RobjectDir))
    system(paste0("mkdir -p ", RobjectDir))
    
    pos_chapman_ctl <- as.character(read.table(file=paste0(ref_dir, 
                                                           "/chapman-roethe_top_DE_RNA_symbols_ids.txt"))[,1])
    neg_chapman_ctl <- as.character(read.table(file=paste0(ref_dir, 
                                                           "/chapman-roethe_non-DE_RNA_symbols_ids.txt"))[,1])
    neg_chapman_ctl2 <- as.character(read.table(file=paste0(ref_dir, 
                                                            "/chapman-roethe_bottom_DE_RNA_symbols_ids.txt"))[,1])
    
    
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
    
    if ( !file.exists(paste0(RobjectDir, "/rp_", Posit, "_", Subset, 
                             "_", as.character(exp_nos[o]), "bp.rds")) ) {
      
      print("Creating expanded repeat annotation...")
      # load repeats annotation:
      rp_annot <- import(paste0(refDir, "custom3rep.hg19.gtf"))
      
      # split annot into GRangesLists by IDs/names:
      rp_annot <- split(rp_annot, rp_annot$ID)
      
      
      
      if ( Subset == "body_no_promoter_subset" ) {
        # subset annotation:
        L1MA9_eg <- rp_annot$L1MA9[start(rp_annot$L1MA9) == 12420804]
        L1M5_eg <- rp_annot$L1M5[start(rp_annot$L1M5) == 197813669]
        L1MC2_eg <- rp_annot$L1MC2[start(rp_annot$L1MC2) == 12416123]
        LSAU_eg <- rp_annot$LSAU[start(rp_annot$LSAU) == 75719517]
        rp_annot <- GRangesList(L1MA9_eg, L1M5_eg, L1MC2_eg, LSAU_eg)
        names(rp_annot) <- c("L1MA9", "L1M5", "L1MC2", "LSAU")
        
      } else if ( Subset == "upstream_subset" ) {
        
        # subset annotation:
        Kanga1_eg <- rp_annot$Kanga1[start(rp_annot$Kanga1) == 54245008]
        AluJr4_eg <- rp_annot$AluJr4[start(rp_annot$AluJr4) == 53778053]
        L1M1_eg <- rp_annot$L1M1[start(rp_annot$L1M1) == 53687498]
        L1MCa_eg <- rp_annot$L1MCa[start(rp_annot$L1MCa) == 56531469]
        rp_annot <- GRangesList(Kanga1_eg, AluJr4_eg, L1M1_eg, L1MCa_eg)
        names(rp_annot) <- c("Kanga1", "AluJr4", "L1M1", "L1MCa")
      }
      
      
      # select subsection of ranges:
      rp_annot <- lapply(rp_annot, exp_annot, Length = exp_nos[o],
                         Posit = Posit)
      
      
      saveRDS(rp_annot, paste0(RobjectDir, "/rp_", Posit, "_", Subset, 
                               "_", exp_nos[o], "bp.rds"))
      
    } else {
      print("Loading expanded repeat annotation...")
      rp_annot <- readRDS(paste0(RobjectDir, "/rp_", Posit, "_", Subset, 
                                 "_", exp_nos[o], "bp.rds"))
    }
    
    if ( !file.exists(paste0(RobjectDir, "/ctl_rp_", Posit, "_", Subset, 
                             "_", as.character(exp_nos[o]), "bp.rds")) ) {
      
      print("Creating expanded control repeat annotation...")
      # load repeats annotation:
      ctl_reps <- import(paste0(refDir, "/ctl_reps.gtf"))
      
      # split annot into GRangesLists by IDs/names:
      ctl_reps <- split(ctl_reps, ctl_reps$ID)
      
      # select subsection of ranges:
      ctl_reps <- lapply(ctl_reps, exp_annot, Length = exp_nos[o],
                         Posit = Posit)
      
      
      saveRDS(ctl_reps, paste0(RobjectDir, "/ctl_rp_", Posit, "_", Subset, 
                               "_", as.character(exp_nos[o]), "bp.rds"))
      
    } else {
      
      print("Loading expanded control repeat annotation...")
      ctl_reps <- readRDS(paste0(RobjectDir, "/ctl_rp_", Posit, "_", Subset, 
                                 "_", as.character(exp_nos[o]), "bp.rds"))
      
    }
    
    if ( !file.exists(paste0(RobjectDir, "/ctls_", Posit, "_", 
                             Subset, "_", exp_nos[o], "_bp_", descrip, ".rds")) ) {
      
      print("Creating expanded control annotation...")
      # load gencode annotation:
      gc <- import(paste0(refDir, "gencode_ercc.v19.annotation.gtf"))
      
      # isolate control ranges from gc:
      ctls <- gc[gc$gene_name %in% pos_chapman_ctl|gc$gene_name %in% neg_chapman_ctl|
                   gc$gene_name %in% neg_chapman_ctl2]
      values(ctls) <- subset(values(ctls), select=gene_name)
      colnames(values(ctls)) <- "ID"
      
      # split annot into GRangesLists by IDs/names:
      ctls <- split(ctls, ctls$ID)
      
      # extract regions of interest from ranges:
      ctls <- lapply(ctls, exp_annot, Length = exp_nos[o], 
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
    
    save.image(file = paste0(RobjectDir, "/annot_expanded_", Subset, 
                             "_", exp_nos[o], "_bp_", descrip, ".rds"))
    
    
    ##########################################################################
    ### 3. Quantitate histone mark enrichment in repeat and control 
    # regions ###
    ##########################################################################
    
    if ( 
      # !file.exists(paste0(RobjectDir,
      #                     "/repeats_peak_overlap_data_", Subset, "_", Posit, "_", exp_nos[o], "bp", 
      #                     descrip, ".rds")) | 
      # !file.exists(paste0(RobjectDir,
      #                     "/ctl_repeats_peak_overlap_data_", Subset, "_", Posit, "_", exp_nos[o], 
      #                     "bp", descrip, ".rds")) | 
      !file.exists(paste0(RobjectDir, 
                          "/ctl_peak_overlap_data_", Subset, "_", Posit, "_", exp_nos[o], "bp", 
                          descrip, ".rds")) ) {
      
      count_repeat_peaks <- function(database, region, is_ctl = F){
        # library(rtracklayer)
        # library(GenomicRanges)
        # library("BSgenome.Hsapiens.UCSC.hg38")
        print(n)
        
        r_region <- reduce(region)
        r_database <- reduce(database)
        
        # count number of ranges in each gene annotation to normalise 
        # results:
        seq_no <- length(r_database)
        
        # count region hits with database:
        hits <- findOverlaps(r_region, r_database)
        
        # create list of results with hit count, hits in region, hits in
        # annotation, total number of ranges in annotation:
        results <- list(length(hits), r_region[queryHits(hits)], 
                        r_database[subjectHits(hits)], seq_no)
        names(results) <- c("hit_count", "region_hits", "database_hits",
                            "total_annot_ranges")
        n <<- n+1
        return(results)
      }
      
      # apply odds ratio function to each element with parLapply):
      for ( j in 1:length(peak_gr) ) {
        if (j==1) {
          
          # print(j)
          # n=1
          # rp_peak_data <- list(lapply(rp_annot, count_repeat_peaks,
          #                                peak_gr[[j]]))
          # names(rp_peak_data)[j] <- names(peak_gr)[j]
          
          # n=1
          # ctl_rp_peak_data <- list(lapply(ctl_reps, count_repeat_peaks,
          #                             peak_gr[[j]]))
          # names(ctl_rp_peak_data)[j] <- names(peak_gr)[j]
          
          n=1
          ctl_peak_data <- list(lapply(ctls, count_repeat_peaks, peak_gr[[j]], 
                                       is_ctl = T))
          names(ctl_peak_data)[j] <- names(peak_gr)[j]
          
        } else {
          
          # print(j)  
          # rp_peak_data[[j]] <- lapply(rp_annot, count_repeat_peaks, peak_gr[[j]])
          # names(rp_peak_data)[j] <- names(peak_gr)[j]
          # 
          # ctl_rp_peak_data[[j]] <- lapply(ctl_reps, count_repeat_peaks, peak_gr[[j]])
          # names(ctl_rp_peak_data)[j] <- names(peak_gr)[j]
          
          ctl_peak_data[[j]] <- lapply(ctls, count_repeat_peaks, peak_gr[[j]], 
                                       is_ctl = T)
          names(ctl_peak_data)[j] <- names(peak_gr)[j]
          
        }
      }
      # print("Saving control and repeat peak data")
      # saveRDS(rp_peak_data, file = paste0(RobjectDir,
      #   "/repeats_peak_overlap_data_", Subset, "_", Posit, "_", exp_nos[o], "bp", descrip, ".rds"))
      # 
      # saveRDS(ctl_rp_peak_data, file = paste0(RobjectDir,
      #    "/ctl_repeats_peak_overlap_data_", Subset, "_", Posit, "_", exp_nos[o], "bp", descrip, ".rds"))
      # 
      saveRDS(ctl_peak_data, file = paste0(RobjectDir, 
                                           "/ctl_peak_overlap_data_", Subset, "_", Posit, "_", exp_nos[o], "bp", descrip, ".rds"))
      
    } else {
      
      # rp_peak_data <- readRDS(file = paste0(RobjectDir,
      #   "/repeats_peak_overlap_data_", Subset, "_", Posit, "_", exp_nos[o], "bp", descrip, ".rds"))
      # 
      # ctl_rp_peak_data <- readRDS(file = paste0(RobjectDir,
      #    "/ctl_repeats_peak_overlap_data_", Subset, "_", Posit, "_", exp_nos[o], "bp", descrip, ".rds"))
      # 
      ctl_peak_data <- readRDS(file = paste0(RobjectDir, 
                                             "/ctl_peak_overlap_data_", Subset, "_", Posit, "_", exp_nos[o], "bp", descrip, ".rds"))
      
    }
    
    
    ########################################################################
    ### 3. Plot counts of peak enrichment in repeat regions #
    ########################################################################
    
    # create function for formatting and converting counts to percentage
    # per sequence:
    format_counts <- function(peak_data) {
      
      # fetch counts only and collapse into data frame:
      Count <- lapply(peak_data, function(x) {
        peak_count <- lapply(x, function(y) {
          return( y[[1]] )
        })
        vec <- unlist(peak_count)
        result <- data.frame(names(vec), vec)
        colnames(result) <- c("id", "count")
        return(result)
      })
      
      # collapse data frames for all peak types into one df:
      df <- do.call("cbind", Count)
      
      # get rid of all but one id column:
      ind <- grep("id", colnames(df), value = T)
      return( df[!(colnames(df) %in% ind)] )
      
    }
    
    # format and melt data for barplot:
    rp_df <- format_counts(rp_peak_data)
    rp_df$group <- "repeat"
    ctl_rp_df <- format_counts(rp_peak_data)
    ctl_rp_df$group <- "ctl_repeat"
    ctl_df <- format_counts(ctl_peak_data)
    ctl_df$group <- "control"
    # only one mark should be recorded per gene if Posit = 
    # upstream/up_and_downstream:
    if ( Posit == "upstream" | Posit == "up_and_downstream" ) {
      ctl_df$HGSOC_H3K4me3_SRR600956.count[ctl_df$HGSOC_H3K4me3_SRR600956.count > 1] <- 1
      ctl_df$HGSOC_H3K27me3_SRR600956.count[ctl_df$HGSOC_H3K27me3_SRR600956.count > 1] <- 1
    }
    
    plot_df <- rbind(ctl_df, ctl_rp_df)
    plot_df <- rbind(plot_df, rp_df)
    in_files <- list.files(inDir, pattern = ".bed", full.names = T, 
                           recursive = T)
    s_ids <- gsub("\\.bed", "", basename(in_files))
    ids <- gsub("\\_SRR.*$", "", s_ids)
    plot_df$id <- rownames(plot_df)
    colnames(plot_df) <- c(ids[1], ids[2], "group", "id")
    
    # fetch DE repeats from RNA-seq with FDR < 0.05 and 0.1
    RNA_DE_0.1 <- read.table(file = paste0(DE_dir, "/sig_reps_FDR_0.1.txt"))
    RNA_DE_0.05 <- read.table(file = paste0(DE_dir, "/sig_reps_FDR_0.05.txt"))
    
    DE_list <- list(rownames(RNA_DE_0.1)[RNA_DE_0.1$logFC > 0])
    DE_list[[2]] <- rownames(RNA_DE_0.05)[RNA_DE_0.05$logFC > 0]
    DE_list[[3]] <- rownames(RNA_DE_0.1)[RNA_DE_0.1$logFC < 0]
    DE_list[[4]] <- rownames(RNA_DE_0.05)[RNA_DE_0.05$logFC < 0]
    names(DE_list) <- c("up_0.1", "up_0.05", "down_0.1", "down_0.05")
    
    # adjust groups column to add control and repeat type information:
    plot_df$group[plot_df$id %in% pos_chapman_ctl] <- "chapman_upregulated_positive_control"
    plot_df$group[plot_df$id %in% neg_chapman_ctl] <- "chapman_non-DE_negative_control"
    plot_df$group[plot_df$id %in% neg_chapman_ctl2] <- "chapman_downregulated_negative_control"
    plot_df$group[plot_df$id %in% DE_list$up_0.1] <- "up_repeat_FDR<0.1"
    #plot_df$group[plot_df$id %in% DE_list$up_0.05] <- "up_repeat_FDR<0.05"
    plot_df$group[plot_df$id %in% DE_list$down_0.1] <- "down_repeat_FDR<0.1"
    #plot_df$group[plot_df$id %in% DE_list$down_0.05] <- "down_repeat_FDR<0.05"
    plot_df$group[plot_df$group == "ctl_repeat"] <- "non-DE_repeat"
    plot_df <- plot_df[plot_df$group!="repeat",]
    plot_df$group <- factor(plot_df$group, levels = c("chapman_non-DE_negative_control",
                                                      "chapman_downregulated_negative_control", "chapman_upregulated_positive_control",
                                                      "non-DE_repeat", "down_repeat_FDR<0.1", "down_repeat_FDR<0.05", "up_repeat_FDR<0.1", 
                                                      "up_repeat_FDR<0.05"))
    
    #     # # do stats between different groups:
    #     # CPMR_ctl <- plot_df[grep("CPM", plot_df$group),]
    #     # hist(CPMR_ctl$HGSOC_H3K27me3)
    #     # hist(CPMR_ctl$HGSOC_H3K4me3)
    #     # 
    #     # CPMR_ctl_K27_tbl <- table(CPMR_ctl$HGSOC_H3K27me3, as.character(CPMR_ctl$group))
    #     # CPMR_ctl_K27_pval <- chisq.test(CPMR_ctl_K27_tbl)
    #     # 
    #     # CPMR_ctl_K4_tbl <- table(CPMR_ctl$HGSOC_H3K4me3, as.character(CPMR_ctl$group))
    #     # CPMR_ctl_K4_pval <- chisq.test(CPMR_ctl_K4_tbl)
    #     # 
    #     # paper_ctl <- plot_df[grep("H3", plot_df$group),]
    #     # hist(paper_ctl$HGSOC_H3K27me3)
    #     # hist(paper_ctl$HGSOC_H3K4me3)
    #     # 
    #     # paper_ctl_K27_tbl <- table(paper_ctl$HGSOC_H3K27me3, as.character(paper_ctl$group))
    #     # paper_ctl_K27_pval <- chisq.test(paper_ctl_K27_tbl)
    #     # 
    #     # paper_ctl_K4_tbl <- table(paper_ctl$HGSOC_H3K4me3, as.character(paper_ctl$group))
    #     # paper_ctl_K4_pval <- chisq.test(paper_ctl_K4_tbl)
    #     # 
    #     # repeat_df <- plot_df[grep("repeat", plot_df$group),]
    #     # hist(repeat_df$HGSOC_H3K27me3)
    #     # hist(repeat_df$HGSOC_H3K4me3)
    #     # 
    #     # paper_ctl_K27_pval1 <- 
    #     #   wilcox.test(repeat_df$HGSOC_H3K27me3[repeat_df$group == "down_repeat_FDR<0.1"],
    #     #   repeat_df$HGSOC_H3K27me3[repeat_df$group == "up_repeat_FDR<0.1"])
    #     # paper_ctl_K27_pval2 <- 
    #     #   wilcox.test(repeat_df$HGSOC_H3K27me3[repeat_df$group == "non_DE_repeat"],
    #     #               repeat_df$HGSOC_H3K27me3[repeat_df$group == "down_repeat_FDR<0.1"])
    #     # paper_ctl_K27_pval3 <- 
    #     #   wilcox.test(repeat_df$HGSOC_H3K27me3[repeat_df$group == "non_DE_repeat"],
    #     #               repeat_df$HGSOC_H3K27me3[repeat_df$group == "up_repeat_FDR<0.1"])
    #     # 
    #     # paper_ctl_K4_pval1 <- 
    #     #   wilcox.test(repeat_df$HGSOC_H3K4me3[repeat_df$group == "down_repeat_FDR<0.1"],
    #     #               repeat_df$HGSOC_H3K4me3[repeat_df$group == "up_repeat_FDR<0.1"])
    #     # paper_ctl_K4_pval2 <- 
    #     #   wilcox.test(repeat_df$HGSOC_H3K4me3[repeat_df$group == "down_repeat_FDR<0.1"],
    #     #               repeat_df$HGSOC_H3K4me3[repeat_df$group == "up_repeat_FDR<0.1"])
    #     # paper_ctl_K4_pval3 <- 
    #     #   wilcox.test(repeat_df$HGSOC_H3K4me3[repeat_df$group == "down_repeat_FDR<0.1"],
    #     #               repeat_df$HGSOC_H3K4me3[repeat_df$group == "up_repeat_FDR<0.1"])
    #     
    #     
    #     
    m_df <- melt(plot_df)
    colnames(m_df) <- c("group", "id", "mark", "count")
    saveRDS(m_df, file=paste0(RobjectDir, "/", Posit, "_", exp_nos[o], "_m_df.rds"))
    
    # plot on barplot:
    p <- ggplot(peak_stats, aes(x=group, y=count, fill=mark))
    p <- p + geom_bar(stat = "identity", position=position_dodge(0.9))
    p <- p + geom_errorbar(aes(ymin=count-se, ymax=count+se), width=0.2,
                           position=position_dodge(0.9))
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
    p <- p + ylab("Histone marks per sequence")
    pdf(paste0(plotDir, "/peak_no_per_seq_500bp_upstream_se.pdf"))
    print(p)
    dev.off()
    
    ######
    df <- m_df[m_df$mark == "HGSOC_H3K4me3",]
    cont <- df[grep("control", df$group),]
    
    # plot on scatterplot:
    p <- ggplot(cont, aes(x=group, y=count, fill=mark, group=id))
    p <- p + geom_point(stat = "identity", position=position_dodge(0.9))
    # p <- p + geom_errorbar(aes(ymin=count-se, ymax=count+se), width=0.2,
    #                        position=position_dodge(0.9))
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
    p <- p + ylab("Absolute histone marks")
    print(p)
    
    p <- ggplot(cont, aes(x=group, y=count))
    p <- p + geom_violin()
    # p <- p + geom_errorbar(aes(ymin=count-se, ymax=count+se), width=0.2,
    #                        position=position_dodge(0.9))
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
    p <- p + ylab("Absolute histone marks")
    print(p)
    
    
    p <- ggplot(reps, aes(x=group, y=count, fill=mark, group=id))
    p <- p + geom_point(stat = "identity", position=position_dodge(0.9))
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
    p <- p + ylab("Absolute histone marks")
    print(p)
    
    return(plot_df)
  })
  names(res) <- c("upstream")
  saveRDS(res, paste0(RobjectDir, "/", exp_nos[o], "bp_peak_no_per_seq_plot_df.rds"))
}


