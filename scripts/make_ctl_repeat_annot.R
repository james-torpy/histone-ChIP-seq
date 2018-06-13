library("GenomicRanges")
library("rtracklayer")

#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/jamestorpy/clusterHome/"
project_dir <- paste0(home_dir, 
	"/projects/hgsoc_repeats/histone-ChIP-seq/chapman-rothe_2013_hg19/")
ref_dir <- paste0(project_dir, "/refs/")
Robject_dir <- paste0(project_dir, "/Robjects/")
RNA_dir <- paste0(home_dir, "/projects/hgsoc_repeats/RNA-seq/")
DE_dir <- paste0(RNA_dir,
	"results/R/exp9/plots/DEplots/htseq_EdgeR_primary_HGSOC_vs_FT/")

system(paste0("mkdir -p ", Robject_dir))

reps <- import(paste0(ref_dir, "/custom3rep.hg19.sorted.gtf"))
reps <- split(reps, reps$ID)

RNA_DE_0.1 <- read.table(file = paste0(DE_dir, "sig_reps_FDR_0.1.txt"))
upRNA <- rownames(RNA_DE_0.1)[RNA_DE_0.1$logFC > 0]
downRNA <- rownames(RNA_DE_0.1)[RNA_DE_0.1$logFC < 0]
bothRNA <- c(upRNA, downRNA)

DEreps <- reps[names(reps) %in% bothRNA]

other_reps <- import(paste0(ref_dir, "/hg19.repeats.gff"))
other_reps$group <- gsub("ID=", "", other_reps$group)
other_reps <- split(other_reps, other_reps$group)

find_ctls <- function(DE, other) {

	# set up empty vector to fill:
	ctl_res <- GRangesList()

	# define lengths of GRanges objects:
	DE_len <- unlist(lapply(DE, function(x) length(x)))
	other_len <- unlist(lapply(other, function(x) length(x)))

	# iterate through DE_len and find similar length in other_len
	# for each one:
	for ( i in 1:length(DE_len) ) {

		print(i)

		added <- FALSE
		maxi <- DE_len[i] + 10
		mini <- DE_len[i] - 10

		for ( j in 1:length(other_len) ) {

			if ( other_len[j] <= maxi & other_len[j] >= mini & length(ctl_res) < i & 
				!(names(other[j]) %in% names(ctl_res)) ) {

				ctl_res <- c(ctl_res, other[j])
				print("error range of 20")
				added <- TRUE

			}
		}

		for ( j in 1:length(other_len) ) {

			if ( !added ) {
	
				maxi <- DE_len[i] + 100
				mini <- DE_len[i] - 100
	
				if ( other_len[j] <= maxi & other_len[j] >= mini 
					& length(ctl_res) < i & !(names(other[j]) %in% names(ctl_res)) ) {
	
					ctl_res <- c(ctl_res, other[j])
					print("error range of 200")
					added <- TRUE
	
				}
			}
		}
		
		for ( j in 1:length(other_len) ) {

			if ( !added ) {
	
				maxi <- DE_len[i] + 1000
				mini <- DE_len[i] - 1000
	
				if ( other_len[j] <= maxi & other_len[j] >= mini 
					& length(ctl_res) < i & !(names(other[j]) %in% names(ctl_res)) ) {
	
					ctl_res <- c(ctl_res, other[j])
					print("error range of 2000")
					added <- TRUE
	
				}
			}
		}
	
		for ( j in 1:length(other_len) ) {

			if ( !added ) {
	
				maxi <- DE_len[i] + 10000
				mini <- DE_len[i] - 10000
	
				if ( other_len[j] <= maxi & other_len[j] >= mini 
					& length(ctl_res) < i & !(names(other[j]) %in% names(ctl_res)) ) {
	
					ctl_res <- c(ctl_res, other[j])
					print("error range of 20000")
					added <- TRUE
	
				}
			}
		}
	
		for ( j in 1:length(other_len) ) {

			if ( !added ) {
	
				maxi <- DE_len[i] + 50000
				mini <- DE_len[i] - 50000
	
				if ( other_len[j] <= maxi & other_len[j] >= mini 
					& length(ctl_res) < i & !(names(other[j]) %in% names(ctl_res)) ) {
	
					ctl_res <- c(ctl_res, other[j])
					print("error range of 100000")
					added <- TRUE
	
				}
			}
		}

		for ( j in 1:length(other_len) ) {

			if ( !added ) {
	
				maxi <- DE_len[i] + 100000
				mini <- DE_len[i] - 100000
	
				if ( other_len[j] <= maxi & other_len[j] >= mini 
					& length(ctl_res) < i & !(names(other[j]) %in% names(ctl_res)) ) {
	
					ctl_res <- c(ctl_res, other[j])
					print("error range of 200000")
					added <- TRUE
	
				}
			}
		}

		for ( j in 1:length(other_len) ) {

			if ( !added ) {
	
				maxi <- DE_len[i] + 500000
				mini <- DE_len[i] - 500000
	
				if ( other_len[j] <= maxi & other_len[j] >= mini 
					& length(ctl_res) < i & !(names(other[j]) %in% names(ctl_res)) ) {
	
					ctl_res <- c(ctl_res, other[j])
					print("error range of 1000000")
					added <- TRUE
	
				}
			}
		}

		for ( j in 1:length(other_len) ) {

			if ( !added & length(ctl_res) < i ) {

				print("nope")
		
			}
		}
	}

	return(ctl_res)

}
ctl_reps <- find_ctls(DEreps, other_reps)

saveRDS(ctl_reps, paste0(Robject_dir, "/ctl_reps.rds"))

# check annotation lengths and plot distributions:
ctl_lengths <- unlist(lapply(ctl_reps, length))
names(ctl_lengths) <- NULL
ctl_lengths <- sort(ctl_lengths)
DE_lengths <- unlist(lapply(DEreps, length))
names(DE_lengths) <- NULL
ctl_lengths <- sort(DE_lengths)

pdf(paste0(ref_dir, "/ctl_repeat_dist.pdf"))
hist(ctl_lengths)
dev.off()

pdf(paste0(ref_dir, "/ctl_repeat_dist_below_30000.pdf"))
hist(ctl_lengths[ctl_lengths < 30000])
dev.off()

pdf(paste0(ref_dir, "/DE_repeat_dist.pdf"))
hist(DE_lengths)
dev.off()

pdf(paste0(ref_dir, "/DE_repeat_dist_below_30000.pdf"))
hist(DE_lengths[DE_lengths < 30000])
dev.off()

write.table(names(ctl_reps), file = paste0(ref_dir, 
  "/control_repeats.txt"), sep = "\t", quote = F)





