library(GenomicRanges)
library(rtracklayer)

reps <- import("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/histone-ChIP-seq/chapman-rothe_2013_hg19/refs/custom3rep.hg19.gtf")

repnames <- c("AluYh3a3", "L1PBa1", "Kanga1d")

for ( i in 1:length(repnames) ) {

	res <- reps[reps$ID==repnames[i]]
	export(res, 
		paste0("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/histone-ChIP-seq/chapman-rothe_2013_hg19/refs/", 
			repnames[i], ".hg19.gtf"))

}
