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


#####

library(org.Hs.eg.db)

#DE <- read.table("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/histone-ChIP-seq/chapman-rothe_2013_hg19/refs/chapman_roethe_DE.txt", 
#	sep = "\t", header = T, fill = T)
DE <- read.table("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/histone-ChIP-seq/chapman-rothe_2013_hg19/refs/chapman_roethe_DE.txt", 
	sep = "\t", header = T, fill = T)

nonDE <- DE[DE$logFC > -.002 & DE$logFC < .002,]
nonDE <- nonDE[!(nonDE$Gene.symbol %in% ""),]
nonDE <- nonDE[grep("LOC|LINC|FAM|\\/|GOLG", nonDE$Gene.symbol, invert=T),]

egENSEMBL <- toTable(org.Hs.egENSEMBL)
egSYMBOL <- toTable(org.Hs.egSYMBOL)

m <- match(nonDE$Gene.symbol, egSYMBOL$symbol)
nonDE$gene_id <- egSYMBOL$gene_id[m]
m <- match(nonDE$gene_id, egENSEMBL$gene_id)
nonDE$ensembl_id <- egENSEMBL$ensembl_id[m]

nonDE <- nonDE[!is.na(nonDE$ensembl_id),]
nonDE <- head(nonDE, 50)

nonDEfinal <- subset(nonDE, select=c(Gene.symbol, ensembl_id))

write.table("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/histone-ChIP-seq/chapman-rothe_2013_hg19/refs/chapman-roethe_non-DE_RNA_symbols_ids.txt",
	col.names = F, row.names = F, quote = F, sep = "\t")

DE <- DE[DE$P.Value < 0.05,]
posDE <- DE[DE$logFC >= 1,]
write.table(posDE, 
  "/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/histone-ChIP-seq/chapman-rothe_2013_hg19/refs/chapman-roethe_up_DE_RNA_symbols_ids.txt")
negDE <- DE[DE$logFC <= 1,]
write.table(negDE, 
            "/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/histone-ChIP-seq/chapman-rothe_2013_hg19/refs/chapman-roethe_down_DE_RNA_symbols_ids.txt")



