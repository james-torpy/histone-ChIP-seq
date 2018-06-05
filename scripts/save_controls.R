DE <- readRDS("/Users/jamestorpy/clusterHome//projects/hgsoc_repeats/RNA-seq/Robjects/exp9/htseq_EdgeR_primary_HGSOC_vs_FT_with_chromatin_remodellers/custom3_DEallGenes.rds")
sigDE <- DE[DE$threshold == T,]
sigDE <- sigDE[order(sigDE$FDR),]
sigDE <- sigDE[order(sigDE$logFC, decreasing = T),]
sigDE <- sigDE[!is.na(sigDE$symbol),]
clean_sigDE <- sigDE[grep("^MT|^GOL|^LINC|^MIR", sigDE$symbol, invert=T),]
up50 <- head(clean_sigDE, 50)
down50 <- tail(clean_sigDE, 50)
write.table(up50, file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/results/R/exp9/plots/DEplots/htseq_EdgeR_primary_HGSOC_vs_FT/top_50_up_genes_allsymbols_DE.txt", quote=F, sep="\t")
write.table(down50, file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/results/R/exp9/plots/DEplots/htseq_EdgeR_primary_HGSOC_vs_FT/top_50_down_genes_allsymbols_DE.txt", quote=F, sep="\t")

CPMR <- readRDS("/Users/jamestorpy/clusterHome//projects/hgsoc_repeats/RNA-seq/results/R/exp9/plots/DEplots/htseq_EdgeR_primary_HGSOC_vs_FT/gc_CPMR.rds")
mean_CPMR <- apply(CPMR, 1, mean)
means_ordered <- mean_CPMR[order(mean_CPMR)]
means_ordered <- means_ordered[means_ordered >=1]

library(org.Hs.eg.db)
egENSEMBL <- toTable(org.Hs.egENSEMBL)
egSYMBOL <- toTable(org.Hs.egSYMBOL)

m <- match(names(means_ordered), egENSEMBL$ensembl_id)
ids <- egENSEMBL$gene_id[m]
m <- match(ids, egSYMBOL$gene_id)
sym <- egSYMBOL$symbol[m] 
names(means_ordered) <- sym

meanCPMR <- means_ordered[!is.na(names(means_ordered))]
meanCPMR <- meanCPMR[grep("COX|ND[0-9]|LINC|orf|LOC|MIR", names(meanCPMR), invert=T)]
topCPMR <- names(tail(meanCPMR, 50))
bottomCPMR <- names(head(meanCPMR, 50))

write.table(topCPMR, 
  "/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/results/R/exp9/plots/DEplots/htseq_EdgeR_primary_HGSOC_vs_FT/top_CPMR.txt",
  quote=F, sep="\t")
write.table(bottomCPMR, 
  "/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/results/R/exp9/plots/DEplots/htseq_EdgeR_primary_HGSOC_vs_FT/bottom_CPMR.txt",
  quote=F, sep="\t")



