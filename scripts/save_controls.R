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