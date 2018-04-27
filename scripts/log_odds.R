#This script takes the results of bam Robjects and makes FPKMs

library(GenomicRanges)
#library(rtracklayer)
#library(ggplot2)
library(RColorBrewer)
#library(gwascat)
library(traseR)
library(R.utils)
library("BSgenome.Hsapiens.UCSC.hg19")
library(gwascat)
library(rtracklayer)
args <- R.utils::commandArgs(asValues=TRUE)
type="GWAScapseq"
snpType="NHGRI" 
if (!is.null(args[["type"]])){type = args$type}
if (!is.null(args[["snpType"]])){snpType = args$snpType}


timeStamp <- format(Sys.time(), "%Y_%m_%d")
chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
homedir="/share/ClusterShare/biodata/contrib/nenbar/"
projectname="tissues"
projectDir =paste0(homedir,"projects/melanoma/project_results/")
robjectsDir = paste(projectDir,"gwas_capseq.Robjects/",sep="")
annotDir = paste0(homedir,"projects/melanoma/annotation/")
tracksDir = paste(annotDir,"captureSpace/",sep="")
outDir=paste0(homedir,"projects/melanoma/project_results/SNPs/allregions/")
outRaw=paste0(homedir,"projects/melanoma/project_results/SNPs/raw_data/")
system(paste0("mkdir -p ",outDir))
######### Load chromosomal locations ###############
#load(paste0("../../project_results/transcript_assembly.Rdata/merged_gtf_",projectname,".Rdata"))

######### Load normalized expression data ###############
#load(paste0("../../project_results/transcript_assembly.Rdata/all_rpkms_",projectname,".Rdata"))



#1. Load the capture regions
#load(paste0("../../project_results/transcript_assembly.Rdata/capture_regions_short_",projectname,".Rdata"))
load(paste0(robjectsDir,"capture_transcripts_start2end_iso.Rdata"))
#load(paste0(robjectsDir,"capture_confident_iso.Rdata")) 

##################### Find overlaps with regions
#capture<-capture_regions_short[[1]]


######### This is easy and comes in few steps

#Fetch information on protein coding, pseudogenes, lncRNAs and novel transcripts
#Extract regions like for the novel (super easy)
#Calculate enrichment odds ratios for SNPs - transcript SNP overlaps versus shuffled transcript SNP overlaps

data(gwrngs19)
gwasGR<-as(gwrngs19,"GRanges")
gwasGRSig<-gwasGR[values(gwasGR)$p.Value<=5e-08]
#import haploblocks from https://academic.oup.com/bioinformatics/article/32/2/283/1743626/Approximately-independent-linkage-disequilibrium#supplementary-data
hb<-import("/share/ClusterShare/biodata/contrib/nenbar/projects/melanoma/scripts/figure4/ldetect-data/EUR/fourier_ls-all.bed")
#hb<-import("/share/ClusterShare/biodata/contrib/nenbar/projects/melanoma/scripts/figure3/haploblocks_with_annotation.bed")
mat<-findOverlaps(hb,gwasGRSig)
#instead of the whole set of haploblocks, use only the part that was analysed
hbShort<-hb[unique(queryHits(mat))]

#Find percentages that overlap with the probes
#tiled<-import("/share/ClusterShare/biodata/contrib/nenbar/projects/melanoma/scripts/figure3/all_haploblocks.bed")
#tiled<-reduce(tiled)
#mat<-findOverlaps(tiled,gr)
#tiled_short<-tiled[unique(queryHits(mat))]
#sum(countOverlaps(gr,tiled_short)>0)

if(snpType=="GWAS"){
	data(taSNPDB)
	snpdb=taSNPDB
} else if (snpType=="NHGRI"){
	#try traseR
	library(gwascat)
	data(gwrngs19)
	gwasGR<-as(gwrngs19,"GRanges")
	gwasGRSig<-gwasGR[values(gwasGR)$p.Value<=5e-08]
	gwasGRSig=gwasGRSig[unique(subjectHits(mat))]
	snpdb=gwasGRSig
} else if(snpType=="haploblock") {
	data(taSNPLDDB)
	snpdb=taSNPLDDB
} else if(snpType=="GTeX"){
	load("/share/ClusterShare/biodata/contrib/nenbar/projects/melanoma/scripts/figure4/snpsGR_normCHR.Rdata")
	snpdb=GTEx_snpsGR_normCHR
} else if(snpType=="SNPs"){
	snps<-import("/home/nenbar/local/lib/R-3.2.0/library/LOLA/extdata/hg19/ucsc_example/backup/snp/snps.bed")
	#if(type=="GWAScapseq"){
	#	mat<-findOverlaps(tiled_short,snps)
	#} else {
		mat<-findOverlaps(hbShort,snps)
	#}
	snps=snps[unique(subjectHits(mat))]
	snpdb=snps
}


#chrGR<-GRanges(seqnames=names(chrs),IRanges(1,chrs))
#strand(exclude)<-"*"
#chrGR<-setdiff(chrGR,exclude)

oddsRatio<-function(region,snpDatabase){

	#for a given nucleotide in region, what are the odds it has a SNP
	totalNuc<-sum(width(region))
	reducedsnpDatabase=reduce(snpDatabase)

	#mat<-findOverlaps(region,reducedsnpDatabase)
	#total hits is the length of intersect
	strand(region)="*"
	strand(reducedsnpDatabase)="*"

	totalHits<-sum(width(intersect(region,reducedsnpDatabase)))

	#totalHits<-length(unique(subjectHits(mat)))
	firstRatio<-totalHits/(totalNuc-totalHits)

	#if(type=="GWAScapseq"){
		#nonRegion<-sum(as.numeric(width(tiled_short)))-totalNuc
		#} else {
	nonRegion<-sum(as.numeric(width(hbShort)))-totalNuc	
		#}

	genomeHits<-sum(as.numeric(width(reducedsnpDatabase)))-totalHits
	secondRatio<-genomeHits/(nonRegion-genomeHits)

	odds<-firstRatio/secondRatio

	se<-sqrt(1/totalNuc+1/totalHits+1/nonRegion+1/genomeHits)
	cat(odds)
	cat("\n")
	cat(se)
	cat("\n")
	result<-list()
	result[[1]]=odds
	result[[2]]=se

	#p-value
	SNPenrichment <-
     matrix(c(as.numeric(totalHits), as.numeric(totalNuc), as.numeric(genomeHits), as.numeric(nonRegion)),
            nrow = 2,
            dimnames = list(Regions = c("WithSNP", "Total"),
                            Genome = c("WithSNP", "Total")))
	result[[3]]=chisq.test(SNPenrichment)$p.value

	return(result)

#backup that works
	#totalNuc<-sum(width(region))
	#reducedsnpDatabase=reduce(snpDatabase)
#
	#mat<-findOverlaps(region,reducedsnpDatabase)
	#totalHits<-length(unique(subjectHits(mat)))
	#firstRatio<-totalHits/(totalNuc-totalHits)
#
	#nonRegion<-sum(as.numeric(chrs))-totalNuc
	#genomeHits<-length(reducedsnpDatabase)-totalHits
	#secondRatio<-genomeHits/(nonRegion-genomeHits)
#
	#odds<-firstRatio/secondRatio
	#se<-sqrt(1/totalNuc+1/totalHits+1/nonRegion+1/genomeHits)
	#cat(odds)
	#cat("\n")
	#cat(se)
	#cat("\n")
	#result<-list()
	#result[[1]]=odds
	#result[[2]]=se
#
	##p-value
	#SNPenrichment <-
   #  matrix(c(as.numeric(totalHits), as.numeric(totalNuc), as.numeric(genomeHits), as.numeric(nonRegion)),
   #         nrow = 2,
   #         dimnames = list(Regions = c("WithSNP", "Total"),
   #                         Genome = c("WithSNP", "Total")))
	#result[[3]]=chisq.test(SNPenrichment)$p.value
#
	#return(result)

}

results<-list()
load(paste0(homedir,"projects/melanoma/project_results/SNPs/regions/cleanRegions_",type,".Rdata"))
cat(type)
cat("\n")
temp<-list()
cleanRegions=cleanRegions[c(6,9)]
for(i in 1:length(cleanRegions)){
	cat(names(cleanRegions)[i])
	cat("\n")
	cleanReg=cleanRegions[[i]]
	strand(cleanReg)="*"
	#if(type=="GWAScapseq"){
	#	intersectGR<-intersect(cleanReg,tiled_short)
	#} else {
		intersectGR<-intersect(cleanReg,hbShort)
	#}
	stats<-unlist(oddsRatio(intersectGR,snpdb))
	names(stats)<-c("oddsRatio","SE","p.value")
	results[[names(cleanRegions)[i]]]<-stats		
}

cat("\n")
cat("\n")

df<-do.call("cbind",results)
write.table(df,file=paste0(outDir,"odds_",type,"_",snpType,".txt"),quote=F)
