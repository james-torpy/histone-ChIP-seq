# before running snakemake, run following commands:
#source activate snakemake
#briR

import os
import re

exp_name = 'H3K27ac_cell_lines'

home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + '/projects/hgsoc_repeats/histone-ChIP-seq/'
in_dir = project_dir + '/raw_files/sra/' + exp_name
SAMPLES = [ re.sub('.sra', '', x) for x in list(os.walk(in_dir))[0][2] ]

#SAMPLES = ['SRR600983-input', 'SRR600956-H3K4me3', 
#	'SRR600559-H3K27me3']

print(SAMPLES)

rule all:
	input:
		expand('results/spp/' + exp_name + '/{sample}/{sample}.spp.result', sample=SAMPLES),
		expand('results/bwa/' + exp_name + '/{sample}/{sample}.sorted.bam.bai', sample=SAMPLES),
		expand('genomes/repeats/bwa/hg38_ercc_all_repbase.bwa.amb', sample=SAMPLES),
		expand('results/fastqc/' + exp_name + '/{sample}/{sample}_pass_fastqc.html', sample=SAMPLES)

rule fastq_dump:
    input:
        expand('raw_files/sra/' + exp_name + '/{sample}.sra', sample=SAMPLES)
    output:
        'raw_files/fastq/' + exp_name + '/{sample}_pass.fastq.gz'
    shell:
        'fastq-dump --outdir raw_files/fastq/' + exp_name + '/ ' +
        '--gzip --skip-technical  --readids --read-filter ' +
        'pass --dumpbase --split-3 --clip {input}'

rule fastqc:
	input:
		'raw_files/fastq/' + exp_name + '/{sample}_pass.fastq.gz'
	output:
		'results/fastqc/' + exp_name + '/{sample}/{sample}_pass_fastqc.html',
		'results/fastqc/' + exp_name + '/{sample}/{sample}_pass_fastqc.zip'
	shell:
		'fastqc -t 2 -o results/fastqc/' + exp_name + '/ {input}'

rule bwa_index:
	input:
		'genomes/repeats/hg38_ercc_all_repbase.fa'
	output:
		'genomes/repeats/bwa/hg38_ercc_all_repbase.bwa.amb'
		'genomes/repeats/bwa/hg38_ercc_all_repbase.bwa.ann'
		'genomes/repeats/bwa/hg38_ercc_all_repbase.bwa.bwt'
		'genomes/repeats/bwa/hg38_ercc_all_repbase.bwa.pac'
		'genomes/repeats/bwa/hg38_ercc_all_repbase.bwa.sa'
	shell:
		'bwa index -p genomes/repeats/bwa/hg38_ercc_all_repbase.bwa {input}'

rule bwa_align:
	input:
		fasq = 'raw_files/fastq/' + exp_name + '/{sample}_pass.fastq.gz',
		fqc = 'results/fastqc/' + exp_name + '/{sample}/{sample}_pass_fastqc.html'
	output:
		'results/bwa/' + exp_name + '/{sample}/{sample}.sai'
	shell:
		'bwa aln -f {output} genomes/repeats/bwa/hg38_ercc_all_repbase.bwa {input.fasq}'

rule bwa_samse:
	input:
		sai = 'results/bwa/' + exp_name + '/{sample}/{sample}.sai',
		fq = 'raw_files/fastq/' + exp_name + '/{sample}/{sample}_pass.fastq.gz'
	output:
		'results/bwa/' + exp_name + '/{sample}/{sample}.sam'
	shell:
		'bwa samse -f {output} genomes/repeats/bwa/hg38_ercc_all_repbase.bwa \
			{input.sai} {input.fq}'

rule samtools_filter:
	input:
		'results/bwa/' + exp_name + '/{sample}/{sample}.sam'
	output:
		'results/bwa/' + exp_name + '/{sample}/{sample}.sorted.bam'
	shell:
		'samtools view -u -F4 {input}  | samtools sort -o {output} \
		-O bam -T results/bwa/' + exp_name + '/{wildcards.sample}/{wildcards.sample}.tmp -@ 2 -'


rule samtools_index:
	input:
		'results/bwa/' + exp_name + '/{sample}/{sample}.sorted.bam'
	output:
		'results/bwa/' + exp_name + '/{sample}/{sample}.sorted.bam.bai'
	shell:
		'samtools index {input}'

rule remove_dupes:
	input:
		bam = 'results/bwa/' + exp_name + '/{sample}/{sample}.sorted.bam',
		bai = 'results/bwa/' + exp_name + '/{sample}/{sample}.sorted.bam.bai'
	output:
		rmdup = 'results/bwa/' + exp_name + '/{sample}/{sample}.sorted.rmdup.bam',
		met = 'results/bwa/' + exp_name + '/{sample}/{sample}.sorted.rmdupmet.bam'
	shell:
		'java -jar /home/jamtor/local/lib/picard/build/libs/picard.jar \
		MarkDuplicates I={input.bam} O={output.rmdup} M={output.met} \
		REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT'

rule spp:
	input:
		'results/bwa/' + exp_name + '/{sample}/{sample}.sorted.rmdup.bam'
	output:
		'results/spp/' + exp_name + '/{sample}/{sample}.spp.result'
	shell:
		'Rscript /home/jamtor/local/lib/phantompeakqualtools/run_spp.R \
		-c={input} -i={input} -savp -p=2 -odir=results/spp -out={output}'

#import subprocess


