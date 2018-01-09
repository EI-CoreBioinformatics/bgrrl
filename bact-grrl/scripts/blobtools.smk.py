import sys
import os
import glob
import csv

TIME_CMD = "/usr/bin/time -v"
SOFTWAREPATH = "/tgac/software/testing"

INPUTDIR = ""
BLOBTOOLS_DIR = os.path.join("Analysis", "qa", "blobtools")
BLOBTOOLS_BWA_DIR = os.path.join(BLOBTOOLS_DIR, "bwa")
BLOBTOOLS_BLAST_DIR = os.path.join(BLOBTOOLS_DIR, "blast")
BLOBTOOLS_BLOB_DIR = os.path.join(BLOBTOOLS_DIR, "blob")
BLOB_BLASTDB = "/ei/public/databases/blast/ncbi/nt_20171013/nt"

BGRRL = "/tgac/workarea/group-pb/schudomc_bact/bact-grrl"
BGRRL_WRAPPERS = os.path.join(BGRRL, "scripts", "wrappers")
BLOB_BWA_WRAPPER = os.path.join(BGRRL_WRAPPERS, "blob_bwa_wrapper")
BLOB_BLAST_WRAPPER = os.path.join(BGRRL_WRAPPERS, "blob_blast_wrapper")
BLOB_WRAPPER = os.path.join(BGRRL_WRAPPERS, "blobtools_wrapper")

INPUTFILES = dict()
for row in csv.reader(open('samples.txt'), delimiter=","):
    INPUTFILES[row[0]] = tuple(row[1:])

from os.path import join
TARGETS = list()
TARGETS.extend(map(lambda s:join(BLOBTOOLS_BWA_DIR, s, s + ".blob_bwa.bam"), INPUTFILES))
TARGETS.extend(map(lambda s:join(BLOBTOOLS_BLAST_DIR, s, s + ".blob_blast.tsv"), INPUTFILES))
TARGETS.extend(map(lambda s:join(BLOBTOOLS_BLOB_DIR, s, s + ".blobDB.table.txt"), INPUTFILES))


def get_sample_files(wc):
	return INPUTFILES[wc.sample]


rule all:
	input: TARGETS

rule blob_bwa_mem:
	input:
		get_sample_files
	output:
		bam = os.path.join(BLOBTOOLS_BWA_DIR, "{sample}", "{sample}.blob_bwa.bam")
	log:
		os.path.join(BLOBTOOLS_BWA_DIR, "log", "{sample}.blob_bwa.log")
		# Analysis/qa/blobtools/bwa/PRO1638_Plate5_C7/Analysis/assembly/unicycler/PRO1638_Plate5_C7/assembly.fasta
	params:
		outdir = os.path.join(BLOBTOOLS_BWA_DIR, "{sample}"),
		index = lambda wildcards: os.path.join(BLOBTOOLS_BWA_DIR, "{sample}", "assembly.fasta"),
		outbam = os.path.join(BLOBTOOLS_BWA_DIR, "{sample}", "{sample}.blob_bwa.bam"),
	threads:
		8
	shell:
		BLOB_BWA_WRAPPER + \
		" {input[0]} {params.outbam} {params.index} -t {threads} {params.index} {input[1]} {input[2]} &> {log}"

rule blob_blast:
	input:
		get_sample_files
	output:
		tsv = os.path.join(BLOBTOOLS_BLAST_DIR, "{sample}", "{sample}.blob_blast.tsv")
	log:
		os.path.join(BLOBTOOLS_BLAST_DIR, "log", "{sample}.blob_blast.log")
	params:
		outdir = os.path.join(BLOBTOOLS_BLAST_DIR, "{sample}")
	threads:
		8
	shell:
		BLOB_BLAST_WRAPPER + \
		" -num_threads {threads} -query {input[0]} -db " + BLOB_BLASTDB + " -out {output.tsv} &> {log}"

rule blob_blobtools:
    input:
        bwa = os.path.join(BLOBTOOLS_BWA_DIR, "{sample}", "{sample}.blob_bwa.bam"),
        blast = os.path.join(BLOBTOOLS_BLAST_DIR, "{sample}", "{sample}.blob_blast.tsv"),
        assembly = get_sample_files
    output:
        os.path.join(BLOBTOOLS_BLOB_DIR, "{sample}", "{sample}.blobDB.table.txt")
    log:
        os.path.join(BLOBTOOLS_BLOB_DIR, "log", "{sample}.blob.log")
    params:
        prefix = os.path.join(BLOBTOOLS_BLOB_DIR, "{sample}", "{sample}")
    threads:
        1
    shell:
        TIME_CMD + \
        " " + BLOB_WRAPPER + \
         " {input.bwa} {input.blast} {input.assembly[0]} {params.prefix} &> {log}"
