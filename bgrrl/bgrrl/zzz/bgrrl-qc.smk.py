import sys
import csv
import os
import glob
from os.path import join, basename

from bgrrl.bgrrl import Sample
from bgrrl import loadPreCmd

TIME_CMD = " /usr/bin/time -v"

BGRRL = "/tgac/workarea/group-pb/schudomc_bact/bact-grrl"
READCOUNT = os.path.join(BGRRL, "scripts", "readcount")

SOFTWAREPATH = "/tgac/software/testing"
BBSUITE_DIR = os.path.join(SOFTWAREPATH, "bbmap", "37.24", "bbmap")
ADAPTERS = os.path.join(BBSUITE_DIR, "resources", "adapters.fa")

# tools
BBDUK = os.path.join(BBSUITE_DIR, "bbduk.sh")
BBNORM = os.path.join(BBSUITE_DIR, "bbnorm.sh")
TADPOLE = os.path.join(BBSUITE_DIR, "tadpole.sh")
FASTQC = os.path.join(SOFTWAREPATH, "fastqc", "0.11.5", "x86_64", "bin", "fastqc")

# wrappers
WRAPPERS = os.path.join(BGRRL, "scripts", "wrappers")
KAT_WRAPPER = os.path.join(WRAPPERS, "kat_wrapper")

OUTPUTDIR = config["out_dir"]

QC_OUTDIR = os.path.join(OUTPUTDIR, "qc")
READCOUNT_DIR = os.path.join(QC_OUTDIR, "readcount")
FASTQC_DIR = os.path.join(QC_OUTDIR, "fastqc")
BBDUK_DIR = os.path.join(QC_OUTDIR, "bbduk")
KATHIST_DIR = os.path.join(QC_OUTDIR, "kat", "hist")
KATGCP_DIR = os.path.join(QC_OUTDIR, "kat", "gcp")
BBNORM_DIR = os.path.join(QC_OUTDIR, "bbnorm")
TADPOLE_DIR = os.path.join(QC_OUTDIR, "tadpole")

# Sample = namedtuple("Sample", "sampleID customerSampleID R1 R2 S taxonomyID taxonomyTxt fastqcR1 fastqcR2 fastqcS".split(" "))
INPUTFILES = dict((row[0], Sample(*row)) for row in csv.reader(open(config["samplesheet"]), delimiter=","))
with open("inputfiles.txt", "w") as input_out:
    print(*INPUTFILES.values(), sep="\n", file=input_out)

TARGETS = list()
TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbduk", s, s + "_R1.bbduk_fastqc.html"), INPUTFILES))
TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbduk", s, s + "_R2.bbduk_fastqc.html"), INPUTFILES))
TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbnorm", s, s + "_R1.bbnorm_fastqc.html"), INPUTFILES))
TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbnorm", s, s + "_R2.bbnorm_fastqc.html"), INPUTFILES))
# TARGETS.extend(map(lambda s:join(TADPOLE_DIR, s, s + "_tadpole_contigs.fasta"), INPUTFILES))
TARGETS.extend(map(lambda s:join(TADPOLE_DIR, s, "quast", "quast.log"), INPUTFILES))

with open("targets.txt", "w") as targets_out:
	print(*TARGETS, sep="\n", file=targets_out)

def get_sample_files(wc):
	return INPUTFILES[wc.sample].R1, INPUTFILES[wc.sample].R2


rule all:
	input: TARGETS

rule qc_bbduk:
	input:
		get_sample_files
	output:
		r1 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
		r2 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
	log:
		os.path.join(QC_OUTDIR, "log", "{sample}.qc_bbduk.log")
	threads:
		8
	params:
		load = loadPreCmd(config["load"]["bbmap"])
	shell:
		"{params.load}" + TIME_CMD + " " + BBDUK + \
		" -Xmx30g t={threads} in1={input[0]} in2={input[1]} out1={output.r1} out2={output.r2}" + \
		" ref=" + ADAPTERS + \
		" ktrim=r k=21 mink=11 hdist=2 qtrim=lr trimq=3 minlen=100 maq=20 tpe tbo &> {log}"

rule qc_fastqc_bbduk:
	input:
		os.path.join(BBDUK_DIR, "{sample}", "{fastq}.bbduk.fastq.gz")
	output:
		fqc = os.path.join(FASTQC_DIR, "bbduk", "{sample}", "{fastq}.bbduk_fastqc.html")
	params:
		outdir = os.path.join(FASTQC_DIR, "bbduk", "{sample}"),
                load = loadPreCmd(config["load"]["fastqc"])
	log:
		os.path.join(QC_OUTDIR, "log", "{fastq}.qc_fastqc_bbduk.log")
	threads:
		2
	shell:
		"{params.load}" + TIME_CMD + " " + FASTQC + \
		" --extract --threads={threads} --outdir={params.outdir} {input} &> {log}"

rule qc_bbnorm:
	input:
		r1 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
		r2 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
	output:
		r1 = os.path.join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
		r2 = os.path.join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz"),
		prehist = os.path.join(BBNORM_DIR, "{sample}", "{sample}.bbnorm.pre.hist"),
		posthist = os.path.join(BBNORM_DIR, "{sample}", "{sample}.bbnorm.post.hist")
	params:
		load = loadPreCmd(config["load"]["bbmap"])
	log:
		os.path.join(QC_OUTDIR, "log", "{sample}.qc_bbnorm.log")
	threads:
		8
	shell:
		"{params.load}" + TIME_CMD + " " + BBNORM + \
		" -Xmx30g t={threads} in={input.r1} in2={input.r2} out={output.r1} out2={output.r2}" + \
		" target=100 min=2 prefilter" + \
		" khist={output.prehist} khistout={output.posthist} &> {log}"

rule qc_fastqc_bbnorm:
	input:
		os.path.join(BBNORM_DIR, "{sample}", "{fastq}.bbnorm.fastq.gz")
	output:
		fqc = os.path.join(FASTQC_DIR, "bbnorm", "{sample}", "{fastq}.bbnorm_fastqc.html")
	params:
		outdir = os.path.join(FASTQC_DIR, "bbnorm", "{sample}"),
                load = loadPreCmd(config["load"]["fastqc"])
	log:
		os.path.join(QC_OUTDIR, "log", "{fastq}.qc_fastqc_bbnorm.log")
	threads:
		2
	shell:
		"{params.load}" + TIME_CMD + " " + FASTQC + \
		" --extract --threads={threads} --outdir={params.outdir} {input} &> {log}"

rule qc_tadpole:
	input:
		r1 = os.path.join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
		r2 = os.path.join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz")
	output:
		contigs = os.path.join(TADPOLE_DIR, "{sample}", "{sample}_tadpole_contigs.fasta")
	params:
		load = loadPreCmd(config["load"]["bbmap"])
	log:
		os.path.join(QC_OUTDIR, "log", "{sample}.qc_tadpole.log")
	threads:
		8
	shell:
		"{params.load}" + TIME_CMD + " " + TADPOLE + \
		" -Xmx30g threads={threads} in={input.r1} in2={input.r2} out={output.contigs} &> {log}"


rule qa_quast:
	input:
		contigs = os.path.join(TADPOLE_DIR, "{sample}", "{sample}_tadpole_contigs.fasta")
	output:
		os.path.join(TADPOLE_DIR, "{sample}", "quast", "quast.log")
	params:
		load = loadPreCmd(config["load"]["quast"]),
		outdir = lambda wildcards: os.path.join(TADPOLE_DIR, wildcards.sample, "quast"),
		cmd = loadPreCmd(config["cmd"]["quast"], is_dependency=False)
	log:
		os.path.join(QC_OUTDIR, "log", "{sample}.qc_quast_tadpole.log")
	threads:
		2
	shell:
		"{params.load}" + TIME_CMD + \
		" {params.cmd} -o {params.outdir} -t {threads} -L -s {input.contigs} --min-contig 0 &> {log}"


'''
rule qc_kat_hist:
	input:
		r1 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
		r2 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
	output:
		kat_hist = os.path.join(KATHIST_DIR, "{sample}", "{sample}.kat.hist")
	log:
		os.path.join(QC_OUTDIR, "log", "{sample}.qc_kathist.log")
	params:
		prefix = lambda wildcards: os.path.join(KATHIST_DIR, wildcards.sample, wildcards.sample + ".kat.hist")
	threads:
		8
	shell:
		KAT_WRAPPER + \
		#" hist -o {params.prefix} -t {threads} -v \<\(zcat {input.r1} {input.r2}\) &> {log}"
		" {input.r1} {input.r2} hist -o {params.prefix} -t {threads} -v &> {log}"

rule qc_kat_gcp:
	input:
		r1 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
		r2 = os.path.join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
	output:
		kat_gcp = os.path.join(KATGCP_DIR, "{sample}", "{sample}.kat-gcp.mx")
	log:
		os.path.join(QC_OUTDIR, "log", "{sample}.qc_katgcp.log")
	params:
		prefix = lambda wildcards: os.path.join(KATGCP_DIR, wildcards.sample, wildcards.sample + ".kat-gcp")
	threads:
		8
	shell:
		KAT_WRAPPER + \
		" {input.r1} {input.r2} gcp -o {params.prefix} -t {threads} -v &> {log}"
'''
