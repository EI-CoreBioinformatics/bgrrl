import sys
import csv
import os
from os.path import join

from bgrrl import TIME_CMD
from bgrrl.samplesheet import readSamplesheet
from eicore.external_process.snakemake_helper import loadPreCmd

DEBUG = config.get("debugmode", False)

# set up i/o
OUTPUTDIR = config["out_dir"]
QC_OUTDIR = join(OUTPUTDIR, "qc")
FASTQC_DIR = join(QC_OUTDIR, "fastqc")
BBDUK_DIR = join(QC_OUTDIR, "bbduk")
KAT_DIR = join(QC_OUTDIR, "kat")
BBNORM_DIR = join(QC_OUTDIR, "bbnorm")
TADPOLE_DIR = join(QC_OUTDIR, "tadpole")

INPUTFILES = dict(readSamplesheet(config["samplesheet"]))

TARGETS = list()
if not config.get("no_normalization", False):
	PRIMARY_READDIR = BBNORM_DIR
	SECONDARY_READDIR = BBDUK_DIR
	PRIMARY_READ_ID = "bbnorm"
	SECONDARY_READ_ID = "bbduk"
	TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbnorm", s, s + "_R1.bbnorm_fastqc.html"), INPUTFILES))
	TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbnorm", s, s + "_R2.bbnorm_fastqc.html"), INPUTFILES))
else:
	PRIMARY_READDIR = BBDUK_DIR
	SECONDARY_READDIR = BBDUK_DIR
	PRIMARY_READ_ID = "bbduk"
	SECONDARY_READ_ID = "bbduk"

TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbduk", s, s + "_R1.bbduk_fastqc.html"), INPUTFILES))
TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbduk", s, s + "_R2.bbduk_fastqc.html"), INPUTFILES))
TARGETS.extend(map(lambda s:join(TADPOLE_DIR, s, s + "_tadpole_contigs.fasta"), INPUTFILES))
TARGETS.extend(map(lambda s:join(KAT_DIR, s, s + ".dist_analysis.json"), INPUTFILES))

if DEBUG:
	with open("inputfiles.txt", "w") as input_out:
		print(*INPUTFILES.values(), sep="\n", file=input_out)
	with open("targets.txt", "w") as targets_out:
		print(*TARGETS, sep="\n", file=targets_out)

def get_sample_files(wc):
	return INPUTFILES[wc.sample].R1, INPUTFILES[wc.sample].R2

localrules: all

rule all:
	input: TARGETS

rule qc_bbduk:
	input:
		get_sample_files
	output:
		r1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
		r2 = join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
	log:
		join(QC_OUTDIR, "log", "{sample}", "{sample}.qc_bbduk.log")
	threads:
		8
	params:
		load = loadPreCmd(config["load"]["bbmap"]),
		bbduk = "bbduk.sh", 
                # config["tools"]["bbduk"],
		adapters = config["resources"]["bb_adapters"],
                bbduk_params = config["params"]["bbduk"]
	shell:
		# "{params.load}" + TIME_CMD + " {params.bbduk}" + \
		"source activate bgqc_env && {params.bbduk}" + \
		" -Xmx30g t={threads} in1={input[0]} in2={input[1]} out1={output.r1} out2={output.r2}" + \
		" ref={params.adapters}" + \
		" {params.bbduk_params}" + \
		" &> {log}"
		# " ktrim=r k=21 mink=7 hdist=1 qtrim=lr trimq=3 minlen=100 maq=20 tpe tbo &> {log}"
		#" ktrim=r k=21 mink=11 hdist=2 qtrim=lr trimq=3 minlen=100 maq=20 tpe tbo &> {log}"

rule qc_fastqc_bbduk:
	input:
		# join(BBDUK_DIR, "{sample}", "{fastq}.bbduk.fastq.gz")
		join(BBDUK_DIR, "{sample}", "{sample}_R{mate}.fastq.gz")
	output:
		# fqc = join(FASTQC_DIR, "bbduk", "{sample}", "{fastq}.bbduk_fastqc.html")
		join(FASTQC_DIR, "bbduk", "{sample}", "{sample}_R{mate}.bbduk_fastqc.html")
	params:
		outdir = join(FASTQC_DIR, "bbduk", "{sample}"),
                load = loadPreCmd(config["load"]["fastqc"]),
		# fastqc = config["tools"]["fastqc"]
		fastqc = "fastqc"
#	log:
#		# join(QC_OUTDIR, "log", "{fastq}.qc_fastqc_bbduk.log")
#		join(QC_OUTDIR, "log", "{sample}", "{sample}_R{mate}.qc_fastqc_bbduk.log")
	threads:
		2
	shell:
		# "{params.load} (" + TIME_CMD + " {params.fastqc}" + \
		#" --extract --threads={threads} --outdir={params.outdir} {input} " + \
		#" || mkdir -p {params.outdir} && touch {output.fqc}) &> {log}"
		"source activate bgqc_env &&" + \
		" ({params.fastqc} --extract --threads={threads} --outdir={params.outdir} {input}" + \
		" || mkdir -p {params.outdir} && touch {output})"  + \
		""
#		" &> {log}"

if not config.get("no_normalization", False):
	rule qc_bbnorm:
		input:
			r1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
			r2 = join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
		output:
			r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
			r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz"),
			prehist = join(BBNORM_DIR, "{sample}", "{sample}.bbnorm.pre.hist"),
			posthist = join(BBNORM_DIR, "{sample}", "{sample}.bbnorm.post.hist")
		params:
			load = loadPreCmd(config["load"]["bbmap"]),
			bbnorm = "bbnorm.sh", 
                        # config["tools"]["bbnorm"],
			bbnorm_params = config["params"]["bbnorm"]
		log:
			join(QC_OUTDIR, "log", "{sample}.qc_bbnorm.log")
		threads:
			8
		shell:
			"source activate bgqc_env &&" + \
			" {params.bbnorm}" + \
			# "{params.load}" + TIME_CMD + " {params.bbnorm}" + \
			" -Xmx30g t={threads} in={input.r1} in2={input.r2} out={output.r1} out2={output.r2}" + \
			# " target=100 min=2 prefilter" + \
			" {params.bbnorm_params}" + \
			" khist={output.prehist} khistout={output.posthist} &> {log}"

	rule qc_fastqc_bbnorm:
		input:
			# join(BBNORM_DIR, "{sample}", "{fastq}.bbnorm.fastq.gz")
			join(BBNORM_DIR, "{sample}", "{sample}_{mate}.bbnorm.fastq.gz")
		output:
			fqc = join(FASTQC_DIR, "bbnorm", "{sample}", "{sample}_{mate}.bbnorm_fastqc.html")
		params:
			outdir = join(FASTQC_DIR, "bbnorm", "{sample}"),
                	load = loadPreCmd(config["load"]["fastqc"]),
			# fastqc = config["tools"]["fastqc"]
			fastqc = "fastqc"
		log:
			join(QC_OUTDIR, "log", "{sample}_{mate}.qc_fastqc_bbnorm.log")
		threads:
			2
		shell:
			# "{params.load} ({params.fastqc}" + \
			# " --extract --threads={threads} --outdir={params.outdir} {input}" + \
			# " || mkdir -p {params.outdir} && touch {output.fqc}) &> {log}"
			"source activate bgqc_env &&" + \
			" ({params.fastqc} --extract --threads={threads} --outdir={params.outdir} {input}" + \
			" || mkdir -p {params.outdir} && touch {output.fqc}) &> {log}"

rule qc_tadpole:
	input:
		# r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
		r1 = join(PRIMARY_READDIR, "{sample}", "{sample}_R1." + PRIMARY_READ_ID + ".fastq.gz"),
		# r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz")
		r2 = join(PRIMARY_READDIR, "{sample}", "{sample}_R2." + PRIMARY_READ_ID + ".fastq.gz")
	output:
		contigs = join(TADPOLE_DIR, "{sample}", "{sample}_tadpole_contigs.fasta")
	params:
		load = loadPreCmd(config["load"]["bbmap"]),
		# tadpole = config["tools"]["tadpole"]
                tadpole = "tadpole.sh"
	log:
		join(QC_OUTDIR, "log", "{sample}.qc_tadpole.log")
	threads:
		8
	shell:
		# "{params.load}" + TIME_CMD + " {params.tadpole}" + \
		"source activate bgqc_env && {params.tadpole}" + \
		" -Xmx30g threads={threads} in={input.r1} in2={input.r2} out={output.contigs} &> {log}"

rule qc_katgcp:
	input:
		# r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
		r1 = join(PRIMARY_READDIR, "{sample}", "{sample}_R1." + PRIMARY_READ_ID + ".fastq.gz"),
		# r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz")
		r2 = join(PRIMARY_READDIR, "{sample}", "{sample}_R2." + PRIMARY_READ_ID + ".fastq.gz")
	output:
		# katgcp = join(KAT_DIR, "{sample}", "{sample}.kat-gcp.mx")
		katgcp = join(KAT_DIR, "{sample}", "{sample}.dist_analysis.json")
	log:
		join(KAT_DIR, "{sample}", "{sample}.kat")
	params:
		# prefix = lambda wildcards: join(KAT_DIR, wildcards.sample, wildcards.sample + ".kat-gcp"),
		prefix = lambda wildcards: join(KAT_DIR, wildcards.sample, wildcards.sample), # + ".dist_analysis.json"),
		load = loadPreCmd(config["load"]["kat"])
	threads:
		2
	shell:
		"touch {output.katgcp} &> {log}"
		# "{params.load}" + \
		# " (kat gcp -o {params.prefix} -t {threads} -v {input.r1} {input.r2} || touch {output.katgcp}) &> {log}"

