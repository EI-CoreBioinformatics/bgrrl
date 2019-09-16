import sys
import csv
import os
from os.path import join

from bgrrl.samplesheet import readSamplesheet, Samplesheet 
from bgrrl.snakemake_helper import get_cmd_call

TIME_V = config.get("tools", dict()).get("time", "time")

DEBUG = config.get("debugmode", False)

SINGLE_CELL_MODE = config.get("single_cell_mode", False)
SKIP_NORMALIZATION = config.get("no_normalization", False)

# set up i/o
OUTPUTDIR = config["out_dir"]
QC_OUTDIR = join(OUTPUTDIR, "qc")
QC_LOGDIR = join(QC_OUTDIR, "log")
FASTQC_DIR = join(QC_OUTDIR, "fastqc")
BBDUK_DIR = join(QC_OUTDIR, "bbduk")
KAT_DIR = join(QC_OUTDIR, "kat")
BBNORM_DIR = join(QC_OUTDIR, "bbnorm")
TADPOLE_DIR = join(QC_OUTDIR, "tadpole")

INPUTFILES = Samplesheet(config["samplesheet"])

# generate target list
TARGETS = list()
TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbduk", s, s + "_R1.bbduk_fastqc.html"), INPUTFILES))
if not SINGLE_CELL_MODE:
	TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbduk", s, s + "_R2.bbduk_fastqc.html"), INPUTFILES))
TARGETS.extend(map(lambda s:join(TADPOLE_DIR, s, s + "_tadpole_contigs.fasta"), INPUTFILES))
TARGETS.extend(map(lambda s:join(KAT_DIR, s, s + ".dist_analysis.json"), INPUTFILES))

# define normalization/no-normalization behaviour
if SKIP_NORMALIZATION:
	PRIMARY_READDIR = BBDUK_DIR
	SECONDARY_READDIR = BBDUK_DIR
	PRIMARY_READ_ID = "bbduk"
	SECONDARY_READ_ID = "bbduk"
else:
	PRIMARY_READDIR = BBNORM_DIR
	SECONDARY_READDIR = BBDUK_DIR
	PRIMARY_READ_ID = "bbnorm"
	SECONDARY_READ_ID = "bbduk"
	TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbnorm", s, s + "_R1.bbnorm_fastqc.html"), INPUTFILES))
	if not SINGLE_CELL_MODE:
		TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbnorm", s, s + "_R2.bbnorm_fastqc.html"), INPUTFILES))

if DEBUG:
	with open("inputfiles.txt", "w") as input_out:
		print(*INPUTFILES.values(), sep="\n", file=input_out)
	with open("targets.txt", "w") as targets_out:
		print(*TARGETS, sep="\n", file=targets_out)

# helper
def get_sample_files(wc):
	return INPUTFILES[wc.sample].R1, INPUTFILES[wc.sample].R2

def get_r1(wc):
	return INPUTFILES[wc.sample].R1
def get_r2(wc):
	return INPUTFILES[wc.sample].R2

CMD_CALL = get_cmd_call(config, "bgrrl_container")


### RULES ###

localrules: all

rule all:
	input: TARGETS


bbduk_command = TIME_V + \
	" {{params.cmd}}" + \
	" -Xmx30g t={{threads}}" + \
	" in1={{input[0]}}{}" + \
	" out1={{output[0]}}{}" + \
	" ref={{params.adapters}}" + \
	" {{params.bbduk_params}}" + \
	" &> {{log}}"

tadpole_command = TIME_V + \
	" {{params.cmd}}" + \
	" -Xmx30g threads={{threads}}" + \
	" in={{input.r1}}{}" + \
	" out={{output.contigs}} &> {{log}}"

katgcp_command = " ({{params.cmd}} gcp -o {{params.prefix}} -t {{threads}} -v {{input.r1}}{} || touch {{output.katgcp}}) &> {{log}}"


if SINGLE_CELL_MODE:
	rule qc_bbduk:
		message:                                                           		
			"Preprocessing read data with bbduk..."
		input:
			get_r1
		output:
			r1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
		log:
			join(QC_LOGDIR, "{sample}", "{sample}.qc_bbduk.log")
		threads:
			8
		params:
			cmd = CMD_CALL + "bbduk.sh",
			adapters = config["resources"]["bb_adapters"],
			bbduk_params = config["params"]["bbduk"]
		shell:
			bbduk_command.format("", "")

	rule qc_tadpole:
		message:
			"Generating survey assemblies with tadpole..."
		input:
			r1 = join(PRIMARY_READDIR, "{sample}", "{sample}_R1." + PRIMARY_READ_ID + ".fastq.gz"),
		output:
			contigs = join(TADPOLE_DIR, "{sample}", "{sample}_tadpole_contigs.fasta")
		params:
			cmd = CMD_CALL + "tadpole.sh"
		log:
			join(QC_LOGDIR, "{sample}", "{sample}.qc_tadpole.log")
		threads:
			8
		shell:
			tadpole_command.format("")

	rule qc_katgcp:
		message:
			"Analyzing k-mer distribution, GC content and estimated genome size with kat..."
		input:
			r1 = join(PRIMARY_READDIR, "{sample}", "{sample}_R1." + PRIMARY_READ_ID + ".fastq.gz"),
		output:
			katgcp = join(KAT_DIR, "{sample}", "{sample}.dist_analysis.json")
		log:
			join(KAT_DIR, "{sample}", "{sample}.kat")
		params:
			prefix = lambda wildcards: join(KAT_DIR, wildcards.sample, wildcards.sample), 
			cmd = CMD_CALL + "kat"
		threads:
			2
		shell:
			katgcp_command.format("")

else:
	rule qc_bbduk:
		message:
			"Preprocessing read data with bbduk..."
		input:
			get_sample_files
		output:
			r1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
			r2 = join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
		log:
			join(QC_LOGDIR, "{sample}", "{sample}.qc_bbduk.log")
		threads:
			8
		params:
			cmd = CMD_CALL + "bbduk.sh",
			adapters = config["resources"]["bb_adapters"],
			bbduk_params = config["params"]["bbduk"]
		shell:
			bbduk_command.format(" in2={input[1]}", " out2={output.r2}")

	rule qc_tadpole:
		message:
			"Generating survey assemblies with tadpole..."
		input:
			r1 = join(PRIMARY_READDIR, "{sample}", "{sample}_R1." + PRIMARY_READ_ID + ".fastq.gz"),
			r2 = join(PRIMARY_READDIR, "{sample}", "{sample}_R2." + PRIMARY_READ_ID + ".fastq.gz")
		output:
			contigs = join(TADPOLE_DIR, "{sample}", "{sample}_tadpole_contigs.fasta")
		params:
			cmd = CMD_CALL + "tadpole.sh"
		log:
			join(QC_LOGDIR, "{sample}", "{sample}.qc_tadpole.log")
		threads:
			8
		shell:
			tadpole_command.format(" in2={input.r2}")
	
	rule qc_katgcp:
		message:
			"Analyzing k-mer distribution, GC content and estimated genome size with kat..."
		input:
			r1 = join(PRIMARY_READDIR, "{sample}", "{sample}_R1." + PRIMARY_READ_ID + ".fastq.gz"),
			r2 = join(PRIMARY_READDIR, "{sample}", "{sample}_R2." + PRIMARY_READ_ID + ".fastq.gz")
		output:
			katgcp = join(KAT_DIR, "{sample}", "{sample}.dist_analysis.json")
		log:
			join(KAT_DIR, "{sample}", "{sample}.kat")
		params:
			prefix = lambda wildcards: join(KAT_DIR, wildcards.sample, wildcards.sample), 
			cmd = CMD_CALL + "kat"
		threads:
			2
		shell:
			katgcp_command.format(" {input.r2}")


rule qc_fastqc_bbduk:
	message:
		"Generating post-preprocessing report with FastQC..."
	input:
		join(BBDUK_DIR, "{sample}", "{sample}_{mate}.bbduk.fastq.gz")
	output:
		fqc = join(FASTQC_DIR, "bbduk", "{sample}", "{sample}_{mate}.bbduk_fastqc.html")
	params:
		outdir = join(FASTQC_DIR, "bbduk", "{sample}"),
		cmd = CMD_CALL + "fastqc" 
	log:
		join(QC_LOGDIR, "{sample}", "{sample}_{mate}.qc_fastqc_bbduk.log")
	threads:
		2
	shell:
		"(" + TIME_V + " {params.cmd}" + \
		" --extract --threads={threads} --outdir={params.outdir} {input} " + \
		" || mkdir -p {params.outdir} && touch {output.fqc}) &> {log}"

"""
BB: 13th June 2019
https://twitter.com/BBToolsBio/status/1139272120125386752
In my experience, BBDuk for adapter removal, then Tadpole error correction, then BBMerge (with rsem flag), then Spades with no error correction, enables the fastest, lowest-memory, and most accurate assemblies.  For single cells BBNorm is suggested also, but not for isolates.
"""

#rule qc_tadpole_error_correction:
#	message:
#		"Performing tadpole error correction on reads..."
#	input:
#		r1 = rules.qc_bbduk.output[0],
#		r2 = rules.qc_bbduk.output[1]
#	output:
#		r1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.corr.fastq.gz"),
#		r2 = join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.corr.fastq.gz")
#	log:
#		join(QC_LOGDIR, "{sample}", "{sample}.qc_tadpole_errc.log")
#	threads:
#		8
#	params:
#		cmd = CMD_CALL + "tadpole.sh"
#	shell:
#		TIME_V + " {params.cmd}" + \
#		" -Xmx30g threads={threads} in={input.r1} in2={input.r2} out={output[0]} out2={output[1]} mode=correct &> {log}"
#
#rule qc_bbmerge:
#	message:
#		"Merging reads with bbmerge..."
#	input:
#		r1 = rules.qc_tadpole_error_correction.output[0],
#		r2 = rules.qc_tadpole_error_correction.output[0]
#	output:
#		# "bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>"
#		merged = join(BBDUK_DIR, "{sample}", "{sample}.merged.fastq.gz"),
#		ur1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.corr.unmerged.fastq.gz"),
#		ur2 = join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.corr.unmerged.fastq.gz")
#	log:
#		join(QC_LOGDIR, "{sample}", "{sample}.qc_bbmerge.log")
#	threads:
#		8
#	params:
#		cmd = CMD_CALL + "bbmerge.sh",
#		extend2 = 0 # ????
#	shell:
#		TIME_V + " {params.cmd}" + \
#		" -Xmx30g threads={threads} in={input.r1} in2={input.r2} out={output.merged} outu1={output.ur1} outu2={output.ur2} rsem=t extend2={params.extend2}"



if not SKIP_NORMALIZATION:

	bbnorm_command = TIME_V + " {{params.cmd}}" + \
		" -Xmx30g t={{threads}} in={{input.r1}}{}" + \
		" out={{output.r1}}{}" + \
		" {{params.bbnorm_params}}" + \
		" khist={{output.prehist}} khistout={{output.posthist}} &> {{log}}"

	if SINGLE_CELL_MODE:
		rule qc_bbnorm:
			message:
				"Normalizing read data with bbnorm..."
			input:
				r1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
			output:
				r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
				prehist = join(BBNORM_DIR, "{sample}", "{sample}.bbnorm.pre.hist"),
				posthist = join(BBNORM_DIR, "{sample}", "{sample}.bbnorm.post.hist")
			params:
				cmd = CMD_CALL + "bbnorm.sh",
				bbnorm_params = config["params"]["bbnorm"]
			log:
				join(QC_LOGDIR, "{sample}", "{sample}.qc_bbnorm.log")
			threads:
				8
			shell:
				bbnorm_command.format("", "")

	else:
		rule qc_bbnorm:
			message:
				"Normalizing read data with bbnorm..."
			input:
				r1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
				r2 = join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
			output:
				r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
				r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz"),
				prehist = join(BBNORM_DIR, "{sample}", "{sample}.bbnorm.pre.hist"),
				posthist = join(BBNORM_DIR, "{sample}", "{sample}.bbnorm.post.hist")
			params:
				cmd = CMD_CALL + "bbnorm.sh",
				bbnorm_params = config["params"]["bbnorm"]
			log:
				join(QC_LOGDIR, "{sample}", "{sample}.qc_bbnorm.log")
			threads:
				8
			shell:
				bbnorm_command.format(" in2={input.r2}", " out2={output.r2}")


	rule qc_fastqc_bbnorm:
		message:
			"Generating post-normalization report with FastQC..."
		input:
			join(BBNORM_DIR, "{sample}", "{sample}_{mate}.bbnorm.fastq.gz")
		output:
			fqc = join(FASTQC_DIR, "bbnorm", "{sample}", "{sample}_{mate}.bbnorm_fastqc.html")
		params:
			outdir = join(FASTQC_DIR, "bbnorm", "{sample}"),
			cmd = CMD_CALL + "fastqc"
		log:
			join(QC_LOGDIR, "{sample}", "{sample}_{mate}.qc_fastqc_bbnorm.log")
		threads:
			2
		shell:
			"({params.cmd}" + \
			" --extract --threads={threads} --outdir={params.outdir} {input}" + \
			" || mkdir -p {params.outdir} && touch {output.fqc}) &> {log}"

#rule qc_tadpole:
#	message:
#		"Generating survey assemblies with tadpole..."
#	input:
#		r1 = join(PRIMARY_READDIR, "{sample}", "{sample}_R1." + PRIMARY_READ_ID + ".fastq.gz"),
#		r2 = join(PRIMARY_READDIR, "{sample}", "{sample}_R2." + PRIMARY_READ_ID + ".fastq.gz")
#	output:
#		contigs = join(TADPOLE_DIR, "{sample}", "{sample}_tadpole_contigs.fasta")
#	params:
#		cmd = CMD_CALL + "tadpole.sh"
#	log:
#		join(QC_LOGDIR, "{sample}", "{sample}.qc_tadpole.log")
#	threads:
#		8
#	shell:
#		TIME_V + " {params.cmd}" + \
#		" -Xmx30g threads={threads} in={input.r1} in2={input.r2} out={output.contigs} &> {log}"
#
#rule qc_katgcp:
#	message:
#		"Analyzing k-mer distribution, GC content and estimated genome size with kat..."
#	input:
#		r1 = join(PRIMARY_READDIR, "{sample}", "{sample}_R1." + PRIMARY_READ_ID + ".fastq.gz"),
#		r2 = join(PRIMARY_READDIR, "{sample}", "{sample}_R2." + PRIMARY_READ_ID + ".fastq.gz")
#	output:
#		katgcp = join(KAT_DIR, "{sample}", "{sample}.dist_analysis.json")
#	log:
#		join(KAT_DIR, "{sample}", "{sample}.kat")
#	params:
#		prefix = lambda wildcards: join(KAT_DIR, wildcards.sample, wildcards.sample), 
#		cmd = CMD_CALL + "kat"
#	threads:
#		2
#	shell:
#		" ({params.cmd} gcp -o {params.prefix} -t {threads} -v {input.r1} {input.r2} || touch {output.katgcp}) &> {log}"
#
