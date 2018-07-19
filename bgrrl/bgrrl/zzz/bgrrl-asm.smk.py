import sys
import csv
import os
from os.path import join, basename, dirname

from bgrrl.samplesheet import readSamplesheet, Samplesheet, ASM_Sample
from bgrrl import TIME_CMD
from eicore.external_process.snakemake_helper import loadPreCmd

DEBUG = config.get("debugmode", False)

# tools
UNICYCLER_WRAPPER = join(config["etc"], "wrappers", "unicycler_wrapper")
ASM_WRAPPER = join(config["etc"], "wrappers", "asm_wrapper")

# set up i/o
OUTPUTDIR = config["out_dir"]
ASSEMBLY_DIR = join(OUTPUTDIR, "assembly")
QC_OUTDIR = os.path.join(OUTPUTDIR, "qc")













BBDUK_DIR = os.path.join(QC_OUTDIR, "bbduk")
BBNORM_DIR = os.path.join(QC_OUTDIR, "bbnorm")


if not config.get("no_normalization", False):
	PRIMARY_READDIR = BBNORM_DIR
	SECONDARY_READDIR = BBDUK_DIR
	PRIMARY_READ_ID = "bbnorm"
	SECONDARY_READ_ID = "bbduk"
else:
	PRIMARY_READDIR = BBDUK_DIR
	SECONDARY_READDIR = BBDUK_DIR
	PRIMARY_READ_ID = "bbduk"
	SECONDARY_READ_ID = "bbduk"

# INPUTFILES = dict(readSamplesheet(config["samplesheet"]))
INPUTFILES = Samplesheet(config["samplesheet"], sampletype=ASM_Sample)
TARGETS = list()
TARGETS.extend(map(lambda s:join(config["cwd"], ASSEMBLY_DIR, s, s + ".assembly.fasta"), INPUTFILES))

if DEBUG:
	print("CONFIG")
	print(config)
	with open("asm-inputfiles.txt", "w") as input_out:
		print(*INPUTFILES.values(), sep="\n", file=input_out)
	with open("asm-targets.txt", "w") as targets_out:
		print(*TARGETS, sep="\n", file=targets_out)

localrules: all

rule all:
	input: TARGETS

# TODO: Adjust asm_wrapper so it deals with no_normalization cases in a different way than just running it twice.
def get_sample_files(wc):
	s = INPUTFILES[wc.sample]
	if not config.get("no_normalization", False):
		return s.R1norm, s.R2norm, s.R1trim, s.R2trim
	else:
		return s.R1norm, s.R2norm, s.R1norm, s.R2norm        


rule asm_assembly:
	input:
		get_sample_files
		#r1 = join(PRIMARY_READDIR, "{sample}", "{sample}" + "_R1.{}.fastq.gz".format(PRIMARY_READ_ID)),
		#r2 = join(PRIMARY_READDIR, "{sample}", "{sample}" + "_R2.{}.fastq.gz".format(PRIMARY_READ_ID)),
		#ur1 = join(SECONDARY_READDIR, "{sample}", "{sample}" + "_R1.{}.fastq.gz".format(SECONDARY_READ_ID)),
		#ur2 = join(SECONDARY_READDIR, "{sample}", "{sample}" + "_R2.{}.fastq.gz".format(SECONDARY_READ_ID)),
	output:
		join(ASSEMBLY_DIR, "{sample}", "assembly.fasta")
	log:
		join(ASSEMBLY_DIR, "log", "{sample}.asm_assembly.log")
	params:
		outdir = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample),
		assembly = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample, "assembly.fasta"),
		final_assembly = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample, wildcards.sample + ".assembly.fasta"),
		assembler = config["assembler"]
	threads:
		8
	shell:
		ASM_WRAPPER + \
		" {params.assembler} {input[0]} {input[1]} {threads} {params.outdir} {input[2]} {input[3]} &> {log}"
		# " {params.assembler} {input.r1} {input.r2} {threads} {params.outdir} {input.ur1} {input.ur2} &> {log}"

rule asm_lengthfilter:
	input:
		assembly = os.path.join(ASSEMBLY_DIR, "{sample}", "assembly.fasta")
	output:
		filtered_assembly = os.path.join(config["cwd"], ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	log:
		join(ASSEMBLY_DIR, "log", "{sample}.asm_lengthfilter.log")
	params:
		minlen = int(config["asm_lengthfilter_contig_minlen"])
	threads:
		1
	run:
		import shutil
		from ktio.ktio import readFasta
		with open(output[0], "w") as seqout:
			for _id, _seq in readFasta(input.assembly):
				if len(_seq) >= params.minlen:
					print(_id, _seq, sep="\n", file=seqout)

if config["reapr_correction"]:

	rule asm_reapr_readlen:
		input:
			assembly = join(ASSEMBLY_DIR, "{sample}", "assembly.fasta"),     		
			r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
			r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz")
		output:
			r1 = join(ASSEMBLY_DIR, "{sample}", "reapr_R1.fastq.gz"),
			r2 = join(ASSEMBLY_DIR, "{sample}", "reapr_R2.fastq.gz") 			
		log:
			join(ASSEMBLY_DIR, "log", "{sample}.asm_reapr_readlen.log")
		params:
			outdir = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample),
			load = loadPreCmd(config["load"]["bbmap"])
		threads:
			2
		shell:
			""

	rule asm_reapr_isize:
		input:
			assembly = join(ASSEMBLY_DIR, "{sample}", "assembly.fasta"),
			r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
                        r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz")
		output:
			isize = join(ASSEMBLY_DIR, "{sample}", "INSERT_SIZE")
		log:
			join(ASSEMBLY_DIR, "log", "{sample}.asm_reapr_isize.log")
		params:
			outdir = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample),
			load = loadPreCmd(config["load"]["bwa"])
		threads:
			8
		shell:
			""

	rule asm_reapr:
		input:
			assembly = join(ASSEMBLY_DIR, "{sample}", "assembly.fasta")
		output:
			reapr_assembly = join(ASSEMBLY_DIR, "{sample}", "assembly.pre_reapr.fasta")
		log:
			join(ASSEMBLY_DIR, "log", "{sample}.asm_reapr.log")
		params:
			load = loadPreCmd(config["load"]["reapr"])
		threads:
			8
		shell:
			"{params.load}" + \
			" (reapr facheck {input.assembly} ||" + \
			"  (reapr facheck {input.assembly} {input.assembly}.reapr_in &&" + \
			"   mv {input.assembly} {input.assembly}.no_reapr &&" + \
			"   mv {input.assembly}.reapr_in {input.assembly})) &&" + \
			" reapr perfectmap {input.assembly} {input.r1} {input.r2} {params.isize} {params.prefix_perfect} &&" + \
			" reapr smaltmap {input.assembly} {input.r1} {input.r2} {params.longbam} &&" + \
			" reapr pipeline {input.assembly} {params.longbam} {params.outdir} {params.prefix_perfect}"
