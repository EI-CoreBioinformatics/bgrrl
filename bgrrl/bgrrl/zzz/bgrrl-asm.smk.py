import sys
import csv
import os
from os.path import join, basename, dirname

from bgrrl.bgrrl import readSamplesheet
from bgrrl import loadPreCmd, TIME_CMD

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

INPUTFILES = dict(readSamplesheet(config["samplesheet"]))
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

rule asm_assembly:
	input:
		r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
		r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz"),
		ur1 = join(BBDUK_DIR, "{sample}", "{sample}_R1.bbduk.fastq.gz"),
		ur2 = join(BBDUK_DIR, "{sample}", "{sample}_R2.bbduk.fastq.gz")
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
		" {params.assembler} {input.r1} {input.r2} {threads} {params.outdir} {input.ur1} {input.ur2} {log}"

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
