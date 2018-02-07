import sys
import csv
import os
from os.path import join, basename, dirname

from bgrrl.bgrrl import readSamplesheet
from bgrrl import loadPreCmd, TIME_CMD


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
with open("asm-inputfiles.txt", "w") as input_out:
    print(*INPUTFILES.values(), sep="\n", file=input_out)

TARGETS = list()
TARGETS.extend(map(lambda s:join(config["cwd"], ASSEMBLY_DIR, s, s + ".assembly.fasta"), INPUTFILES))
# ????
#if config.get("use_asm_lengthfilter", ""):
# TARGETS.extend(map(lambda s:join(ASSEMBLY_DIR, s, s + ".assembly.fasta"), INPUTFILES))

print("CONFIG")
print(config)

print(config.get("use_asm_length_filter", "NA"))

with open("asm-targets.txt", "w") as targets_out:
	print(*TARGETS, sep="\n", file=targets_out)

localrules: all, asm_link_assembly

rule all:
	input: TARGETS

if True:
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
			# load = loadPreCmd(config["load"]["unicycler"])
		threads:
			8
		shell:
			ASM_WRAPPER + \
			" {params.assembler} {input.r1} {input.r2} {threads} {params.outdir} {input.ur1} {input.ur2} {log}"
			#" &> {log}"
			# UNICYCLER_WRAPPER + \
			# " -1 {input.r1} -2 {input.r2} -t {threads} -o {params.outdir} &> {log}" 
			# + \
			# " && ln -s {params.assembly} {params.final_assembly} &> {log}"
"""
if False:
	rule asm_assembly:
		input:
			r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
			r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz")
		output:
			join(ASSEMBLY_DIR, "{sample}", "assembly.fasta")
		log:
			join(ASSEMBLY_DIR, "log", "{sample}.asm_assembly.log")
"""

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
		# original_assembly = join(os.path.dirname(input.assembly), "assembly.fasta")
		# shutil.copy2(original_assembly, original_assembly + ".full")
		# shutil.copy2(output[0], wildcards.sample + "." + original_assembly)
