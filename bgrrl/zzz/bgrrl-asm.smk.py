import sys
import csv
import os
from os.path import join, basename, dirname

from bgrrl import TIME_CMD
from bgrrl.samplesheet import readSamplesheet, Samplesheet, ASM_Sample
from bgrrl.snakemake_helper import loadPreCmd

DEBUG = config.get("debugmode", False)

SKIP_NORMALIZATION = config.get("no_normalization", False)

# set up i/o
OUTPUTDIR = config["out_dir"]
ASSEMBLY_DIR = join(OUTPUTDIR, "assembly")
QC_OUTDIR = join(OUTPUTDIR, "qc")
BBDUK_DIR = join(QC_OUTDIR, "bbduk")
BBNORM_DIR = join(QC_OUTDIR, "bbnorm")

INPUTFILES = Samplesheet(config["samplesheet"], sampletype=ASM_Sample)

# generate target list
TARGETS = list()
TARGETS.extend(map(lambda s:join(config["cwd"], ASSEMBLY_DIR, s, s + ".assembly.fasta"), INPUTFILES))

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

if DEBUG:
	print("CONFIG")
	print(config)
	with open("asm-inputfiles.txt", "w") as input_out:
		print(*INPUTFILES.values(), sep="\n", file=input_out)
	with open("asm-targets.txt", "w") as targets_out:
		print(*TARGETS, sep="\n", file=targets_out)

# helper
def get_sample_files(wc):
	s = INPUTFILES[wc.sample]
	# default case is with normalization, i.e. if --no-normalization option isn't present, it should be "False"
	# skip_normalization = config.get("no_normalization", False)
	if SKIP_NORMALIZATION:
		print("skipping normalization")
		return s.R1trim, s.R2trim
	else:
		return s.R1norm, s.R1trim, s.R2norm, s.R2trim


### RULES ###

localrules: all

rule all:
	input: TARGETS

rule asm_assembly:
	message:
		"Generating assembly with " + config["assembler"] + "..."
	input:
		get_sample_files
	output:
		join(ASSEMBLY_DIR, "{sample}", "assembly.fasta")
	log:
		join(ASSEMBLY_DIR, "log", "{sample}.asm_assembly.log")
	params:
		outdir = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample),
		assembly = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample, "assembly.fasta"),
		assembler = config["assembler"],
		r1 = lambda wildcards: get_sample_files(wildcards)[0] if SKIP_NORMALIZATION else ",".join(get_sample_files(wildcards)[:2]),
		r2 = lambda wildcards: get_sample_files(wildcards)[1] if SKIP_NORMALIZATION else ",".join(get_sample_files(wildcards)[2:])
	threads:
		8
	shell:
		"set +u && source activate bgasm_env" + \
		" && asm_wrapper --threads {threads} {params.assembler} {params.r1} {params.r2} {params.outdir} &> {log}"

rule asm_postprocess:
	message:
		"Postprocessing assembly (minimum contig length={})...".format(config["asm_lengthfilter_contig_minlen"]) 
	input:
		assembly = join(ASSEMBLY_DIR, "{sample}", "assembly.fasta")
	output:
		filtered_assembly = join(config["cwd"], ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	log:
		join(ASSEMBLY_DIR, "log", "{sample}.asm_postprocess.log")
	params:
		minlen = int(config["asm_lengthfilter_contig_minlen"]),
		load = loadPreCmd(config["load"]["bbmap"]),
		reformat = "reformat.sh"
	threads:
		1
	shell:
		"{params.load}" + TIME_CMD + " {params.reformat}" + \
		" in={input.assembly} out={output[0]} minlength={params.minlen}"


#### ! UNDER DEVELOPMENT ####

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
