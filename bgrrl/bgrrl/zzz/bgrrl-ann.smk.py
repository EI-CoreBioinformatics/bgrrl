import sys
import csv
import os
from os.path import join, basename, dirname

from bgrrl.bgrrl import readSamplesheet
from bgrrl import loadPreCmd, TIME_CMD

OUTPUTDIR = config["out_dir"]
ASSEMBLY_DIR = join(OUTPUTDIR, "assembly")
ANNOTATION_DIR = join(OUTPUTDIR, "annotation")
PROKKA_DIR = join(ANNOTATION_DIR, "prokka")
RATT_DIR = join(ANNOTATION_DIR, "ratt")

PROKKA_WRAPPER = join(config["etc"], "wrappers", "prokka_wrapper")

INPUTFILES = dict(readSamplesheet(config["samplesheet"]))
with open("ann-inputfiles.txt", "w") as input_out:
    print(*INPUTFILES.values(), sep="\n", file=input_out)

TARGETS = list()
if config["run_prokka"]:
	TARGETS.extend(map(lambda s:join(PROKKA_DIR, s, s + ".log"), INPUTFILES))
if config["run_ratt"]:
	TARGETS.extend(map(lambda s:join(RATT_DIR, s, s + ".done"), INPUTFILES))

CWD = os.getcwd()

print("CONFIG")
print(config)

with open("ann-targets.txt", "w") as targets_out:
	print(*TARGETS, sep="\n", file=targets_out)

localrules: all

rule all:
	input: TARGETS

if config["run_prokka"]:
	rule ann_prokka:
		input:
			contigs = join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
		output:
			log = join(PROKKA_DIR, "{sample}", "{sample}.log"),
			faa = join(join(config["cwd"]), PROKKA_DIR, "{sample}", "{sample}.faa"),
			ffn = join(join(config["cwd"]), PROKKA_DIR, "{sample}", "{sample}.ffn")
		log:
			join(ANNOTATION_DIR, "log", "{sample}_ann_prokka.log")
		params:
			outdir = lambda wildcards: join(PROKKA_DIR, wildcards.sample),
			prefix = lambda wildcards: wildcards.sample
		threads:
			8
		shell:
			PROKKA_WRAPPER + \
			" {params.outdir} {params.prefix} {input.contigs} {log} {threads}"

if config["run_ratt"]:
	rule ann_ratt:
		input:
			contigs = join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
		output:
			done = join(RATT_DIR, "{sample}", "{sample}.done")
		log:
			join(config["cwd"], ANNOTATION_DIR, "log", "{sample}_ann_ratt.log")
		params:	
			outdir = lambda wildcards: join(RATT_DIR, wildcards.sample),
			load = loadPreCmd(config["load"]["ratt"]),
			prefix = lambda wildcards: wildcards.sample,
			reference = config["ratt_reference"]
		threads:
			1
		shell:
			"{params.load}" + \
			" (cd {params.outdir} && " + TIME_CMD + \
			" $RATT_HOME/start.ratt.sh {params.reference} " + \
			join(config["cwd"], "{input.contigs}") + " {params.prefix} Strain" + \
			" && touch " + join(config["cwd"], "{output.done}") + \
			" && extractfeat -sequence *.final.embl -outseq {params.prefix}.fasta -type CDS" + \
			" && transeq -sequence {params.prefix}.fasta {params.prefix}.pep.fasta" + \
			" && mv {params.prefix}.fasta {params.prefix}.ffn" + \
			" && mv {params.prefix}.pep.fasta {params.prefix}.faa" + \
			" && cd " + CWD + ") &> {log}"
			

