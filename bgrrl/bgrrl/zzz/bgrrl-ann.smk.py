import sys
import csv
import os
from os.path import join, basename, dirname

from bgrrl.bgrrl import readSamplesheet
from bgrrl import loadPreCmd, TIME_CMD

OUTPUTDIR = config["out_dir"]
ASSEMBLY_DIR = join(OUTPUTDIR, "assembly")
ANNOTATION_DIR = join(OUTPUTDIR, "annotation")

PROKKA_WRAPPER = join(config["etc"], "wrappers", "prokka_wrapper")

INPUTFILES = dict(readSamplesheet(config["samplesheet"]))
with open("ann-inputfiles.txt", "w") as input_out:
    print(*INPUTFILES.values(), sep="\n", file=input_out)

TARGETS = list()
TARGETS.extend(map(lambda s:join(ANNOTATION_DIR, s, s + ".log"), INPUTFILES))

print("CONFIG")
print(config)

with open("ann-targets.txt", "w") as targets_out:
	print(*TARGETS, sep="\n", file=targets_out)

localrules: all

rule all:
	input: TARGETS

rule ann_prokka:
	input:
		contigs = join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	output:
		log = join(ANNOTATION_DIR, "{sample}", "{sample}.log"),
		faa = join(join(config["cwd"]), ANNOTATION_DIR, "{sample}", "{sample}.faa"),
		ffn = join(join(config["cwd"]), ANNOTATION_DIR, "{sample}", "{sample}.ffn")
	log:
		join(ANNOTATION_DIR, "log", "{sample}_ann_prokka.log")
	params:
		outdir = lambda wildcards: join(ANNOTATION_DIR, wildcards.sample),
		prefix = lambda wildcards: wildcards.sample
	threads:
		8
	shell:
		PROKKA_WRAPPER + \
		" {params.outdir} {params.prefix} {input.contigs} {log} {threads}"

