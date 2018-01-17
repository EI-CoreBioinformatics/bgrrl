import sys
import csv
import os
import glob
from os.path import join, basename, dirname

from bgrrl.bgrrl import readSamplesheet
from bgrrl import loadPreCmd, TIME_CMD

SOFTWAREPATH = "/tgac/software/testing"
BBSUITE_DIR = os.path.join(SOFTWAREPATH, "bbmap", "37.24", "bbmap")
ADAPTERS = os.path.join(BBSUITE_DIR, "resources", "adapters.fa")

# tools
# print("THIS:", __file__)
UNICYCLER_WRAPPER = join(config["etc"], "wrappers", "unicycler_wrapper")
BUSCO_INIT_DIR = join(config["etc"], "util", "busco_init_dir")
BUSCO_DATA = "/tgac/workarea/group-pb/schudomc_bact/bact-grrl/data/busco/bacteria_odb9"


OUTPUTDIR = config["out_dir"]
ASSEMBLY_DIR = join(OUTPUTDIR, "assembly")
QC_OUTDIR = os.path.join(OUTPUTDIR, "qc")
FASTQC_DIR = os.path.join(QC_OUTDIR, "fastqc")
BBDUK_DIR = os.path.join(QC_OUTDIR, "bbduk")
KAT_DIR = os.path.join(QC_OUTDIR, "kat")
BBNORM_DIR = os.path.join(QC_OUTDIR, "bbnorm")
TADPOLE_DIR = os.path.join(QC_OUTDIR, "tadpole")

INPUTFILES = dict(readSamplesheet(config["samplesheet"]))
with open("asm-inputfiles.txt", "w") as input_out:
    print(*INPUTFILES.values(), sep="\n", file=input_out)

TARGETS = list()
TARGETS.extend(map(lambda s:join(ASSEMBLY_DIR, s, s + ".assembly.fasta"), INPUTFILES))
TARGETS.extend(map(lambda s:join(ASSEMBLY_DIR, s, "quast", "quast.log"), INPUTFILES))
# TARGETS.extend(map(lambda s:join(ASSEMBLY_DIR, s, "busco", "run_geno", "short_summary_geno.txt"), INPUTFILES))

with open("asm-targets.txt", "w") as targets_out:
	print(*TARGETS, sep="\n", file=targets_out)

rule all:
	input: TARGETS

rule asm_assembly:
	input:
		r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
		r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz")
	output:
		join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	log:
		join(ASSEMBLY_DIR, "log", "{sample}.asm_assembly.log")
	params:
		outdir = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample),
		assembly = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample, "assembly.fasta"),
		final_assembly = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample, wildcards.sample + ".assembly.fasta")
		# load = loadPreCmd(config["load"]["unicycler"])
	threads:
		8
	shell:
		UNICYCLER_WRAPPER + \
		" -1 {input.r1} -2 {input.r2} -t {threads} -o {params.outdir}" + \
		" && cp {params.assembly} {params.final_assembly} &> {log}"
"""
rule qa_busco_geno:
	input:
		assembly = join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	output:
		join(ASSEMBLY_DIR, "{sample}", "busco", "run_geno", "short_summary_geno.txt")
	log:	
		join(ASSEMBLY_DIR, "log", "{sample}_busco_geno.log")
	params:
		outdir = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample, "busco", "run_geno"),
		tmp = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample, "busco", "tmp"),
		load = loadPreCmd(config["load"]["busco"])
	threads:
		8
	shell:
		BUSCO_INIT_DIR + " {params.outdir} &&" + \
		" {params.load}" + TIME_CMD + \
		" run_BUSCO.py -i {input.assembly} -c {threads} -m geno" + \
		" --force -t {params.tmp} -l " + BUSCO_DATA + " -o geno &> {log}"
"""
rule qa_quast:
        input:
                assembly = os.path.join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
        output:
                os.path.join(ASSEMBLY_DIR, "{sample}", "quast", "quast.log")
        params:
                load = loadPreCmd(config["load"]["quast"]),
                outdir = lambda wildcards: os.path.join(ASSEMBLY_DIR, wildcards.sample, "quast"),
                cmd = loadPreCmd(config["cmd"]["quast"], is_dependency=False)
        log:
                os.path.join(ASSEMBLY_DIR, "log", "{sample}.asm_quast_assembly.log")
        threads:
                2
        shell:
                "{params.load}" + TIME_CMD + \
                " {params.cmd} -o {params.outdir} -t {threads} -L -s {input.assembly} --min-contig 1000 &> {log}"

