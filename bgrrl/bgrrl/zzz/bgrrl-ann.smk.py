import sys
import csv
import os
import glob
from os.path import join, basename, dirname

from bgrrl.bgrrl import readSamplesheet
from bgrrl import loadPreCmd, TIME_CMD

SOFTWAREPATH = "/tgac/software/testing"

# tools
BUSCO_INIT_DIR = join(config["etc"], "util", "busco_init_dir")
BUSCO_DATA = "/tgac/workarea/group-pb/schudomc_bact/bact-grrl/data/busco/bacteria_odb9"
# needed for BUSCO
CWD = os.getcwd()

OUTPUTDIR = config["out_dir"]
ASSEMBLY_DIR = join(OUTPUTDIR, "assembly")
ANNOTATION_DIR = join(OUTPUTDIR, "annotation")
QA_DIR = join(OUTPUTDIR, "qa")

PROKKA_WRAPPER = "/tgac/workarea/group-pb/schudomc_bact/bact-grrl/scripts/wrappers/prokka_wrapper"

INPUTFILES = dict(readSamplesheet(config["samplesheet"]))
with open("ann-inputfiles.txt", "w") as input_out:
    print(*INPUTFILES.values(), sep="\n", file=input_out)

TARGETS = list()
# TARGETS.extend(map(lambda s:join(config["cwd"], ASSEMBLY_DIR, s, s + ".assembly.fasta"), INPUTFILES))
TARGETS.extend(map(lambda s:join(config["cwd"], QA_DIR, "busco", "tran", s, "short_summary_{}.txt".format(s)), INPUTFILES))
TARGETS.extend(map(lambda s:join(config["cwd"], QA_DIR, "busco", "prot", s, "short_summary_{}.txt".format(s)), INPUTFILES))

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

rule qa_busco_tran:
	input:
		transcripts = join(config["cwd"], ANNOTATION_DIR, "{sample}", "{sample}.ffn")
	output:
		join(config["cwd"], QA_DIR, "busco", "tran", "{sample}", "short_summary_{sample}.txt")
	log:
		join(config["cwd"], QA_DIR, "log", "{sample}_busco_tran.log")
	params:
		outdir = lambda wildcards: join(config["cwd"], QA_DIR, "busco", "tran", "run_" + wildcards.sample),
		final_outdir = lambda wildcards: join(config["cwd"], QA_DIR, "busco", wildcards.sample),
		tmp = lambda wildcards: join(config["cwd"], QA_DIR, "busco", "tmp", "tran_" + wildcards.sample),
		load = loadPreCmd(config["load"]["busco"])
	threads:
		8
	shell:
		BUSCO_INIT_DIR + " {params.outdir} && cd {params.outdir} && cd .." + \
		" && {params.load}" + TIME_CMD + \
		" run_BUSCO.py -i {input.transcripts} -c {threads} -m tran" + \
		" --force -t {params.tmp} -l " + BUSCO_DATA + " -o {wildcards.sample}" + \
		" && cd " + CWD + \
		" && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
		" && rm -rf {params.outdir} &> {log}"

rule qa_busco_prot:
        input:
                proteins = join(config["cwd"], ANNOTATION_DIR, "{sample}", "{sample}.faa")
        output:
                join(config["cwd"], QA_DIR, "busco", "prot", "{sample}", "short_summary_{sample}.txt")
        log:
                join(config["cwd"], QA_DIR, "log", "{sample}_busco_prot.log")
        params:
                outdir = lambda wildcards: join(config["cwd"], QA_DIR, "busco", "prot", "run_" + wildcards.sample),
                final_outdir = lambda wildcards: join(config["cwd"], QA_DIR, "busco", "prot", wildcards.sample),
                tmp = lambda wildcards: join(config["cwd"], QA_DIR, "busco", "tmp", "prot_" + wildcards.sample),
                load = loadPreCmd(config["load"]["busco"])
        threads:
                8
        shell:
                BUSCO_INIT_DIR + " {params.outdir} && cd {params.outdir} && cd .." + \
                " && {params.load}" + TIME_CMD + \
                " run_BUSCO.py -i {input.proteins} -c {threads} -m prot" + \
                " --force -t {params.tmp} -l " + BUSCO_DATA + " -o {wildcards.sample}" + \
                " && cd " + CWD + \
                " && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
                " && rm -rf {params.outdir}" + \
                " &> {log}"
