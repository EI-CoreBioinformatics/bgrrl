import sys
import csv
import os
from os.path import join

from bgrrl.bgrrl import readSamplesheet
from bgrrl import loadPreCmd, TIME_CMD

# set up i/o
OUTPUTDIR = config["out_dir"]
QC_OUTDIR = join(OUTPUTDIR, "qc")
FASTQC_DIR = join(QC_OUTDIR, "fastqc")
BBDUK_DIR = join(QC_OUTDIR, "bbduk")
KAT_DIR = join(QC_OUTDIR, "kat")
BBNORM_DIR = join(QC_OUTDIR, "bbnorm")
TADPOLE_DIR = join(QC_OUTDIR, "tadpole")

INPUTFILES = dict(readSamplesheet(config["samplesheet"]))
with open("inputfiles.txt", "w") as input_out:
    print(*INPUTFILES.values(), sep="\n", file=input_out)

TARGETS = list()
TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbduk", s, s + "_R1.bbduk_fastqc.html"), INPUTFILES))
TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbduk", s, s + "_R2.bbduk_fastqc.html"), INPUTFILES))
TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbnorm", s, s + "_R1.bbnorm_fastqc.html"), INPUTFILES))
TARGETS.extend(map(lambda s:join(FASTQC_DIR, "bbnorm", s, s + "_R2.bbnorm_fastqc.html"), INPUTFILES))
TARGETS.extend(map(lambda s:join(TADPOLE_DIR, s, s + "_tadpole_contigs.fasta"), INPUTFILES))
TARGETS.extend(map(lambda s:join(KAT_DIR, s, s + ".kat-gcp.mx"), INPUTFILES))

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
		join(QC_OUTDIR, "log", "{sample}.qc_bbduk.log")
	threads:
		8
	params:
		load = loadPreCmd(config["load"]["bbmap"]),
		bbduk = config["tools"]["bbduk"],
		adapters = config["resources"]["bb_adapters"]
	shell:
		"{params.load}" + TIME_CMD + " {params.bbduk}" + \
		" -Xmx30g t={threads} in1={input[0]} in2={input[1]} out1={output.r1} out2={output.r2}" + \
		" ref={params.adapters}" + \
		" ktrim=r k=21 mink=7 hdist=1 qtrim=lr trimq=3 minlen=100 maq=20 tpe tbo &> {log}"
		#" ktrim=r k=21 mink=11 hdist=2 qtrim=lr trimq=3 minlen=100 maq=20 tpe tbo &> {log}"

rule qc_fastqc_bbduk:
	input:
		join(BBDUK_DIR, "{sample}", "{fastq}.bbduk.fastq.gz")
	output:
		fqc = join(FASTQC_DIR, "bbduk", "{sample}", "{fastq}.bbduk_fastqc.html")
	params:
		outdir = join(FASTQC_DIR, "bbduk", "{sample}"),
                load = loadPreCmd(config["load"]["fastqc"]),
		fastqc = config["tools"]["fastqc"]
	log:
		join(QC_OUTDIR, "log", "{fastq}.qc_fastqc_bbduk.log")
	threads:
		2
	shell:
		"{params.load} " + TIME_CMD + " {params.fastqc}" + \
		" --extract --threads={threads} --outdir={params.outdir} {input} " + \
		" || mkdir -p {params.outdir} && touch {output.fqc} &> {log}"
		# " --extract --threads={threads} --outdir={params.outdir} {input}) || (mkdir -p {params.outdir} && touch {output.fqc}) &> {log}"

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
		bbnorm = config["tools"]["bbnorm"]
	log:
		join(QC_OUTDIR, "log", "{sample}.qc_bbnorm.log")
	threads:
		8
	shell:
		"{params.load}" + TIME_CMD + " {params.bbnorm}" + \
		" -Xmx30g t={threads} in={input.r1} in2={input.r2} out={output.r1} out2={output.r2}" + \
		" target=100 min=2 prefilter" + \
		" khist={output.prehist} khistout={output.posthist} &> {log}"

rule qc_fastqc_bbnorm:
	input:
		join(BBNORM_DIR, "{sample}", "{fastq}.bbnorm.fastq.gz")
	output:
		fqc = join(FASTQC_DIR, "bbnorm", "{sample}", "{fastq}.bbnorm_fastqc.html")
	params:
		outdir = join(FASTQC_DIR, "bbnorm", "{sample}"),
                load = loadPreCmd(config["load"]["fastqc"]),
		fastqc = config["tools"]["fastqc"]
	log:
		join(QC_OUTDIR, "log", "{fastq}.qc_fastqc_bbnorm.log")
	threads:
		2
	shell:
		"{params.load} ({params.fastqc}" + \
		" --extract --threads={threads} --outdir={params.outdir} {input}" + \
		" || mkdir -p {params.outdir} && touch {output.fqc}) &> {log}"
		# "{params.load} (" + TIME_CMD + " {params.fastqc}" + \
		# " --extract --threads={threads} --outdir={params.outdir} {input} && touch {output.fqc}) &> {log}"
		# " --extract --threads={threads} --outdir={params.outdir} {input}) || (mkdir -p {params.outdir} && touch {output.fqc}) &> {log}"

rule qc_tadpole:
	input:
		r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
		r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz")
	output:
		contigs = join(TADPOLE_DIR, "{sample}", "{sample}_tadpole_contigs.fasta")
	params:
		load = loadPreCmd(config["load"]["bbmap"]),
		tadpole = config["tools"]["tadpole"]
	log:
		join(QC_OUTDIR, "log", "{sample}.qc_tadpole.log")
	threads:
		8
	shell:
		"{params.load}" + TIME_CMD + " {params.tadpole}" + \
		" -Xmx30g threads={threads} in={input.r1} in2={input.r2} out={output.contigs} &> {log}"

"""
rule qa_quast:
	input:
		contigs = join(TADPOLE_DIR, "{sample}", "{sample}_tadpole_contigs.fasta")
	output:
		join(TADPOLE_DIR, "{sample}", "quast", "quast.log")
	params:
		load = loadPreCmd(config["load"]["quast"]),
		outdir = lambda wildcards: join(TADPOLE_DIR, wildcards.sample, "quast"),
		cmd = loadPreCmd(config["cmd"]["quast"], is_dependency=False)
	log:
		join(QC_OUTDIR, "log", "{sample}.qc_quast_tadpole.log")
	threads:
		2
	shell:
		"{params.load}" + " (" + TIME_CMD + \
		" {params.cmd} -o {params.outdir} -t {threads} -L -s {input.contigs} --min-contig 0) || (mkdir -p {params.outdir} && touch {output[0]}) &> {log}"
"""

rule qc_katgcp:
	input:
		r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
		r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz")
	output:
		katgcp = join(KAT_DIR, "{sample}", "{sample}.kat-gcp.mx")
	log:
		# join(QC_OUTDIR, "log", "{sample}.qc_katgcp.log")
		join(KAT_DIR, "{sample}", "{sample}.kat")
	params:
		prefix = lambda wildcards: join(KAT_DIR, wildcards.sample, wildcards.sample + ".kat-gcp"),
		load = loadPreCmd(config["load"]["kat"])
	threads:
		2
	shell:
		"{params.load}" + \
		" (kat gcp -o {params.prefix} -t {threads} -v {input.r1} {input.r2} || touch {output.katgcp}) &> {log}"

"""
rule qc_katda:
	input:
		katgcp = join(KAT_DIR, "{sample}", "{sample}.kat-gcp.mx")
	output:
		katda = join(KAT_DIR, "{sample}", "{sample}.kat-gcp.mx.dist")
	log:
		join(QC_OUTDIR, "log", "{sample}.qc_katda.log")
	params:
		load = loadPreCmd(config["load"]["kat"])
	threads:
		1
	shell:
		"{params.load}" + \                                                                    		
		" (kat_distanalysis {input.katgcp} || touch {output.katda}) > {output.katda} 2> {log}"
"""
