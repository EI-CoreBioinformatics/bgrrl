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
UNICYCLER_WRAPPER = join(config["etc"], "wrappers", "unicycler_wrapper")
BUSCO_INIT_DIR = join(config["etc"], "util", "busco_init_dir")
BUSCO_DATA = "/tgac/workarea/group-pb/schudomc_bact/bact-grrl/data/busco/bacteria_odb9"
# needed for BUSCO
CWD = os.getcwd()

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
TARGETS.extend(map(lambda s:join(config["cwd"], ASSEMBLY_DIR, s, "busco", "run_geno", "short_summary_geno.txt"), INPUTFILES))
TARGETS.extend(map(lambda s:join(ASSEMBLY_DIR, s, "blobtools", "bwa", s + ".blob_bwa.bam"), INPUTFILES))
TARGETS.extend(map(lambda s:join(ASSEMBLY_DIR, s, "blobtools", "blast", s + ".blob_blast.tsv"), INPUTFILES))


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

rule qa_busco_geno:
	input:
		assembly = join(config["cwd"], ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	output:
		join(config["cwd"], ASSEMBLY_DIR, "{sample}", "busco", "run_geno", "short_summary_geno.txt")
	log:	
		join(config["cwd"], ASSEMBLY_DIR, "log", "{sample}_busco_geno.log")
	params:
		outdir = lambda wildcards: join(config["cwd"], ASSEMBLY_DIR, wildcards.sample, "busco", "run_geno"),
		tmp = lambda wildcards: join(config["cwd"], ASSEMBLY_DIR, wildcards.sample, "busco", "tmp"),
		load = loadPreCmd(config["load"]["busco"])
	threads:
		8
	shell:
		BUSCO_INIT_DIR + " {params.outdir} && cd {params.outdir} && cd .. && " + \
		" {params.load}" + TIME_CMD + \
		" run_BUSCO.py -i {input.assembly} -c {threads} -m geno" + \
		" --force -t {params.tmp} -l " + BUSCO_DATA + " -o geno &> {log} && cd " + CWD + " &> {log}"

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
                join(ASSEMBLY_DIR, "log", "{sample}.asm_quast_assembly.log")
        threads:
                2
        shell:
                "{params.load}" + TIME_CMD + \
                " {params.cmd} -o {params.outdir} -t {threads} -L -s {input.assembly} --min-contig 1000 &> {log}"

rule qa_blob_bwa_mem:
	input:
		r1 = join(BBNORM_DIR, "{sample}", "{sample}_R1.bbnorm.fastq.gz"),
		r2 = join(BBNORM_DIR, "{sample}", "{sample}_R2.bbnorm.fastq.gz"),
                assembly = os.path.join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	output:
		bam = join(ASSEMBLY_DIR, "{sample}", "blobtools", "bwa", "{sample}.blob_bwa.bam")
	log:
		join(ASSEMBLY_DIR, "log", "{sample}.blob_bwa.log")
	params:
		outdir = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample, "blobtools", "bwa"),
		index = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample, "blobtools", "bwa", wildcards.sample + ".assembly.fasta"),
		outbam = lambda wildcards: join(ASSEMBLY_DIR, wildcards.sample, "blobtools", "bwa", wildcards.sample + ".blob_bwa.bam"),
		load = loadPreCmd(config["load"]["blob_bwa"])
	threads:
		8
	shell:
		"{params.load}" + TIME_CMD + \
		" bwa index -p {params.index} {input.assembly} &&" + \
		" bwa mem -t {threads} {params.index} {input.r1} {input.r2} | samtools view -buSh - > {output.bam} 2> {log}"

rule qa_blob_blast:
	input:
                assembly = os.path.join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	output:
		tsv = join(ASSEMBLY_DIR, "{sample}", "blobtools", "blast", "{sample}.blob_blast.tsv")
	log:
		join(ASSEMBLY_DIR, "log", "{sample}.blob_blast.log")
	params:
		outdir = join(ASSEMBLY_DIR, "{sample}", "blobtools", "blast"),
		load = loadPreCmd(config["load"]["blob_blast"])
	threads:
		8
	shell:
		"{params.load}" + TIME_CMD + \
		" blastn -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 " + \
		" -num_threads {threads} -query {input.assembly} -db " + config["blob_blastdb"] + " -out {output.tsv} &> {log}"

rule qa_blobtools:
	input:
		bwa = join(ASSEMBLY_DIR, "{sample}", "blobtools", "bwa", "{sample}.blob_bwa.bam"),
		blast = join(ASSEMBLY_DIR, "{sample}", "blobtools", "blast", "{sample}.blob_blast.tsv"),
		assembly = os.path.join(ASSEMBLY_DIR, "{sample}", "{sample}.assembly.fasta")
	output:
		blobtable = join(ASSEMBLY_DIR, "{sample}", "blobtools", "blob", "{sample}.blobDB.table.txt")
	log:
		join(ASSEMBLY_DIR, "log", "{sample}.blobtools.log")
	params:
		prefix = lambda wildcards: join(ASSEMBLY_DIR, "{sample}", "blobtools", wildcards.sample, wildcards.sample),
		load = loadPreCmd(config["load"]["blobtools"])
	threads:
		1
	shell:
		"{params.load}" + TIME_CMD + \
		" blobtools create -t {input.blast} -b {input.bwa} -i {input.assembly} -o {params.prefix} -x bestsumorder &&" + \
		"blobtools view -i {params.prefix}.blobDB.json -o $(dirname {params.prefix}) -x bestsumorder -r species &&" + \
		"blobtools blobplot -r species -l 1000 -i {params.prefix}.blobDB.json -o $(dirname {params.prefix}) &> {log}"
